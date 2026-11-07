# Module D - stepped-wedge design diagnostics (Halkiewicz and Alexeev,
# "Concurrent controls and heterogeneity tolerance in stepped-wedge trials",
# arXiv:2610.35853).
#
# A design is a numeric vector of length T: cluster counts (or shares) for the
# start periods 2, ..., T followed by the wait-list (never treated during
# follow-up). Exposure time in treated cell (g, t) is e = t - g + 1, with
# L = T - 1 exposures. The target is the equal-weight TATE u = 1/L unless a
# target weight vector is supplied.
#
# Equal cluster-period sizes: period effects are removed by subtracting the
# share-weighted cross-sectional mean (Web Appendix A.2), so the information for
# the treatment parameters is
#   S = sum_g n_g X_g' Omega X_g - G Xbar' Omega Xbar,  Omega = V^{-1}
# (TWFE: Omega = I - J/T). Unequal sizes: cluster-level GLS with explicit
# period effects and V_c = diag(sigma^2/m_ct) + tau_g^2 I + tau_a^2 J (unit total
# variance), or cell-size-weighted TWFE.

.sw_check_design <- function(design) {
  if (!is.numeric(design) || length(design) < 3)
    stop("`design` must be numeric of length T >= 3: counts for start periods 2..T, then the wait-list")
  if (any(!is.finite(design)) || any(design < 0))
    stop("`design` entries must be finite and non-negative")
  if (design[1] <= 0)
    stop("the first switching sequence (start period 2) must be non-empty, so that every exposure 1..T-1 is observed")
  as.numeric(design)
}

.sw_X <- function(T, g, model) {
  t <- seq_len(T)
  if (model == "IT") return(matrix(as.numeric(t >= g), T, 1))
  L <- T - 1
  Z <- matrix(0, T, L)
  for (k in seq_len(T)) {
    e <- t[k] - g + 1
    if (e >= 1 && e <= L) Z[k, e] <- 1
  }
  Z
}

.sw_omega <- function(T, phi = NULL, icc = NULL, cac = 1, m = NULL,
                      estimator = "GLS") {
  Q <- diag(T) - matrix(1 / T, T, T)
  P <- matrix(1 / T, T, T)
  if (estimator == "TWFE") return(Q)
  if (!is.null(phi) && !is.null(icc))
    stop("supply either `phi` or (`icc`, `cac`, `m`), not both")
  if (!is.null(phi)) {
    if (!(phi >= 0)) stop("`phi` must be non-negative")
    lam <- if (is.infinite(phi)) 0 else 1 / (1 + T * phi)
    return(Q + lam * P)
  }
  if (is.null(icc) || is.null(m))
    stop("GLS needs `phi`, or `icc`, `cac` and `m`")
  php <- sw_phi_nested(icc, cac, m)
  c0 <- (1 - cac) * icc + (1 - icc) / m
  (Q + P / (1 + T * php)) / c0
}

# Equal-size engine. Returns S (information after removing period and cluster
# effects) and the per-sequence period-demeaned regressors.
.sw_equal_info <- function(design, Omega, model) {
  T <- length(design)
  G <- sum(design)
  gs <- 2:(T + 1)
  Xs <- lapply(gs, function(g) .sw_X(T, g, model))
  Xbar <- Reduce(`+`, Map(function(n, X) n * X, design, Xs)) / G
  k <- ncol(Xs[[1]])
  S <- matrix(0, k, k)
  for (i in seq_along(gs))
    if (design[i] > 0) S <- S + design[i] * crossprod(Xs[[i]], Omega %*% Xs[[i]])
  S <- S - G * crossprod(Xbar, Omega %*% Xbar)
  list(S = S, Xs = Xs, Xbar = Xbar, gs = gs)
}

.sw_solve <- function(S, what = "model") {
  rc <- tryCatch(rcond(S), error = function(e) 0)
  if (!is.finite(rc) || rc < 1e-12)
    stop("the design does not identify the ", what, " (singular information)")
  solve(S)
}

# Cluster-level engine for unequal cluster-period sizes.
.sw_cluster_fit <- function(design, sizes, icc, cac, estimator, model) {
  T <- length(design)
  if (any(abs(design - round(design)) > 0))
    stop("with `sizes`, `design` must contain integer cluster counts")
  starts <- rep(2:(T + 1), times = design)
  sizes <- as.matrix(sizes)
  if (nrow(sizes) != length(starts) || ncol(sizes) != T)
    stop("`sizes` must be a ", length(starts), " x ", T,
         " matrix (clusters in design order: sequence 2 first, wait-list last)")
  if (any(!is.finite(sizes)) || any(sizes <= 0)) stop("`sizes` must be positive")
  twfe <- estimator == "TWFE"
  P <- diag(T)
  if (twfe) P <- P[, -1, drop = FALSE]
  Om <- vector("list", length(starts))
  for (i in seq_along(starts)) {
    m <- sizes[i, ]
    if (twfe) {
      Om[[i]] <- diag(m) - tcrossprod(m) / sum(m)
    } else {
      if (is.null(icc)) stop("GLS with `sizes` needs `icc` (and `cac`)")
      V <- diag((1 - icc) / m, T) + (1 - cac) * icc * diag(T) +
        cac * icc * matrix(1, T, T)
      Om[[i]] <- solve(V)
    }
  }
  Xs <- lapply(starts, function(g) cbind(P, .sw_X(T, g, model)))
  p <- ncol(P)
  M <- Reduce(`+`, Map(function(X, O) crossprod(X, O %*% X), Xs, Om))
  Mi <- .sw_solve(M, model)
  idx <- (p + 1):ncol(M)
  list(Mi = Mi, Xs = Xs, Om = Om, starts = starts, idx = idx, T = T)
}

#' Effective Hussey--Hughes parameter under nested-exchangeable correlation
#'
#' For cross-sectional sampling with \code{m} individuals in every
#' cluster-period, within-period intracluster correlation \code{icc} and
#' cluster autocorrelation \code{cac}, the covariance of a cluster's period
#' means is proportional to \eqn{I + \phi' J} with
#' \deqn{\phi' = m\,r\,\rho / \{m(1-r)\rho + 1 - \rho\},}
#' so every Hussey--Hughes result holds with \eqn{\phi} replaced by
#' \eqn{\phi'} (Corollary 1 of Halkiewicz and Alexeev). At \code{cac = 1},
#' \eqn{\phi' = m\rho/(1-\rho)}.
#' @param icc within-period intracluster correlation \eqn{\rho \in [0, 1)}
#' @param cac cluster autocorrelation \eqn{r \in [0, 1]}
#' @param m individuals per cluster-period
#' @return the effective parameter \eqn{\phi'} (vectorised)
#' @examples
#' sw_phi_nested(icc = 0.05, cac = 0.5, m = 20)
#' sw_phi_nested(icc = 0.05, cac = 1, m = 20)   # = 20 * 0.05 / 0.95
#' @export
sw_phi_nested <- function(icc, cac, m) {
  if (any(icc < 0 | icc >= 1)) stop("`icc` must lie in [0, 1)")
  if (any(cac < 0 | cac > 1)) stop("`cac` must lie in [0, 1]")
  if (any(m <= 0)) stop("`m` must be positive")
  m * cac * icc / (m * (1 - cac) * icc + 1 - icc)
}

#' Exposure weights of the immediate-treatment estimator in a stepped-wedge design
#'
#' Computes the implicit exposure-time weights \eqn{W_e}, \eqn{e = 1, \ldots, T-1},
#' of the immediate-treatment (IT) estimator, so that
#' \eqn{E(\hat\beta) = \sum_e W_e \theta_e} under exposure-time heterogeneity.
#' The weights sum to one.
#'
#' @param design numeric vector of length \code{T}: cluster counts (or shares)
#'   for start periods \code{2, ..., T}, followed by the number wait-listed.
#'   The first entry must be positive.
#' @param phi Hussey--Hughes parameter \eqn{\phi = m\rho/(1-\rho)} (GLS);
#'   \code{Inf} gives the fixed-effects limit.
#' @param icc,cac,m alternative nested-exchangeable specification: within-period
#'   ICC, cluster autocorrelation and individuals per cluster-period
#'   (see [sw_phi_nested()]). With \code{sizes}, \code{m} is not needed.
#' @param estimator \code{"GLS"} (default) or \code{"TWFE"} (two-way fixed
#'   effects, which does not use the covariance).
#' @param sizes optional \code{G x T} matrix of cluster-period sizes, one row
#'   per cluster in design order (sequence 2 first, wait-list last). Requires
#'   integer counts in \code{design} and, for GLS, \code{icc} and \code{cac}.
#' @return named numeric vector \code{W1, ..., W(T-1)}
#' @examples
#' sw_exposure_weights(c(6, 6, 6, 6, 0), phi = 2)          # complete wedge
#' sw_exposure_weights(c(4, 4, 4, 4, 8), phi = 2)          # 8 wait-listed
#' sw_exposure_weights(c(6, 6, 6, 6, 0), estimator = "TWFE")
#' @export
sw_exposure_weights <- function(design, phi = NULL, icc = NULL, cac = 1, m = NULL,
                                estimator = c("GLS", "TWFE"), sizes = NULL) {
  design <- .sw_check_design(design)
  estimator <- match.arg(estimator)
  T <- length(design); L <- T - 1
  W <- numeric(L)
  if (is.null(sizes)) {
    Om <- .sw_omega(T, phi, icc, cac, m, estimator)
    E <- .sw_equal_info(design, Om, "IT")
    h <- .sw_solve(E$S, "IT model")[1, 1]
    for (i in seq_along(E$gs)) {
      if (design[i] == 0) next
      a <- design[i] * h * as.vector(Om %*% (E$Xs[[i]] - E$Xbar))
      g <- E$gs[i]
      for (t in seq_len(T)) {
        e <- t - g + 1
        if (e >= 1 && e <= L) W[e] <- W[e] + a[t]
      }
    }
  } else {
    F <- .sw_cluster_fit(design, sizes, icc, cac, estimator, "IT")
    h <- F$Mi[F$idx, , drop = FALSE]
    for (i in seq_along(F$starts)) {
      a <- as.vector(h %*% t(F$Xs[[i]]) %*% F$Om[[i]])
      g <- F$starts[i]
      for (t in seq_len(T)) {
        e <- t - g + 1
        if (e >= 1 && e <= L) W[e] <- W[e] + a[t]
      }
    }
  }
  stats::setNames(W, paste0("W", seq_len(L)))
}

.sw_resolve_W <- function(W, design, ...) {
  if (is.null(W) == is.null(design)) stop("supply exactly one of `W` or `design`")
  if (!is.null(design)) W <- sw_exposure_weights(design, ...)
  as.numeric(W)
}

#' Worst-case estimand mismatch of the immediate-treatment estimator
#'
#' \eqn{\Gamma = \sqrt{L}\,\lVert W - v \rVert_2}: the worst-case bias of the IT
#' estimator for the target \eqn{v^\top\theta} per unit RMS exposure-time
#' heterogeneity (Proposition 1 of Halkiewicz and Alexeev).
#' @param W exposure weights (as from [sw_exposure_weights()]); alternatively
#'   supply \code{design} and the covariance arguments in \code{...}
#' @param design optional design, passed with \code{...} to [sw_exposure_weights()]
#' @param target target weights \eqn{v} summing to one; default TATE
#'   (\eqn{v = 1/L})
#' @param ... covariance and estimator arguments for [sw_exposure_weights()]
#' @return \eqn{\Gamma}
#' @examples
#' sw_gamma(design = c(6, 6, 6, 6, 0), phi = 2)
#' @export
sw_gamma <- function(W = NULL, design = NULL, target = NULL, ...) {
  W <- .sw_resolve_W(W, design, ...)
  L <- length(W)
  v <- if (is.null(target)) rep(1 / L, L) else as.numeric(target)
  if (length(v) != L) stop("`target` must have length ", L)
  if (abs(sum(v) - 1) > 1e-10) stop("`target` must sum to one")
  sqrt(L) * sqrt(sum((W - v)^2))
}

#' Adequacy radius (heterogeneity tolerance) of a stepped-wedge design
#'
#' \eqn{\kappa^* = \eta_{\alpha,\epsilon} / \{(z_{1-\alpha/2} + z_{power})\Gamma\}}:
#' the largest RMS exposure-time heterogeneity, relative to the planned effect,
#' at which the worst-case coverage of the nominal \eqn{1-\alpha} IT interval
#' falls by at most \eqn{\epsilon}, for a trial planned with the stated power
#' against a constant effect.
#' @inheritParams sw_gamma
#' @param alpha,epsilon test level and tolerated coverage loss (defaults 0.05)
#' @param power planned power (default 0.80)
#' @return \eqn{\kappa^*} (0.2329 / \eqn{\Gamma} at the defaults)
#' @examples
#' sw_radius(design = c(6, 6, 6, 6, 0), phi = 2)    # 0.133
#' sw_radius(design = c(4, 4, 4, 4, 8), phi = 2)    # 0.214
#' @export
sw_radius <- function(W = NULL, design = NULL, target = NULL, alpha = 0.05,
                      epsilon = 0.05, power = 0.80, ...) {
  G <- sw_gamma(W, design, target, ...)
  .eta_dagger(alpha, epsilon) /
    ((stats::qnorm(1 - alpha / 2) + stats::qnorm(power)) * G)
}

#' Cost of the estimand-aligned (ETI) analysis of a stepped-wedge design
#'
#' Variance of the exposure-time-indicator (ETI) estimator of the target
#' (default TATE), the IT variance of the same design, and both relative to
#' the IT variance of the balanced complete wedge with the same number of
#' clusters, periods and covariance.
#' @inheritParams sw_exposure_weights
#' @param target ETI contrast weights (default TATE, \code{1/L})
#' @return list with \code{var_eti}, \code{var_it}, \code{var_ref} (balanced
#'   complete wedge, IT, equal sizes at the same covariance; with \code{sizes},
#'   the mean size), \code{eti_cost = var_eti / var_ref} and
#'   \code{it_ratio = var_it / var_ref}. Variances are for the stated number of
#'   clusters on the covariance scale of the model (GLS: \eqn{V = I + \phi J};
#'   TWFE: i.i.d. unit errors).
#' @examples
#' sw_eti_cost(c(6, 6, 6, 6, 0), phi = 2)$eti_cost    # 2.27
#' sw_eti_cost(c(6, 6, 1, 0, 11), phi = 2)$eti_cost   # 1.23
#' @export
sw_eti_cost <- function(design, phi = NULL, icc = NULL, cac = 1, m = NULL,
                        estimator = c("GLS", "TWFE"), sizes = NULL, target = NULL) {
  design <- .sw_check_design(design)
  estimator <- match.arg(estimator)
  T <- length(design); L <- T - 1
  cvec <- if (is.null(target)) rep(1 / L, L) else as.numeric(target)
  if (length(cvec) != L) stop("`target` must have length ", L)
  G <- sum(design)
  if (is.null(sizes)) {
    Om <- .sw_omega(T, phi, icc, cac, m, estimator)
    Se <- .sw_equal_info(design, Om, "ETI")$S
    Si <- .sw_equal_info(design, Om, "IT")$S
    var_eti <- drop(crossprod(cvec, .sw_solve(Se, "ETI model") %*% cvec))
    var_it <- 1 / Si[1, 1]
    if (!(var_it > 0) || !is.finite(var_it)) stop("the design does not identify the IT model")
    mref <- m
  } else {
    Fe <- .sw_cluster_fit(design, sizes, icc, cac, estimator, "ETI")
    Fi <- .sw_cluster_fit(design, sizes, icc, cac, estimator, "IT")
    var_eti <- drop(crossprod(cvec, Fe$Mi[Fe$idx, Fe$idx] %*% cvec))
    var_it <- Fi$Mi[Fi$idx, Fi$idx]
    mref <- mean(sizes)
    if (estimator == "GLS") phi <- NULL
  }
  ref <- c(rep(G / (T - 1), T - 1), 0)
  if (estimator == "TWFE" && !is.null(sizes)) {
    Oref <- .sw_omega(T, estimator = "TWFE") * mref
  } else {
    Oref <- .sw_omega(T, phi, icc, cac, mref, estimator)
  }
  var_ref <- 1 / .sw_equal_info(ref, Oref, "IT")$S[1, 1]
  list(var_eti = var_eti, var_it = var_it, var_ref = var_ref,
       eti_cost = var_eti / var_ref, it_ratio = var_it / var_ref)
}

.sw_compositions <- function(G, k) {
  # all non-negative integer vectors of length k summing to G (lexicographic)
  if (k == 1) return(matrix(G, 1, 1))
  out <- NULL
  for (first in 0:G) {
    rest <- .sw_compositions(G - first, k - 1)
    out <- rbind(out, cbind(first, rest, deparse.level = 0))
  }
  out
}

#' Exhaustive search over stepped-wedge and wait-list allocations
#'
#' Enumerates every integer allocation of \code{G} clusters to start periods
#' \code{2, ..., T} and a wait-list (at least one cluster starting in period 2),
#' subject to a rollout cap and a maximum wait-list size, and evaluates the ETI
#' cost, IT variance ratio and adequacy radius of each (equal cluster-period
#' sizes). The optimum is over this finite set only.
#' @param G number of clusters
#' @param T number of periods
#' @param cap maximum number of clusters starting in any one period
#' @param max_wait maximum number of wait-listed clusters
#' @inheritParams sw_exposure_weights
#' @param objective \code{"eti_cost"} (minimise the ETI cost) or
#'   \code{"radius"} (maximise \eqn{\kappa^*} subject to
#'   \code{it_ratio <= max_it_ratio})
#' @param max_it_ratio IT-variance constraint for \code{objective = "radius"}
#' @param top number of best designs returned
#' @param alpha,epsilon,power radius calibration (see [sw_radius()])
#' @return list with \code{best} (data frame of the \code{top} designs),
#'   \code{frontier} (best objective for each wait-list size bound
#'   \code{0, ..., max_wait}) and \code{designs} (all evaluated designs)
#' @examples
#' \donttest{
#' s <- sw_design_search(G = 24, T = 5, cap = 6, max_wait = 12, phi = 2)
#' s$frontier[s$frontier$max_wait %in% c(0, 4, 6, 8, 12), ]
#' }
#' @export
sw_design_search <- function(G, T, cap = Inf, max_wait = G - 1, phi = NULL,
                             icc = NULL, cac = 1, m = NULL,
                             estimator = c("GLS", "TWFE"),
                             objective = c("eti_cost", "radius"),
                             max_it_ratio = Inf, top = 1,
                             alpha = 0.05, epsilon = 0.05, power = 0.80) {
  estimator <- match.arg(estimator); objective <- match.arg(objective)
  if (T < 3 || G < 2) stop("need T >= 3 and G >= 2")
  ncomb <- choose(G + T - 1, T - 1)
  if (ncomb > 5e5) stop("too many allocations (", ncomb, "); reduce G or T")
  C <- .sw_compositions(G, T)
  keep <- C[, 1] >= 1 & apply(C[, 1:(T - 1), drop = FALSE], 1, max) <= cap &
    C[, T] <= max_wait
  C <- C[keep, , drop = FALSE]
  Om <- .sw_omega(T, phi, icc, cac, m, estimator)
  ref <- c(rep(G / (T - 1), T - 1), 0)
  var_ref <- 1 / .sw_equal_info(ref, Om, "IT")$S[1, 1]
  L <- T - 1; u <- rep(1 / L, L)
  kf <- .eta_dagger(alpha, epsilon) / (stats::qnorm(1 - alpha / 2) + stats::qnorm(power))
  # precompute per-sequence blocks
  gs <- 2:(T + 1)
  XE <- lapply(gs, function(g) .sw_X(T, g, "ETI"))
  XI <- lapply(gs, function(g) .sw_X(T, g, "IT"))
  AE <- lapply(XE, function(X) crossprod(X, Om %*% X))
  AI <- vapply(XI, function(X) drop(crossprod(X, Om %*% X)), 0)
  expo <- lapply(gs, function(g) { t <- seq_len(T); e <- t - g + 1; ifelse(e >= 1 & e <= L, e, 0) })
  res <- matrix(NA_real_, nrow(C), 3)
  for (r in seq_len(nrow(C))) {
    n <- C[r, ]
    Xb <- Reduce(`+`, Map(`*`, n, XE)) / G
    Se <- Reduce(`+`, Map(`*`, n, AE)) - G * crossprod(Xb, Om %*% Xb)
    rc <- rcond(Se)
    if (!is.finite(rc) || rc < 1e-12) next
    ve <- drop(crossprod(u, solve(Se, u)))
    db <- Reduce(`+`, Map(`*`, n, XI)) / G
    si <- sum(n * AI) - G * drop(crossprod(db, Om %*% db))
    W <- numeric(L)
    for (i in seq_along(gs)) {
      if (n[i] == 0) next
      a <- n[i] * as.vector(Om %*% (XI[[i]] - db)) / si
      ex <- expo[[i]]
      W[ex[ex > 0]] <- W[ex[ex > 0]] + a[ex > 0]
    }
    res[r, ] <- c(ve / var_ref, (1 / si) / var_ref, kf / (sqrt(L) * sqrt(sum((W - u)^2))))
  }
  ok <- !is.na(res[, 1])
  D <- data.frame(C[ok, , drop = FALSE])
  names(D) <- c(paste0("start", 2:T), "wait")
  D$eti_cost <- res[ok, 1]; D$it_ratio <- res[ok, 2]; D$radius <- res[ok, 3]
  score <- if (objective == "eti_cost") D$eti_cost else
    ifelse(D$it_ratio <= max_it_ratio + 1e-12, -D$radius, Inf)
  ord <- order(score, seq_len(nrow(D)))
  best <- D[ord[seq_len(min(top, sum(is.finite(score))))], , drop = FALSE]
  fr <- lapply(0:max_wait, function(k) {
    idx <- which(D$wait <= k & is.finite(score))
    if (!length(idx)) return(NULL)
    j <- idx[order(score[idx], idx)][1]
    cbind(data.frame(max_wait = k), D[j, , drop = FALSE])
  })
  frontier <- do.call(rbind, fr)
  rownames(best) <- NULL; rownames(frontier) <- NULL; rownames(D) <- NULL
  list(best = best, frontier = frontier, designs = D,
       settings = list(G = G, T = T, cap = cap, max_wait = max_wait,
                       estimator = estimator, objective = objective,
                       var_ref = var_ref))
}

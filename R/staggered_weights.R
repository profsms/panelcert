# Module C - staggered-DiD / TWFE-heterogeneity adequacy (Paper C).
#
# Design statistics (pre-outcome, from the adoption pattern ALONE):
#   Gamma   = sqrt(N1)||w - u|| (Def. def-gamma); Gamma = 0 iff block design.
#   Gamma_S = sqrt(N1)||Pi_S(w - u)|| restricted to a heterogeneity subspace S
#             (Prop. prop-restricted): cohort, event-time, additive combined,
#             and saturated group-time; in balanced panels Gamma_gt = Gamma.
#   worst-case |eta| = (c/sigma) Gamma_S (Cor. cor-gamma); threshold
#             (c/sigma) Gamma_{S,CR} <= eta-dagger (Cor. cor-cv).
# Inference layer:
#   cluster   Gamma_{S,CR} = Gamma_S/sqrt(psi), with direct CR1 as the
#             headline normalization and AR(1)/iid/user-supplied alternatives.
#   point     COVARIANCE-AWARE (eq-pilot): c_S^2 = n_w max{0, (1/N1)||Pi_S d||^2
#             - (1/N1) tr(Pi_S~ A Omega A' Pi_S~)}, Omega = Cov(hat Delta_{g,t})
#             from a fixed-design wild cluster bootstrap. This is descriptive.
#   bounds    the lower test uses an HC2 multiplier radius and reverse triangle
#             inequality.  The upper certificate uses covariance-aware
#             noncentral-chi-square inversion and requires full projected rank.
#             HC3 is a sensitivity.  The full confidence-ball construction is
#             retained as an optional, more general but less powerful fallback.
#   direction eta_dir,S = sqrt(n_w) (w-u)'P_S Delta_hat / q_hat, a signed
#             companion to (not a replacement for) the worst-case envelope.
# Always-treated units (adopted before the sample) are DROPPED (setup g >= 2).

.treatment_indicator <- function(time, first_treat) {
  as.numeric(!is.na(first_treat) & time >= first_treat)
}

.staggered_codes <- function(unit, time, first_treat) {
  n <- length(unit)
  if (length(time) != n || length(first_treat) != n)
    stop("unit, time, first_treat must have equal length")
  uid <- match(unit, unique(unit))
  stimes <- sort(unique(time))
  tid <- match(time, stimes)
  N <- max(uid); T <- length(stimes)
  ft_of <- rep(NA_real_, N); seen <- logical(N)
  for (k in seq_len(n)) {
    f <- if (is.na(first_treat[k])) NA_real_ else as.numeric(first_treat[k])
    if (!seen[uid[k]]) { ft_of[uid[k]] <- f; seen[uid[k]] <- TRUE }
    else if (!identical(is.na(ft_of[uid[k]]), is.na(f)) ||
             (!is.na(f) && ft_of[uid[k]] != f))
      stop("first_treat varies within unit ", unit[k])
  }
  ftc <- numeric(N)
  for (i in seq_len(N)) {
    f <- ft_of[i]
    if (is.na(f) || !is.finite(f) || f > stimes[T]) ftc[i] <- Inf
    else if (f <= stimes[1]) ftc[i] <- 0
    else {
      pos <- which(as.numeric(stimes) == f)
      if (length(pos) == 0) stop("first_treat value ", f, " is not an observed period")
      ftc[i] <- pos
    }
  }
  list(uid = uid, tid = tid, N = N, T = T, ftc = ftc)
}

# Drop always-treated units (cohort code 0); re-index survivors to 1:N'.
.drop_always_treated <- function(uid, tid, ftc) {
  N <- max(uid)
  drop <- which(ftc == 0)
  if (length(drop) == 0) return(list(uid = uid, tid = tid, ftc = ftc, n_drop = 0L))
  keep <- !(uid %in% drop)
  survivors <- setdiff(seq_len(N), drop)
  remap <- stats::setNames(seq_along(survivors), as.character(survivors))
  list(uid = unname(remap[as.character(uid[keep])]), tid = tid[keep],
       ftc = ftc[survivors], keep = keep, n_drop = length(drop))
}

.design_stats <- function(D, Dt) {
  treated <- which(D == 1)
  N1 <- length(treated)
  if (N1 < 2) stop("fewer than 2 treated cells - no staggered design")
  n_w <- sum(Dt^2)
  if (n_w <= 1e-12) stop("treatment has no within variation")
  w <- Dt[treated] / sum(Dt[treated])
  list(treated = treated, N1 = N1, n_w = n_w, w = w,
       Gamma = sqrt(N1) * sqrt(sum((w - 1 / N1)^2)), neg_share = sum(w < 0) / N1)
}

# ---- restricted-profile machinery (Prop. prop-restricted) --------------------
.pinv <- function(M, tol = 1e-10) {
  s <- svd(M)
  d <- ifelse(s$d > tol * max(s$d, 1), 1 / s$d, 0)
  s$v %*% (d * t(s$u))
}
.indicator_basis <- function(labels) {
  labs <- sort(unique(labels))
  B <- matrix(0, length(labels), length(labs))
  B[cbind(seq_along(labels), match(labels, labs))] <- 1
  B
}
.proj <- function(B, x) B %*% (.pinv(crossprod(B)) %*% (t(B) %*% x))
.hat  <- function(B) B %*% .pinv(crossprod(B)) %*% t(B)

.centered_basis <- function(B, tol = 1e-10) {
  Bc <- sweep(B, 2L, colMeans(B), "-")
  s <- svd(Bc, nu = min(dim(Bc)), nv = 0L)
  cutoff <- tol * max(c(s$d, 1))
  keep <- s$d > cutoff
  if (!any(keep)) return(matrix(0, nrow(B), 0L))
  s$u[, keep, drop = FALSE]
}

.noncentrality_upper <- function(wald, rank, error_prob) {
  if (stats::pchisq(wald, df = rank, ncp = 0) <= error_prob) return(0)
  upper <- max(1, wald + rank)
  while (stats::pchisq(wald, df = rank, ncp = upper) > error_prob) {
    upper <- 2 * upper
    if (!is.finite(upper)) stop("failed to bracket noncentrality upper bound")
  }
  stats::uniroot(
    function(value) stats::pchisq(wald, df = rank, ncp = value) - error_prob,
    c(0, upper), tol = 1e-12
  )$root
}

.equivalence_upper <- function(delta_value, Q, A, score, Gamma, se,
                               error_prob) {
  expected_rank <- ncol(Q)
  if (expected_rank == 0L || abs(Gamma) <= 1e-14) {
    return(list(expected_rank = expected_rank, covariance_rank = 0L,
                supported = TRUE, wald = 0, lambda_max = 0,
                ncp_upper = 0, K_upper = 0))
  }
  N1 <- length(delta_value)
  factor <- Gamma / se
  effect <- as.numeric(crossprod(Q, delta_value)) / sqrt(N1)
  score_coordinates <- crossprod(Q, A) %*% score / sqrt(N1)
  covariance <- factor^2 * tcrossprod(score_coordinates)
  covariance <- (covariance + t(covariance)) / 2
  eig <- eigen(covariance, symmetric = TRUE)
  largest <- if (length(eig$values)) max(eig$values) else 0
  tolerance <- max(largest, .Machine$double.xmin) *
    max(dim(covariance)) * .Machine$double.eps * 100
  keep <- eig$values > tolerance
  covariance_rank <- sum(keep)
  answer <- list(expected_rank = expected_rank,
                 covariance_rank = covariance_rank,
                 supported = covariance_rank == expected_rank,
                 wald = NA_real_, lambda_max = NA_real_,
                 ncp_upper = NA_real_, K_upper = NA_real_)
  if (covariance_rank == 0L) return(answer)
  values <- eig$values[keep]
  vectors <- eig$vectors[, keep, drop = FALSE]
  coordinates <- as.numeric(crossprod(vectors, factor * effect))
  answer$wald <- sum(coordinates^2 / values)
  answer$lambda_max <- max(values)
  if (!answer$supported) return(answer)
  answer$ncp_upper <- .noncentrality_upper(answer$wald, covariance_rank,
                                           error_prob)
  answer$K_upper <- sqrt(answer$lambda_max * answer$ncp_upper)
  answer
}

.restricted_gammas <- function(w, g_cell, t_cell, e_cell, N1) {
  d <- w - 1 / N1
  Bcoh <- .indicator_basis(g_cell); Bevt <- .indicator_basis(e_cell)
  Bcmb <- cbind(Bcoh, Bevt)
  Bgt <- .indicator_basis(paste(g_cell, t_cell, sep = "_"))
  gam <- function(B) sqrt(N1) * sqrt(sum(.proj(B, d)^2))
  out <- list(unr = sqrt(N1) * sqrt(sum(d^2)), coh = gam(Bcoh),
              evt = gam(Bevt), cmb = gam(Bcmb), gt = gam(Bgt),
              Bcoh = Bcoh, Bevt = Bevt, Bcmb = Bcmb, Bgt = Bgt)
  if (!isTRUE(all.equal(out$gt, out$unr, tolerance = 1e-9)))
    stop("internal group-time closure check failed")
  out
}

# ---- group-time ATTs (not-yet-treated difference-in-means at (g,t)) ----------
.group_time_atts <- function(uid, tid, ftc, y, N, T,
                             controls = c("not_yet", "never")) {
  controls <- match.arg(controls)
  Ymat <- matrix(NA_real_, N, T)
  Ymat[cbind(uid, tid)] <- y
  cohorts <- sort(unique(ftc[is.finite(ftc) & ftc > 1]))
  keys <- list(); vals <- numeric(0)
  for (g in cohorts) {
    gunits <- which(ftc == g); base <- as.integer(g) - 1L
    for (t in as.integer(g):T) {
      ctrl <- if (controls == "never") which(!is.finite(ftc)) else which(ftc > t)
      gd <- Ymat[gunits, t] - Ymat[gunits, base]
      cd <- Ymat[ctrl, t]  - Ymat[ctrl, base]
      gd <- gd[!is.na(gd)]; cd <- cd[!is.na(cd)]
      if (length(gd) == 0 || length(cd) == 0) next
      keys[[length(keys) + 1L]] <- c(as.integer(g), t)
      vals <- c(vals, mean(gd) - mean(cd))
    }
  }
  list(keys = keys, vals = vals, cohorts = cohorts)
}
.gt_lookup <- function(gt) {
  if (length(gt$keys) == 0) return(function(g, t) NA_real_)
  nm <- vapply(gt$keys, function(k) paste(k, collapse = "_"), "")
  v <- stats::setNames(gt$vals, nm)
  function(g, t) { z <- v[paste(as.integer(g), t, sep = "_")]; if (is.na(z)) NA_real_ else unname(z) }
}
.cohort_means <- function(gt) {
  if (length(gt$keys) == 0) return(numeric(0))
  gvec <- vapply(gt$keys, `[`, integer(1), 1L)
  tapply(gt$vals, gvec, mean)
}

# ---- covariance-aware pilot (eq-pilot) ---------------------------------------
.delta_cell <- function(look, g_cell, t_cell, N1)
  vapply(seq_len(N1), function(k) look(g_cell[k], t_cell[k]), numeric(1))

.cov_pilots <- function(look, g_cell, t_cell, N1, bases, traces, n_w, sigma) {
  delta <- .delta_cell(look, g_cell, t_cell, N1)
  cov <- !is.na(delta); dc <- delta - mean(delta[cov])
  pil <- function(B, tr) {
    comp <- .proj(B[cov, , drop = FALSE], dc[cov])
    raw <- sum(comp^2) / sum(cov)
    sqrt(n_w * max(0, raw - tr) / sigma^2)
  }
  list(coh = pil(bases$Bcoh, traces$coh), evt = pil(bases$Bevt, traces$evt),
       cmb = pil(bases$Bcmb, traces$cmb), gt = pil(bases$Bgt, traces$gt))
}

.directional_profile <- function(look, g_cell, t_cell, N1, B, w,
                                 n_w, sigma, psi) {
  delta <- .delta_cell(look, g_cell, t_cell, N1)
  covered <- is.finite(delta)
  if (!all(covered) || !is.finite(sigma) || sigma <= 0 ||
      !is.finite(psi) || psi <= 0)
    return(list(eta = NA_real_, alignment = NA_real_))
  profile <- .proj(B, delta - mean(delta))
  exposure <- .proj(B, w - 1 / N1)
  bias <- sum(exposure * profile)
  den <- sqrt(sum(exposure^2) * sum(profile^2))
  list(eta = sqrt(n_w) * bias / (sigma * sqrt(psi)),
       alignment = if (den > 0) bias / den else 0)
}

.group_time_map <- function(uid, tid, ftc, N, T,
                            controls = c("not_yet", "never")) {
  controls <- match.arg(controls)
  n <- length(uid)
  obs <- matrix(NA_integer_, N, T)
  obs[cbind(uid, tid)] <- seq_len(n)
  cohorts <- sort(unique(ftc[is.finite(ftc) & ftc > 1]))
  keys <- list()
  for (g in cohorts) for (tt in as.integer(g):T) {
    ctrl <- if (controls == "never") which(!is.finite(ftc)) else which(ftc > tt)
    gunits <- which(ftc == g)
    gunits <- gunits[!is.na(obs[gunits, tt]) & !is.na(obs[gunits, g - 1L])]
    ctrl <- ctrl[!is.na(obs[ctrl, tt]) & !is.na(obs[ctrl, g - 1L])]
    if (length(gunits) && length(ctrl))
      keys[[length(keys) + 1L]] <- c(as.integer(g), tt)
  }
  key_names <- vapply(keys, function(k) paste(k, collapse = "_"), "")
  L <- matrix(0, n, length(keys))
  for (j in seq_along(keys)) {
    g <- keys[[j]][1]; tt <- keys[[j]][2]; base <- g - 1L
    gunits <- which(ftc == g)
    ctrl <- if (controls == "never") which(!is.finite(ftc)) else which(ftc > tt)
    gunits <- gunits[!is.na(obs[gunits, tt]) & !is.na(obs[gunits, base])]
    ctrl <- ctrl[!is.na(obs[ctrl, tt]) & !is.na(obs[ctrl, base])]
    L[obs[gunits, tt], j] <- 1 / length(gunits)
    L[obs[gunits, base], j] <- -1 / length(gunits)
    L[obs[ctrl, tt], j] <- -1 / length(ctrl)
    L[obs[ctrl, base], j] <- 1 / length(ctrl)
  }
  list(keys = keys, names = key_names, L = L)
}

# ---- fixed-design wild cluster bootstrap -------------------------------------
.wild_bootstrap <- function(uid, tid, ftc, D, Dt, y, N, T, n_w, treated, N1,
                            g_cell, t_cell, d_K, gt0, B, seed, controls) {
  n <- length(uid)
  # Absorb unit and time effects before fitting the modest set of treated
  # group-time indicators. This is algebraically the same saturated surface as
  # [unit FE, time FE, treated (g,t)] without an n x N dummy matrix.
  gtmap <- .group_time_map(uid, tid, ftc, N, T, controls)
  keys_nm <- gtmap$names
  m <- length(keys_nm)
  Z <- matrix(0, n, m)
  treated_obs <- which(D == 1)
  zcol <- match(paste(ftc[uid[treated_obs]], tid[treated_obs], sep = "_"), keys_nm)
  keep_z <- !is.na(zcol)
  Z[cbind(treated_obs[keep_z], zcol[keep_z])] <- 1
  Zt <- vapply(seq_len(m), function(j)
    .twoway_demean_codes(Z[, j], uid, tid, N, T), numeric(n))
  yt <- .twoway_demean_codes(y, uid, tid, N, T)
  fit <- as.numeric(Zt %*% (.pinv(crossprod(Zt)) %*% (t(Zt) %*% yt)))
  e <- yt - fit
  yhat <- y - e
  att_center <- as.numeric(crossprod(gtmap$L, yhat))
  score <- matrix(0, m, N)
  # HC2/HC3 adjustments apply to the low-dimensional treated group-time
  # surface after unit and time effects have been removed. They quantify pilot
  # uncertainty; they are not alternative SEs for the TWFE test under audit.
  Hcoef <- Zt %*% .pinv(crossprod(Zt))
  hgt <- rowSums(Hcoef * Zt)
  one_minus_h <- pmax(1 - hgt, 1e-10)
  e_hc2 <- e / sqrt(one_minus_h)
  e_hc3 <- e / one_minus_h
  score_hc2 <- matrix(0, m, N)
  score_hc3 <- matrix(0, m, N)
  for (j in seq_len(m)) {
    score[j, ] <- rowsum(gtmap$L[, j] * e, uid, reorder = FALSE)[, 1L]
    score_hc2[j, ] <- rowsum(gtmap$L[, j] * e_hc2, uid,
                             reorder = FALSE)[, 1L]
    score_hc3[j, ] <- rowsum(gtmap$L[, j] * e_hc3, uid,
                             reorder = FALSE)[, 1L]
  }
  dof <- n - d_K - 1
  set.seed(seed)
  sig <- numeric(0); psi_ar1 <- numeric(0); psi_direct <- numeric(0); gtm <- list()
  error_hc2 <- matrix(NA_real_, m, B)
  error_hc3 <- matrix(NA_real_, m, B)
  for (b in seq_len(B)) {
    v <- sample(c(-1, 1), N, replace = TRUE)
    error_hc2[, b] <- as.numeric(score_hc2 %*% v)
    error_hc3[, b] <- as.numeric(score_hc3 %*% v)
    ystar <- yhat + v[uid] * e
    yt <- .twoway_demean_codes(ystar, uid, tid, N, T)
    beta <- sum(Dt * yt) / n_w
    resid <- yt - beta * Dt
    s <- sqrt(sum(resid^2) / dof)
    rho <- .rho_ar1(resid, uid, tid)
    pa <- .psi_parametric(Dt, uid, tid, rho, kind = "ar1")
    pd <- .psi_direct(Dt, resid, uid, n_w, s^2, N)
    if (!is.finite(pa) || pa <= 0 || !is.finite(pd) || pd <= 0) next
    vec <- stats::setNames(as.numeric(att_center + score %*% v), keys_nm)
    sig <- c(sig, s); psi_ar1 <- c(psi_ar1, pa)
    psi_direct <- c(psi_direct, pd); gtm[[length(gtm) + 1L]] <- vec
  }
  # Exact covariance under the fitted Rademacher wild-bootstrap distribution.
  Omega <- score %*% t(score)
  list(Omega = Omega, sigma = sig, psi_ar1 = psi_ar1,
       psi_direct = psi_direct, gtm = gtm, keys = keys_nm,
       score_hc2 = score_hc2, score_hc3 = score_hc3,
       error_hc2 = error_hc2, error_hc3 = error_hc3)
}

#' Pre-outcome TWFE design-statistic ladder
#'
#' The design statistic \code{Gamma = sqrt(N1)||w - u||} and its restricted
#' variants \code{Gamma_coh}, \code{Gamma_evt}, \code{Gamma_c+e}, and
#' \code{Gamma_gt} (Prop. prop-restricted), the negative-weight share, and the
#' breakdown ratio from the adoption pattern alone. In a balanced panel,
#' \code{Gamma_gt = Gamma}. Always-treated units are dropped.
#' @param unit,time raw identifiers (time numeric)
#' @param first_treat adoption time per observation (constant within unit; NA = never)
#' @param alpha,delta level and size tolerance
#' @return an \code{AdequacyReport}
#' @export
twfe_design <- function(unit, time, first_treat, alpha = 0.05, delta = 0.05) {
  sc <- .staggered_codes(unit, time, first_treat)
  dr <- .drop_always_treated(sc$uid, sc$tid, sc$ftc)
  uid <- dr$uid; tid <- dr$tid; ftc <- dr$ftc; N <- max(uid); T <- sc$T
  n <- length(uid)
  D <- as.numeric(is.finite(ftc[uid]) & tid >= ftc[uid])
  Dt <- .twoway_demean_codes(D, uid, tid, N, T)
  ds <- .design_stats(D, Dt)
  g_cell <- ftc[uid[ds$treated]]; t_cell <- tid[ds$treated]
  e_cell <- t_cell - g_cell
  G <- .restricted_gammas(ds$w, g_cell, t_cell, e_cell, ds$N1)
  fd <- fe_dimension(uid, tid, N, T)
  design <- .design_summary_codes(uid, tid, N, T, xt = Dt)
  eta_dag <- .eta_dagger(alpha, delta)
  notes <- character(0)
  if (dr$n_drop > 0)
    notes <- c(notes, sprintf("%d always-treated unit(s) dropped (setup g >= 2)", dr$n_drop))
  statistic <- list(Gamma = ds$Gamma, Gamma_coh = G$coh, Gamma_evt = G$evt,
                    Gamma_cmb = G$cmb, Gamma_gt = G$gt,
                    neg_share = ds$neg_share, N1 = ds$N1, n_w = ds$n_w)
  if (ds$Gamma <= 1e-8) {
    notes <- c(notes, "block design: within-transformed treatment is uniform (Prop. prop-gamma0)")
    return(.new_AdequacyReport("twfe_heterogeneity", design, statistic, NULL,
                               eta_dag, Inf, NULL, "CERTIFIED", alpha, delta, notes))
  }
  breakdown <- eta_dag / G$gt
  notes <- c(notes, sprintf("pre-outcome: the saturated group-time class requires c/sigma <= %.3f (= eta-dagger/Gamma_gt); supply the outcome (twfe_adequacy) to calibrate c/sigma and obtain a boundary-robust certificate", breakdown))
  .new_AdequacyReport("twfe_heterogeneity", design, statistic, NULL, eta_dag,
                      breakdown, NULL, "INCONCLUSIVE", alpha, delta, notes)
}

#' TWFE design-statistic ladder
#'
#' Convenience accessor returning the unrestricted, cohort, event-time,
#' combined additive, and saturated group-time design statistics from the
#' adoption pattern.
#'
#' @inheritParams twfe_design
#' @return A list with `unr`, `coh`, `evt`, `cmb`, `gt`, `neg_share`, `N1`,
#'   and `n_w`.
#' @export
twfe_gammas <- function(unit, time, first_treat) {
  s <- twfe_design(unit, time, first_treat)$statistic
  list(unr = s$Gamma, coh = s$Gamma_coh, evt = s$Gamma_evt,
       cmb = s$Gamma_cmb, gt = s$Gamma_gt,
       neg_share = s$neg_share, N1 = s$N1,
       n_w = s$n_w)
}

#' TWFE-heterogeneity adequacy (covariance-aware, wild-bootstrap)
#'
#' Restricted design-statistic ladder, direct CR1 rescaling
#' \code{Gamma_{S,CR} = Gamma_S/sqrt(psi)}, covariance-aware pilots
#' \code{c_S/sigma}, saturated group-time and restricted point envelopes, and
#' the signed directional plug-in. A fixed-design cluster-multiplier radius
#' gives the regular one-sided lower bound. The upper certificate inverts the
#' noncentral chi-square law of the projected Wald statistic and is issued only
#' when the score covariance has full rank in the prespecified class.
#' Point-envelope bootstrap percentiles remain descriptive. Always-treated
#' units are dropped.
#' @param object outcome vector
#' @param unit,time,first_treat as in [twfe_design()]
#' @param alpha,delta level and size tolerance
#' @param cluster normalization: \code{"direct"} (default, the reported CR1
#'   score scale), \code{"ar1"}, or \code{"iid"}
#' @param psi optional positive user-supplied variance-inflation factor; when
#'   supplied it overrides \code{cluster}
#' @param controls comparison group for group-time effects: not-yet-treated
#'   (including never-treated) or never-treated only
#' @param heterogeneity_class prespecified class used for the report verdict:
#'   saturated group-time, additive cohort-plus-event-time, cohort, or event-time
#' @param bootstrap number of wild-cluster draws (>0 enables the covariance
#'   correction, descriptive point summaries, and one-sided procedures)
#' @param gamma error probability for each reported one-sided procedure
#' @param q_band optional relative half-width for the denominator band in the
#'   more general full confidence-ball construction. It does not affect the
#'   main decision-specific verdict.
#' @param seed RNG seed for the wild bootstrap
#' @param ... passed between methods
#' @return an \code{AdequacyReport}
#' @export
twfe_adequacy <- function(object, ...) UseMethod("twfe_adequacy")

#' @rdname twfe_adequacy
#' @export
twfe_adequacy.default <- function(object, unit, time, first_treat,
                                  alpha = 0.05, delta = 0.05,
                                   cluster = c("direct", "ar1", "iid"), psi = NULL,
                                   controls = c("not_yet", "never"),
                                   heterogeneity_class = c("group_time", "additive",
                                                           "cohort", "event"),
                                   bootstrap = 999L, seed = 20260715L,
                                   gamma = 0.05, q_band = NULL, ...) {
  y0 <- object; cluster <- match.arg(cluster); controls <- match.arg(controls)
  heterogeneity_class <- match.arg(heterogeneity_class)
  if (!is.null(psi) && (!is.numeric(psi) || length(psi) != 1L ||
                        !is.finite(psi) || psi <= 0))
    stop("psi must be a single positive finite number")
  if (!is.numeric(gamma) || length(gamma) != 1L || !is.finite(gamma) ||
      gamma <= 0 || gamma >= 0.5)
    stop("gamma must be a single number in (0, 0.5)")
  if (!is.null(q_band) && (!is.numeric(q_band) || length(q_band) != 1L ||
                           !is.finite(q_band) || q_band < 0 || q_band >= 1))
    stop("q_band must be NULL or a single number in [0, 1)")
  sc <- .staggered_codes(unit, time, first_treat)
  if (length(y0) != length(sc$uid)) stop("y must have length n = ", length(sc$uid))
  dr <- .drop_always_treated(sc$uid, sc$tid, sc$ftc)
  uid <- dr$uid; tid <- dr$tid; ftc <- dr$ftc; N <- max(uid); T <- sc$T
  y <- if (dr$n_drop > 0) y0[dr$keep] else y0
  n <- length(uid)
  D <- as.numeric(is.finite(ftc[uid]) & tid >= ftc[uid])
  Dt <- .twoway_demean_codes(D, uid, tid, N, T)
  ds <- .design_stats(D, Dt)
  g_cell <- ftc[uid[ds$treated]]; t_cell <- tid[ds$treated]
  e_cell <- t_cell - g_cell
  G <- .restricted_gammas(ds$w, g_cell, t_cell, e_cell, ds$N1)
  fd <- fe_dimension(uid, tid, N, T); dof <- n - fd$d_K - 1
  if (dof <= 0) stop("no residual degrees of freedom")
  yt <- .twoway_demean_codes(y, uid, tid, N, T)
  beta <- sum(Dt * yt) / ds$n_w
  resid <- yt - beta * Dt
  sigma <- sqrt(sum(resid^2) / dof)

  notes <- character(0)
  if (dr$n_drop > 0)
    notes <- c(notes, sprintf("%d always-treated unit(s) dropped (setup g >= 2)", dr$n_drop))

  rho_ar1 <- .rho_ar1(resid, uid, tid)
  psi_ar1 <- .psi_parametric(Dt, uid, tid, rho_ar1, "ar1")
  psi_direct <- .psi_direct(Dt, resid, uid, ds$n_w, sigma^2, N)
  normalization <- if (!is.null(psi)) "user-supplied" else cluster
  psi_hat <- if (!is.null(psi)) as.numeric(psi) else switch(
    cluster, direct = psi_direct, ar1 = psi_ar1, iid = 1)
  CR <- list(unr = G$unr / sqrt(psi_hat), coh = G$coh / sqrt(psi_hat),
             evt = G$evt / sqrt(psi_hat), cmb = G$cmb / sqrt(psi_hat),
             gt = G$gt / sqrt(psi_hat))

  gt0 <- .group_time_atts(uid, tid, ftc, y, N, T, controls)
  look0 <- .gt_lookup(gt0)
  bases <- list(Bcoh = G$Bcoh, Bevt = G$Bevt, Bcmb = G$Bcmb,
                Bgt = G$Bgt)
  directional <- .directional_profile(
    look0, g_cell, t_cell, ds$N1, G$Bcmb, ds$w,
    ds$n_w, sigma, psi_hat
  )
  boot <- NULL
  formal <- NULL
  if (bootstrap > 0 && length(gt0$keys) > 0) {
    wb <- .wild_bootstrap(uid, tid, ftc, D, Dt, y, N, T, ds$n_w, ds$treated,
                          ds$N1, g_cell, t_cell, fd$d_K, gt0, as.integer(bootstrap),
                          as.integer(seed), controls)
    keys_nm <- wb$keys
    A <- matrix(0, ds$N1, length(keys_nm))
    cellkey <- paste(as.integer(g_cell), t_cell, sep = "_")
    for (k in seq_len(ds$N1)) { j <- match(cellkey[k], keys_nm); if (!is.na(j)) A[k, j] <- 1 }
    bnames <- c("coh", "evt", "cmb", "gt")
    Blist <- list(coh = G$Bcoh, evt = G$Bevt, cmb = G$Bcmb, gt = G$Bgt)
    Qlist <- lapply(Blist, .centered_basis)
    # tr(Q' A Omega A' Q) is evaluated in group-time coordinates. This avoids
    # materializing an N1-by-N1 treated-cell covariance matrix.
    traces <- lapply(Qlist, function(Q) {
      coordinates <- crossprod(Q, A)
      sum((coordinates %*% wb$Omega) * coordinates) / ds$N1
    })
    pilots <- .cov_pilots(look0, g_cell, t_cell, ds$N1, bases, traces, ds$n_w, sigma)
    size_pt <- list(coh = .noncentral_size(pilots$coh * CR$coh, alpha),
                    evt = .noncentral_size(pilots$evt * CR$evt, alpha),
                    cmb = .noncentral_size(pilots$cmb * CR$cmb, alpha),
                    gt = .noncentral_size(pilots$gt * CR$gt, alpha))
    draw_envelope <- numeric(0); draw_directional <- numeric(0)
    draw_directional_eta <- numeric(0); draw_alignment <- numeric(0)
    draw_psi <- if (!is.null(psi)) rep(as.numeric(psi), length(wb$sigma)) else
      switch(cluster, direct = wb$psi_direct, ar1 = wb$psi_ar1,
             iid = rep(1, length(wb$sigma)))
    for (i in seq_along(wb$sigma)) {
      vec <- wb$gtm[[i]]
      look <- function(g, t) { z <- vec[paste(as.integer(g), t, sep = "_")]; if (is.na(z)) NA_real_ else unname(z) }
      p <- .cov_pilots(look, g_cell, t_cell, ds$N1, bases, traces, ds$n_w, wb$sigma[i])
      draw_envelope <- c(draw_envelope,
        .noncentral_size(p$cmb * G$cmb / sqrt(draw_psi[i]), alpha))
      dp <- .directional_profile(look, g_cell, t_cell, ds$N1, G$Bcmb,
                                 ds$w, ds$n_w, wb$sigma[i], draw_psi[i])
      draw_directional_eta <- c(draw_directional_eta, dp$eta)
      draw_directional <- c(draw_directional,
                            .noncentral_size(dp$eta, alpha))
      draw_alignment <- c(draw_alignment, dp$alignment)
    }
    q <- function(v, p) stats::quantile(v[is.finite(v)], p, names = FALSE)

    # The trace-debiased quadratic above remains a descriptive point
    # calibration.  The lower test applies reverse triangle inequality to the
    # multiplier norm radius.  The upper certificate uses covariance-aware
    # noncentral-chi-square inversion and therefore checks projected rank.
    delta0 <- .delta_cell(look0, g_cell, t_cell, ds$N1)
    if (all(is.finite(delta0)) && all(rowSums(A) == 1)) {
      projected_scale <- function(Q, delta_value) {
        dc <- delta_value - mean(delta_value)
        sqrt(ds$n_w / ds$N1) * sqrt(sum(crossprod(Q, dc)^2))
      }
      error_norms <- function(Q, E) {
        projected_coordinates <- crossprod(Q, A) %*% E
        sqrt(ds$n_w / ds$N1) * sqrt(colSums(projected_coordinates^2))
      }
      scales <- stats::setNames(vapply(Qlist, projected_scale, numeric(1),
                                      delta_value = delta0), bnames)
      radii <- stats::setNames(vapply(Qlist, function(Q)
        q(error_norms(Q, wb$error_hc2), 1 - gamma), numeric(1)), bnames)
      radii_hc3 <- stats::setNames(vapply(Qlist, function(Q)
        q(error_norms(Q, wb$error_hc3), 1 - gamma), numeric(1)), bnames)
      gamma_values <- c(coh = G$coh, evt = G$evt, cmb = G$cmb, gt = G$gt)
      q_scale <- sigma * sqrt(psi_hat)
      se_scale <- q_scale / sqrt(ds$n_w)
      lower <- gamma_values * pmax(0, scales - radii) / q_scale
      lower_hc3 <- gamma_values * pmax(0, scales - radii_hc3) / q_scale
      equivalence <- lapply(bnames, function(name)
        .equivalence_upper(delta0, Qlist[[name]], A, wb$score_hc2,
                           gamma_values[[name]], se_scale, gamma))
      names(equivalence) <- bnames
      equivalence_hc3 <- lapply(bnames, function(name)
        .equivalence_upper(delta0, Qlist[[name]], A, wb$score_hc3,
                           gamma_values[[name]], se_scale, gamma))
      names(equivalence_hc3) <- bnames
      upper <- stats::setNames(vapply(equivalence, `[[`, numeric(1), "K_upper"),
                               bnames)
      upper_hc3 <- stats::setNames(vapply(equivalence_hc3, `[[`, numeric(1),
                                           "K_upper"), bnames)
      ball_upper_raw <- gamma_values * (scales + radii) / q_scale
      ball_upper_raw_hc3 <- gamma_values * (scales + radii_hc3) / q_scale
      ball_lower <- if (is.null(q_band)) rep(NA_real_, length(lower)) else
        (1 - q_band) * lower
      ball_upper <- if (is.null(q_band)) rep(NA_real_, length(lower)) else
        (1 + q_band) * ball_upper_raw
      names(ball_lower) <- names(ball_upper) <- bnames
      formal <- list(confidence = 1 - gamma, scale = scales,
                     radius = radii, radius_hc3 = radii_hc3,
                     lower = lower, upper = upper,
                     lower_hc3 = lower_hc3, upper_hc3 = upper_hc3,
                     equivalence = equivalence,
                     equivalence_hc3 = equivalence_hc3,
                     ball_upper_raw = ball_upper_raw,
                     ball_upper_raw_hc3 = ball_upper_raw_hc3,
                     ball_lower = ball_lower, ball_upper = ball_upper,
                     q_band = q_band)
    } else {
      notes <- c(notes, paste0(
        "one-sided procedures unavailable: the group-time estimator ",
        "does not cover every treated cell"))
    }

    boot <- list(n = length(wb$sigma),
                 envelope_med = q(draw_envelope, 0.5),
                 envelope_lo = q(draw_envelope, 0.025),
                 envelope_hi = q(draw_envelope, 0.975),
                 envelope_p95 = q(draw_envelope, 0.95),
                 directional_med = q(draw_directional, 0.5),
                 directional_lo = q(draw_directional, 0.025),
                 directional_hi = q(draw_directional, 0.975),
                 directional_eta_med = q(draw_directional_eta, 0.5),
                 directional_eta_lo = q(draw_directional_eta, 0.025),
                 directional_eta_hi = q(draw_directional_eta, 0.975),
                 alignment_med = q(draw_alignment, 0.5),
                 psi_lo = q(draw_psi, 0.025), psi_hi = q(draw_psi, 0.975),
                 psi_direct_lo = q(wb$psi_direct, 0.025),
                 psi_direct_hi = q(wb$psi_direct, 0.975),
                 psi_ar1_lo = q(wb$psi_ar1, 0.025),
                 psi_ar1_hi = q(wb$psi_ar1, 0.975),
                 Omega_trace_cmb = traces$cmb,
                 Omega_trace_gt = traces$gt,
                 norm = formal)
    # v0.6 compatibility aliases; these are envelope, never realized-size, fields.
    boot$cmb_med <- boot$envelope_med
    boot$cmb_lo <- boot$envelope_lo
    boot$cmb_hi <- boot$envelope_hi
    notes <- c(notes, sprintf(paste0(
      "point pilot: covariance trace computed from the fitted Rademacher ",
      "score covariance; %d multiplier draws calibrate the %.1f%% HC2 ",
      "lower norm radius; the upper procedure uses projected-Wald ",
      "noncentrality inversion"), boot$n, 100 * (1 - gamma)))
  } else {
    rawpil <- function(B) {
      delta <- .delta_cell(look0, g_cell, t_cell, ds$N1)
      cov <- !is.na(delta); dc <- delta - mean(delta[cov])
      comp <- .proj(B[cov, , drop = FALSE], dc[cov])
      sqrt(ds$n_w * sum(comp^2) / sum(cov) / sigma^2)
    }
    pilots <- list(coh = rawpil(G$Bcoh), evt = rawpil(G$Bevt),
                   cmb = rawpil(G$Bcmb), gt = rawpil(G$Bgt))
    size_pt <- list(coh = .noncentral_size(pilots$coh * CR$coh, alpha),
                    evt = .noncentral_size(pilots$evt * CR$evt, alpha),
                    cmb = .noncentral_size(pilots$cmb * CR$cmb, alpha),
                    gt = .noncentral_size(pilots$gt * CR$gt, alpha))
    notes <- c(notes, paste0(
      "bootstrap disabled: pilots are raw projected dispersion and no ",
      "boundary-robust certificate is available; set bootstrap > 0"))
  }

  notes <- c(notes, sprintf(
    "normalization = %s: psi_hat = %.3f; direct CR1 psi = %.3f; AR(1) psi = %.3f (rho = %.3f)",
    normalization, psi_hat, psi_direct, psi_ar1, rho_ar1))
  if (!is.finite(directional$eta))
    notes <- c(notes, "directional plug-in unavailable because the group-time profile does not cover every treated cell")
  if (N < 40)
    notes <- c(notes, sprintf("few clusters (G = %d): wild bootstrap refinements advisable (Rem. sec-clusters)", N))
  if (T / N > 0.25)
    notes <- c(notes, sprintf("fixed-T, many-cluster approximation is strained (G = %d, T = %d)", N, T))

  eta_dag <- .eta_dagger(alpha, delta)
  class_keys <- c(group_time = "gt", additive = "cmb", cohort = "coh",
                  event = "evt")
  selected_key <- unname(class_keys[[heterogeneity_class]])
  class_names <- c(coh = "cohort", evt = "event-time", cmb = "additive",
                   gt = "saturated group-time")
  decide <- function(key) {
    if (is.null(formal) || !is.finite(formal$lower[[key]]))
      return("INCONCLUSIVE")
    if (formal$lower[[key]] > eta_dag) return("FLAGGED")
    if (isTRUE(formal$equivalence[[key]]$supported) &&
        is.finite(formal$upper[[key]]) && formal$upper[[key]] <= eta_dag)
      return("CERTIFIED")
    "INCONCLUSIVE"
  }
  class_verdicts <- stats::setNames(vapply(c("coh", "evt", "cmb", "gt"),
                                            decide, character(1)),
                                      c("coh", "evt", "cmb", "gt"))
  verdict <- unname(class_verdicts[[selected_key]])
  eta_selected <- pilots[[selected_key]] * CR[[selected_key]]
  if (!is.null(formal)) {
    eq_selected <- formal$equivalence[[selected_key]]
    if (verdict == "CERTIFIED") {
      notes <- c(notes, sprintf(
        "%s upper bound %.3f is below eta-dagger %.3f",
        class_names[[selected_key]], formal$upper[[selected_key]], eta_dag))
    } else if (verdict == "FLAGGED") {
      notes <- c(notes, sprintf(paste0(
        "%s lower bound %.3f exceeds eta-dagger %.3f; uniform ",
        "certification is withheld, not the published inference invalidated"),
        class_names[[selected_key]], formal$lower[[selected_key]], eta_dag))
    } else if (!isTRUE(eq_selected$supported)) {
      notes <- c(notes, sprintf(paste0(
        "%s upper certificate unavailable: projected score covariance rank ",
        "%d/%d; result is inconclusive for this structural reason"),
        class_names[[selected_key]], eq_selected$covariance_rank,
        eq_selected$expected_rank))
    } else {
      notes <- c(notes, sprintf(paste0(
        "%s one-sided lower %.3f and upper %.3f do not determine a verdict ",
        "at eta-dagger %.3f"), class_names[[selected_key]],
        formal$lower[[selected_key]], formal$upper[[selected_key]], eta_dag))
    }
    if (!is.null(q_band))
      notes <- c(notes, sprintf(paste0(
        "general confidence-ball fallback reported with denominator band ",
        "q_band = %.4f; it does not determine the main verdict"), q_band))
  }
  design <- .design_summary_codes(uid, tid, N, T, xt = Dt)
  statistic <- list(Gamma = ds$Gamma, Gamma_coh = G$coh, Gamma_evt = G$evt,
                    Gamma_cmb = G$cmb, Gamma_gt = G$gt,
                    neg_share = ds$neg_share, psi_hat = psi_hat,
                    normalization = normalization, psi_direct = psi_direct,
                    psi_ar1 = psi_ar1, psi_driven = psi_direct,
                    rho_ar1 = rho_ar1, Gamma_CR = CR$unr,
                    Gamma_coh_CR = CR$coh, Gamma_evt_CR = CR$evt,
                    Gamma_cmb_CR = CR$cmb, Gamma_gt_CR = CR$gt,
                    beta = beta, sigma = sigma,
                    se_cr1 = sigma * sqrt(psi_direct / ds$n_w),
                    q_hat = sigma * sqrt(psi_direct),
                    sign_reversal_rms = if (ds$Gamma > 0) abs(beta) / ds$Gamma else Inf,
                    N1 = ds$N1,
                    n_w = ds$n_w, n_cohorts = length(gt0$cohorts),
                    att_target = {
                      d0 <- .delta_cell(look0, g_cell, t_cell, ds$N1)
                      if (all(is.finite(d0))) mean(d0) else NA_real_
                    },
                    pilot_coh = pilots$coh, pilot_evt = pilots$evt,
                    pilot_cmb = pilots$cmb, pilot_gt = pilots$gt,
                     size_coh = size_pt$coh, size_evt = size_pt$evt,
                     size_cmb = size_pt$cmb, size_gt = size_pt$gt,
                     selected_class = heterogeneity_class,
                     verdict_coh = unname(class_verdicts["coh"]),
                     verdict_evt = unname(class_verdicts["evt"]),
                     verdict_cmb = unname(class_verdicts["cmb"]),
                     verdict_gt = unname(class_verdicts["gt"]),
                     K_lower_coh = if (is.null(formal)) NA_real_ else unname(formal$lower["coh"]),
                     K_upper_coh = if (is.null(formal)) NA_real_ else unname(formal$upper["coh"]),
                     K_lower_evt = if (is.null(formal)) NA_real_ else unname(formal$lower["evt"]),
                     K_upper_evt = if (is.null(formal)) NA_real_ else unname(formal$upper["evt"]),
                     K_lower_cmb = if (is.null(formal)) NA_real_ else unname(formal$lower["cmb"]),
                     K_upper_cmb = if (is.null(formal)) NA_real_ else unname(formal$upper["cmb"]),
                     K_lower_gt = if (is.null(formal)) NA_real_ else unname(formal$lower["gt"]),
                     K_upper_gt = if (is.null(formal)) NA_real_ else unname(formal$upper["gt"]),
                    K_lower_gt_hc3 = if (is.null(formal)) NA_real_ else unname(formal$lower_hc3["gt"]),
                     K_upper_gt_hc3 = if (is.null(formal)) NA_real_ else unname(formal$upper_hc3["gt"]),
                     expected_rank_coh = if (is.null(formal)) NA_integer_ else formal$equivalence$coh$expected_rank,
                     covariance_rank_coh = if (is.null(formal)) NA_integer_ else formal$equivalence$coh$covariance_rank,
                     expected_rank_evt = if (is.null(formal)) NA_integer_ else formal$equivalence$evt$expected_rank,
                     covariance_rank_evt = if (is.null(formal)) NA_integer_ else formal$equivalence$evt$covariance_rank,
                     expected_rank_cmb = if (is.null(formal)) NA_integer_ else formal$equivalence$cmb$expected_rank,
                     covariance_rank_cmb = if (is.null(formal)) NA_integer_ else formal$equivalence$cmb$covariance_rank,
                     expected_rank_gt = if (is.null(formal)) NA_integer_ else formal$equivalence$gt$expected_rank,
                     covariance_rank_gt = if (is.null(formal)) NA_integer_ else formal$equivalence$gt$covariance_rank,
                     q_band = if (is.null(q_band)) NA_real_ else q_band,
                     K_ball_lower_gt = if (is.null(formal)) NA_real_ else unname(formal$ball_lower["gt"]),
                     K_ball_upper_gt = if (is.null(formal)) NA_real_ else unname(formal$ball_upper["gt"]),
                     gamma = gamma,
                    eta_directional = directional$eta,
                    size_directional = .noncentral_size(directional$eta, alpha),
                    directional_alignment = directional$alignment,
                    # Deprecated v0.6 aliases retained for code compatibility.
                    eta_real_cr = directional$eta,
                    size_realized = .noncentral_size(directional$eta, alpha),
                    controls = controls, boot = boot)
  .new_AdequacyReport("twfe_heterogeneity", design, statistic, eta_selected,
                       eta_dag, eta_dag / CR[[selected_key]],
                       size_pt[[selected_key]], verdict,
                       alpha, delta, notes)
}

# ---- psi machinery (Paper C Def. def-psi) ------------------------------------
.rho_ar1 <- function(resid, uid, tid) {
  o <- order(uid, tid); u <- uid[o]; t <- tid[o]; e <- resid[o]
  j <- 2:length(o)
  keep <- u[j] == u[j - 1] & t[j] == t[j - 1] + 1
  den <- sum(e[j - 1][keep]^2)
  if (den > 0) sum(e[j][keep] * e[j - 1][keep]) / den else 0
}
.psi_parametric <- function(Dt, uid, tid, rho, kind = c("ar1", "exchangeable")) {
  kind <- match.arg(kind); n_w <- sum(Dt^2); nwcr <- 0
  for (i in unique(uid)) {
    sel <- uid == i; d <- Dt[sel]; t <- tid[sel]
    if (kind == "exchangeable") nwcr <- nwcr + (1 - rho) * sum(d^2) + rho * sum(d)^2
    else { R <- rho^abs(outer(t, t, "-")); nwcr <- nwcr + as.numeric(t(d) %*% R %*% d) }
  }
  nwcr / n_w
}
.psi_direct <- function(Dt, resid, uid, n_w, sigma2, G) {
  meat <- rowsum(Dt * resid, uid)[, 1L]
  (G / (G - 1)) * sum(meat^2) / (n_w * sigma2)
}
.psi_driven <- .psi_direct

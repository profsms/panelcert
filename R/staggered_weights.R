# Module C - staggered-DiD / TWFE-heterogeneity adequacy (Paper C).
#
# Design statistics (pre-outcome, from the adoption pattern ALONE):
#   Gamma   = sqrt(N1)||w - u|| (Def. def-gamma); Gamma = 0 iff block design.
#   Gamma_S = sqrt(N1)||Pi_S(w - u)|| restricted to a heterogeneity subspace S
#             (Prop. prop-restricted): cohort, event-time, additive combined;
#             Gamma_coh, Gamma_evt <= Gamma_c+e <= Gamma.
#   worst-case |eta| = (c/sigma) Gamma_S (Cor. cor-gamma); threshold
#             (c/sigma) Gamma_{S,CR} <= eta-dagger (Cor. cor-cv).
# Inference layer:
#   cluster   Gamma_{S,CR} = Gamma_S/sqrt(psi), psi on the realized design.
#   pilot     COVARIANCE-AWARE (eq-pilot): c_S^2 = n_w max{0, (1/N1)||Pi_S d||^2
#             - (1/N1) tr(Pi_S~ A Omega A' Pi_S~)}, Omega = Cov(hat Delta_{g,t})
#             from a fixed-design wild cluster bootstrap.
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

.restricted_gammas <- function(w, g_cell, e_cell, N1) {
  d <- w - 1 / N1
  Bcoh <- .indicator_basis(g_cell); Bevt <- .indicator_basis(e_cell)
  Bcmb <- cbind(Bcoh, Bevt)
  gam <- function(B) sqrt(N1) * sqrt(sum(.proj(B, d)^2))
  list(unr = sqrt(N1) * sqrt(sum(d^2)), coh = gam(Bcoh), evt = gam(Bevt),
       cmb = gam(Bcmb), Bcoh = Bcoh, Bevt = Bevt, Bcmb = Bcmb)
}

# ---- group-time ATTs (not-yet-treated difference-in-means at (g,t)) ----------
.group_time_atts <- function(uid, tid, ftc, y, N, T) {
  Ymat <- matrix(NA_real_, N, T)
  Ymat[cbind(uid, tid)] <- y
  cohorts <- sort(unique(ftc[is.finite(ftc) & ftc > 1]))
  keys <- list(); vals <- numeric(0)
  for (g in cohorts) {
    gunits <- which(ftc == g); base <- as.integer(g) - 1L
    for (t in as.integer(g):T) {
      ctrl <- which(ftc > t)                       # not-yet + never
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
       cmb = pil(bases$Bcmb, traces$cmb))
}

# ---- fixed-design wild cluster bootstrap -------------------------------------
.wild_bootstrap <- function(uid, tid, ftc, D, Dt, y, N, T, n_w, treated, N1,
                            g_cell, t_cell, d_K, gt0, B, seed) {
  n <- length(uid)
  Umat <- matrix(0, n, N); Umat[cbind(seq_len(n), uid)] <- 1
  Tmat <- matrix(0, n, T); Tmat[cbind(seq_len(n), tid)] <- 1
  gtlab <- rep("none", n)
  treated_obs <- which(D == 1)
  gtlab[treated_obs] <- paste(ftc[uid[treated_obs]], tid[treated_obs], sep = "_")
  GTmat <- .indicator_basis(gtlab)
  X <- cbind(Umat, Tmat, GTmat)
  yhat <- as.numeric(X %*% (.pinv(crossprod(X)) %*% (t(X) %*% y)))
  e <- y - yhat
  keys_nm <- if (length(gt0$keys)) sort(vapply(gt0$keys, function(k) paste(k, collapse = "_"), "")) else character(0)
  flat <- function(look_keys) look_keys
  dof <- n - d_K - 1
  set.seed(seed)
  sig <- numeric(0); psi <- numeric(0); gtm <- list()
  for (b in seq_len(B)) {
    v <- sample(c(-1, 1), N, replace = TRUE)
    ystar <- yhat + v[uid] * e
    yt <- .twoway_demean_codes(ystar, uid, tid, N, T)
    beta <- sum(Dt * yt) / n_w
    resid <- yt - beta * Dt
    s <- sqrt(sum(resid^2) / dof)
    rho <- .rho_ar1(resid, uid, tid)
    p <- .psi_parametric(Dt, uid, tid, rho, kind = "ar1")
    if (!is.finite(p) || p <= 0) next
    gt <- .group_time_atts(uid, tid, ftc, ystar, N, T)
    nm <- if (length(gt$keys)) vapply(gt$keys, function(k) paste(k, collapse = "_"), "") else character(0)
    vec <- stats::setNames(rep(NA_real_, length(keys_nm)), keys_nm)
    vec[nm] <- gt$vals
    sig <- c(sig, s); psi <- c(psi, p); gtm[[length(gtm) + 1L]] <- vec
  }
  G <- do.call(rbind, gtm)                # ndraw x m
  Omega <- stats::cov(G)
  list(Omega = Omega, sigma = sig, psi = psi, gtm = gtm, keys = keys_nm)
}

#' Pre-outcome design vetting (Paper C design-statistic ladder)
#'
#' The design statistic \code{Gamma = sqrt(N1)||w - u||} and its restricted
#' variants \code{Gamma_coh}, \code{Gamma_evt}, \code{Gamma_c+e}
#' (Prop. prop-restricted), the negative-weight share, and the breakdown ratio,
#' from the adoption pattern alone. Always-treated units are dropped.
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
  G <- .restricted_gammas(ds$w, g_cell, e_cell, ds$N1)
  fd <- fe_dimension(uid, tid, N, T)
  design <- structure(list(n = n, N = N, T = T, d_K = fd$d_K, rho = fd$d_K / n,
                           ncomponents = fd$ncomponents, tau_star2 = ds$n_w),
                      class = "DesignSummary")
  eta_dag <- .eta_dagger(alpha, delta)
  notes <- character(0)
  if (dr$n_drop > 0)
    notes <- c(notes, sprintf("%d always-treated unit(s) dropped (setup g >= 2)", dr$n_drop))
  statistic <- list(Gamma = ds$Gamma, Gamma_coh = G$coh, Gamma_evt = G$evt,
                    Gamma_cmb = G$cmb, neg_share = ds$neg_share, N1 = ds$N1, n_w = ds$n_w)
  if (ds$Gamma <= 1e-8) {
    notes <- c(notes, "block design: within-transformed treatment is uniform (Prop. prop-gamma0)")
    return(.new_AdequacyReport("twfe_heterogeneity", design, statistic, NULL,
                               eta_dag, Inf, NULL, "CERTIFIED", alpha, delta, notes))
  }
  breakdown <- eta_dag / G$cmb
  notes <- c(notes, sprintf("pre-outcome: naive TWFE inference is size-controlled iff c/sigma <= %.3f (= eta-dagger/Gamma_c+e); supply the outcome (twfe_adequacy) to pilot c/sigma", breakdown))
  .new_AdequacyReport("twfe_heterogeneity", design, statistic, NULL, eta_dag,
                      breakdown, NULL, "INCONCLUSIVE", alpha, delta, notes)
}

#' TWFE-heterogeneity adequacy (covariance-aware, wild-bootstrap)
#'
#' Restricted design-statistic ladder, cluster-robust rescaling
#' \code{Gamma_{S,CR} = Gamma_S/sqrt(psi)}, covariance-aware pilots
#' \code{c_S/sigma} (eq-pilot; Omega from a fixed-design wild cluster bootstrap),
#' the combined-class worst-case implied size (headline) with its bootstrap
#' median and interval, and the verdict. Always-treated units are dropped.
#' @param object outcome vector
#' @param unit,time,first_treat as in [twfe_design()]
#' @param alpha,delta level and size tolerance
#' @param cluster "ar1" (default) or "iid"
#' @param psi optional user-supplied variance-inflation factor
#' @param bootstrap number of wild-cluster draws (>0 enables the covariance
#'   correction and the size intervals)
#' @param seed RNG seed for the wild bootstrap
#' @param ... passed between methods
#' @return an \code{AdequacyReport}
#' @export
twfe_adequacy <- function(object, ...) UseMethod("twfe_adequacy")

#' @rdname twfe_adequacy
#' @export
twfe_adequacy.default <- function(object, unit, time, first_treat,
                                  alpha = 0.05, delta = 0.05,
                                  cluster = c("ar1", "iid"), psi = NULL,
                                  bootstrap = 999L, seed = 20260715L, ...) {
  y0 <- object; cluster <- match.arg(cluster)
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
  G <- .restricted_gammas(ds$w, g_cell, e_cell, ds$N1)
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
  psi_driven <- .psi_driven(Dt, resid, uid, ds$n_w, sigma^2, N)
  psi_hat <- if (!is.null(psi)) as.numeric(psi)
             else if (cluster == "iid") 1 else .psi_parametric(Dt, uid, tid, rho_ar1, "ar1")
  CR <- list(unr = G$unr / sqrt(psi_hat), coh = G$coh / sqrt(psi_hat),
             evt = G$evt / sqrt(psi_hat), cmb = G$cmb / sqrt(psi_hat))

  gt0 <- .group_time_atts(uid, tid, ftc, y, N, T)
  look0 <- .gt_lookup(gt0)
  cm0 <- .cohort_means(gt0)
  cell0 <- unname(cm0[as.character(g_cell)])
  good0 <- !is.na(cell0)
  att_bar <- if (sum(good0) >= 1) mean(cell0[good0]) else NA_real_
  eta_real_iid <- if (sum(good0) >= 2)
    (sum(ds$w[good0] * cell0[good0]) / sum(ds$w[good0]) - att_bar) * sqrt(ds$n_w) / sigma else NA_real_
  eta_real_cr <- eta_real_iid / sqrt(psi_hat)

  bases <- list(Bcoh = G$Bcoh, Bevt = G$Bevt, Bcmb = G$Bcmb)
  boot <- NULL
  if (bootstrap > 0 && length(gt0$keys) > 0) {
    wb <- .wild_bootstrap(uid, tid, ftc, D, Dt, y, N, T, ds$n_w, ds$treated,
                          ds$N1, g_cell, t_cell, fd$d_K, gt0, as.integer(bootstrap),
                          as.integer(seed))
    keys_nm <- wb$keys
    kmat <- do.call(rbind, strsplit(keys_nm, "_"))
    A <- matrix(0, ds$N1, length(keys_nm))
    cellkey <- paste(as.integer(g_cell), t_cell, sep = "_")
    for (k in seq_len(ds$N1)) { j <- match(cellkey[k], keys_nm); if (!is.na(j)) A[k, j] <- 1 }
    Jm <- matrix(1 / ds$N1, ds$N1, ds$N1)
    Pt <- list(coh = .hat(G$Bcoh) - Jm, evt = .hat(G$Bevt) - Jm, cmb = .hat(G$Bcmb) - Jm)
    M <- A %*% wb$Omega %*% t(A)
    traces <- list(coh = sum(Pt$coh * M) / ds$N1, evt = sum(Pt$evt * M) / ds$N1,
                   cmb = sum(Pt$cmb * M) / ds$N1)
    pilots <- .cov_pilots(look0, g_cell, t_cell, ds$N1, bases, traces, ds$n_w, sigma)
    size_pt <- list(coh = .noncentral_size(pilots$coh * CR$coh, alpha),
                    evt = .noncentral_size(pilots$evt * CR$evt, alpha),
                    cmb = .noncentral_size(pilots$cmb * CR$cmb, alpha))
    draw_cmb <- numeric(0)
    for (i in seq_along(wb$sigma)) {
      vec <- wb$gtm[[i]]
      look <- function(g, t) { z <- vec[paste(as.integer(g), t, sep = "_")]; if (is.na(z)) NA_real_ else unname(z) }
      p <- .cov_pilots(look, g_cell, t_cell, ds$N1, bases, traces, ds$n_w, wb$sigma[i])
      draw_cmb <- c(draw_cmb, .noncentral_size(p$cmb * G$cmb / sqrt(wb$psi[i]), alpha))
    }
    q <- function(v, p) stats::quantile(v[is.finite(v)], p, names = FALSE)
    boot <- list(n = length(wb$sigma), cmb_med = q(draw_cmb, 0.5),
                 cmb_lo = q(draw_cmb, 0.025), cmb_hi = q(draw_cmb, 0.975),
                 psi_lo = q(wb$psi, 0.025), psi_hi = q(wb$psi, 0.975),
                 Omega_trace_cmb = traces$cmb)
    notes <- c(notes, sprintf("covariance-aware pilot (eq-pilot): Omega from %d wild-cluster draws", boot$n))
  } else {
    rawpil <- function(B) {
      delta <- .delta_cell(look0, g_cell, t_cell, ds$N1)
      cov <- !is.na(delta); dc <- delta - mean(delta[cov])
      comp <- .proj(B[cov, , drop = FALSE], dc[cov])
      sqrt(ds$n_w * sum(comp^2) / sum(cov) / sigma^2)
    }
    pilots <- list(coh = rawpil(G$Bcoh), evt = rawpil(G$Bevt), cmb = rawpil(G$Bcmb))
    size_pt <- list(coh = .noncentral_size(pilots$coh * CR$coh, alpha),
                    evt = .noncentral_size(pilots$evt * CR$evt, alpha),
                    cmb = .noncentral_size(pilots$cmb * CR$cmb, alpha))
    notes <- c(notes, "bootstrap disabled: pilots are RAW projected dispersion, upward-biased (set bootstrap>0)")
  }

  if (psi_hat != 1) {
    dir <- if (psi_hat > 1) "clustering shrinks the non-centrality here; iid size is an upper bound"
           else "clustering WORSENS the distortion here (psi < 1)"
    notes <- c(notes, sprintf("psi_hat = %.3f (AR(1) rho = %.3f): %s", psi_hat, rho_ar1, dir))
  }
  if (psi_hat > 0 && max(psi_driven / psi_hat, psi_hat / psi_driven) > 1.5)
    notes <- c(notes, sprintf("estimator-driven cross-check psi = %.2f vs parametric %.2f", psi_driven, psi_hat))
  if (N < 40)
    notes <- c(notes, sprintf("few clusters (G = %d): wild bootstrap refinements advisable (Rem. sec-clusters)", N))

  eta_dag <- .eta_dagger(alpha, delta)
  verdict <- if (size_pt$cmb <= alpha + delta) "CERTIFIED" else "FLAGGED"
  design <- structure(list(n = n, N = N, T = T, d_K = fd$d_K, rho = fd$d_K / n,
                           ncomponents = fd$ncomponents, tau_star2 = ds$n_w),
                      class = "DesignSummary")
  statistic <- list(Gamma = ds$Gamma, Gamma_coh = G$coh, Gamma_evt = G$evt,
                    Gamma_cmb = G$cmb, neg_share = ds$neg_share, psi_hat = psi_hat,
                    psi_driven = psi_driven, rho_ar1 = rho_ar1, Gamma_CR = CR$unr,
                    Gamma_coh_CR = CR$coh, Gamma_evt_CR = CR$evt, Gamma_cmb_CR = CR$cmb,
                    beta = beta, sigma = sigma, att_bar = att_bar, N1 = ds$N1,
                    n_w = ds$n_w, n_cohorts = length(gt0$cohorts),
                    pilot_coh = pilots$coh, pilot_evt = pilots$evt, pilot_cmb = pilots$cmb,
                    size_coh = size_pt$coh, size_evt = size_pt$evt, size_cmb = size_pt$cmb,
                    size_realized = .noncentral_size(eta_real_cr, alpha),
                    eta_real_cr = eta_real_cr, boot = boot)
  .new_AdequacyReport("twfe_heterogeneity", design, statistic, pilots$cmb * CR$cmb,
                      eta_dag, eta_dag / CR$cmb, size_pt$cmb, verdict, alpha, delta, notes)
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
.psi_driven <- function(Dt, resid, uid, n_w, sigma2, G) {
  meat <- rowsum(Dt * resid, uid)[, 1L]
  (G / (G - 1)) * sum(meat^2) / (n_w * sigma2)
}

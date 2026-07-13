# Module C \u2014 staggered-DiD / TWFE-heterogeneity adequacy (spec section 5;
# Paper C). Formula sources: Gamma = sqrt(N1)||w - u|| (Def. def-gamma);
# worst-case |eta| = (c/sigma) Gamma (Cor. cor-gamma); adequacy threshold
# (c/sigma) Gamma_CR <= eta\u2020 (Cor. cor-cv); SHRINKAGE-CORRECTED dispersion
# var_shrunk = max(0, var(ATT_g) - mean(se_g^2)) (eq-shrink); cluster layer
# psi = sum_i d_i' R_i d_i / n_w, Gamma_CR = Gamma/sqrt(psi) (Thm. thm-cluster)
# with psi COMPUTED on the realized design, direction never asserted.
#
# API-FREEZE CAVEAT (spec section 5.3): the cluster-robust inference layer
# follows the paper_c_twfe_v1 DRAFT; its inputs stay flexible until Paper C is
# final. Everything downstream of psi_hat carries a note in the report.

.treatment_indicator <- function(time, first_treat) {
  as.numeric(!is.na(first_treat) & time >= first_treat)
}

# Sorted time codes plus per-unit cohort codes: 0 = always-treated (adopted
# before the sample; D = 1 but not a cohort), Inf = never treated in sample,
# else the 1-based code of the adoption period.
.staggered_codes <- function(unit, time, first_treat) {
  n <- length(unit)
  if (length(time) != n || length(first_treat) != n)
    stop("unit, time, first_treat must have equal length")
  uid <- match(unit, unique(unit))
  stimes <- sort(unique(time))
  tid <- match(time, stimes)
  N <- max(uid); T <- length(stimes)

  ft_of <- rep(NA_real_, N)
  seen <- logical(N)
  for (k in seq_len(n)) {
    f <- if (is.na(first_treat[k])) NA_real_ else as.numeric(first_treat[k])
    if (!seen[uid[k]]) {
      ft_of[uid[k]] <- f
      seen[uid[k]] <- TRUE
    } else if (!identical(is.na(ft_of[uid[k]]), is.na(f)) ||
               (!is.na(f) && ft_of[uid[k]] != f)) {
      stop("first_treat varies within unit ", unit[k])
    }
  }
  ftc <- numeric(N)
  for (i in seq_len(N)) {
    f <- ft_of[i]
    if (is.na(f) || !is.finite(f) || f > stimes[T]) {
      ftc[i] <- Inf
    } else if (f <= stimes[1]) {
      ftc[i] <- 0
    } else {
      pos <- which(as.numeric(stimes) == f)
      if (length(pos) == 0)
        stop("first_treat value ", f, " is not an observed period")
      ftc[i] <- pos
    }
  }
  list(uid = uid, tid = tid, N = N, T = T, ftc = ftc)
}

.design_stats <- function(D, Dt) {
  treated <- which(D == 1)
  N1 <- length(treated)
  if (N1 < 2) stop("fewer than 2 treated cells \u2014 no staggered design")
  n_w <- sum(Dt^2)
  if (n_w <= 1e-12)
    stop("treatment has no within variation (single common adoption date?)")
  w <- Dt[treated] / sum(Dt[treated])
  Gamma <- sqrt(N1) * sqrt(sum((w - 1 / N1)^2))
  list(treated = treated, N1 = N1, n_w = n_w, w = w, Gamma = Gamma,
       neg_share = sum(w < 0) / N1)
}

#' Module C pre-outcome design vetting (Paper C's design statistic)
#'
#' The design statistic `Gamma = sqrt(N1)||w - u||`, the negative-weight share,
#' and the breakdown heterogeneity-to-noise ratio `(c/sigma)-dagger = eta-dagger/Gamma`,
#' computed from the adoption pattern ALONE --- before any outcome exists.
#' `first_treat` is the unit's adoption time on the same scale as `time`
#' (`NA` = never treated; values before the sample = always-treated).
#'
#' @param unit,time raw identifiers (time numeric)
#' @param first_treat adoption time per observation (constant within unit)
#' @param alpha,delta level and size tolerance
#' @return an object of class \code{AdequacyReport} (CERTIFIED for block
#'   designs, else INCONCLUSIVE with the breakdown heterogeneity-to-noise
#'   ratio)
#' @examples
#' # a staggered three-cohort design with no never-treated reservoir
#' N <- 48; T <- 12
#' unit <- rep(1:N, each = T); time <- rep(1:T, times = N)
#' ft <- ifelse(unit <= 16, 3, ifelse(unit <= 32, 7, 11))
#' twfe_design(unit, time, ft)   # Gamma = 1.54, 11% negative weights
#' @export
twfe_design <- function(unit, time, first_treat, alpha = 0.05, delta = 0.05) {
  sc <- .staggered_codes(unit, time, first_treat)
  n <- length(sc$uid)
  D <- .treatment_indicator(time, first_treat)
  Dt <- .twoway_demean_codes(D, sc$uid, sc$tid, sc$N, sc$T)
  ds <- .design_stats(D, Dt)
  fd <- fe_dimension(sc$uid, sc$tid, sc$N, sc$T)
  design <- structure(list(n = n, N = sc$N, T = sc$T, d_K = fd$d_K,
                           rho = fd$d_K / n, ncomponents = fd$ncomponents,
                           tau_star2 = ds$n_w), class = "DesignSummary")
  eta_dag <- .eta_dagger(alpha, delta)
  statistic <- list(Gamma = ds$Gamma, neg_share = ds$neg_share, N1 = ds$N1,
                    n_w = ds$n_w)
  if (ds$Gamma <= 1e-8) {
    notes <- "block design: within-transformed treatment is uniform on treated cells, so NO heterogeneity profile can distort the t-test (eta = 0 for every profile; Paper C, Prop. prop-gamma0)"
    return(.new_AdequacyReport("twfe_heterogeneity", design, statistic, NULL,
                               eta_dag, Inf, NULL, "CERTIFIED", alpha, delta,
                               notes))
  }
  breakdown <- eta_dag / ds$Gamma
  notes <- sprintf(
    "pre-outcome design statistic: naive TWFE inference is size-controlled iff the heterogeneity-to-noise ratio c/sigma <= %.3f (= eta\u2020/Gamma; Paper C, Cor. cor-cv) \u2014 supply the outcome (twfe_adequacy) to pilot c/sigma",
    breakdown)
  .new_AdequacyReport("twfe_heterogeneity", design, statistic, NULL, eta_dag,
                      breakdown, NULL, "INCONCLUSIVE", alpha, delta, notes)
}

# Simplified not-yet-treated difference-in-means cohort pilot. Deliberately
# NOT a Callaway-Sant'Anna implementation (spec section 0 scope discipline).
.cohort_pilot <- function(uid, tid, ftc, y, N, T) {
  Ymat <- matrix(NA_real_, N, T)
  Ymat[cbind(uid, tid)] <- y
  cohorts <- sort(unique(ftc[is.finite(ftc) & ftc > 1]))
  atts <- numeric(0); ses <- numeric(0); gs <- numeric(0)
  for (g in cohorts) {
    gunits <- which(ftc == g)
    base <- as.integer(g) - 1L
    diffs <- numeric(0)
    for (t in as.integer(g):T) {
      ctrl <- which(ftc > t | ftc <= 0)  # not-yet + never + always-treated
      gt <- Ymat[gunits, t] - Ymat[gunits, base]
      ct <- Ymat[ctrl, t] - Ymat[ctrl, base]
      gt <- gt[!is.na(gt)]; ct <- ct[!is.na(ct)]
      if (length(gt) == 0 || length(ct) == 0) next
      diffs <- c(diffs, mean(gt) - mean(ct))
    }
    if (length(diffs) == 0) next
    gs <- c(gs, g)
    atts <- c(atts, mean(diffs))
    ses <- c(ses, if (length(diffs) > 1) stats::sd(diffs) / sqrt(length(diffs)) else NA_real_)
  }
  list(gs = gs, atts = atts, ses = ses)
}

#' Module C inference layer: TWFE-heterogeneity adequacy
#'
#' Realized and worst-case non-centrality of the TWFE t-test under
#' treatment-effect heterogeneity, with the SHRINKAGE-CORRECTED
#' cohort-dispersion pilot (Paper C eq-shrink) and the cluster-robust rescaling
#' `Gamma_CR = Gamma/sqrt(psi_hat)`, `psi_hat` computed on the realized design
#' (no direction assumed).
#'
#' @param object outcome vector (default method), or a fitted \code{fixest} /
#'   \code{plm} / \code{lm} TWFE model of the outcome on the treatment dummy
#' @param unit,time,first_treat as in [twfe_design()]; for the model methods,
#'   \code{first_treat} must align with the ESTIMATION sample (same length and
#'   order as the data the model was fitted on, after any dropped rows)
#' @param alpha nominal test level
#' @param delta size-distortion tolerance
#' @param cluster \code{"ar1"} (default; Paper C's parametric route) or
#'   \code{"iid"}
#' @param psi optional user-supplied variance-inflation factor (overrides
#'   \code{cluster})
#' @param cohort_effects,cohort_ses optional cohort-level effect estimates
#'   (aligned with the SORTED adoption times) from a heterogeneity-robust
#'   estimator (did / csdid / fixest::sunab). If omitted, an internal
#'   order-of-magnitude pilot is used and labelled as such.
#' @param ... passed between methods
#' @return an object of class \code{AdequacyReport}
#' @references Halkiewicz, S. M. S. Stock--Yogo critical values for the
#'   two-way fixed-effects t-test under treatment-effect heterogeneity
#'   (Paper C).
#' @examples
#' N <- 30; T <- 8
#' unit <- rep(1:N, each = T); time <- rep(1:T, times = N)
#' ft <- ifelse(unit <= 8, 3, ifelse(unit <= 16, 6, NA))
#' y <- rnorm(N * T) + 0.5 * as.numeric(!is.na(ft) & time >= ft)
#' twfe_adequacy(y, unit, time, ft)
#' @export
twfe_adequacy <- function(object, ...) UseMethod("twfe_adequacy")

#' @rdname twfe_adequacy
#' @export
twfe_adequacy.default <- function(object, unit, time, first_treat,
                                  alpha = 0.05, delta = 0.05,
                                  cluster = c("ar1", "iid"), psi = NULL,
                                  cohort_effects = NULL, cohort_ses = NULL,
                                  ...) {
  y <- object
  cluster <- match.arg(cluster)
  sc <- .staggered_codes(unit, time, first_treat)
  uid <- sc$uid; tid <- sc$tid; N <- sc$N; T <- sc$T; ftc <- sc$ftc
  n <- length(uid)
  if (length(y) != n) stop("y must have length n = ", n)
  D <- .treatment_indicator(time, first_treat)
  Dt <- .twoway_demean_codes(D, uid, tid, N, T)
  ds <- .design_stats(D, Dt)
  fd <- fe_dimension(uid, tid, N, T)
  dof <- n - fd$d_K - 1
  if (dof <= 0) stop("no residual degrees of freedom")

  yt <- .twoway_demean_codes(y, uid, tid, N, T)
  beta <- sum(Dt * yt) / ds$n_w
  resid <- yt - beta * Dt
  sigma <- sqrt(sum(resid^2) / dof)

  notes <- character(0)

  # ---- cohort-effect pilot ----
  cohorts <- sort(unique(ftc[is.finite(ftc) & ftc > 1]))
  if (is.null(cohort_effects)) {
    cp <- .cohort_pilot(uid, tid, ftc, y, N, T)
    atts <- cp$atts; ses <- cp$ses; gs <- cp$gs
    notes <- c(notes, "cohort-effect pilot is a simplified not-yet-treated difference-in-means \u2014 an order-of-magnitude pilot, NOT a Callaway\u2013Sant'Anna implementation; pass cohort_effects from csdid/did/sunab for a paper-grade pilot (spec \u00a70 scope discipline)")
  } else {
    if (length(cohort_effects) != length(cohorts))
      stop("cohort_effects must have one entry per adoption cohort (",
           length(cohorts), ", sorted by adoption time)")
    atts <- as.numeric(cohort_effects)
    if (is.null(cohort_ses)) {
      ses <- rep(0, length(atts))
      notes <- c(notes, "no cohort SEs supplied: cohort dispersion NOT shrinkage-corrected (assumed noiseless)")
    } else {
      if (length(cohort_ses) != length(atts))
        stop("cohort_ses must match cohort_effects")
      ses <- as.numeric(cohort_ses)
    }
    gs <- cohorts
  }
  good <- !is.na(ses)
  att_of <- stats::setNames(atts, as.character(gs))

  # realized eta: cell-level cohort deviations under the dCdH weights
  cell_att <- unname(att_of[as.character(ftc[uid[ds$treated]])])
  gcell <- !is.na(cell_att)
  att_bar <- mean(cell_att[gcell])
  eta_real_iid <- (sum(ds$w[gcell] * cell_att[gcell]) / sum(ds$w[gcell]) - att_bar) *
    sqrt(ds$n_w) / sigma

  # shrinkage-corrected dispersion (eq-shrink) \u2014 REQUIRED default (spec 5.2)
  if (sum(good) >= 2) {
    raw_var <- stats::var(atts[good])
    noise <- mean(ses[good]^2)
    shr_var <- max(0, raw_var - noise)
    csd_raw <- sqrt(raw_var); csd_shr <- sqrt(shr_var)
    shrink_factor <- if (raw_var > 0) sqrt(shr_var / raw_var) else 1
  } else {
    csd_raw <- 0; csd_shr <- 0; shrink_factor <- 1
    notes <- c(notes, "fewer than two cohorts with usable SEs: cohort dispersion not estimable; worst-case eta unavailable")
  }
  if (csd_shr == 0 && csd_raw > 0) {
    notes <- c(notes, sprintf(
      "cohort dispersion (raw sd %.4g) is entirely attributable to sampling noise \u2014 shrunk sd = 0; worst-case eta at the shrunk pilot is 0 (raw-pilot worst-case reported for comparison; Paper C eq-shrink)",
      csd_raw))
  } else if (csd_raw > 0) {
    notes <- c(notes, sprintf(
      "shrinkage-corrected cohort dispersion: raw sd %.4g -> shrunk %.4g (factor %.3f) \u2014 heterogeneity %s (Paper C eq-shrink)",
      csd_raw, csd_shr, shrink_factor,
      if (shrink_factor > 0.9) "genuine, not a small-cohort artifact" else "partly sampling noise"))
  }

  # ---- cluster layer (spec 5.3 \u2014 API flexible until Paper C final) ----
  rho_ar1 <- .rho_ar1(resid, uid, tid)
  psi_driven <- .psi_driven(Dt, resid, uid, ds$n_w, sigma^2, N)
  psi_hat <- if (!is.null(psi)) as.numeric(psi)
             else if (cluster == "iid") 1
             else .psi_parametric(Dt, uid, tid, rho_ar1, kind = "ar1")
  Gamma_CR <- ds$Gamma / sqrt(psi_hat)

  eta <- eta_real_iid / sqrt(psi_hat)
  eta_worst <- (csd_shr * sqrt(ds$n_w) / sigma) * Gamma_CR
  eta_worst_raw <- (csd_raw * sqrt(ds$n_w) / sigma) * Gamma_CR

  if (psi_hat != 1) {
    dir <- if (psi_hat > 1)
      "clustering shrinks the non-centrality on THIS design; the iid implied size is an upper bound here"
    else "clustering WORSENS the distortion on THIS design (psi < 1)"
    notes <- c(notes, sprintf(
      "realized variance-inflation psi_hat = %.3f (AR(1) fit rho = %.3f): %s \u2014 direction computed from the design, never assumed (Paper C \u00a75.3)",
      psi_hat, rho_ar1, dir))
  }
  if (psi_hat > 0 && max(psi_driven / psi_hat, psi_hat / psi_driven) > 1.5)
    notes <- c(notes, sprintf(
      "estimator-driven cross-check psi = %.2f diverges from the parametric %.2f; route (ii) is only valid when the homoskedastic sigma_hat is consistent (Paper C \u00a7sec-feasible) \u2014 parametric route reported",
      psi_driven, psi_hat))
  if (N < 40)
    notes <- c(notes, sprintf(
      "few clusters (G = %d): feasible CR1 has its own finite-cluster error; CR3/jackknife or wild cluster bootstrap refinements are advisable (Paper C, Rem. sec-clusters)",
      N))
  notes <- c(notes, "cluster-robust layer follows the paper_c_twfe_v1 DRAFT; feasible CR1-t validity requires moderate-to-large (G, T), so clustered sizes are design-computed approximations \u2014 this layer's API stays flexible until Paper C is final (spec \u00a75.3)")

  eta_dag <- .eta_dagger(alpha, delta)
  verdict <- if (max(abs(eta), eta_worst) <= eta_dag) "CERTIFIED" else "FLAGGED"
  design <- structure(list(n = n, N = N, T = T, d_K = fd$d_K, rho = fd$d_K / n,
                           ncomponents = fd$ncomponents, tau_star2 = ds$n_w),
                      class = "DesignSummary")
  statistic <- list(Gamma = ds$Gamma, neg_share = ds$neg_share,
                    Gamma_CR = Gamma_CR, psi_hat = psi_hat,
                    psi_driven = psi_driven, rho_ar1 = rho_ar1, beta = beta,
                    sigma = sigma, att_bar = att_bar, cohort_sd_raw = csd_raw,
                    cohort_sd_shrunk = csd_shr, shrink_factor = shrink_factor,
                    eta_worst = eta_worst, eta_worst_raw = eta_worst_raw,
                    N1 = ds$N1, n_w = ds$n_w, n_cohorts = length(cohorts))
  .new_AdequacyReport("twfe_heterogeneity", design, statistic, eta, eta_dag,
                      eta_dag / Gamma_CR, .noncentral_size(eta, alpha), verdict,
                      alpha, delta, notes)
}

# ---- psi machinery (Paper C Def. def-psi, section sec-feasible) ----

# Pooled within-unit lag-1 residual autocorrelation (consecutive periods only).
.rho_ar1 <- function(resid, uid, tid) {
  o <- order(uid, tid)
  u <- uid[o]; t <- tid[o]; e <- resid[o]
  j <- 2:length(o)
  keep <- u[j] == u[j - 1] & t[j] == t[j - 1] + 1
  num <- sum(e[j][keep] * e[j - 1][keep])
  den <- sum(e[j - 1][keep]^2)
  if (den > 0) num / den else 0
}

# psi = sum_i d_i' R_i d_i / n_w for parametric R.
.psi_parametric <- function(Dt, uid, tid, rho, kind = c("ar1", "exchangeable")) {
  kind <- match.arg(kind)
  n_w <- sum(Dt^2)
  nwcr <- 0
  for (i in unique(uid)) {
    sel <- uid == i
    d <- Dt[sel]; t <- tid[sel]
    if (kind == "exchangeable") {
      # d'Rd = (1-rho)||d||^2 + rho (d'1)^2; d'1 = 0 after the within
      # transform (Paper C, Lem. lem-center), so this is exact
      nwcr <- nwcr + (1 - rho) * sum(d^2) + rho * sum(d)^2
    } else {
      R <- rho^abs(outer(t, t, "-"))
      nwcr <- nwcr + as.numeric(t(d) %*% R %*% d)
    }
  }
  nwcr / n_w
}

# Estimator-driven cross-check: CR1 sandwich over homoskedastic sigma^2.
.psi_driven <- function(Dt, resid, uid, n_w, sigma2, G) {
  meat <- rowsum(Dt * resid, uid)[, 1L]
  (G / (G - 1)) * sum(meat^2) / (n_w * sigma2)
}

# Measurement-error adequacy for "Breakdown Reliability for Saturated Fixed-Effect Inference".
# Formula sources: feasible non-centrality (Cor. cor-feasible); corrected pilot
# beta*/lambda (Prop. prop-pilot); pilot se s_CR/(lambda sqrt(tau*2)) with
# s_CR = sigma sqrt(psi) \u2014 the CLUSTER-ROBUST scale, not the i.i.d. sigma
# (Rem. rem-plugin, Prop. prop-certificate); exact-inversion threshold eta\u2020
# with the quadratic closed form as companion \u2014 NEVER the discarded linear
# surrogate (Rem. rem-exact-cv); FIXED-POINT breakdown reliability
# lambda\u2020 = t*/(t* + eta\u2020), t* = |beta*|sqrt(tau*2)/sigma
# (Def. def-breakdown); cluster layer psi_hat deflating |eta| and t* by
# sqrt(psi) (Rem. rem-cluster); exact non-central implied size.
#
# NOTE on the pilot se: Cor. cor-slope gives the i.i.d. law
# beta0_corr => beta0 + N(0, sigma^2/(lambda^2 tau*2)). Under clustering the
# score variance is Psi_n ->p psi sigma^2 tau*2, so the correct limit carries
# psi: N(0, psi sigma^2/(lambda^2 tau*2)). This file has always used s_CR and
# is therefore right; the manuscript's Remark rem-plugin cited the i.i.d.
# variance until the 2026-07-31 audit and was corrected to match this code.

#' Measurement-error SDs from published credible-interval bounds
#'
#' V-Dem convention: the interval brackets one posterior SD, so
#' `sigma_nu = (codehigh - codelow)/2` (lead application in the accompanying article).
#'
#' @param codelow,codehigh interval bounds, one per observation
#' @return numeric vector of per-observation measurement-error SDs
#' @examples
#' reliability_from_interval(c(0.20, 0.31), c(0.24, 0.35))
#' @export
reliability_from_interval <- function(codelow, codehigh) {
  if (length(codelow) != length(codehigh))
    stop("codelow and codehigh must have equal length")
  if (any(!is.finite(codelow)) || any(!is.finite(codehigh)))
    stop("interval bounds must be finite")
  if (any(codehigh < codelow))
    stop("codehigh must be at least codelow for every observation")
  (codehigh - codelow) / 2
}

#' Measurement-error SD implied by an external reliability ratio
#'
#' With the documented default `scale = "observed"`, `within_sd` is the SD of
#' the observed regressor and `sigma_nu = sqrt(1-r) * within_sd`. Set
#' `scale = "signal"` only when `within_sd` is the latent-signal SD; then
#' `sigma_nu = sqrt((1-r)/r) * within_sd`.
#'
#' @param r reliability ratio in (0, 1]
#' @param within_sd within-SD on the scale selected by `scale`
#' @param scale `"observed"` (default) or `"signal"`
#' @return scalar measurement-error SD
#' @examples
#' reliability_from_ratio(0.75, within_sd = 0.15)
#' @export
reliability_from_ratio <- function(r, within_sd,
                                   scale = c("observed", "signal")) {
  scale <- match.arg(scale)
  if (!(r > 0 && r <= 1)) stop("reliability ratio must be in (0, 1]")
  if (!is.finite(within_sd) || within_sd < 0)
    stop("within_sd must be finite and non-negative")
  factor <- if (scale == "observed") sqrt(1 - r) else sqrt((1 - r) / r)
  factor * within_sd
}

#' Reliability from two measurements of the same regressor
#'
#' Apply the identical projection and complete-case sample to both measurements
#' before calling this helper. The default covariance estimator is
#' `cov(first_measure, second_measure) / var(first_measure)` and identifies the
#' first measurement's reliability when the two classical reporting errors are
#' uncorrelated; it does not require equal error variances. The
#' `"equal_variance"` estimator is
#' `1 - var(first_measure - second_measure) / (2 * var(first_measure))` and
#' additionally imposes equal reporting-error variances.
#'
#' Correlation across reporting errors invalidates both formulas and must be
#' modeled or examined as a sensitivity. The estimate is returned without
#' truncation so assumption or sampling failures remain visible;
#' [eiv_adequacy()] requires a reliability in `(0, 1]`.
#'
#' @param first_measure,second_measure numeric vectors containing two projected
#'   measurements on the same complete-case sample.
#' @param method `"covariance"` (default) or `"equal_variance"`.
#' @return A scalar estimate of the first measurement's reliability.
#' @examples
#' x <- c(-2, -1, 1, 2)
#' z <- c(-1.8, -1.2, 0.9, 2.1)
#' reliability_from_repeats(x, z)
#' reliability_from_repeats(x, z, method = "equal_variance")
#' @export
reliability_from_repeats <- function(first_measure, second_measure,
                                     method = c("covariance", "equal_variance")) {
  method <- match.arg(method)
  if (length(first_measure) != length(second_measure))
    stop("first_measure and second_measure must have equal length")
  if (length(first_measure) < 2L)
    stop("repeated measurements must contain at least two observations")
  if (!is.numeric(first_measure) || !is.numeric(second_measure) ||
      any(!is.finite(first_measure)) || any(!is.finite(second_measure)))
    stop("repeated measurements must be finite numeric vectors on a common complete-case sample")
  x <- first_measure - mean(first_measure)
  z <- second_measure - mean(second_measure)
  xx <- sum(x^2)
  if (xx <= 0) stop("first_measure has zero variance")
  if (method == "covariance") return(sum(x * z) / xx)
  1 - sum((x - z)^2) / (2 * xx)
}

#' Self-consistent breakdown reliability
#'
#' The fixed point \eqn{\lambda = \lambda^\dagger(\lambda)} when the corrected
#' pilot \eqn{\beta^*/\lambda} is evaluated at the reliability being solved
#' for (Definition def-breakdown in the accompanying article). Closed form
#' \eqn{\lambda^\dagger = t^*/(t^* + \eta^\dagger)} with
#' \eqn{t^* = |\beta^*|\sqrt{\tau^{*2}}/\sigma} --- the specification's
#' conventional t-statistic, so no noise input is needed. Under cluster-robust
#' standardization pass the variance-inflation factor \code{psi} (Remark
#' rem-cluster): \eqn{t^*} is deflated by \eqn{\sqrt{\psi}}. Uses the
#' exact-inversion root \eqn{\eta^\dagger}.
#'
#' @param beta_star attenuated slope from the FE regression
#' @param sigma residual standard deviation
#' @param tau_star2 observed within variation of the regressor
#' @param alpha nominal test level
#' @param delta size-distortion tolerance
#' @param psi cluster variance-inflation factor (1 = i.i.d. standardization)
#' @return the breakdown reliability in [0, 1]
#' @examples
#' breakdown_reliability(0.0908759, 0.5095233, 537.2959) # twins: ~0.864
#' @export
breakdown_reliability <- function(beta_star, sigma, tau_star2,
                                  alpha = 0.05, delta = 0.05, psi = 1) {
  if (!is.finite(sigma) || sigma <= 0) stop("sigma must be positive")
  if (!is.finite(tau_star2) || tau_star2 <= 0)
    stop("tau_star2 must be positive")
  if (beta_star == 0) return(0)
  if (!is.finite(psi) || psi <= 0) stop("psi must be positive")
  eta_dag <- .eta_dagger(alpha, delta)
  t_star <- abs(beta_star) * sqrt(tau_star2) / (sigma * sqrt(psi))
  t_star / (t_star + eta_dag)
}

#' Certified breakdown reliability
#'
#' Replaces the reported absolute t-statistic by
#' \eqn{|t| + z_{1-\gamma_\beta}} in the point-breakdown formula:
#' \deqn{\lambda^\dagger_{\gamma_\beta}=
#' (|t|+z_{1-\gamma_\beta})/(|t|+z_{1-\gamma_\beta}+\eta^\dagger).}
#' Compare the result with a known/consistent within reliability or a lower
#' confidence bound. If that lower bound has coverage error
#' \eqn{\gamma_\lambda}, the total false-certification bound is
#' \eqn{\gamma_\beta+\gamma_\lambda}; the second budget belongs to the bound,
#' not to this threshold formula.
#'
#' @inheritParams breakdown_reliability
#' @param gamma_beta coefficient-uncertainty error budget in (0, 0.5]
#' @return the certified breakdown reliability in [0, 1]
#' @examples
#' certified_breakdown_reliability(0.0908759, 0.5095233, 537.2959)
#' @export
certified_breakdown_reliability <- function(beta_star, sigma, tau_star2,
                                             alpha = 0.05, delta = 0.05,
                                             gamma_beta = 0.05, psi = 1) {
  if (!is.finite(sigma) || sigma <= 0) stop("sigma must be positive")
  if (!is.finite(tau_star2) || tau_star2 <= 0)
    stop("tau_star2 must be positive")
  if (!is.finite(psi) || psi <= 0) stop("psi must be positive")
  if (!(is.numeric(gamma_beta) && length(gamma_beta) == 1L &&
        is.finite(gamma_beta) && gamma_beta > 0 && gamma_beta <= 0.5))
    stop("gamma_beta must be a single number in (0, 0.5]")
  eta_dag <- .eta_dagger(alpha, delta)
  t_reported <- abs(beta_star) * sqrt(tau_star2) / (sigma * sqrt(psi))
  tz <- t_reported + stats::qnorm(1 - gamma_beta)
  tz / (tz + eta_dag)
}

#' Measurement-error adequacy diagnostic
#'
#' Certifies whether naive inference on an already-estimated FE regression is
#' size-controlled under classical measurement error in the regressor. The
#' regression is reproduced via Frisch--Waugh from raw data (default method) or
#' consumed directly from a fitted \pkg{fixest}, \pkg{plm}, or \code{lm} object
#' --- the model is never re-specified. Supply exactly one noise input:
#' \code{sigma_nu} (per-observation or scalar measurement-error SD),
#' \code{codelow} + \code{codehigh} (V-Dem-style posterior interval bounds), or
#' \code{reliability} (the within reliability \eqn{\hat\lambda} directly).
#'
#' Correct-by-default honesty machinery: \code{pilot = "conservative"}
#' (default) compares a reliability lower bound with the certified breakdown
#' obtained by replacing \eqn{|t|} by \eqn{|t|+z_{1-\gamma_\beta}};
#' \code{"point"} is
#' the descriptive corrected-pilot verdict; \code{"naive"} plugs in the
#' attenuated \eqn{\hat\beta^*} and is anti-conservative --- exposed for
#' comparison only, and labelled as such in the report.
#'
#' @param object outcome vector (default method), or a fitted \code{fixest} /
#'   \code{plm} / \code{lm} model with one regressor and two-way fixed effects
#' @param x observed (mismeasured) regressor; for the \code{lm} method, the
#'   NAME of the regressor in the model frame
#' @param unit unit identifiers (any type); for the \code{lm} method, the name
#'   of the unit variable in the model frame
#' @param time period identifiers (any type); for the \code{lm} method, the
#'   name of the time variable in the model frame
#' @param sigma_nu per-observation (or scalar) measurement-error SD
#' @param codelow lower posterior-interval bound (with \code{codehigh})
#' @param codehigh upper posterior-interval bound (with \code{codelow})
#' @param reliability the within reliability \eqn{\hat\lambda} in (0, 1]
#' @param reliability_lower optional lower confidence bound for within
#'   reliability. If omitted, the computed reliability is treated as
#'   known/consistent, so certification is conditional on that treatment.
#' @param alpha nominal test level
#' @param delta size-distortion tolerance
#' @param gamma coefficient-uncertainty budget \eqn{\gamma_\beta}; the name is
#'   retained for backward compatibility
#' @param gamma_lambda coverage-error budget for \code{reliability_lower}
#' @param pilot \code{"conservative"} (default), \code{"point"}, or
#'   \code{"naive"}
#' @param cluster standardization of the t-test (cluster-robust extension):
#'   \code{"iid"} (default), \code{"crve"} (by-unit Arellano CR1
#'   variance-inflation), or \code{"ar1"} (parametric within-unit AR(1));
#'   \eqn{|\eta|} and the breakdown are deflated by \eqn{\sqrt{\hat\psi}}
#' @param psi supply your own variance-inflation factor (overrides
#'   \code{cluster})
#' @param ... passed between methods
#' @return an object of class \code{AdequacyReport}
#' @references Halkiewicz, S. M. S. Breakdown Reliability for Saturated
#'   Fixed-Effect Inference.
#' @examples
#' n <- 200; unit <- rep(1:20, each = 10); time <- rep(1:10, times = 20)
#' x <- rnorm(n) + 0.3 * unit; y <- 0.5 * x + rnorm(n)
#' eiv_adequacy(y, x, unit, time, reliability = 0.9)
#' @export
eiv_adequacy <- function(object, ...) UseMethod("eiv_adequacy")

#' @rdname eiv_adequacy
#' @export
eiv_adequacy.default <- function(object, x, unit, time, sigma_nu = NULL,
                                 codelow = NULL, codehigh = NULL,
                                 reliability = NULL, alpha = 0.05,
                                 delta = 0.05, gamma = 0.05,
                                 reliability_lower = NULL,
                                 gamma_lambda = 0,
                                 pilot = c("conservative", "point", "naive"),
                                 cluster = c("iid", "crve", "ar1"),
                                 psi = NULL, ...) {
  y <- object
  pilot <- match.arg(pilot)
  cluster <- match.arg(cluster)
  cc <- .integer_codes(unit, time)
  uid <- cc$uid; tid <- cc$tid; N <- cc$N; T <- cc$T
  n <- length(uid)
  if (length(y) != n || length(x) != n)
    stop("y, x, unit, time must have equal length")
  fd <- fe_dimension(uid, tid, N, T)
  d_K <- fd$d_K
  dof <- n - d_K - 1
  if (dof <= 0) stop("no residual degrees of freedom")

  xt <- .twoway_demean_codes(x, uid, tid, N, T)
  yt <- .twoway_demean_codes(y, uid, tid, N, T)
  tau_star2 <- sum(xt^2)
  if (tau_star2 <= 1e-12 * max(sum(as.numeric(x)^2), 1))
    stop("regressor has no within variation")
  beta_star <- sum(xt * yt) / tau_star2
  u <- yt - beta_star * xt
  sigma <- sqrt(sum(u^2) / dof)

  # noise pilot: exactly one source
  if (!is.null(codelow) || !is.null(codehigh)) {
    if (is.null(codelow) || is.null(codehigh))
      stop("supply both codelow and codehigh")
    if (!is.null(sigma_nu)) stop("multiple noise inputs")
    sigma_nu <- reliability_from_interval(codelow, codehigh)
  }
  if ((!is.null(sigma_nu)) + (!is.null(reliability)) != 1L)
    stop("supply exactly one noise input: sigma_nu, (codelow, codehigh), or reliability")
  if (!is.null(reliability)) {
    if (!(reliability > 0 && reliability <= 1))
      stop("reliability must be in (0, 1]")
    lambda <- reliability
  } else {
    if (length(sigma_nu) != 1L && length(sigma_nu) != n)
      stop("sigma_nu must be scalar or have one value per observation")
    if (any(!is.finite(sigma_nu)) || any(sigma_nu < 0))
      stop("sigma_nu must be finite and non-negative")
    s2 <- mean(sigma_nu^2)
    a_hat <- s2 * (n - d_K)
    lambda <- 1 - a_hat / tau_star2
  }

  extra <- character(0)
  if (length(unique(x)) <= 2)
    extra <- c(extra, "binary treatment detected: errors in binary treatments are MISCLASSIFICATION (nonclassical); this classical-EIV threshold does not apply (scope section of the measurement-error article)")

  # cluster variance-inflation factor psi_hat (cluster-robust extension)
  rho_ar1 <- NULL
  if (!is.null(psi)) {
    psi_hat <- psi
  } else if (cluster == "crve") {
    psi_hat <- .psi_driven(xt, u, uid, tau_star2, sigma^2, N)
  } else if (cluster == "ar1") {
    rho_ar1 <- .rho_ar1(u, uid, tid)
    psi_hat <- .psi_parametric(xt, uid, tid, rho_ar1, kind = "ar1")
  } else {
    psi_hat <- 1
  }

  # Assumption ass-cluster and the checkable form of projection compatibility
  # (lem-nest(b)). Computed whenever a cluster layer is actually in force.
  cdiag <- NULL
  if (psi_hat != 1) {
    cdiag <- cluster_diagnostics(xt, uid, list(unit = uid, time = tid),
                                 tau_star2 = tau_star2)
    pcomp <- projection_compatibility(xt, uid, unit, time,
                                      tau_star2 = tau_star2)
    cdiag$projection_ratio <- pcomp$ratio
    cdiag$projection_cells <- pcomp$cells
    if (cdiag$ratio_ne > 0.20)
      extra <- c(extra, sprintf(
        "PROJECTION-COMPATIBILITY SHORTCUT STRAINED: %d fixed effects cut across only %d clusters (d_ne/G = %.3f). This shortcut needs the spectral and energy balance conditions (N1)--(N2); it is not the assumption itself. Nested dimensions here: %s.",
        cdiag$d_ne, cdiag$G, cdiag$ratio_ne,
        if (length(cdiag$nested)) paste(cdiag$nested, collapse = ", ") else "none"))
    if (is.finite(cdiag$projection_ratio))
      extra <- c(extra, sprintf(
        "direct projection-compatibility diagnostic chi_proj = %.5f = sum_g ||M a^(g)-a^(g)||^2/tau*2 (measurement-error protocol). This finite-panel number evaluates the named sample quantity but does not itself prove the asymptotic sequence condition.",
        cdiag$projection_ratio))
    else
      extra <- c(extra, sprintf(
        "direct projection-compatibility diagnostic skipped because n*G = %.0f exceeds its allocation guard; call projection_compatibility(..., max_cells=...) deliberately to compute it.",
        cdiag$projection_cells))
    if (cdiag$max_energy > 0.10)
      extra <- c(extra, sprintf(
        "one cluster carries %.0f%% of the residualized signal (max_g A_g/tau*2): the no-dominant-cluster condition ass-cluster(ii) is strained and the cluster CLT may not apply",
        100 * cdiag$max_energy))
    if (cdiag$G < 30)
      extra <- c(extra, sprintf(
        "only G = %d clusters: ass-cluster(ii) is a many-clusters condition and the cluster-robust reading is unreliable at this G",
        cdiag$G))
  }

  design <- .design_summary_codes(uid, tid, N, T, xt = xt)
  .eiv_core(design, beta_star, sigma, tau_star2, lambda, alpha = alpha,
            delta = delta, gamma = gamma, pilot = pilot, extra_notes = extra,
            reliability_lower = reliability_lower,
            gamma_lambda = gamma_lambda,
            psi_hat = psi_hat, rho_ar1 = rho_ar1, cluster_diag = cdiag)
}

#' Module B diagnostic from regression summary output
#'
#' For users who have only fitted-model output (no raw data): supply the
#' attenuated slope, residual SD, observed within variation, n and d_K, plus
#' exactly one of `reliability` or `sigma_nu2` (the mean squared
#' measurement-error SD).
#'
#' @param beta_star,sigma,tau_star2 regression output: attenuated slope,
#'   residual SD, observed within variation of the regressor
#' @param n,d_K sample size and fixed-effect dimension
#' @param reliability,sigma_nu2 noise input (exactly one): within reliability,
#'   or the mean squared measurement-error SD
#' @param N,T optional design shape (0 = unknown)
#' @param alpha,delta,gamma,gamma_lambda,reliability_lower,pilot as in
#'   [eiv_adequacy()]
#' @param psi cluster variance-inflation factor (Remark rem-cluster in the
#'   accompanying measurement-error article);
#'   the \code{cluster} presets are unavailable without raw data
#' @return an object of class \code{AdequacyReport}
#' @examples
#' # Repeated-report twins application: both independent-report estimates flag.
#' eiv_adequacy_summary(0.0908759, 1.0, (0.0219815^-2), n = 147, d_K = 4,
#'                      reliability = 0.574904, pilot = "point")
#' @export
eiv_adequacy_summary <- function(beta_star, sigma, tau_star2, n, d_K,
                                 reliability = NULL, sigma_nu2 = NULL,
                                 N = 0L, T = 0L, alpha = 0.05, delta = 0.05,
                                 gamma = 0.05,
                                 reliability_lower = NULL,
                                 gamma_lambda = 0,
                                 pilot = c("conservative", "point", "naive"),
                                 psi = NULL) {
  pilot <- match.arg(pilot)
  if ((!is.null(reliability)) + (!is.null(sigma_nu2)) != 1L)
    stop("supply exactly one of reliability or sigma_nu2")
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) || n <= 0)
    stop("n must be positive")
  if (!is.numeric(d_K) || length(d_K) != 1L || !is.finite(d_K) ||
      d_K < 0 || d_K >= n)
    stop("d_K must satisfy 0 <= d_K < n")
  if (!is.finite(sigma) || sigma <= 0) stop("sigma must be positive")
  if (!is.finite(tau_star2) || tau_star2 <= 0)
    stop("tau_star2 must be positive")
  lambda <- if (!is.null(reliability)) {
    if (!(reliability > 0 && reliability <= 1))
      stop("reliability must be in (0, 1]")
    reliability
  } else {
    if (!is.finite(sigma_nu2) || sigma_nu2 < 0)
      stop("sigma_nu2 must be finite and non-negative")
    1 - sigma_nu2 * (n - d_K) / tau_star2
  }
  design <- structure(list(n = n, N = N, T = T, d_K = d_K, rho = d_K / n,
                           ncomponents = 1L, tau_star2 = tau_star2,
                           lambda_n = NULL, n_eff = NULL),
                      class = "DesignSummary")
  .eiv_core(design, beta_star, sigma, tau_star2, lambda, alpha = alpha,
            delta = delta, gamma = gamma, pilot = pilot,
            reliability_lower = reliability_lower,
            gamma_lambda = gamma_lambda,
            extra_notes = character(0),
            psi_hat = if (is.null(psi)) 1 else psi)
}

.eiv_core <- function(design, beta_star, sigma, tau_star2, lambda, alpha,
                      delta, gamma, pilot, extra_notes,
                      reliability_lower = NULL, gamma_lambda = 0,
                      psi_hat = 1, rho_ar1 = NULL, cluster_diag = NULL) {
  if (!(alpha > 0 && alpha < 1)) stop("alpha must lie in (0, 1)")
  if (!(delta > 0 && delta < 1 - alpha))
    stop("delta must lie in (0, 1-alpha)")
  if (!is.finite(sigma) || sigma <= 0) stop("sigma must be positive")
  if (!is.finite(tau_star2) || tau_star2 <= 0)
    stop("tau_star2 must be positive")
  if (!is.finite(lambda) || lambda > 1)
    stop("reliability cannot exceed one and must be finite")
  if (!is.finite(psi_hat) || psi_hat <= 0) stop("psi must be positive")
  # gamma is gamma_beta and must leave z_{1-gamma} >= 0. At
  # gamma > 1/2 the "upper" bound U_n = |beta0_corr| + z_{1-gamma} se would
  # DEFLATE the pilot, making the conservative verdict weaker than the point
  # verdict while still being labelled CERTIFIED.
  if (!(is.numeric(gamma) && length(gamma) == 1L && is.finite(gamma) &&
        gamma > 0 && gamma <= 0.5))
    stop("gamma must be a single number in (0, 0.5]: at gamma > 0.5 the certificate's upper confidence bound deflates rather than inflates the pilot and the verdict is no longer conservative")
  if (!(is.numeric(gamma_lambda) && length(gamma_lambda) == 1L &&
        is.finite(gamma_lambda) && gamma_lambda >= 0 && gamma_lambda < 1))
    stop("gamma_lambda must be a single number in [0, 1)")
  if (gamma + gamma_lambda >= 1)
    stop("gamma + gamma_lambda must be less than one")
  eta_dag <- .eta_dagger(alpha, delta)
  sqpsi <- sqrt(psi_hat)
  notes <- extra_notes

  # fixed-point breakdown lambda-dagger = t*/(t* + eta-dagger),
  # t* = |beta*| sqrt(tau*2)/sigma (Definition def-breakdown); under
  # clustering t* is deflated by sqrt(psi)
  t_star <- abs(beta_star) * sqrt(tau_star2) / sigma
  t_reported <- t_star / sqpsi
  breakdown <- t_reported / (t_reported + eta_dag)
  z_beta <- stats::qnorm(1 - gamma)
  certified_breakdown <- (t_reported + z_beta) /
    (t_reported + z_beta + eta_dag)

  if (lambda <= 0) {
    notes <- c(notes, sprintf(
      "implied noise exceeds ALL residual within variation (lambda_hat = %.3f <= 0): the noise pilot may be misscaled, or attenuation is total; exact size = 1",
      lambda))
    statistic <- list(lambda_hat = lambda, noise_ratio = Inf,
                      beta_star = beta_star, sigma = sigma,
                      t_star = t_star, psi_hat = psi_hat)
    return(.new_AdequacyReport("measurement_error", design, statistic, Inf,
                               eta_dag, breakdown, 1, "FLAGGED", alpha, delta,
                               notes))
  }

  ell <- if (is.null(reliability_lower)) lambda else reliability_lower
  if (!(is.numeric(ell) && length(ell) == 1L && is.finite(ell) &&
        ell > 0 && ell <= 1))
    stop("reliability_lower must be a single number in (0, 1]")

  # Cluster-robust scale (cor-cluster-feasible): s_CR^2 = sigma_CJN^2 * psi =
  # V^sc_CR / tau*2, an algebraic identity. EVERY scale below is s_CR --
  # including the pilot standard error. Using the i.i.d. sigma there (as the
  # pre-correction code did) understates the certificate by sqrt(psi) and is
  # anti-conservative whenever psi > 1.
  s_CR <- sigma * sqpsi

  beta_corr <- beta_star / lambda
  se_corr <- s_CR / (lambda * sqrt(tau_star2))
  b_pilot <- if (pilot == "naive") abs(beta_star) else abs(beta_corr)
  eta_point <- (b_pilot / s_CR) * (1 - lambda) * sqrt(tau_star2)
  eta_upper <- NULL
  if (pilot == "conservative")
    eta_upper <- (t_reported + z_beta) * (1 - ell) / ell
  eta_used <- if (pilot == "conservative") eta_upper else eta_point

  # rem-plugin: only the certificate is size-controlled. A passing point pilot
  # is a POINT PASS, not a certificate.
  passes <- if (pilot == "conservative") ell >= certified_breakdown else
    eta_used <= eta_dag
  verdict <- if (passes) {
    if (pilot == "conservative") "CERTIFIED" else "POINT_PASS"
  } else "FLAGGED"
  implied_size <- .noncentral_size(eta_point, alpha)

  if (pilot == "conservative") {
    notes <- c(notes, sprintf(
      "formal certificate (prop-certificate): certified breakdown lambda_dagger_gamma = %.3f is obtained by replacing |t| with |t| + z_(1-gamma_beta); the comparison uses reliability lower bound ell = %.3f. False certification is at most gamma_beta + gamma_lambda = %.3g + %.3g = %.3g, without requiring independence. Implied size shown is at the point pilot.",
      certified_breakdown, ell, gamma, gamma_lambda, gamma + gamma_lambda))
    if (is.null(reliability_lower))
      notes <- c(notes, "CONDITIONAL RELIABILITY TREATMENT: no reliability_lower was supplied, so lambda_hat is treated as known/consistent and gamma_lambda = 0. A noisy finite-sample reliability estimate requires a lower confidence bound and its coverage-error budget.")
  } else if (pilot == "point") {
    notes <- c(notes, "POINT PASS, not a certificate (rem-plugin): the corrected pilot at its point estimate is descriptive. Under weak information eta_hat converges to a nondegenerate random multiple (|B|/|beta0|)|eta| of the target -- median close to it, but no concentration. Its sampling variability is a first-order feature of the regime, not a vanishing approximation error. For a size-controlled statement use pilot = \"conservative\".")
  } else {
    notes <- c(notes, "ANTI-CONSERVATIVE naive pilot (attenuated beta*) \u2014 for comparison only; understates |eta| by the factor lambda (Proposition prop-pilot(i))")
  }
  if (pilot != "naive")
    notes <- c(notes, "corrected pilot beta*/lambda_hat is not a consistent point estimate under weak information; reported with its sampling band (Corollary cor-slope)")
  if (eta_point > 1)
    notes <- c(notes, "far from the threshold (|eta| > 1): the local quadratic approximation is uninformative here; the verdict uses exact inversion (Remark rem-exact-cv)")
  if (lambda < 1)
    notes <- c(notes, sprintf(
      "power tax: local power slope attenuated by sqrt(lambda) = %.2f (Proposition prop-power)",
      sqrt(lambda)))
  if (psi_hat != 1) {
    dir <- if (psi_hat > 1)
      "on THIS design clustering deflates the measured distortion"
    else
      "on THIS design clustering WORSENS the distortion (psi < 1) -- the i.i.d. reading is anti-conservative here"
    notes <- c(notes, sprintf(
      "cluster-robust standardization: psi_hat = %.3f, so the scale is s_CR = sigma*sqrt(psi) = %.4g and |eta|, t* are deflated by 1/sqrt(psi) = %.3f; %s. The sign of psi-1 is not free and not guessable (rem-psi-sign): an equicorrelated component at a level that is itself a fixed effect is annihilated exactly, and a serially dependent error with a within-cluster serially independent regressor gives psi < 1. What pushes psi above one is persistence in the regressor and the error together. Read the direction off psi_hat.",
      psi_hat, s_CR, 1 / sqpsi, dir))
    notes <- c(notes, sprintf(
      "lambda_dagger_CR = %.3f is def-breakdown evaluated at the reported cluster-robust t-statistic t^CR = t*/sqrt(psi) = %.3f -- an algebraic identity (cor-cluster-feasible(c)) requiring no limit theory: the cluster-robust diagnostic IS the i.i.d. diagnostic run on the reported cluster-robust t.",
      breakdown, t_star / sqpsi))
    notes <- c(notes, "the CRVE here omits the conventional (n-1)/(n-K) small-sample factor, which converges to 1/(1-rho) rather than 1 and so OVER-corrects when rho is non-negligible (lem-crve). Default software applies it: reproducing psi_hat with such defaults will inflate it by roughly 1/(1-rho).")
  }

  statistic <- list(lambda_hat = lambda, reliability_lower = ell,
                    noise_ratio = (1 - lambda) / lambda,
                    beta_star = beta_star, beta_corr = beta_corr,
                    se_beta_corr = se_corr, sigma = sigma, s_CR = s_CR,
                    t_star = t_star, t_CR = t_star / sqpsi,
                    psi_hat = psi_hat,
                    breakdown_certified = certified_breakdown,
                    gamma = gamma, gamma_beta = gamma,
                    gamma_lambda = gamma_lambda,
                    false_certification_bound = gamma + gamma_lambda,
                    pilot = pilot,
                    eta_quad_threshold = .eta_quad(alpha, delta))
  if (!is.null(rho_ar1)) statistic$rho_ar1 <- rho_ar1
  if (!is.null(cluster_diag)) statistic$cluster <- cluster_diag
  if (!is.null(eta_upper)) statistic$eta_upper <- eta_upper

  .new_AdequacyReport("measurement_error", design, statistic, eta_point,
                      eta_dag, breakdown, implied_size, verdict, alpha, delta,
                      notes)
}


#' Checkable conditions behind the cluster-robust layer
#'
#' Assumption ass-cluster and Lemma lem-nest of the accompanying article. All are computable from
#' the design alone, before any outcome is examined.
#'
#' Condition (iv) of ass-cluster -- projection compatibility -- is the one with
#' no counterpart in the i.i.d. theory, and it fails silently when the fixed
#' effects cut across clusters in a high-dimensional way.
#'
#' Lemma lem-nest(b) bounds it by `varpi_n * max_g A_g`, which is
#' unconditional. The familiar reduction to `d_ne / G -> 0` is *not*: it holds
#' only under two balance conditions, which this function reports rather than
#' assumes.
#' \itemize{
#'   \item (N1) spectral: `lambda_min^+(Lambda_n) ~ n / d_ne`. Exact in a
#'     balanced two-way panel, where `Lambda_n = G (I_T - 11'/T)`, but equal
#'     cell counts alone do not control the smallest non-zero eigenvalue of the
#'     residualized Gram matrix. Not computed here.
#'   \item (N2) energy: `max_g A_g = O(tau*2 / G)`, reported as
#'     \code{max_energy}. If the residualized signal concentrates in
#'     `sqrt(G)` clusters then `max_energy ~ G^-1/2` and the requirement
#'     becomes `d_ne / sqrt(G) -> 0` instead.
#' }
#' So read \code{ratio_ne} and \code{max_energy} together: a small
#' \code{ratio_ne} carries no warrant on its own. When the two disagree, or
#' when (N1) is in doubt, use [projection_compatibility()], which evaluates
#' `sum_g ||M a^(g) - a^(g)||^2 / tau*2` directly and needs neither condition.
#'
#' In the V-Dem application country effects nest in country clusters while the
#' 59 year effects do not, against `G = 163`, so `ratio_ne ~ 0.36` and the
#' direct diagnostic is load-bearing. The repeated-report twins application
#' samples independent pairs and therefore uses the baseline i.i.d. theorem;
#' it does not invoke this cluster condition or a `psi` rescaling.
#'
#' @param xt residualized regressor.
#' @param cluster cluster identifier, one per observation.
#' @param fe_levels named list of fixed-effect id vectors, one per FE dimension.
#' @param tau_star2 the within variation `x*'M x*`.
#' @return A list with `G`, `max_size`, `min_size`, `max_energy`, `d_ne`,
#'   `ratio_ne` and `nested`.
#' @export
cluster_diagnostics <- function(xt, cluster, fe_levels, tau_star2) {
  n <- length(xt)
  if (length(cluster) != n) stop("cluster must have the same length as the data")
  cid <- as.integer(factor(cluster, levels = unique(cluster)))
  G <- max(cid)
  sizes <- tabulate(cid, nbins = G)
  energy <- vapply(seq_len(G), function(g) sum(xt[cid == g]^2), numeric(1))

  nested <- character(0)
  d_ne <- 0L
  for (nm in names(fe_levels)) {
    ids <- fe_levels[[nm]]
    if (length(ids) != n)
      stop(sprintf("fixed-effect ids for %s must have length n", nm))
    f <- as.integer(factor(ids, levels = unique(ids)))
    # nested iff every level of this FE dimension sits in exactly one cluster
    per <- tapply(cid, f, function(v) length(unique(v)))
    if (all(per == 1L)) {
      nested <- c(nested, nm)
    } else {
      d_ne <- d_ne + length(unique(f))
    }
  }

  list(G = G, max_size = max(sizes), min_size = min(sizes),
       max_energy = max(energy) / tau_star2,
       d_ne = d_ne, ratio_ne = d_ne / G, nested = nested)
}

#' Projection compatibility, evaluated directly
#'
#' Assumption ass-cluster(iv) of the accompanying article asks that
#' \eqn{\sum_g \|M a^{(g)} - a^{(g)}\|^2 = o_p(\tau^{*2})}, where
#' \eqn{a^{(g)}} is the residualized regressor restricted to cluster \eqn{g}.
#' This function computes that ratio as it stands, for a two-way (unit and
#' time) fixed-effect design.
#'
#' Unlike the \code{ratio_ne} shortcut in [cluster_diagnostics()], this needs
#' neither of the balance conditions (N1) and (N2) of Lemma lem-nest(b): it is
#' the quantity the assumption actually names. Use it when \code{ratio_ne} and
#' \code{max_energy} point in different directions, when the panel is far from
#' balanced, or whenever the cluster-robust column is load-bearing.
#'
#' Under Lemma lem-nest(a) the ratio is exactly zero when every fixed-effect
#' cell sits inside one cluster, so a non-zero value is entirely the work of
#' the fixed effects that cut across clusters.
#'
#' Cost is one alternating-projection sweep over an \code{n * G} matrix. The
#' \code{max_cells} guard returns \code{NA} rather than allocating a matrix
#' larger than that; raise it deliberately if you want the number anyway.
#'
#' @param xt residualized regressor (the same \code{M x*} used elsewhere).
#' @param cluster cluster identifier, one per observation.
#' @param unit,time raw fixed-effect identifiers, one per observation.
#' @param tau_star2 the within variation \code{sum(xt^2)}; defaults to
#'   \code{sum(xt^2)}.
#' @param tol,maxit convergence controls for the alternating projections.
#' @param max_cells refuse to allocate more than this many matrix cells.
#' @return A list with `ratio` (the quantity above; `NA` if the guard tripped),
#'   `G`, `n` and `cells`.
#' @examples
#' n <- 200; unit <- rep(1:20, each = 10); time <- rep(1:10, times = 20)
#' x <- rnorm(n)
#' xt <- twoway_demean(x, unit, time)
#' # unit FE nest in unit clusters, the 10 time effects do not
#' projection_compatibility(xt, unit, unit, time)$ratio
#' @export
projection_compatibility <- function(xt, cluster, unit, time,
                                     tau_star2 = sum(xt^2),
                                     tol = 1e-10, maxit = 10000L,
                                     max_cells = 5e6) {
  xt <- as.numeric(xt)
  n <- length(xt)
  if (length(cluster) != n || length(unit) != n || length(time) != n)
    stop("xt, cluster, unit and time must have equal length")
  if (!(tau_star2 > 0)) stop("tau_star2 must be positive")

  cid <- as.integer(factor(cluster, levels = unique(cluster)))
  G <- max(cid)
  cells <- as.numeric(n) * G
  if (cells > max_cells)
    return(list(ratio = NA_real_, G = G, n = n, cells = cells))

  cc <- .integer_codes(unit, time)
  uid <- cc$uid; tid <- cc$tid; N <- cc$N; T <- cc$T
  ucnt <- pmax(tabulate(uid, N), 1L)
  tcnt <- pmax(tabulate(tid, T), 1L)

  # columns of A are the cluster loading vectors a^(g)
  A <- matrix(0, n, G)
  A[cbind(seq_len(n), cid)] <- xt

  W <- A
  converged <- FALSE
  for (it in seq_len(maxit)) {
    um <- rowsum(W, uid) / ucnt
    W <- W - um[uid, , drop = FALSE]
    tm <- rowsum(W, tid) / tcnt
    delta <- max(abs(tm))
    W <- W - tm[tid, , drop = FALSE]
    if (delta < tol) { converged <- TRUE; break }
  }
  if (!converged) warning("two-way demeaning did not converge within maxit")

  # W is M A, so A - W is P A and sum((A - W)^2) = sum_g ||M a^(g) - a^(g)||^2
  list(ratio = sum((A - W)^2) / tau_star2, G = G, n = n, cells = cells)
}

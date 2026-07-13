# Module B \u2014 measurement-error adequacy (spec section 4; Paper B).
# Formula sources: feasible non-centrality (Cor. cor-feasible); corrected pilot
# beta*/lambda (Prop. prop-pilot); pilot se sigma/(lambda sqrt(tau*2))
# (Cor. cor-slope); exact-inversion threshold eta\u2020 with the quadratic closed
# form as companion \u2014 NEVER the discarded linear surrogate (Rem. rem-exact-cv);
# breakdown reliability (Def. def-breakdown); exact non-central implied size.

#' Measurement-error SDs from published credible-interval bounds
#'
#' V-Dem convention: the interval brackets one posterior SD, so
#' `sigma_nu = (codehigh - codelow)/2` (Paper B, lead application).
#'
#' @param codelow,codehigh interval bounds, one per observation
#' @return numeric vector of per-observation measurement-error SDs
#' @examples
#' reliability_from_interval(c(0.20, 0.31), c(0.24, 0.35))
#' @export
reliability_from_interval <- function(codelow, codehigh) {
  if (length(codelow) != length(codehigh))
    stop("codelow and codehigh must have equal length")
  (codehigh - codelow) / 2
}

#' Measurement-error SD implied by an external reliability ratio
#'
#' PSID-style validation studies: `sigma_nu = sqrt((1-r)/r) * within_sd`.
#'
#' @param r reliability ratio in (0, 1]
#' @param within_sd within-SD of the observed regressor
#' @return scalar measurement-error SD
#' @examples
#' reliability_from_ratio(0.65, within_sd = 0.15)  # PSID within reliability
#' @export
reliability_from_ratio <- function(r, within_sd) {
  if (!(r > 0 && r <= 1)) stop("reliability ratio must be in (0, 1]")
  sqrt((1 - r) / r) * within_sd
}

#' Self-consistent breakdown reliability
#'
#' The fixed point \eqn{\lambda = \lambda^\dagger(\lambda)} when the corrected
#' pilot \eqn{\beta^*/\lambda} is evaluated at the reliability being solved
#' for: \eqn{\lambda^\dagger = 1/(1 + \eta^\dagger \sigma/(|\beta^*|
#' \sqrt{\tau^{*2}}))} (Paper B, PSID application). Needs no external noise
#' input; uses the exact-inversion root \eqn{\eta^\dagger}.
#'
#' @param beta_star attenuated slope from the FE regression
#' @param sigma residual standard deviation
#' @param tau_star2 observed within variation of the regressor
#' @param alpha nominal test level
#' @param delta size-distortion tolerance
#' @return the breakdown reliability in [0, 1]
#' @examples
#' breakdown_reliability(0.733, 4.25, 83.2)  # PSID application: ~0.71
#' @export
breakdown_reliability <- function(beta_star, sigma, tau_star2,
                                  alpha = 0.05, delta = 0.05) {
  if (beta_star == 0) return(0)
  eta_dag <- .eta_dagger(alpha, delta)
  1 / (1 + eta_dag * sigma / (abs(beta_star) * sqrt(tau_star2)))
}

#' Module B diagnostic: measurement-error adequacy (Paper B, flagship)
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
#' Correct-by-default honesty machinery (Paper B, protocol steps 4--5):
#' \code{pilot = "conservative"} (default) evaluates the verdict at the upper
#' \code{1 - gamma} confidence bound of the attenuation-corrected pilot
#' \eqn{\hat\beta^*/\hat\lambda} --- the formal certificate; \code{"point"} is
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
#' @param alpha nominal test level
#' @param delta size-distortion tolerance
#' @param gamma level for the conservative certificate's upper bound
#' @param pilot \code{"conservative"} (default), \code{"point"}, or
#'   \code{"naive"}
#' @param ... passed between methods
#' @return an object of class \code{AdequacyReport}
#' @references Halkiewicz, S. M. S. Stock--Yogo critical values for
#'   fixed-effect saturation under measurement error (Paper B).
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
                                 pilot = c("conservative", "point", "naive"),
                                 ...) {
  y <- object
  pilot <- match.arg(pilot)
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
    s2 <- mean(sigma_nu^2)
    a_hat <- s2 * (n - d_K)
    lambda <- 1 - a_hat / tau_star2
  }

  extra <- character(0)
  if (length(unique(x)) <= 2)
    extra <- c(extra, "binary treatment detected: errors in binary treatments are MISCLASSIFICATION (nonclassical); this classical-EIV threshold does not apply (Paper B \u00a75.3)")

  design <- structure(list(n = n, N = N, T = T, d_K = d_K, rho = d_K / n,
                           ncomponents = fd$ncomponents, tau_star2 = tau_star2),
                      class = "DesignSummary")
  .eiv_core(design, beta_star, sigma, tau_star2, lambda, alpha = alpha,
            delta = delta, gamma = gamma, pilot = pilot, extra_notes = extra)
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
#' @param alpha,delta,gamma,pilot as in [eiv_adequacy()]
#' @return an object of class \code{AdequacyReport}
#' @examples
#' # PSID application (Paper B): flagged at the within reliability 0.65
#' eiv_adequacy_summary(0.733, 4.25, 83.2, n = 4165, d_K = 601,
#'                      reliability = 0.65, pilot = "point")
#' @export
eiv_adequacy_summary <- function(beta_star, sigma, tau_star2, n, d_K,
                                 reliability = NULL, sigma_nu2 = NULL,
                                 N = 0L, T = 0L, alpha = 0.05, delta = 0.05,
                                 gamma = 0.05,
                                 pilot = c("conservative", "point", "naive")) {
  pilot <- match.arg(pilot)
  if ((!is.null(reliability)) + (!is.null(sigma_nu2)) != 1L)
    stop("supply exactly one of reliability or sigma_nu2")
  lambda <- if (!is.null(reliability)) {
    if (!(reliability > 0 && reliability <= 1))
      stop("reliability must be in (0, 1]")
    reliability
  } else {
    1 - sigma_nu2 * (n - d_K) / tau_star2
  }
  design <- structure(list(n = n, N = N, T = T, d_K = d_K, rho = d_K / n,
                           ncomponents = 1L, tau_star2 = tau_star2),
                      class = "DesignSummary")
  .eiv_core(design, beta_star, sigma, tau_star2, lambda, alpha = alpha,
            delta = delta, gamma = gamma, pilot = pilot,
            extra_notes = character(0))
}

.eiv_core <- function(design, beta_star, sigma, tau_star2, lambda, alpha,
                      delta, gamma, pilot, extra_notes) {
  eta_dag <- .eta_dagger(alpha, delta)
  notes <- extra_notes

  if (lambda <= 0) {
    notes <- c(notes, sprintf(
      "implied noise exceeds ALL residual within variation (lambda_hat = %.3f <= 0): the noise pilot may be misscaled, or attenuation is total; exact size = 1",
      lambda))
    statistic <- list(lambda_hat = lambda, noise_ratio = Inf,
                      beta_star = beta_star, sigma = sigma)
    return(.new_AdequacyReport("measurement_error", design, statistic, Inf,
                               eta_dag, 1, 1, "FLAGGED", alpha, delta, notes))
  }

  beta_corr <- beta_star / lambda
  se_corr <- sigma / (lambda * sqrt(tau_star2))
  b_pilot <- if (pilot == "naive") abs(beta_star) else abs(beta_corr)
  eta_point <- (b_pilot / sigma) * (1 - lambda) * sqrt(tau_star2)
  eta_upper <- NULL
  if (pilot == "conservative")
    eta_upper <- ((abs(beta_corr) + stats::qnorm(1 - gamma) * se_corr) / sigma) *
      (1 - lambda) * sqrt(tau_star2)
  eta_used <- if (pilot == "conservative") eta_upper else eta_point

  verdict <- if (eta_used <= eta_dag) "CERTIFIED" else "FLAGGED"
  implied_size <- .noncentral_size(eta_point, alpha)
  breakdown <- if (b_pilot == 0) 0 else
    min(max(1 - eta_dag * sigma / (b_pilot * sqrt(tau_star2)), 0), 1)

  if (pilot == "conservative") {
    notes <- c(notes, sprintf(
      "formal certificate: verdict uses the upper %.0f%% confidence bound of the corrected pilot, |eta|_ub = %.3f (Paper B, protocol step 5); implied size shown is at the point pilot",
      100 * (1 - gamma), eta_upper))
  } else if (pilot == "point") {
    notes <- c(notes, "point diagnostic (descriptive): corrected pilot at its point estimate \u2014 not a formally size-controlled certificate (Paper B, Remark rem-corr-noise)")
  } else {
    notes <- c(notes, "ANTI-CONSERVATIVE naive pilot (attenuated beta*) \u2014 for comparison only; understates |eta| by the factor lambda (Paper B, Prop. prop-pilot(i))")
  }
  if (pilot != "naive")
    notes <- c(notes, "corrected pilot beta*/lambda_hat is not a consistent point estimate under weak information; reported with its sampling band (Paper B, Cor. cor-slope)")
  if (eta_point > 1)
    notes <- c(notes, "far from the threshold (|eta| > 1): the local quadratic approximation is uninformative here; the verdict uses exact inversion (Paper B, Remark rem-exact-cv)")
  if (lambda < 1)
    notes <- c(notes, sprintf(
      "power tax: local power slope attenuated by sqrt(lambda) = %.2f (Paper B, Prop. prop-power)",
      sqrt(lambda)))

  statistic <- list(lambda_hat = lambda, noise_ratio = (1 - lambda) / lambda,
                    beta_star = beta_star, beta_corr = beta_corr,
                    se_beta_corr = se_corr, sigma = sigma,
                    eta_quad_threshold = .eta_quad(alpha, delta))
  if (!is.null(eta_upper)) statistic$eta_upper <- eta_upper

  .new_AdequacyReport("measurement_error", design, statistic, eta_point,
                      eta_dag, breakdown, implied_size, verdict, alpha, delta,
                      notes)
}

# Normal / non-central size machinery shared by all modules. These constants
# are cross-language parity targets: the Julia reference engine implements the
# identical bisection, so the two packages agree to ~1e-12.
#
# The size is EVEN in eta; its leading distortion is QUADRATIC z*phi(z)*eta^2.
# The discarded LINEAR surrogate delta/(2*phi) must never be used
# (Paper B, Remark rem-exact-cv; Paper C, Remark on the linear heuristic).

# Exact two-sided size of the nominal-alpha test when T ~ N(eta, 1).
.noncentral_size <- function(eta, alpha) {
  z <- stats::qnorm(1 - alpha / 2)
  stats::pnorm(-z - abs(eta)) + 1 - stats::pnorm(z - abs(eta))
}

# Exact-inversion threshold eta\u2020(alpha, delta): unique positive root of
# size(eta) = alpha + delta. ~0.652 at alpha = delta = 0.05.
.eta_dagger <- function(alpha, delta) {
  if (!(delta > 0 && delta < 1 - alpha)) stop("need 0 < delta < 1 - alpha")
  target <- alpha + delta
  lo <- 0; hi <- 1
  while (.noncentral_size(hi, alpha) < target) {
    hi <- hi * 2
    if (hi > 1e6) stop("eta_dagger bracket failed")
  }
  for (i in 1:200) {
    mid <- (lo + hi) / 2
    if (.noncentral_size(mid, alpha) < target) lo <- mid else hi <- mid
  }
  (lo + hi) / 2
}

# Quadratic closed-form companion threshold sqrt(delta / (z phi(z))).
.eta_quad <- function(alpha, delta) {
  z <- stats::qnorm(1 - alpha / 2)
  sqrt(delta / (z * stats::dnorm(z)))
}

# Diffuse-companion asymptotic size maps (homoskedastic / design-balanced limits).
.size_naive <- function(rho, alpha) {
  z <- stats::qnorm(1 - alpha / 2)
  2 * (1 - stats::pnorm(z * sqrt(1 - rho)))
}

.size_hc3 <- function(rho, alpha) {
  z <- stats::qnorm(1 - alpha / 2)
  2 * (1 - stats::pnorm(z / sqrt(1 - rho)))
}

# Breakdown saturation rho\u2020 at which the HC0/naive size reaches alpha + delta.
.rho_dagger <- function(alpha, delta) {
  1 - (stats::qnorm(1 - (alpha + delta) / 2) / stats::qnorm(1 - alpha / 2))^2
}

# Retained name: HC0 shares the naive limit under conditional homoskedasticity.
.size_hc0 <- function(rho, alpha) .size_naive(rho, alpha)

#' Stock-Yogo critical value for the residual treatment variance
#'
#' The default is Paper B's exact-normal inversion:
#' `tau2_crit = beta0^2 c^4 (1-rho)/(sigma^2 eta_dagger^2) - c^2`.
#' Set `method = "quadratic"` for the closed-form companion
#' `beta0^2 c^4 (1-rho) z phi(z)/(sigma^2 delta) - c^2`. The quadratic
#' boundary is mildly anti-conservative, so it is not the default.
#' Two features trace to the corrected Hessian limit and are easy to get wrong:
#' the leading term is LINEAR in `(1-rho)`, not quadratic, and the subtracted
#' term is `c^2`, not `c^2 (1-rho)`.
#'
#' @param rho fixed-effect saturation `d_K/n`, in `[0,1)`.
#' @param c2 the drift constant `c^2` in `sigma_nu^2 = c^2/n`.
#' @param beta0,sigma coefficient and residual scale.
#' @param alpha,delta level and size tolerance.
#' @param method `"exact"` (default) or `"quadratic"`.
#' @return The critical value; negative when the design is adequate at every `tau^2`.
#' @export
tau2_crit <- function(rho, c2, beta0, sigma, alpha = 0.05, delta = 0.05,
                      method = c("exact", "quadratic")) {
  method <- match.arg(method)
  if (rho < 0 || rho >= 1) stop("need 0 <= rho < 1")
  if (c2 < 0) stop("c2 must be non-negative")
  if (sigma <= 0) stop("sigma must be positive")
  if (!(alpha > 0 && alpha < 1)) stop("alpha must lie in (0, 1)")
  if (!(delta > 0 && delta < 1 - alpha))
    stop("delta must lie in (0, 1-alpha)")
  factor <- if (method == "exact") {
    1 / .eta_dagger(alpha, delta)^2
  } else {
    z <- stats::qnorm(1 - alpha / 2)
    z * stats::dnorm(z) / delta
  }
  beta0^2 * c2^2 * (1 - rho) * factor / sigma^2 - c2
}

#' Finite-n non-centrality mapping
#'
#' Paper B eq-eta-n under the local drift `sigma_nu^2 = c^2/n`:
#' `eta_n = -beta0 c^2 sqrt(1-rho) / (sigma sqrt(tau^2 + c^2))`. The
#' `sqrt(1-rho)` factor is the common-scale shrinkage implied by
#' `X'M_K X/(nQ_K) ->p 1-rho`; it was absent from the pre-correction formula.
#'
#' @param beta0,c2,sigma,tau2,rho model primitives.
#' @return The non-centrality `eta_n`.
#' @export
eta_finite_n <- function(beta0, c2, sigma, tau2, rho) {
  if (rho < 0 || rho >= 1) stop("need 0 <= rho < 1")
  if (c2 < 0 || tau2 < 0) stop("c2 and tau2 must be non-negative")
  if (sigma <= 0) stop("sigma must be positive")
  -beta0 * c2 * sqrt(1 - rho) / (sigma * sqrt(tau2 + c2))
}

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

# Paper A asymptotic size maps (homoskedastic / design-balanced limits).
.size_hc0 <- function(rho, alpha) {
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

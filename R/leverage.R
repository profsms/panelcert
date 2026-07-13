# Module A \u2014 leverage / variance diagnostics (spec section 3; Paper A).
# Formula sources: full-regression leverage H_ii = (P_K)_ii + Xt_i^2/(X'MX)
# (eq. 2.3); HC0-HC3 weights (section 3.2); asymptotic sizes (Thm 3.1, Cor 3.1);
# HC3 over-correction regime rho > 0.1 (section 6); HC1-vs-LO divergence under
# non-uniform leverage (Remark 3.4).

#' Diagonal of the two-way fixed-effect projection
#'
#' `(P_K)_ii` for each observation, computed exactly from the (N+T)x(N+T)
#' FE Gram matrix via pseudo-inverse (handles the rank deficiency
#' `d_K = N + T - #components`).
#'
#' @param unit,time raw identifier vectors
#' @return numeric vector of FE leverages, one per observation
#' @export
fe_leverage <- function(unit, time) {
  cc <- .integer_codes(unit, time)
  .fe_leverage_codes(cc$uid, cc$tid, cc$N, cc$T)
}

.fe_leverage_codes <- function(uid, tid, N, T) {
  G <- matrix(0, N + T, N + T)
  for (k in seq_along(uid)) {
    u <- uid[k]; t <- N + tid[k]
    G[u, u] <- G[u, u] + 1
    G[t, t] <- G[t, t] + 1
    G[u, t] <- G[u, t] + 1
    G[t, u] <- G[t, u] + 1
  }
  Gp <- .pinv(G)
  ui <- uid; ti <- N + tid
  Gp[cbind(ui, ui)] + 2 * Gp[cbind(ui, ti)] + Gp[cbind(ti, ti)]
}

#' Module A diagnostic: variance-estimator adequacy under FE saturation
#'
#' Reproduces the user's FE regression of `y` on `x` with unit and time fixed
#' effects via Frisch-Waugh (never re-specified), then reports the naive / CJN /
#' HC0-HC3 variance hierarchy, the leverage diagnostics, and whether the
#' variance-estimator choice materially changes inference at tolerance `delta`.
#'
#' Verdict: FLAGGED when the asymptotic HC0/naive over-rejection exceeds
#' `alpha + delta` at this design's saturation `rho`, or when significance at
#' level `alpha` flips across the estimators on this data.
#'
#' @param object outcome vector (default method), or a fitted \code{fixest} /
#'   \code{plm} / \code{lm} model with one regressor and two-way fixed effects
#' @param x regressor of interest; for the \code{lm} method, its NAME in the
#'   model frame
#' @param unit unit identifiers (any type); for the \code{lm} method, the name
#'   of the unit variable in the model frame
#' @param time period identifiers (any type); for the \code{lm} method, the
#'   name of the time variable in the model frame
#' @param alpha nominal test level
#' @param delta size-distortion tolerance
#' @param ... passed between methods
#' @return an object of class \code{AdequacyReport}
#' @references Halkiewicz, S. M. S. Variance estimation for saturated
#'   fixed-effect specifications (Paper A).
#' @examples
#' n <- 200; unit <- rep(1:20, each = 10); time <- rep(1:10, times = 20)
#' x <- rnorm(n); y <- 0.5 * x + rnorm(n)
#' leverage_report(y, x, unit, time)
#' @export
leverage_report <- function(object, ...) UseMethod("leverage_report")

#' @rdname leverage_report
#' @export
leverage_report.default <- function(object, x, unit, time, alpha = 0.05,
                                    delta = 0.05, ...) {
  y <- object
  cc <- .integer_codes(unit, time)
  uid <- cc$uid; tid <- cc$tid; N <- cc$N; T <- cc$T
  n <- length(uid)
  if (length(y) != n || length(x) != n)
    stop("y, x, unit, time must have equal length")
  fd <- fe_dimension(uid, tid, N, T)
  d_K <- fd$d_K
  dof <- n - d_K - 1
  if (dof <= 0) stop("no residual degrees of freedom (n - d_K - 1 = ", dof, " <= 0)")

  xt <- .twoway_demean_codes(x, uid, tid, N, T)
  yt <- .twoway_demean_codes(y, uid, tid, N, T)
  tau_star2 <- sum(xt^2)
  if (tau_star2 <= 1e-12 * max(sum(as.numeric(x)^2), 1))
    stop("regressor has no within variation (collinear with the fixed effects)")

  beta <- sum(xt * yt) / tau_star2
  u <- yt - beta * xt
  rss <- sum(u^2)

  p_fe <- .fe_leverage_codes(uid, tid, N, T)
  H <- p_fe + xt^2 / tau_star2
  maxH <- max(H)
  if (maxH >= 1 - 1e-10)
    stop("an observation has full leverage H_ii = 1; HC2/HC3 are undefined ",
         "(degenerate cell \u2014 Paper A's bounded-leverage condition fails)")

  hc0 <- sum(xt^2 * u^2)
  hc1 <- n / dof * hc0
  hc2 <- sum(xt^2 * u^2 / (1 - H))
  hc3 <- sum(xt^2 * u^2 / (1 - H)^2)

  se_naive <- sqrt(rss / n / tau_star2)
  se_cjn   <- sqrt(rss / dof / tau_star2)
  se_hc0   <- sqrt(hc0) / tau_star2
  se_hc1   <- sqrt(hc1) / tau_star2
  se_hc2   <- sqrt(hc2) / tau_star2
  se_hc3   <- sqrt(hc3) / tau_star2

  rho <- d_K / n
  z <- stats::qnorm(1 - alpha / 2)
  size_hc0 <- .size_hc0(rho, alpha)
  size_hc3 <- .size_hc3(rho, alpha)
  rho_dag <- .rho_dagger(alpha, delta)

  sig <- abs(beta / c(se_cjn, se_hc0, se_hc1, se_hc2, se_hc3)) > z
  flip <- any(sig != sig[1])

  notes <- c(
    "leave-one-out (HC2) is the recommended default for saturated FE (Paper A, Remark 3.6)",
    "implied sizes use the homoskedastic / design-balanced limits of Theorem 3.3 (sigma_eff^2 = omega^2)"
  )
  if (rho > 0.1)
    notes <- c(notes, sprintf(
      "HC3 over-correcting regime (rho = %.3f > 0.1): HC3 intervals are artificially conservative (implied size %.1f%%) \u2014 re-run with HC2/LO (Paper A \u00a76)",
      rho, 100 * size_hc3))
  spread <- maxH / min(H)
  if (spread > 2)
    notes <- c(notes, sprintf(
      "leverage non-uniform (hmax/hmin = %.2f): HC1 is unreliable here; prefer HC2/LO (Paper A, Remark 3.4)",
      spread))
  if (flip)
    notes <- c(notes, "significance at level alpha flips across variance estimators \u2014 inference is estimator-dependent; use HC2/LO")
  max_x_share <- max(xt^2) / tau_star2
  if (max_x_share > 0.1)
    notes <- c(notes, sprintf(
      "treatment leverage-ratio condition strained: one observation carries %.0f%% of the within variation",
      100 * max_x_share))

  verdict <- if (size_hc0 - alpha > delta || flip) "FLAGGED" else "CERTIFIED"

  design <- structure(list(n = n, N = N, T = T, d_K = d_K, rho = rho,
                           ncomponents = fd$ncomponents, tau_star2 = tau_star2),
                      class = "DesignSummary")
  statistic <- list(beta = beta, se_naive = se_naive, se_cjn = se_cjn,
                    se_hc0 = se_hc0, se_hc1 = se_hc1, se_hc2 = se_hc2,
                    se_hc3 = se_hc3, t_hc2 = beta / se_hc2,
                    max_leverage = maxH, max_fe_leverage = max(p_fe),
                    leverage_spread = spread, max_x_share = max_x_share,
                    implied_size_hc3 = size_hc3)
  .new_AdequacyReport("leverage", design, statistic, NULL, NULL, rho_dag,
                      size_hc0, verdict, alpha, delta, notes)
}

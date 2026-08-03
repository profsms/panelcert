# Leverage / variance diagnostics for the diffuse-regime companion paper.
# Formula sources: full-regression leverage H_ii = (P_K)_ii + Xt_i^2/(X'MX)
# (eq. 2.3); HC0-HC3 weights (section 3.2); asymptotic sizes (Thm 3.1, Cor 3.1);
# HC3 over-correction regime rho > 0.1 (section 6); and the failure of scalar
# degrees-of-freedom corrections under non-uniform leverage.

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

.score_concentration <- function(xt, u) {
  if (length(xt) != length(u)) stop("xt and residuals must have equal length")
  mass <- sum(xt^2 * u^2)
  if (!(mass > 0))
    return(list(lambda_score = NA_real_, H_score = NA_real_,
                n_eff_score = NA_real_))
  shares <- xt^2 * u^2 / mass
  H <- sum(shares^2)
  list(lambda_score = max(shares), H_score = H, n_eff_score = 1 / H)
}

#' Realized score concentration
#'
#' Computes the largest residual-score variance share, its Herfindahl index,
#' and inverse-Herfindahl effective support size for a two-way fixed-effect
#' regression. This is a one-outcome warning diagnostic, not by itself a
#' consistent estimator under unrestricted heteroskedasticity.
#'
#' @param y,x outcome and regressor vectors.
#' @param unit,time fixed-effect identifiers.
#' @param controls optional numeric nuisance-covariate vector or matrix.
#' @return A list with `lambda_score`, `H_score`, and `n_eff_score`.
#' @export
score_concentration <- function(y, x, unit, time, controls = NULL) {
  cc <- .integer_codes(unit, time)
  n <- length(cc$uid)
  if (length(y) != n || length(x) != n)
    stop("y, x, unit, time must have equal length")
  partial <- .partial_within_codes(x, cc$uid, cc$tid, cc$N, cc$T,
                                   controls = controls)
  xt <- partial$xt
  yt <- .partial_outcome_codes(y, cc$uid, cc$tid, cc$N, cc$T, partial$Q)
  Vn <- sum(xt^2)
  if (Vn <= 1e-12 * max(sum(as.numeric(x)^2), 1))
    stop("regressor has no within variation (collinear with the fixed effects)")
  beta <- sum(xt * yt) / Vn
  .score_concentration(xt, yt - beta * xt)
}

#' Module A diagnostic: variance-estimator adequacy under FE saturation
#'
#' Reproduces the user's FE regression of `y` on `x` with unit and time fixed
#' effects via Frisch-Waugh (never re-specified), then reports the naive / df-corrected /
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
#' @param controls optional numeric nuisance-covariate vector or matrix.
#' @param ... passed between methods
#' @return an object of class \code{AdequacyReport}
#' @references Halkiewicz, S. M. S. Corrected diffuse-regime variance
#'   framework for saturated fixed-effect specifications.
#' @examples
#' n <- 200; unit <- rep(1:20, each = 10); time <- rep(1:10, times = 20)
#' x <- rnorm(n); y <- 0.5 * x + rnorm(n)
#' leverage_report(y, x, unit, time)
#' @export
leverage_report <- function(object, ...) UseMethod("leverage_report")

#' @rdname leverage_report
#' @export
leverage_report.default <- function(object, x, unit, time, alpha = 0.05,
                                    delta = 0.05, controls = NULL, ...) {
  y <- object
  cc <- .integer_codes(unit, time)
  uid <- cc$uid; tid <- cc$tid; N <- cc$N; T <- cc$T
  n <- length(uid)
  if (length(y) != n || length(x) != n)
    stop("y, x, unit, time must have equal length")
  partial <- .partial_within_codes(x, uid, tid, N, T, controls = controls)
  xt <- partial$xt
  design <- .design_summary_codes(uid, tid, N, T, xt = xt)
  d_K <- design$d_K
  dof <- n - d_K - partial$rank - 1
  if (dof <= 0) stop("no residual degrees of freedom after fixed effects, controls, and the target regressor (", dof, " <= 0)")
  yt <- .partial_outcome_codes(y, uid, tid, N, T, partial$Q)
  tau_star2 <- design$tau_star2
  if (tau_star2 <= 1e-12 * max(sum(as.numeric(x)^2), 1))
    stop("regressor has no within variation (collinear with the fixed effects)")

  beta <- sum(xt * yt) / tau_star2
  u <- yt - beta * xt
  rss <- sum(u^2)

  p_fe <- .fe_leverage_codes(uid, tid, N, T)
  h_controls <- if (partial$rank) rowSums(partial$Q^2) else rep(0, n)
  H <- p_fe + h_controls + xt^2 / tau_star2
  maxH <- max(H)
  if (maxH >= 1 - 1e-10)
    stop("an observation has full leverage H_ii = 1; HC2/HC3 are undefined ",
         "(degenerate cell -- the diffuse framework's bounded-leverage condition fails)")

  hc0 <- sum(xt^2 * u^2)
  hc1 <- n / dof * hc0
  hc2 <- sum(xt^2 * u^2 / (1 - H))
  hc3 <- sum(xt^2 * u^2 / (1 - H)^2)

  se_naive <- sqrt(rss / n / tau_star2)
  se_df    <- sqrt(rss / dof / tau_star2)
  se_hc0   <- sqrt(hc0) / tau_star2
  se_hc1   <- sqrt(hc1) / tau_star2
  se_hc2   <- sqrt(hc2) / tau_star2
  se_hc3   <- sqrt(hc3) / tau_star2

  rho <- d_K / n
  z <- stats::qnorm(1 - alpha / 2)
  size_naive <- .size_naive(rho, alpha)
  size_hc3 <- .size_hc3(rho, alpha)
  rho_dag <- .rho_dagger(alpha, delta)

  sig <- abs(beta / c(se_df, se_hc0, se_hc1, se_hc2, se_hc3)) > z
  flip <- any(sig != sig[1])

  notes <- c(
    "HC2 is the recommended default among the HC0--HC3 estimators reported here; it is leverage-adjusted but is not the Kline--Saggio--Solvsten leave-out estimator",
    "implied sizes are the conditional-homoskedastic limits of thm:hc (omega_eff^2 = sigma^2); under heteroskedasticity the HCc limits carry the extra factor omega^2/omega_eff^2, omega_eff^2 = (1-rho) omega^2 + rho mu"
  )
  if (rho > 0.1)
    notes <- c(notes, sprintf(
      "HC3 over-correcting regime (rho = %.3f > 0.1): HC3 intervals are artificially conservative (implied size %.1f%%); HC2 is the preferred member of the HC0--HC3 family reported here",
      rho, 100 * size_hc3))
  spread <- maxH / min(H)
  if (spread > 2)
    notes <- c(notes, sprintf(
      "leverage non-uniform (hmax/hmin = %.2f): HC1 is unreliable here; prefer a leverage-adjusted estimator such as HC2, or a separately implemented leave-out estimator",
      spread))
  if (flip)
    notes <- c(notes, "significance at level alpha flips across variance estimators; inference is estimator-dependent, and HC2 is the preferred member of the estimators reported here")
  # the two design conditions the corrected theory actually uses
  lambda_n <- design$lambda_n                       # ass:des(ii)
  n_eff <- design$n_eff
  unif_gap <- max(abs(H - rho))                     # ass:hc(ii)
  Q_hat <- tau_star2 / (n * (1 - rho))              # lem:hess(b) deflation
  notes <- c(notes, sprintf(
    "sample Hessian V_n = %.6g estimates (1-rho) n Qbar, not n Qbar (lem:hess(b)); the deflation-corrected primitive is Qhat = V_n/(n(1-rho)) = %.6g",
    tau_star2, Q_hat))
  lam_bad <- lambda_n > 0.10
  if (lam_bad)
    notes <- c(notes, sprintf(
      "CONCENTRATION WARNING: lambda_n = max_i Xt_i^2/V_n = %.3f does not look negligible (N_eff = %.1f). The Gaussian limit requires lambda_n -> 0. Along persistently concentrated sequences the limit is not fixed across error distributions, and at the fully concentrated boundary no fixed distribution-free critical value is uniformly valid; use the exact contrast test (cycle_report, Paper A) instead.",
      lambda_n, n_eff))
  unif_bad <- unif_gap > 0.25
  if (unif_bad)
    notes <- c(notes, sprintf(
      "uniform-leverage condition strained: max_i |H_ii - rho| = %.3f (rho = %.3f). The HC limits of thm:hc are derived under approximate balance; unbalanced bipartite designs are the leave-out literature's territory, not this module's.",
      unif_gap, rho))
  verdict <- if (size_naive - alpha > delta || flip || lam_bad || unif_bad)
    "FLAGGED" else "CERTIFIED"

  score <- .score_concentration(xt, u)
  if (is.finite(score$lambda_score)) {
    notes <- c(notes, sprintf(
      "realized score diagnostic (Paper A): lambda_score = %.4f, N_eff,score = %.1f. This is a one-realization warning statistic, not by itself a consistent population concentration estimate.",
      score$lambda_score, score$n_eff_score))
  }
  statistic <- list(beta = beta, se_naive = se_naive, se_df = se_df,
                    se_hc0 = se_hc0, se_hc1 = se_hc1, se_hc2 = se_hc2,
                    se_hc3 = se_hc3, t_hc2 = beta / se_hc2,
                    max_leverage = maxH, max_fe_leverage = max(p_fe),
                    leverage_spread = spread, lambda_n = lambda_n,
                    n_eff = n_eff, uniform_leverage_gap = unif_gap,
                    lambda_score = score$lambda_score,
                    H_score = score$H_score,
                    n_eff_score = score$n_eff_score,
                    score_lambda_n = score$lambda_score,
                    score_H_n = score$H_score,
                    score_n_eff = score$n_eff_score,
                    control_rank = partial$rank,
                    V_n = tau_star2, Q_hat = Q_hat,
                    implied_size_hc3 = size_hc3)
  .new_AdequacyReport("leverage", design, statistic, NULL, NULL, rho_dag,
                      size_naive, verdict, alpha, delta, notes)
}

#' Pre-flight applicability check for Paper A
#'
#' For a binary treatment with `n1` treated and `n0` untreated observations, no
#' more than `min(n1, n0)` disjoint contrasts can carry a nonzero loading. A
#' two-sided level-`alpha` test therefore requires
#' `min(n1, n0) >= 1 + log2(1/alpha)`. Other designs pass this
#' structural check; packing quality is assessed separately by
#' [cycle_report()].
#'
#' @param x treatment or regressor vector.
#' @param unit,time fixed-effect identifiers.
#' @param alpha nominal test level.
#' @param controls optional numeric nuisance-covariate vector or matrix.
#' @return A named list with logical `ok` and a plain-language `reason`.
#' @export
applicable <- function(x, unit, time, alpha = 0.05, controls = NULL) {
  if (!(alpha > 0 && alpha < 1)) stop("alpha must lie in (0, 1)")
  design_summary(unit, time, x = x, controls = controls)
  info <- .binary_floor_info(x, alpha)
  list(ok = info$binary_floor_ok, reason = info$reason)
}

#' One-row panel adequacy screen
#'
#' Composes design and cycle diagnostics into one flat named record. The base
#' record requires no outcome. Supplying `y` appends realized score
#' concentration and the Paper A verdict.
#'
#' @param x treatment or regressor vector.
#' @param unit,time fixed-effect identifiers.
#' @param y optional outcome vector.
#' @param method cycle packing method: `"structured"`, `"sparse"`, or
#'   `"greedy"`.
#' @param alpha nominal level used for the Lei--Bickel feasibility condition
#'   and, when `y` is supplied, the cycle-report verdict.
#' @param controls optional numeric nuisance-covariate vector or matrix.
#' @return A flat named list suitable for conversion to a one-row data frame.
#' @export
adequacy_row <- function(x, unit, time, y = NULL,
                         method = c("structured", "sparse", "greedy"),
                         alpha = 0.05, controls = NULL) {
  method <- match.arg(method)
  if (!(alpha > 0 && alpha < 1)) stop("alpha must lie in (0, 1)")
  d <- design_summary(unit, time, x = x, controls = controls)
  if (!(d$tau_star2 > 0)) stop("regressor has no within variation")
  greedy <- cycle_contrasts(x, unit, time, method = "greedy", controls = controls)
  designed <- if (method == "greedy") greedy else
    cycle_contrasts(x, unit, time, method = method, controls = controls)
  binary <- .binary_floor_info(x, alpha)
  out <- list(n = d$n, N = d$N, T = d$T, d_K = d$d_K, rho = d$rho,
              treated_obs = sum(as.numeric(x) != 0), Vn = d$tau_star2,
              n_treated = binary$n_treated,
              binary_floor_ok = binary$binary_floor_ok,
              lambda_n = d$lambda_n, n_eff = d$n_eff,
              lei_bickel_feasible = d$n / (d$d_K + 1) >= 1 / alpha - 1,
              kappa_greedy = greedy$kappa,
              kappa_designed = designed$kappa,
              se_price = if (designed$kappa > 0) 1 / sqrt(designed$kappa) else Inf)
  if (is.null(y)) return(out)
  if (length(y) != d$n) stop("y, x, unit, time must have equal length")
  score <- score_concentration(y, x, unit, time, controls = controls)
  report <- cycle_report(y, x, unit, time, alpha = alpha, method = method,
                         controls = controls, interval = FALSE)
  c(out, score, list(verdict = report$verdict))
}

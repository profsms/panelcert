# The unified result object (spec section 2.2). One S3 class returned by every
# diagnostic; the print method is the adoption-critical legible verdict.

.PATHOLOGY_TITLES <- c(
  leverage            = "Leverage / Variance (Paper A)",
  measurement_error   = "Measurement Error (Paper B)",
  twfe_heterogeneity  = "TWFE Heterogeneity (Paper C)"
)

.new_AdequacyReport <- function(pathology, design, statistic, eta, threshold,
                                breakdown, implied_size, verdict, alpha, delta,
                                notes) {
  if (!pathology %in% names(.PATHOLOGY_TITLES))
    stop("unknown pathology '", pathology, "'")
  if (!verdict %in% c("CERTIFIED", "FLAGGED", "INCONCLUSIVE"))
    stop("unknown verdict '", verdict, "'")
  structure(list(pathology = pathology, design = design, statistic = statistic,
                 eta = eta, threshold = threshold, breakdown = breakdown,
                 implied_size = implied_size, verdict = verdict,
                 alpha = alpha, delta = delta, notes = notes),
            class = "AdequacyReport")
}

# Module-specific rendering of the statistic line(s); mirrors the Julia engine.
.statistic_lines <- function(pathology, s) {
  lines <- character(0)
  has <- function(k) !is.null(s[[k]])
  if (pathology == "measurement_error" && has("lambda_hat")) {
    line <- sprintf("Within reliability lambda_hat = %.3f", s$lambda_hat)
    if (has("noise_ratio") && is.finite(s$noise_ratio))
      line <- paste0(line, sprintf("   ((1-lambda)/lambda = %.3f)", s$noise_ratio))
    lines <- c(lines, line)
    if (has("beta_corr"))
      lines <- c(lines, sprintf("Pilot: beta* = %.4g -> corrected beta0 = %.4g (se %.3g)",
                                s$beta_star, s$beta_corr, s$se_beta_corr))
    if (has("eta_upper"))
      lines <- c(lines, sprintf("Conservative |eta| (upper-bound pilot) = %.3f",
                                s$eta_upper))
  } else if (pathology == "twfe_heterogeneity" && has("Gamma")) {
    line <- sprintf("Design statistic Gamma = %.3f", s$Gamma)
    if (has("neg_share"))
      line <- paste0(line, sprintf("   negative-weight share = %.1f%%", 100 * s$neg_share))
    lines <- c(lines, line)
    if (has("Gamma_CR"))
      lines <- c(lines, sprintf("Cluster-robust Gamma_CR = %.3f (psi_hat = %.3f)",
                                s$Gamma_CR, s$psi_hat))
    if (has("beta"))
      lines <- c(lines, sprintf("TWFE beta_hat = %.4g   sigma = %.4g", s$beta, s$sigma))
    if (has("cohort_sd_raw"))
      lines <- c(lines, sprintf("Cohort dispersion: raw sd %.4g -> shrunk %.4g (factor %.3f)",
                                s$cohort_sd_raw, s$cohort_sd_shrunk, s$shrink_factor))
    if (has("eta_worst"))
      lines <- c(lines, sprintf("Worst-case |eta| (shrunk pilot) = %.3g", s$eta_worst))
  } else if (pathology == "leverage" && has("max_leverage")) {
    line <- sprintf("Max leverage max_i H_ii = %.3f", s$max_leverage)
    if (has("leverage_spread"))
      line <- paste0(line, sprintf(" | spread hmax/hmin = %.2f", s$leverage_spread))
    lines <- c(lines, line)
    if (has("se_cjn")) {
      lines <- c(lines, sprintf("SE(beta): CJN %.4g | HC0 %.4g | HC2/LO %.4g | HC3 %.4g",
                                s$se_cjn, s$se_hc0, s$se_hc2, s$se_hc3))
      lines <- c(lines, sprintf("beta_hat = %.4g   t (HC2/LO) = %.2f", s$beta, s$t_hc2))
    }
  } else {
    for (k in names(s)) lines <- c(lines, paste0(k, " = ", format(s[[k]])))
  }
  lines
}

#' Print an adequacy report
#'
#' Renders the plain-language verdict block (spec section 2.2): design line,
#' pathology-specific statistics, non-centrality vs threshold, implied size,
#' and the VERDICT with any honesty caveats. The format is identical across
#' the Julia, R, and Stata implementations.
#'
#' @param x an \code{AdequacyReport}
#' @param ... unused
#' @return \code{x}, invisibly
#' @export
print.AdequacyReport <- function(x, ...) {
  cat("Panel Adequacy Report \u2014 ", .PATHOLOGY_TITLES[[x$pathology]], "\n", sep = "")
  d <- x$design
  if (d$N > 0) {
    cat(sprintf("Design: n=%d, N=%d, T=%d, d_K=%d, rho=%.4f\n",
                d$n, d$N, d$T, d$d_K, d$rho))
  } else {  # summary-form input: unit/time structure not supplied
    cat(sprintf("Design: n=%d, d_K=%d, rho=%.4f (from summary input)\n",
                d$n, d$d_K, d$rho))
  }
  for (line in .statistic_lines(x$pathology, x$statistic)) cat(line, "\n", sep = "")
  if (!is.null(x$eta)) {
    cat(sprintf("Non-centrality |eta| = %.3f", abs(x$eta)))
    if (!is.null(x$threshold))
      cat(sprintf("   Threshold (delta=%.2g) = %.3f", x$delta, x$threshold))
    cat("\n")
  }
  if (!is.null(x$breakdown))
    cat(sprintf("Breakdown threshold = %.3f\n", x$breakdown))
  if (!is.null(x$implied_size))
    cat(sprintf("Implied size of nominal %.0f%% test: %.1f%%\n",
                100 * x$alpha, 100 * x$implied_size))
  if (x$verdict == "INCONCLUSIVE") {
    cat("VERDICT: INCONCLUSIVE")
  } else {
    cat(sprintf("VERDICT: %s at delta=%.2g", x$verdict, x$delta))
  }
  for (note in x$notes) cat("\nNote: ", note, sep = "")
  cat("\n")
  invisible(x)
}

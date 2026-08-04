# The unified result object (spec section 2.2). One S3 class returned by every
# diagnostic; the print method is the adoption-critical legible verdict.

.PATHOLOGY_TITLES <- c(
  leverage            = "Leverage / Variance (diffuse-regime companion)",
  measurement_error   = "Measurement Error (Paper B)",
  twfe_heterogeneity  = "TWFE Heterogeneity (Paper C)",
  cycle_inference     = "Concentrated Identifying Variation (Paper A)"
)

.new_AdequacyReport <- function(pathology, design, statistic, eta, threshold,
                                breakdown, implied_size, verdict, alpha, delta,
                                notes) {
  if (!pathology %in% names(.PATHOLOGY_TITLES))
    stop("unknown pathology '", pathology, "'")
  if (!verdict %in% c("CERTIFIED", "POINT_PASS", "FLAGGED", "INCONCLUSIVE"))
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
    if (has("Gamma_cmb"))
      lines <- c(lines, sprintf("  restricted ladder: Gamma_c+e = %.3f | Gamma_evt = %.3f | Gamma_coh = %.3f",
                                s$Gamma_cmb, s$Gamma_evt, s$Gamma_coh))
    if (has("Gamma_CR"))
      lines <- c(lines, sprintf("Cluster-robust (psi_hat = %.3f): Gamma_c+e,CR = %.3f | Gamma_CR = %.3f",
                                s$psi_hat, s$Gamma_cmb_CR, s$Gamma_CR))
    if (has("beta"))
      lines <- c(lines, sprintf("TWFE beta_hat = %.4g   sigma = %.4g", s$beta, s$sigma))
    if (has("pilot_cmb"))
      lines <- c(lines, sprintf("Covariance-corrected pilots c_S/sigma: c+e = %.3g | evt = %.3g | coh = %.3g",
                                s$pilot_cmb, s$pilot_evt, s$pilot_coh))
    if (has("size_cmb"))
      lines <- c(lines, sprintf("Worst-case size: combined-class = %.1f%% (headline) | cohort %.1f%% | event %.1f%%",
                                100 * s$size_cmb, 100 * s$size_coh, 100 * s$size_evt))
    if (has("boot") && !is.null(s$boot))
      lines <- c(lines, sprintf("  wild bootstrap (B=%d): combined median %.1f%%, 95%% [%.1f, %.1f]; psi in [%.2f, %.2f]",
                                s$boot$n, 100 * s$boot$cmb_med, 100 * s$boot$cmb_lo,
                                100 * s$boot$cmb_hi, s$boot$psi_lo, s$boot$psi_hi))
    if (has("size_realized"))
      lines <- c(lines, sprintf("Realized-profile size (CR) = %.1f%%", 100 * s$size_realized))
  } else if (pathology == "cycle_inference" && has("kappa")) {
    lines <- c(lines, sprintf("Concentration: lambda_n = %.4f (N_eff = %.1f)",
                              s$lambda_n, s$n_eff))
    if (has("score_lambda_n") && is.finite(s$score_lambda_n))
      lines <- c(lines, sprintf("Realized score concentration: lambda_score = %.4f (N_eff,score = %.1f)",
                                s$score_lambda_n, s$score_n_eff))
    ctext <- if (has("effective_C") && s$effective_C != s$C)
      sprintf("%d supports (%d treatment-loaded)", s$C, s$effective_C)
    else sprintf("%d supports", s$C)
    lines <- c(lines, sprintf(
      "Capture kappa_C = %.4f over %s (cycle-space dim %d) | capture-implied SE ratio %.3fx | max share %.3f",
      s$kappa, ctext, s$cycle_dim, s$se_price, s$max_share))
    if (has("reason") && nzchar(s$reason))
      lines <- c(lines, paste0("Reason: ", s$reason))
    if (!is.null(s$beta_tilde)) {
      l <- sprintf("Contrast estimate beta~ = %.4g", s$beta_tilde)
      if (!is.null(s$ci_lo) && !is.na(s$ci_lo)) {
        l <- paste0(l, sprintf("   exact %.0f%% set: [%.4g, %.4g]%s",
                    100 * (if (has("ci_level")) s$ci_level else 0.95),
                    s$ci_lo, s$ci_hi,
                    if (isTRUE(s$ci_grid_truncated)) " (conservative: grid boundary reached)" else ""))
      } else {
        l <- paste0(l, "   exact set: EMPTY at this level")
      }
      lines <- c(lines, l)
    }
  } else if (pathology == "leverage" && has("max_leverage")) {
    line <- sprintf("Max leverage max_i H_ii = %.3f", s$max_leverage)
    if (has("leverage_spread"))
      line <- paste0(line, sprintf(" | spread hmax/hmin = %.2f", s$leverage_spread))
    lines <- c(lines, line)
    if (has("lambda_n"))
      lines <- c(lines, sprintf("Design conditions: lambda_n = %.4f (N_eff = %.1f) | max|H_ii - rho| = %.3f",
                                s$lambda_n, s$n_eff, s$uniform_leverage_gap))
    if (has("score_lambda_n") && is.finite(s$score_lambda_n))
      lines <- c(lines, sprintf("Realized score concentration: lambda_score = %.4f (N_eff,score = %.1f)",
                                s$score_lambda_n, s$score_n_eff))
    if (has("se_df")) {
      lines <- c(lines, sprintf("SE(beta): df-corrected %.4g | HC0 %.4g | HC2 %.4g | HC3 %.4g",
                                s$se_df, s$se_hc0, s$se_hc2, s$se_hc3))
      lines <- c(lines, sprintf("beta_hat = %.4g   t (HC2) = %.2f", s$beta, s$t_hc2))
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
#' and the VERDICT. Cycle-inference caveats remain in `x$notes` and are shown
#' only when requested.
#'
#' @param x an \code{AdequacyReport}
#' @param notes logical; whether to print detailed diagnostic notes. Defaults
#'   to `FALSE` for cycle-inference reports and `TRUE` for other reports.
#' @param ... unused
#' @return \code{x}, invisibly
#' @export
print.AdequacyReport <- function(
    x, notes = x$pathology != "cycle_inference", ...) {
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
  if (x$verdict == "POINT_PASS") {
    cat(sprintf("VERDICT: POINT PASS at delta=%.2g (descriptive \u2014 not a certificate)",
                x$delta))
  } else if (x$verdict == "INCONCLUSIVE") {
    cat("VERDICT: INCONCLUSIVE")
  } else {
    cat(sprintf("VERDICT: %s at delta=%.2g", x$verdict, x$delta))
  }
  if (isTRUE(notes)) {
    if (x$pathology == "cycle_inference") {
      cat("\n")
      .print_report_notes(x)
    } else {
      for (note in x$notes) cat("\nNote: ", note, sep = "")
    }
  } else if (x$pathology == "cycle_inference" && length(x$notes)) {
    cat(sprintf(
      "\nDiagnostic notes hidden (%d); call show_notes(report) to display them.",
      length(x$notes)))
  }
  cat("\n")
  invisible(x)
}

.print_report_notes <- function(x) {
  if (!length(x$notes)) {
    cat("No diagnostic notes.\n")
    return(invisible(NULL))
  }
  cat(sprintf("Diagnostic notes for %s (%d):\n",
              .PATHOLOGY_TITLES[[x$pathology]], length(x$notes)))
  for (i in seq_along(x$notes)) {
    if (i > 1L) cat("\n")
    cat(i, ". ", x$notes[[i]], "\n", sep = "")
  }
  invisible(NULL)
}

#' Display detailed diagnostic notes
#'
#' Prints the notes stored in an [AdequacyReport] as a numbered list. These
#' notes remain programmatically available even when the default report display
#' is concise.
#'
#' @param x an `AdequacyReport`.
#' @return `x`, invisibly.
#' @export
show_notes <- function(x) {
  if (!inherits(x, "AdequacyReport"))
    stop("x must be an AdequacyReport")
  .print_report_notes(x)
  invisible(x)
}

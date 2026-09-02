#' @keywords internal
#' @details
#' Inference and diagnostics: [cycle_report()] and [contrast_system()] (Paper A,
#' exact contrast inference), [eiv_adequacy()] (Paper B, measurement error),
#' [twfe_design()] / [twfe_adequacy()] (staggered-DiD heterogeneity),
#' and [leverage_report()] (the diffuse-regime variance companion). Shared
#' primitives include [design_summary()], [twoway_demean()],
#' [multiway_demean()], and [fe_leverage()]. Reports use a plain-language
#' verdict; see `vignette("panelcert")` and `vignette("cycle-inference")`.
"_PACKAGE"

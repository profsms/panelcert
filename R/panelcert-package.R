#' @keywords internal
#' @details
#' Diagnostics: [leverage_report()] (Module A, variance estimators under FE
#' saturation), [eiv_adequacy()] (Module B, measurement error), and
#' [twfe_design()] / [twfe_adequacy()] (Module C, staggered-DiD TWFE
#' heterogeneity). Shared primitives: [design_summary()], [twoway_demean()],
#' [fe_leverage()]. Every diagnostic returns an `AdequacyReport` with a
#' plain-language verdict; see `vignette("paneldiagnostics")`.
"_PACKAGE"

# Documentation for the bundled reference datasets (spec section 7 / the paper
# applications). Shipped as lazy-loaded objects: after library(panelcert) they
# are available by name, so the article applications run with no data download.

#' V-Dem democracy--growth panel (Paper B lead application)
#'
#' Country--year panel used in the measurement-error application of Paper B:
#' log GDP per capita against continuous V-Dem institutional indices, each
#' accompanied by the measurement model's posterior standard deviation (so the
#' measurement-error input is data, not a calibration). The two-pole contrast
#' the paper reports -- polyarchy certified, the constraint sub-indices flagged
#' -- is reproduced directly from this object.
#'
#' @format A data frame with 8931 country-year rows and 9 variables:
#' \describe{
#'   \item{iso}{country ISO code (unit id)}
#'   \item{year}{calendar year (time id)}
#'   \item{ly}{log GDP per capita (Maddison Project 2020)}
#'   \item{v2x_polyarchy, v2x_polyarchy_sd}{electoral-democracy index and its posterior SD}
#'   \item{v2xlg_legcon, v2xlg_legcon_sd}{legislative-constraints index and its posterior SD}
#'   \item{v2x_jucon, v2x_jucon_sd}{judicial-constraints index and its posterior SD}
#' }
#' @source V-Dem dataset (measurement-model posterior SDs); Maddison Project
#'   Database 2020.
#' @examples
#' d <- vdem[stats::complete.cases(vdem$ly, vdem$v2x_polyarchy, vdem$v2x_polyarchy_sd), ]
#' eiv_adequacy(d$ly, d$v2x_polyarchy, d$iso, d$year, sigma_nu = d$v2x_polyarchy_sd)
"vdem"

#' Cornwell--Rupert PSID earnings panel (Paper B second application)
#'
#' Person--year panel of weeks worked and the log wage, the canonical noisy
#' regressor of the validation literature. Used with an external reliability
#' ratio (rather than a per-observation posterior) to show the fixed-effect
#' transformation moving the same regressor from certified to flagged.
#'
#' @format A data frame with 4165 person-year rows and 4 variables:
#' \describe{
#'   \item{id}{person id (unit)}
#'   \item{year}{year (time)}
#'   \item{lwage}{log wage (the mismeasured regressor)}
#'   \item{wks}{weeks worked (outcome)}
#' }
#' @source Cornwell--Rupert PSID extract, distributed with the \pkg{plm} package.
#' @examples
#' eiv_adequacy(psid$wks, psid$lwage, psid$id, psid$year,
#'              reliability = 0.65, pilot = "point")
"psid"

#' Castle-doctrine adoption panel (Paper C, certified pole)
#'
#' State--year panel for the Cheng--Hoekstra castle-doctrine design: a large
#' never-treated reservoir, no negative weights, design statistic Gamma = 0.21.
#' The certified pole of Paper C's two-design audit.
#'
#' @format A data frame with 550 state-year rows and 4 variables:
#' \describe{
#'   \item{uid}{state id (unit)}
#'   \item{tid}{year code (time)}
#'   \item{ft}{first-treatment period; \code{NA} for never-treated states}
#'   \item{y}{log homicide rate (outcome)}
#' }
#' @source Cheng and Hoekstra castle-doctrine replication data; analysis panel
#'   derived as in the Paper C audit.
#' @examples
#' twfe_design(castle$uid, castle$tid, castle$ft)
#' twfe_adequacy(castle$y, castle$uid, castle$tid, castle$ft)
"castle"

#' No-fault-divorce adoption panel (Paper C, flagged pole)
#'
#' State--year panel for the Stevenson--Wolfers no-fault-divorce design
#' (Goodman-Bacon's pathology example): near-universal eventual adoption, a
#' 7.3\% negative-weight share, design statistic Gamma = 0.86. The flagged pole
#' of Paper C's audit -- implied size 60\% under i.i.d. errors, 42\% after the
#' cluster correction.
#'
#' @format A data frame with 1377 state-year rows and 4 variables:
#' \describe{
#'   \item{uid}{state id (unit)}
#'   \item{tid}{year code (time)}
#'   \item{ft}{first-treatment period; \code{NA} never-treated, \code{0} always-treated}
#'   \item{y}{female suicide rate per 100k (outcome)}
#' }
#' @source Stevenson and Wolfers divorce data (via the \pkg{bacondecomp}
#'   distribution); analysis panel derived as in the Paper C audit.
#' @examples
#' twfe_adequacy(divorce$y, divorce$uid, divorce$tid, divorce$ft)
"divorce"

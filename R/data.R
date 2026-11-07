# Documentation for the bundled reference datasets (spec section 7 / the paper
# applications). Shipped as lazy-loaded objects: after library(panelcert) they
# are available by name, so the article applications run with no data download.

#' V-Dem democracy--growth measurement-error application
#'
#' Country--year panel used in *Breakdown Reliability for Saturated
#' Fixed-Effect Inference*:
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
#' @source V-Dem dataset at commit
#'   `f4dd26922e658442524dfd954bf14f7ebe622d5d` (measurement-model posterior
#'   SDs); Maddison Project Database 2020. The raw V-Dem artifact is pinned by
#'   SHA-256 in `inst/DATA_SOURCES.md`.
#' @examples
#' d <- vdem[stats::complete.cases(vdem$ly, vdem$v2x_polyarchy, vdem$v2x_polyarchy_sd), ]
#' eiv_adequacy(d$ly, d$v2x_polyarchy, d$iso, d$year, sigma_nu = d$v2x_polyarchy_sd)
"vdem"

#' Ashenfelter--Krueger repeated-report twins extract
#'
#' Public teaching extract from the Ashenfelter--Krueger identical-twins design.
#' Each twin reports both siblings' schooling, producing two measurements of the
#' within-pair schooling difference. The article's controlled calculation uses
#' the 147 complete rows on `DLHRWAGE`, `DEDUC1`, `DEDUC2`, `DTEN`, `DMARRIED`,
#' and `DUNCOV`, applying the identical control projection and sample to both
#' schooling reports.
#'
#' @format A data frame with 183 twin-pair rows and 16 variables:
#' \describe{
#'   \item{DLHRWAGE}{within-pair difference in log hourly wages}
#'   \item{DEDUC1}{schooling difference from each twin's self-report}
#'   \item{AGE}{age of the twin pair}
#'   \item{AGESQ}{squared age}
#'   \item{HRWAGEH}{hourly wage of the H-labelled twin}
#'   \item{WHITEH}{white indicator for the H-labelled twin}
#'   \item{MALEH}{male indicator for the H-labelled twin}
#'   \item{EDUCH}{years of schooling of the H-labelled twin}
#'   \item{HRWAGEL}{hourly wage of the L-labelled twin}
#'   \item{WHITEL}{white indicator for the L-labelled twin}
#'   \item{MALEL}{male indicator for the L-labelled twin}
#'   \item{EDUCL}{years of schooling of the L-labelled twin}
#'   \item{DEDUC2}{schooling difference from the co-twin reports}
#'   \item{DTEN}{within-pair tenure difference}
#'   \item{DMARRIED}{within-pair married-status difference}
#'   \item{DUNCOV}{within-pair union-coverage difference}
#' }
#' @source `RbyExample::twins` version 0.0.100, attributed to Ashenfelter and
#'   Krueger (1994), *American Economic Review* 84(5), 1157--1173. Upstream
#'   package license: GPL (>= 2).
#' @examples
#' keep <- stats::complete.cases(twins[c("DLHRWAGE", "DEDUC1", "DEDUC2",
#'                                       "DTEN", "DMARRIED", "DUNCOV")])
#' W <- stats::model.matrix(~ DTEN + DMARRIED + DUNCOV, data = twins[keep, ])
#' x <- qr.resid(qr(W), twins$DEDUC1[keep])
#' z <- qr.resid(qr(W), twins$DEDUC2[keep])
#' reliability_from_repeats(x, z)
"twins"

#' Cornwell--Rupert PSID earnings panel
#'
#' Person--year panel of weeks worked and the log wage, the canonical noisy
#' regressor of the validation literature. Retained as a general panel example
#' and for backward compatibility; it is not an application in the current
#' measurement-error article.
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

#' Castle-doctrine adoption panel (additive-class certificate)
#'
#' State--year panel for the Cheng--Hoekstra castle-doctrine design: a large
#' never-treated reservoir, no negative weights, design statistic Gamma = 0.21.
#' The homicide application certifies in the prespecified additive class. The
#' saturated group-time class is structurally inconclusive because its
#' projected score covariance has rank 18 rather than 19.
#'
#' @format A data frame with 550 state-year rows and 4 variables:
#' \describe{
#'   \item{uid}{state id (unit)}
#'   \item{tid}{year code (time)}
#'   \item{ft}{first-treatment period; \code{NA} for never-treated states}
#'   \item{y}{log homicide rate (outcome)}
#' }
#' @source Cheng and Hoekstra castle-doctrine replication data; analysis panel
#'   derived as in the TWFE-heterogeneity audit.
#' @examples
#' twfe_design(castle$uid, castle$tid, castle$ft)
#' twfe_adequacy(castle$y, castle$uid, castle$tid, castle$ft,
#'               heterogeneity_class = "additive")
"castle"

#' No-fault-divorce adoption panel (rank-inconclusive application)
#'
#' State--year panel for the Stevenson--Wolfers no-fault-divorce design
#' (Goodman-Bacon's pathology example): near-universal eventual adoption.
#' After the two always-treated states are dropped, the analysis design has a
#' 1.1\% negative-weight share and Gamma = 0.64. Its direct-CR1 point envelope
#' is 41.5\%, but bootstrap uncertainty is wide and the fixed-\eqn{T} warning
#' applies.
#'
#' @format A data frame with 1377 state-year rows and 4 variables:
#' \describe{
#'   \item{uid}{state id (unit)}
#'   \item{tid}{year code (time)}
#'   \item{ft}{first-treatment period; \code{NA} never-treated, \code{0} always-treated}
#'   \item{y}{female suicide rate per 100k (outcome)}
#' }
#' @source Stevenson and Wolfers divorce data (via the \pkg{bacondecomp}
#'   distribution); analysis panel derived as in the TWFE-heterogeneity audit.
#' @examples
#' \donttest{
#' twfe_adequacy(divorce$y, divorce$uid, divorce$tid, divorce$ft)
#' }
"divorce"

#' Minimum-wage county panel (headline TWFE application)
#'
#' Exact balanced county--year panel for the unconditional specification in
#' Callaway and Sant'Anna (2021). The current diagnostic uses never-treated counties as the
#' comparison group. All treated-cell TWFE weights are positive, yet the
#' additive cohort-plus-event class fails uniform inference certification.
#'
#' @format A data frame with 15,988 county-year rows and 4 variables:
#' \describe{
#'   \item{uid}{county id (unit)}
#'   \item{tid}{year code, 1--7 (time)}
#'   \item{ft}{first-treatment period; \code{NA} for never-treated counties}
#'   \item{y}{log teen employment (outcome)}
#' }
#' @source Callaway and Sant'Anna (2021) public minimum-wage replication data;
#'   exact published-sample extract used in the TWFE-heterogeneity article.
#' @examples
#' twfe_design(minimum_wage$uid, minimum_wage$tid, minimum_wage$ft)
#' \donttest{
#' twfe_adequacy(minimum_wage$y, minimum_wage$uid, minimum_wage$tid,
#'               minimum_wage$ft, controls = "never", bootstrap = 19)
#' }
"minimum_wage"

#' Brazilian property-tax panel (large-panel TWFE application)
#'
#' Balanced 2004--2015 municipal panel for the Christensen--Garfias property-
#' tax specification. Municipalities treated in 2004 are excluded so every
#' treated group-time effect has an in-sample untreated baseline. This is the
#' balanced-subset reanalysis in Paper C, rather than the full published
#' estimation sample.
#'
#' @format A data frame with 34,080 municipality-year rows and 4 variables:
#' \describe{
#'   \item{uid}{consecutive municipality id}
#'   \item{tid}{period code, 1--12 for 2004--2015}
#'   \item{ft}{first-treatment period; \code{NA} for never-treated municipalities}
#'   \item{y}{log property-tax revenue (\code{logiptu})}
#' }
#' @source Christensen and Garfias (2021) replication data distributed in the
#'   public Chiu, Lan, Liu, and Xu causal-panel reanalysis archive,
#'   doi:10.7910/DVN/9RJFZF.
#' @examples
#' twfe_design(brazil$uid, brazil$tid, brazil$ft)
#' \donttest{
#' twfe_adequacy(brazil$y, brazil$uid, brazil$tid, brazil$ft,
#'               bootstrap = 19)
#' }
"brazil"

#' Piotroski F-Score / Visegrad firm panel (Paper A application)
#'
#' Hand-collected firm--year panel of Piotroski F-Scores and one-year-ahead
#' returns for firms on the Warsaw, Budapest and Prague exchanges, 2010--2024,
#' consolidated from the three exchange production files. Paper A uses the
#' firm and country--year specification for a complete exact-inference workflow.
#' The diffuse-regime variance companion also uses seven outcome/subsample/FE
#' variants to compare HC0--HC3 behavior.
#'
#' @format A data frame with 217 firm-year rows and 6 variables:
#' \describe{
#'   \item{uid}{firm ticker (unit id)}
#'   \item{year}{fiscal year (time id)}
#'   \item{country}{Poland, Hungary or Czech Republic}
#'   \item{status}{\code{"active"} or \code{"delisted"}}
#'   \item{fscore}{Piotroski F-Score, 0--9 integer composite (the regressor)}
#'   \item{ret}{one-year-ahead return (outcome)}
#' }
#' @source Hand-collected from Warsaw (WSE), Budapest (BSE) and Prague (PSE)
#'   exchange filings; the panel of Paper A's empirical application.
#' @examples
#' # Paper A specification: log return on F-Score, firm + country-year effects
#' cy <- interaction(fscore$country, fscore$year, drop = TRUE)
#' cycle_capture(fscore$fscore, fscore$uid, cy)
"fscore"

#' Canonical 11-firm Grunfeld investment panel
#'
#' The complete corrected 1935--1954 panel used for Paper A's concentrated-
#' variation showcase. The canonical regression is investment on capital and
#' market value with firm and year fixed effects. Several 10- and five-firm
#' variants circulate; Kleiber and Zeileis (2010) document their omissions and
#' transcription errors.
#'
#' @format A data frame with 220 firm-year rows and 5 variables:
#' \describe{
#'   \item{invest}{gross investment in 1947 dollars}
#'   \item{value}{market value at year end in 1947 dollars}
#'   \item{capital}{stock of plant and equipment in 1947 dollars}
#'   \item{firm}{firm name}
#'   \item{year}{calendar year, 1935--1954}
#' }
#' @source Public-domain statsmodels Grunfeld data at commit
#'   \code{57169f2cfc7089141513ed6103f79fa56ab213ac}, accessed 2026-08-02;
#'   reconstructed from Grunfeld (1958) by Kleiber and Zeileis (2010).
#' @examples
#' applicable(grunfeld$capital, grunfeld$firm, grunfeld$year,
#'            controls = grunfeld$value)
#' \donttest{
#' cycle_report(grunfeld$invest, grunfeld$capital, grunfeld$firm,
#'              grunfeld$year, controls = grunfeld$value, interval = FALSE)
#' }
"grunfeld"

#' Weekend allied-health stepped-wedge trial (Haines et al. 2017, Trial 1)
#'
#' Analysis extract of Trial 1 of Haines et al. (2017): withdrawal of weekend
#' allied-health services from 12 acute wards in two hospitals, a complete
#' stepped wedge with 7 periods and 2 wards starting in each of periods 2--7.
#' One row per patient episode. Used by Halkiewicz and Alexeev
#' (arXiv:2610.35853) as the application with visible exposure-time
#' heterogeneity; the exposure weights of the IT estimator at the trial's actual
#' cell sizes give an adequacy radius of 0.124.
#'
#' @format A data frame with 14,834 rows and 6 variables:
#' \describe{
#'   \item{ward}{ward identifier (cluster)}
#'   \item{site}{hospital}
#'   \item{period}{study period 1--7}
#'   \item{start}{first period in which the ward was treated (sequence)}
#'   \item{trt}{1 if the weekend service had been withdrawn (intervention)}
#'   \item{acute_los}{acute length of stay in days (analysed on the log scale)}
#' }
#' @source Haines TP, Bowles K-A, Mitchell D, et al. (2017). Impact of
#'   disinvestment from weekend allied health services across acute medical
#'   and surgical wards: 2 stepped-wedge cluster randomised controlled trials.
#'   \emph{PLOS Medicine} 14(10): e1002412, \doi{10.1371/journal.pmed.1002412};
#'   S2 Data, \doi{10.1371/journal.pmed.1002412.s002}. Licence: CC BY 4.0 (not
#'   covered by the package's MIT licence). Only the six variables above are
#'   retained; see \code{data-raw/make_sw_trials.py} and
#'   \code{inst/DATA_SOURCES.md}.
#' @examples
#' n <- with(haines2017, tapply(rep(1, nrow(haines2017)), list(ward, period), sum))
#' starts <- tapply(haines2017$start, haines2017$ward, `[`, 1)
#' ord <- order(starts, as.numeric(names(starts)))
#' sw_radius(design = c(2, 2, 2, 2, 2, 2, 0), icc = 0.0709, cac = 1,
#'           sizes = n[ord, ])
"haines2017"

#' 99DOTS digital adherence stepped-wedge trial, Uganda (Cattamanchi et al. 2021)
#'
#' Patient-level public data of the 99DOTS stepped-wedge trial: 18 health
#' facilities in six randomisation blocks, eight trial months, treatment
#' success as the primary outcome. Used by Halkiewicz and Alexeev
#' (arXiv:2610.35853) as the application with no detectable heterogeneity but
#' low power to rule out heterogeneity that would damage IT inference.
#'
#' @format A data frame with 2,224 rows and 13 variables:
#' \describe{
#'   \item{block}{intervention group (randomisation block 1--6)}
#'   \item{facility}{health facility 1--18 (cluster)}
#'   \item{facility_type}{facility type}
#'   \item{month}{trial month 0--7}
#'   \item{period}{study period: Control, Buffer, or Intervention}
#'   \item{sex, age, hiv, disease_class, retreatment}{patient characteristics}
#'   \item{treated_successfully, persisting, not_ltfu}{outcomes (Yes/No)}
#' }
#' @source Cattamanchi A, Crowder R, Kityamuwesi A, et al. (2021). Digital
#'   adherence technology for tuberculosis treatment supervision: a
#'   stepped-wedge cluster-randomized trial in Uganda. \emph{PLOS Medicine}
#'   18(5): e1003628, \doi{10.1371/journal.pmed.1003628}; S1 Data,
#'   \doi{10.1371/journal.pmed.1003628.s002}. Licence: CC BY 4.0 (not covered
#'   by the package's MIT licence). Column names are snake_case; commas inside
#'   text labels are replaced by semicolons. See \code{inst/DATA_SOURCES.md}.
#' @examples
#' table(dots99$block, dots99$period)
"dots99"

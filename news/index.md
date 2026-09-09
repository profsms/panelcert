# Changelog

## panelcert 0.8.0

- Replaced the TWFE confidence-ball upper decision with covariance-aware
  projected-Wald noncentral-chi-square inversion. An upper certificate
  is now issued only when the estimated score covariance spans the full
  prespecified heterogeneity class; rank-deficient classes are reported
  as structurally `INCONCLUSIVE` unless the one-sided lower test already
  flags them.
- Added `heterogeneity_class` so the report verdict can be based on a
  prespecified saturated group-time, additive, cohort, or event-time
  class. Every report exposes the separate lower bound, upper bound,
  rank, and verdict for all four classes.
- Retained the regular multiplier/reverse-triangle lower test without
  the old denominator buffer under the bounded local experiment. The
  more general full confidence-ball construction is available through
  `q_band`; it is reported separately and never determines the main
  verdict.
- Updated the castle-doctrine lock: the additive homicide class
  certifies with an upper bound near 0.236, while the saturated class is
  inconclusive because its projected covariance has rank 18/19. The
  divorce saturated class is likewise reported as rank-deficient
  (49/167).
- Reworked the covariance algebra in group-time coordinates, eliminating
  the quadratic treated-cell covariance allocation. Bundled the
  reproducible 2,840-municipality Brazil property-tax extract and locked
  the finding that every class is flagged by its lower test.

## panelcert 0.7.1

- [`eiv_adequacy()`](https://profsms.github.io/panelcert/reference/eiv_adequacy.md)
  now uses the exact leverage-weighted noise trace
  `sum((1 - h_ii) * sigma_nu_i^2)` for observation-specific
  measurement-error standard deviations. Both this trace and the
  common-SD formula now use the full nuisance projection (fixed effects
  plus supplied controls).
- Added nuisance-control support to the raw-data method and the `lm`
  adapter, including correct residual degrees of freedom and
  control-aware projection- compatibility diagnostics.
- Replaced the bundled V-Dem interval-half-width proxy with V-Dem’s
  direct posterior-standard-deviation fields and refreshed the article
  reproduction locks.
  [`reliability_from_interval()`](https://profsms.github.io/panelcert/reference/reliability_from_interval.md)
  now documents that its conversion is valid only when an interval
  half-width is substantively calibrated as one error standard
  deviation.

## panelcert 0.7.0

- Added
  [`certified_breakdown_reliability()`](https://profsms.github.io/panelcert/reference/certified_breakdown_reliability.md)
  and exposed the paper’s closed-form certified threshold in every
  measurement-error report.
- Separated coefficient and reliability uncertainty: `reliability_lower`
  and `gamma_lambda` implement the lower-bound rule, while the existing
  `gamma` argument remains the coefficient budget. Reports state when
  certification is conditional on treating a reliability input as
  known/consistent.
- Aligned the TWFE module with the current *Positive Weights Do Not
  Certify TWFE Inference* theory: direct CR1 cluster-score normalization
  is now the default, while AR(1), i.i.d., and user-supplied scales
  remain explicit sensitivity routes.
- Added the saturated post-treatment group-time class and verify its
  balanced- panel identity `Gamma_gt = Gamma`. Reports now distinguish
  that full class from the cohort, event-time, and additive
  cohort-plus-event-time ladder.
- Replaced point-envelope certification by a one-sided confidence
  construction for the projected group-time vector. The resulting HC2
  projected-norm bounds remain valid at zero heterogeneity; HC3 bounds
  are reported as a sensitivity check. Trace-debiased quadratic pilots
  and percentile summaries remain explicitly descriptive.
- Separated all worst-case size envelopes from the signed directional
  plug-in. Reports label the former as uniform upper bounds, never as
  realized rejection probabilities, and expose directional alignment.
- Replaced the dense fixed-effect bootstrap regression with an absorbed,
  cluster-score implementation. This makes large panels practical while
  retaining the covariance-aware group-time pilot.
- Added the exact 15,988-row, 2,284-county Callaway–Sant’Anna
  minimum-wage analysis extract. The reported target is now the equally
  weighted treated- cell ATT (about -5.2%), matching the theory rather
  than the differently weighted group aggregate. The saturated-class
  point envelope is about 31.6%, and its projected-norm lower bound
  formally withholds certification.
- Added never-treated versus not-yet-treated comparison-group selection,
  direct and AR(1) scale outputs, CR1 standard errors, sign-reversal
  RMS, and a fixed-panel-length warning when `T` is large relative to
  the cluster count.
- Distinguished Gaussian mean-shift calculations from procedure-specific
  finite-cluster inference: alternative standard errors or bootstrap
  tests are not treated as drop-in changes to the Gaussian rejection
  map.
- Retained `eta_real_cr`, `size_realized`, and the old bootstrap
  interval names as deprecated compatibility aliases; they now point to
  explicitly named directional or envelope fields.

## panelcert 0.6.0

- Added
  [`reliability_from_repeats()`](https://profsms.github.io/panelcert/reference/reliability_from_repeats.md)
  with covariance and equal-variance methods.
- Bundled the public Ashenfelter–Krueger twins extract and reproduced
  the current twins table, including the Rouse correlated-report
  sensitivity.
- Replaced obsolete PSID article locks while retaining the dataset for
  backward compatibility; locked the Design 4/5 exact-normal endpoints
  at 3.58 and 3.67.
- Corrected stale cluster-direction guidance and refreshed documentation
  and data provenance.

## panelcert 0.5.1

- [`cycle_report()`](https://profsms.github.io/panelcert/reference/cycle_report.md)
  keeps all detailed caveats in `report$notes` but no longer prints them
  by default.
- Added `show_notes(report)` and `print(report, notes = TRUE)` for an
  explicit, numbered display; structural `INCONCLUSIVE` reasons remain
  visible in the concise report.

## panelcert 0.5.0

- Bundle the canonical public-domain 11-firm Grunfeld panel used by
  Paper A’s concentrated-regime showcase, with pinned source and version
  provenance.
- [`design_summary()`](https://profsms.github.io/panelcert/reference/design_summary.md)
  now returns `lambda_n` and `n_eff`;
  [`leverage_report()`](https://profsms.github.io/panelcert/reference/leverage_report.md)
  uses the same calculation.
- Added
  [`twfe_gammas()`](https://profsms.github.io/panelcert/reference/twfe_gammas.md)
  for Julia parity.
- Added
  [`score_concentration()`](https://profsms.github.io/panelcert/reference/score_concentration.md),
  multiway `design_summary(fe_levels = ...)`, and
  [`adequacy_row()`](https://profsms.github.io/panelcert/reference/adequacy_row.md).
- Added
  [`applicable()`](https://profsms.github.io/panelcert/reference/applicable.md)
  and surfaced the structural binary-treatment granularity floor in
  [`cycle_report()`](https://profsms.github.io/panelcert/reference/cycle_report.md)
  and
  [`adequacy_row()`](https://profsms.github.io/panelcert/reference/adequacy_row.md).
- Canonicalized cycle packing and added deterministic multi-start
  traversal so input row order cannot change the selected design.
- Locked the stable v0.5.0 captures for the public fixtures: KSS match
  `0.5110`, KSS wage `0.6112`, F-score `0.8274`, and calibrated dense
  `0.8341`; each meets or improves on the article’s former lower-bound
  construction.
- Added cross-language fixtures, parity checks, property tests, and
  published-result regression locks.
- Added the external-use README, concentrated-design vignette, CI,
  coverage, pkgdown, and release metadata.
- Added explicit numeric nuisance-control support to design, score,
  leverage, screen, and exact-cycle APIs. For one continuous control,
  dense exact inference uses locally projected 2-by-3 supports; model
  adapters no longer diagnose a different univariate regression
  silently.
- Locked the canonical 11-firm Grunfeld capital specification at
  `lambda_score = 0.7388` and valid controlled capture
  `kappa_C = 0.6270` over 32 supports.

# panelcert 0.7.0

- Aligned the TWFE module with the current *Is Bias Correction Enough?*
  theory: direct CR1 cluster-score normalization is now the default, while
  AR(1), i.i.d., and user-supplied scales remain explicit sensitivity routes.
- Separated the combined-class worst-case size envelope from the signed
  directional plug-in. Reports now label the former as a uniform upper bound,
  never as realized rejection probability, and expose directional alignment.
- Replaced the dense fixed-effect bootstrap regression with an absorbed,
  cluster-score implementation. This makes large panels practical while
  retaining the covariance-aware group-time pilot.
- Added the exact 15,988-row, 2,284-county Callaway--Sant'Anna minimum-wage
  analysis extract and locked its 29.3% envelope, 27.8% directional diagnostic,
  and sign-reversal calibration.
- Added never-treated versus not-yet-treated comparison-group selection,
  direct and AR(1) scale outputs, CR1 standard errors, sign-reversal RMS, and a
  fixed-panel-length warning when `T` is large relative to the cluster count.
- Retained `eta_real_cr`, `size_realized`, and the old bootstrap interval names
  as deprecated compatibility aliases; they now point to explicitly named
  directional or envelope fields.

# panelcert 0.6.0

- Added `reliability_from_repeats()` with covariance and equal-variance methods.
- Bundled the public Ashenfelter--Krueger twins extract and reproduced the
  current twins table, including the Rouse correlated-report sensitivity.
- Replaced obsolete PSID article locks while retaining the dataset for backward
  compatibility; locked the Design 4/5 exact-normal endpoints at 3.58 and 3.67.
- Corrected stale cluster-direction guidance and refreshed documentation and
  data provenance.

# panelcert 0.5.1

- `cycle_report()` keeps all detailed caveats in `report$notes` but no longer prints them by default.
- Added `show_notes(report)` and `print(report, notes = TRUE)` for an explicit, numbered display; structural `INCONCLUSIVE` reasons remain visible in the concise report.

# panelcert 0.5.0

* Bundle the canonical public-domain 11-firm Grunfeld panel used by Paper A's concentrated-regime showcase, with pinned source and version provenance.
- `design_summary()` now returns `lambda_n` and `n_eff`; `leverage_report()` uses the same calculation.
- Added `twfe_gammas()` for Julia parity.
- Added `score_concentration()`, multiway `design_summary(fe_levels = ...)`, and `adequacy_row()`.
- Added `applicable()` and surfaced the structural binary-treatment granularity floor in `cycle_report()` and `adequacy_row()`.
- Canonicalized cycle packing and added deterministic multi-start traversal so input row order cannot change the selected design.
- Locked the stable v0.5.0 captures for the public fixtures: KSS match `0.5110`, KSS wage `0.6112`, F-score `0.8274`, and calibrated dense `0.8341`; each meets or improves on the article's former lower-bound construction.
- Added cross-language fixtures, parity checks, property tests, and published-result regression locks.
- Added the external-use README, concentrated-design vignette, CI, coverage, pkgdown, and release metadata.
- Added explicit numeric nuisance-control support to design, score, leverage, screen, and exact-cycle APIs. For one continuous control, dense exact inference uses locally projected 2-by-3 supports; model adapters no longer diagnose a different univariate regression silently.
- Locked the canonical 11-firm Grunfeld capital specification at `lambda_score = 0.7388` and valid controlled capture `kappa_C = 0.6270` over 32 supports.

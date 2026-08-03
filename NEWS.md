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

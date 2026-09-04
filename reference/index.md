# Package index

## Design and screening

Outcome-free panel geometry and one-row prevalence screens.

- [`design_summary()`](https://profsms.github.io/panelcert/reference/design_summary.md)
  : Design-summary primitive (spec section 2.1)
- [`applicable()`](https://profsms.github.io/panelcert/reference/applicable.md)
  : Pre-flight applicability check for Paper A
- [`adequacy_row()`](https://profsms.github.io/panelcert/reference/adequacy_row.md)
  : One-row panel adequacy screen
- [`print(`*`<DesignSummary>`*`)`](https://profsms.github.io/panelcert/reference/print.DesignSummary.md)
  : Print a design summary
- [`twoway_demean()`](https://profsms.github.io/panelcert/reference/twoway_demean.md)
  : Two-way within transformation by alternating projections
  (unbalanced-safe)
- [`multiway_demean()`](https://profsms.github.io/panelcert/reference/multiway_demean.md)
  : Multiway within transformation by alternating projections
- [`fe_dimension()`](https://profsms.github.io/panelcert/reference/fe_dimension.md)
  : Fixed-effect dimension of a two-way design

## Leverage and score concentration

Diffuse-regime leverage and realized-score diagnostics.

- [`leverage_report()`](https://profsms.github.io/panelcert/reference/leverage_report.md)
  : Module A diagnostic: variance-estimator adequacy under FE saturation
- [`fe_leverage()`](https://profsms.github.io/panelcert/reference/fe_leverage.md)
  : Diagonal of the two-way fixed-effect projection
- [`score_concentration()`](https://profsms.github.io/panelcert/reference/score_concentration.md)
  : Realized score concentration
- [`print(`*`<AdequacyReport>`*`)`](https://profsms.github.io/panelcert/reference/print.AdequacyReport.md)
  : Print an adequacy report

## Concentrated identifying variation

Paper A cycle packing, capture, and exact sign-flip inference.

- [`cycle_report()`](https://profsms.github.io/panelcert/reference/cycle_report.md)
  : Paper A adequacy report
- [`show_notes()`](https://profsms.github.io/panelcert/reference/show_notes.md)
  : Display detailed diagnostic notes
- [`cycle_capture()`](https://profsms.github.io/panelcert/reference/cycle_capture.md)
  : Cycle capture ratio
- [`cycle_contrasts()`](https://profsms.github.io/panelcert/reference/cycle_contrasts.md)
  : Nuisance-annihilating contrast system for a two-way design
- [`contrast_system()`](https://profsms.github.io/panelcert/reference/contrast_system.md)
  : Validated design-only annihilating contrast system
- [`support_compatibility()`](https://profsms.github.io/panelcert/reference/support_compatibility.md)
  : Check dependence-block compatibility of contrast supports
- [`signflip_test()`](https://profsms.github.io/panelcert/reference/signflip_test.md)
  : Exact sign-flip randomization test
- [`signflip_interval()`](https://profsms.github.io/panelcert/reference/signflip_interval.md)
  : Exact confidence set by test inversion

## Measurement error

Reliability, exact-normal thresholds, clustering, and certificates.

- [`eiv_adequacy()`](https://profsms.github.io/panelcert/reference/eiv_adequacy.md)
  : Measurement-error adequacy diagnostic
- [`eiv_adequacy_summary()`](https://profsms.github.io/panelcert/reference/eiv_adequacy_summary.md)
  : Module B diagnostic from regression summary output
- [`reliability_from_interval()`](https://profsms.github.io/panelcert/reference/reliability_from_interval.md)
  : Measurement-error SDs from published credible-interval bounds
- [`reliability_from_ratio()`](https://profsms.github.io/panelcert/reference/reliability_from_ratio.md)
  : Measurement-error SD implied by an external reliability ratio
- [`reliability_from_repeats()`](https://profsms.github.io/panelcert/reference/reliability_from_repeats.md)
  : Reliability from two measurements of the same regressor
- [`breakdown_reliability()`](https://profsms.github.io/panelcert/reference/breakdown_reliability.md)
  : Self-consistent breakdown reliability
- [`certified_breakdown_reliability()`](https://profsms.github.io/panelcert/reference/certified_breakdown_reliability.md)
  : Certified breakdown reliability
- [`tau2_crit()`](https://profsms.github.io/panelcert/reference/tau2_crit.md)
  : Stock-Yogo critical value for the residual treatment variance
- [`eta_finite_n()`](https://profsms.github.io/panelcert/reference/eta_finite_n.md)
  : Finite-n non-centrality mapping
- [`cluster_diagnostics()`](https://profsms.github.io/panelcert/reference/cluster_diagnostics.md)
  : Checkable conditions behind the cluster-robust layer
- [`projection_compatibility()`](https://profsms.github.io/panelcert/reference/projection_compatibility.md)
  : Projection compatibility, evaluated directly

## TWFE heterogeneity

Staggered-adoption design weights and TWFE adequacy diagnostics.

- [`twfe_design()`](https://profsms.github.io/panelcert/reference/twfe_design.md)
  : Pre-outcome TWFE design-statistic ladder
- [`twfe_gammas()`](https://profsms.github.io/panelcert/reference/twfe_gammas.md)
  : TWFE design-statistic ladder
- [`twfe_adequacy()`](https://profsms.github.io/panelcert/reference/twfe_adequacy.md)
  : TWFE-heterogeneity adequacy (covariance-aware, wild-bootstrap)

## Bundled applications

- [`vdem`](https://profsms.github.io/panelcert/reference/vdem.md) :
  V-Dem democracy–growth measurement-error application
- [`psid`](https://profsms.github.io/panelcert/reference/psid.md) :
  Cornwell–Rupert PSID earnings panel
- [`twins`](https://profsms.github.io/panelcert/reference/twins.md) :
  Ashenfelter–Krueger repeated-report twins extract
- [`castle`](https://profsms.github.io/panelcert/reference/castle.md) :
  Castle-doctrine adoption panel (certified application)
- [`divorce`](https://profsms.github.io/panelcert/reference/divorce.md)
  : No-fault-divorce adoption panel (flagged application)
- [`fscore`](https://profsms.github.io/panelcert/reference/fscore.md) :
  Piotroski F-Score / Visegrad firm panel (Paper A application)
- [`grunfeld`](https://profsms.github.io/panelcert/reference/grunfeld.md)
  : Canonical 11-firm Grunfeld investment panel
- [`minimum_wage`](https://profsms.github.io/panelcert/reference/minimum_wage.md)
  : Minimum-wage county panel (headline TWFE application)

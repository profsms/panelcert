
# panelcert

panelcert screens fixed-effect panel designs for concentrated
identifying variation, leverage, measurement error, and TWFE
heterogeneity, and supplies exact contrast inference where Gaussian
approximations are fragile.

[![R-CMD-check](https://github.com/profsms/panelcert/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/profsms/panelcert/actions/workflows/R-CMD-check.yaml)
[![Cross-language
parity](https://github.com/profsms/panelcert/actions/workflows/parity.yaml/badge.svg)](https://github.com/profsms/panelcert/actions/workflows/parity.yaml)
[![Codecov](https://codecov.io/gh/profsms/panelcert/branch/main/graph/badge.svg)](https://app.codecov.io/gh/profsms/panelcert)

## Install

Install the current release from GitHub:

``` r
install.packages("remotes")
remotes::install_github("profsms/panelcert@v0.5.1")
```

The r-universe build will be installable with
`install.packages("panelcert", repos = c("https://profsms.r-universe.dev", "https://cloud.r-project.org"))`.
CRAN installation will be documented after acceptance.

## Quickstart

``` r
library(panelcert); data(grunfeld)
applicable(grunfeld$capital, grunfeld$firm, grunfeld$year, controls = grunfeld$value)
#> $ok
#> [1] TRUE
#>
#> $reason
#> [1] "not a binary treatment; the binary granularity floor does not apply."
report <- cycle_report(grunfeld$invest, grunfeld$capital, grunfeld$firm, grunfeld$year, controls = grunfeld$value, interval = FALSE)
c(report$verdict, lambda_score = report$statistic$lambda_score, kappa = report$statistic$kappa)
#>                            lambda_score               kappa
#>           "FLAGGED" "0.738792811202239" "0.626959467924497"
```

## Entry Points

- `leverage_report(y, x, unit, time)` is the diffuse-regime diagnostic.
  It checks treatment concentration and leverage balance and compares
  df-corrected and HC0-HC3 inference.
- `cycle_report(y, x, unit, time)` is the concentrated-regime workflow
  from Paper A. It reports capture, granularity, an exact sign-flip
  test, and an exact confidence set.
- `eiv_adequacy(y, x, unit, time, ...)` screens continuous-regressor
  measurement error under Paper B’s exact-normal mapping and
  conservative certificate.
- `twfe_adequacy(y, unit, time, first_treat)` screens staggered-DiD/TWFE
  designs for heterogeneous-effect exposure.

Default `cycle_report()` printing is concise. Detailed caveats remain
available in `report$notes`, `show_notes(report)`, or
`print(report, notes = TRUE)`.

`adequacy_row(x, unit, time, y = NULL)` composes the design and capture
diagnostics into one flat record for prevalence screens.

The regression diagnostics accept `controls =` for numeric nuisance
covariates. Exact cycle packing supports one linearly independent
control and constructs locally annihilating 2-by-3 supports; it never
substitutes globally residualized outcomes for valid disjoint contrasts.

`applicable(x, unit, time)` is the outcome-free pre-flight check; it
reports structural binary-treatment granularity failures before packing.

## Packing Methods

- `"structured"` exploits four-cycle structure in dense panels and
  2-by-3 local projections when one continuous control is supplied.
- `"sparse"` is the scalable method for AKM-style mobility networks.
- `"greedy"` is a cheap lower bound that keeps the best of four
  deterministic DFS traversals.

## Verdicts

- `POINT_PASS` is a descriptive pass at the requested finite-design
  threshold; it is not a theorem-level certificate.
- `FLAGGED` marks a concentration or adequacy failure for which
  conventional inference is not licensed.
- `INCONCLUSIVE` is returned when the enumeration floor
  `2^(1-effective_C)` exceeds `alpha`. Its `reason` distinguishes the
  coding-invariant binary-treatment floor based on
  `min(n_treated, n_untreated)`, which repacking cannot fix, from a
  selected packing that may be improved.

Measurement-error reports can additionally return `CERTIFIED` when the
formal upper-bound pilot passes.

## Design-Only Use

`design_summary`, `cycle_capture`, and `adequacy_row(..., y = NULL)`
require only `(x, D, controls)`, not an outcome. In particular, `rho`,
`V_n`, `lambda_n`, `n_eff`, and `kappa_C` are available before
estimation.

## Julia Package And Data

The sibling Julia package is
[PanelAdequacy.jl](https://github.com/profsms/PanelAdequacy.jl). Fixed
CSV designs and CI enforce agreement to `1e-12` on numerical outputs and
exact agreement on integer fields, verdicts, and packed supports.

Dataset access follows language conventions: R exposes bundled panels
through `LazyData`; Julia exposes `datasets()`, `datapath()`, and
`load_dataset()`. Redistributable paper panels, including the
public-domain canonical 11-firm Grunfeld showcase, are bundled for
offline replication. `inst/DATA_SOURCES.md` records pinned provenance
and explains why the KSS test extract uses a checksum-pinned direct
download instead.

## Citation

Please cite both working papers when the corresponding diagnostics are
used:

- Halkiewicz, Stanislaw M. S. (2026). *Exact Inference in Fixed-Effect
  Regressions with Concentrated Identifying Variation*.
- Halkiewicz, Stanislaw M. S. (2026). *Fixed-Effect Saturation Is Not
  Weak Identification: Certifying Inference under Measurement Error*.

Run `citation("panelcert")` for machine-readable package metadata.

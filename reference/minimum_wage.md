# Minimum-wage county panel (headline TWFE application)

Exact balanced county–year panel for the unconditional specification in
Callaway and Sant'Anna (2021). The current diagnostic uses never-treated
counties as the comparison group. All treated-cell TWFE weights are
positive, yet the additive cohort-plus-event class fails uniform
inference certification.

## Usage

``` r
minimum_wage
```

## Format

A data frame with 15,988 county-year rows and 4 variables:

- uid:

  county id (unit)

- tid:

  year code, 1–7 (time)

- ft:

  first-treatment period; `NA` for never-treated counties

- y:

  log teen employment (outcome)

## Source

Callaway and Sant'Anna (2021) public minimum-wage replication data;
exact published-sample extract used in the TWFE-heterogeneity article.

## Examples

``` r
twfe_design(minimum_wage$uid, minimum_wage$tid, minimum_wage$ft)
#> Panel Adequacy Report — TWFE Heterogeneity
#> Design: n=15988, N=2284, T=7, d_K=2290, rho=0.1432
#> Design statistic Gamma = 0.270   negative-weight share = 0.0%
#>   restricted ladder: Gamma_gt = 0.270 | Gamma_c+e = 0.257 | Gamma_evt = 0.251 | Gamma_coh = 0.166
#> Breakdown threshold = 2.419
#> VERDICT: INCONCLUSIVE
#> Note: pre-outcome: the saturated group-time class requires c/sigma <= 2.419 (= eta-dagger/Gamma_gt); supply the outcome (twfe_adequacy) to calibrate c/sigma and obtain a boundary-robust certificate
# \donttest{
twfe_adequacy(minimum_wage$y, minimum_wage$uid, minimum_wage$tid,
              minimum_wage$ft, controls = "never", bootstrap = 19)
#> Panel Adequacy Report — TWFE Heterogeneity
#> Design: n=15988, N=2284, T=7, d_K=2290, rho=0.1432
#> Design statistic Gamma = 0.270   negative-weight share = 0.0%
#>   restricted ladder: Gamma_gt = 0.270 | Gamma_c+e = 0.257 | Gamma_evt = 0.251 | Gamma_coh = 0.166
#> Cluster normalization (direct; psi_hat = 1.385): Gamma_gt,CR = 0.229 | Gamma_c+e,CR = 0.219
#> TWFE beta_hat = -0.03661   sigma = 0.1432
#> Target-matched robust ATT = -0.05158
#> Trace-debiased point pilots c_S/sigma: gt = 6.46 | c+e = 6.47 | evt = 6.62 | coh = 4.48
#> Point worst-case envelopes: group-time = 31.6% | c+e = 29.3% | cohort 9.7% | event 29.2%
#> Selected heterogeneity class: group_time (each bound is a separate one-sided 95.0% statement)
#>   cohort: lower 0.276; upper 1.110; INCONCLUSIVE
#>   event-time: lower 0.963; upper 1.994; FLAGGED
#>   additive: lower 0.825; upper 3.622; FLAGGED
#>   group-time: lower 0.855; upper 3.961; FLAGGED
#>   descriptive point-envelope bootstrap (B=19): median 32.3%, 95% [22.3, 47.4]; psi in [1.36, 1.40]
#> Directional plug-in: eta = +1.370 (alignment +0.925), size 27.8%
#>   wild directional size: median 23.1%, 95% [12.5, 38.7]
#> Sign-reversal RMS threshold = 0.1358
#> Worst-case |eta| envelope = 1.480   Threshold (delta=0.05) = 0.652
#> Breakdown threshold = 2.847
#> Worst-case size envelope for nominal 5% test: 31.6%
#> VERDICT: FLAGGED at delta=0.05 (uniform certificate withheld)
#> Note: point pilot: covariance trace computed from the fitted Rademacher score covariance; 19 multiplier draws calibrate the 95.0% HC2 lower norm radius; the upper procedure uses projected-Wald noncentrality inversion
#> Note: normalization = direct: psi_hat = 1.385; direct CR1 psi = 1.385; AR(1) psi = 1.091 (rho = 0.266)
#> Note: saturated group-time lower bound 0.855 exceeds eta-dagger 0.652; uniform certification is withheld, not the published inference invalidated
# }
```

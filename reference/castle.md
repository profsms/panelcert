# Castle-doctrine adoption panel (additive-class certificate)

State–year panel for the Cheng–Hoekstra castle-doctrine design: a large
never-treated reservoir, no negative weights, design statistic Gamma =
0.21. The homicide application certifies in the prespecified additive
class. The saturated group-time class is structurally inconclusive
because its projected score covariance has rank 18 rather than 19.

## Usage

``` r
castle
```

## Format

A data frame with 550 state-year rows and 4 variables:

- uid:

  state id (unit)

- tid:

  year code (time)

- ft:

  first-treatment period; `NA` for never-treated states

- y:

  log homicide rate (outcome)

## Source

Cheng and Hoekstra castle-doctrine replication data; analysis panel
derived as in the TWFE-heterogeneity audit.

## Examples

``` r
twfe_design(castle$uid, castle$tid, castle$ft)
#> Panel Adequacy Report — TWFE Heterogeneity
#> Design: n=550, N=50, T=11, d_K=60, rho=0.1091
#> Design statistic Gamma = 0.211   negative-weight share = 0.0%
#>   restricted ladder: Gamma_gt = 0.211 | Gamma_c+e = 0.198 | Gamma_evt = 0.168 | Gamma_coh = 0.142
#> Breakdown threshold = 3.086
#> VERDICT: INCONCLUSIVE
#> Note: pre-outcome: the saturated group-time class requires c/sigma <= 3.086 (= eta-dagger/Gamma_gt); supply the outcome (twfe_adequacy) to calibrate c/sigma and obtain a boundary-robust certificate
twfe_adequacy(castle$y, castle$uid, castle$tid, castle$ft,
              heterogeneity_class = "additive")
#> Panel Adequacy Report — TWFE Heterogeneity
#> Design: n=550, N=50, T=11, d_K=60, rho=0.1091
#> Design statistic Gamma = 0.211   negative-weight share = 0.0%
#>   restricted ladder: Gamma_gt = 0.211 | Gamma_c+e = 0.198 | Gamma_evt = 0.168 | Gamma_coh = 0.142
#> Cluster normalization (direct; psi_hat = 3.372): Gamma_gt,CR = 0.115 | Gamma_c+e,CR = 0.108
#> TWFE beta_hat = 0.08181   sigma = 0.187
#> Target-matched robust ATT = 0.1094
#> Trace-debiased point pilots c_S/sigma: gt = 0.748 | c+e = 0 | evt = 0 | coh = 0
#> Point worst-case envelopes: group-time = 5.1% | c+e = 5.0% | cohort 5.0% | event 5.0%
#> Selected heterogeneity class: additive (each bound is a separate one-sided 95.0% statement)
#>   cohort: lower 0.000; upper 0.000; CERTIFIED
#>   event-time: lower 0.000; upper 0.083; CERTIFIED
#>   additive: lower 0.000; upper 0.236; CERTIFIED
#>   group-time: lower 0.000; upper unavailable (covariance rank 18/19); INCONCLUSIVE
#>   descriptive point-envelope bootstrap (B=999): median 5.1%, 95% [5.0, 5.9]; psi in [3.09, 3.51]
#> Directional plug-in: eta = -0.013 (alignment -0.164), size 5.0%
#>   wild directional size: median 5.0%, 95% [5.0, 5.4]
#> Sign-reversal RMS threshold = 0.387
#> Worst-case |eta| envelope = 0.000   Threshold (delta=0.05) = 0.652
#> Breakdown threshold = 6.045
#> Worst-case size envelope for nominal 5% test: 5.0%
#> VERDICT: FORMALLY CERTIFIED at (alpha, delta, gamma) = (0.05, 0.05, 0.05)
#> Note: point pilot: covariance trace computed from the fitted Rademacher score covariance; 999 multiplier draws calibrate the 95.0% HC2 lower norm radius; the upper procedure uses projected-Wald noncentrality inversion
#> Note: normalization = direct: psi_hat = 3.372; direct CR1 psi = 3.372; AR(1) psi = 1.336 (rho = 0.226)
#> Note: additive upper bound 0.236 is below eta-dagger 0.652
```

# Brazilian property-tax panel (large-panel TWFE application)

Balanced 2004–2015 municipal panel for the Christensen–Garfias property-
tax specification. Municipalities treated in 2004 are excluded so every
treated group-time effect has an in-sample untreated baseline. This is
the balanced-subset reanalysis in Paper C, rather than the full
published estimation sample.

## Usage

``` r
brazil
```

## Format

A data frame with 34,080 municipality-year rows and 4 variables:

- uid:

  consecutive municipality id

- tid:

  period code, 1–12 for 2004–2015

- ft:

  first-treatment period; `NA` for never-treated municipalities

- y:

  log property-tax revenue (`logiptu`)

## Source

Christensen and Garfias (2021) replication data distributed in the
public Chiu, Lan, Liu, and Xu causal-panel reanalysis archive,
doi:10.7910/DVN/9RJFZF.

## Examples

``` r
twfe_design(brazil$uid, brazil$tid, brazil$ft)
#> Panel Adequacy Report — TWFE Heterogeneity
#> Design: n=34080, N=2840, T=12, d_K=2851, rho=0.0837
#> Design statistic Gamma = 1.138   negative-weight share = 17.3%
#>   restricted ladder: Gamma_gt = 1.138 | Gamma_c+e = 1.065 | Gamma_evt = 1.025 | Gamma_coh = 0.463
#> Breakdown threshold = 0.573
#> VERDICT: INCONCLUSIVE
#> Note: pre-outcome: the saturated group-time class requires c/sigma <= 0.573 (= eta-dagger/Gamma_gt); supply the outcome (twfe_adequacy) to calibrate c/sigma and obtain a boundary-robust certificate
# \donttest{
twfe_adequacy(brazil$y, brazil$uid, brazil$tid, brazil$ft,
              bootstrap = 19)
#> Panel Adequacy Report — TWFE Heterogeneity
#> Design: n=34080, N=2840, T=12, d_K=2851, rho=0.0837
#> Design statistic Gamma = 1.138   negative-weight share = 17.3%
#>   restricted ladder: Gamma_gt = 1.138 | Gamma_c+e = 1.065 | Gamma_evt = 1.025 | Gamma_coh = 0.463
#> Cluster normalization (direct; psi_hat = 1.463): Gamma_gt,CR = 0.941 | Gamma_c+e,CR = 0.880
#> TWFE beta_hat = 0.05935   sigma = 1.106
#> Target-matched robust ATT = 0.02846
#> Trace-debiased point pilots c_S/sigma: gt = 10.3 | c+e = 9.79 | evt = 7.1 | coh = 7.98
#> Point worst-case envelopes: group-time = 100.0% | c+e = 100.0% | cohort 86.3% | event 100.0%
#> Selected heterogeneity class: group_time (each bound is a separate one-sided 95.0% statement)
#>   cohort: lower 0.684; upper 6.100; FLAGGED
#>   event-time: lower 2.704; upper 9.610; FLAGGED
#>   additive: lower 3.028; upper 16.042; FLAGGED
#>   group-time: lower 3.245; upper 26.224; FLAGGED
#>   descriptive point-envelope bootstrap (B=19): median 98.5%, 95% [26.1, 100.0]; psi in [1.44, 1.46]
#> Directional plug-in: eta = +5.145 (alignment +0.518), size 99.9%
#>   wild directional size: median 49.8%, 95% [6.9, 97.6]
#> Sign-reversal RMS threshold = 0.05215
#> Worst-case |eta| envelope = 9.688   Threshold (delta=0.05) = 0.652
#> Breakdown threshold = 0.693
#> Worst-case size envelope for nominal 5% test: 100.0%
#> VERDICT: FLAGGED at delta=0.05 (uniform certificate withheld)
#> Note: point pilot: covariance trace computed from the fitted Rademacher score covariance; 19 multiplier draws calibrate the 95.0% HC2 lower norm radius; the upper procedure uses projected-Wald noncentrality inversion
#> Note: normalization = direct: psi_hat = 1.463; direct CR1 psi = 1.463; AR(1) psi = 1.164 (rho = 0.181)
#> Note: saturated group-time lower bound 3.245 exceeds eta-dagger 0.652; uniform certification is withheld, not the published inference invalidated
# }
```

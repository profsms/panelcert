# No-fault-divorce adoption panel (flagged application)

State–year panel for the Stevenson–Wolfers no-fault-divorce design
(Goodman-Bacon's pathology example): near-universal eventual adoption.
After the two always-treated states are dropped, the analysis design has
a 1.1% negative-weight share and Gamma = 0.64. Its direct-CR1 point
envelope is 41.5%, but bootstrap uncertainty is wide and the fixed-\\T\\
warning applies.

## Usage

``` r
divorce
```

## Format

A data frame with 1377 state-year rows and 4 variables:

- uid:

  state id (unit)

- tid:

  year code (time)

- ft:

  first-treatment period; `NA` never-treated, `0` always-treated

- y:

  female suicide rate per 100k (outcome)

## Source

Stevenson and Wolfers divorce data (via the bacondecomp distribution);
analysis panel derived as in the TWFE-heterogeneity audit.

## Examples

``` r
# \donttest{
twfe_adequacy(divorce$y, divorce$uid, divorce$tid, divorce$ft)
#> Panel Adequacy Report — TWFE Heterogeneity
#> Design: n=1323, N=49, T=27, d_K=75, rho=0.0567
#> Design statistic Gamma = 0.643   negative-weight share = 1.1%
#>   restricted ladder: Gamma_gt = 0.643 | Gamma_c+e = 0.562 | Gamma_evt = 0.468 | Gamma_coh = 0.381
#> Cluster normalization (direct; psi_hat = 2.917): Gamma_gt,CR = 0.377 | Gamma_c+e,CR = 0.329
#> TWFE beta_hat = -0.01578   sigma = 0.1961
#> Target-matched robust ATT = -0.07176
#> Trace-debiased point pilots c_S/sigma: gt = 6.6 | c+e = 5.3 | evt = 2.95 | coh = 5.08
#> Point worst-case envelopes: group-time = 70.0% | c+e = 41.5% | cohort 20.5% | event 12.8%
#> Boundary-robust group-time K interval (95.0%, HC2): [0.120, 7.523]; HC3 [0.000, 7.961]
#>   descriptive point-envelope bootstrap (B=999): median 37.8%, 95% [5.0, 86.0]; psi in [2.32, 3.29]
#> Directional plug-in: eta = +1.629 (alignment +0.568), size 37.1%
#>   wild directional size: median 16.3%, 95% [5.0, 48.2]
#> Sign-reversal RMS threshold = 0.02454
#> Worst-case |eta| envelope = 2.485   Threshold (delta=0.05) = 0.652
#> Breakdown threshold = 1.732
#> Worst-case size envelope for nominal 5% test: 70.0%
#> VERDICT: INCONCLUSIVE
#> Note: 2 always-treated unit(s) dropped (setup g >= 2)
#> Note: point pilot: covariance trace computed from the fitted Rademacher score covariance; 999 multiplier draws calibrate the 95.0% HC2 projected-vector radius
#> Note: normalization = direct: psi_hat = 2.917; direct CR1 psi = 2.917; AR(1) psi = 1.620 (rho = 0.296)
#> Note: fixed-T, many-cluster approximation is strained (G = 49, T = 27)
#> Note: saturated group-time interval [0.120, 7.523] crosses eta-dagger 0.652
# }
```

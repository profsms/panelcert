# V-Dem democracy–growth measurement-error application

Country–year panel used in \*Breakdown Reliability for Saturated
Fixed-Effect Inference\*: log GDP per capita against continuous V-Dem
institutional indices, each accompanied by the measurement model's
posterior standard deviation (so the measurement-error input is data,
not a calibration). The two-pole contrast the paper reports – polyarchy
certified, the constraint sub-indices flagged – is reproduced directly
from this object.

## Usage

``` r
vdem
```

## Format

A data frame with 8931 country-year rows and 9 variables:

- iso:

  country ISO code (unit id)

- year:

  calendar year (time id)

- ly:

  log GDP per capita (Maddison Project 2020)

- v2x_polyarchy, v2x_polyarchy_sd:

  electoral-democracy index and its posterior SD

- v2xlg_legcon, v2xlg_legcon_sd:

  legislative-constraints index and its posterior SD

- v2x_jucon, v2x_jucon_sd:

  judicial-constraints index and its posterior SD

## Source

V-Dem dataset at commit \`f4dd26922e658442524dfd954bf14f7ebe622d5d\`
(measurement-model posterior SDs); Maddison Project Database 2020. The
raw V-Dem artifact is pinned by SHA-256 in \`inst/DATA_SOURCES.md\`.

## Examples

``` r
d <- vdem[stats::complete.cases(vdem$ly, vdem$v2x_polyarchy, vdem$v2x_polyarchy_sd), ]
eiv_adequacy(d$ly, d$v2x_polyarchy, d$iso, d$year, sigma_nu = d$v2x_polyarchy_sd)
#> Panel Adequacy Report — Measurement Error
#> Design: n=8930, N=163, T=59, d_K=221, rho=0.0247
#> Within reliability lambda_hat = 0.894   ((1-lambda)/lambda = 0.119)
#> Pilot: beta* = 0.06096 -> corrected beta0 = 0.06821 (se 0.0326)
#> Certified breakdown = 0.851   reliability lower bound = 0.894   |eta| upper bound = 0.445
#> Non-centrality |eta| = 0.249   Threshold (delta=0.05) = 0.652
#> Breakdown threshold = 0.762
#> Implied size of nominal 5% test: 5.7%
#> VERDICT: FORMALLY CERTIFIED at (alpha, delta, gamma) = (0.05, 0.05, 0.05)
#> Note: formal certificate (prop-certificate): certified breakdown lambda_dagger_gamma = 0.851 is obtained by replacing |t| with |t| + z_(1-gamma_beta); the comparison uses reliability lower bound ell = 0.894. False certification is at most gamma_beta + gamma_lambda = 0.05 + 0 = 0.05, without requiring independence. Implied size shown is at the point pilot.
#> Note: CONDITIONAL RELIABILITY TREATMENT: no reliability_lower was supplied, so lambda_hat is treated as known/consistent and gamma_lambda = 0. A noisy finite-sample reliability estimate requires a lower confidence bound and its coverage-error budget.
#> Note: corrected pilot beta*/lambda_hat is not a consistent point estimate under weak information; reported with its sampling band (Corollary cor-slope)
#> Note: power tax: local power slope attenuated by sqrt(lambda) = 0.95 (Proposition prop-power)
```

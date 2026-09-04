# Cornwell–Rupert PSID earnings panel

Person–year panel of weeks worked and the log wage, the canonical noisy
regressor of the validation literature. Retained as a general panel
example and for backward compatibility; it is not an application in the
current measurement-error article.

## Usage

``` r
psid
```

## Format

A data frame with 4165 person-year rows and 4 variables:

- id:

  person id (unit)

- year:

  year (time)

- lwage:

  log wage (the mismeasured regressor)

- wks:

  weeks worked (outcome)

## Source

Cornwell–Rupert PSID extract, distributed with the plm package.

## Examples

``` r
eiv_adequacy(psid$wks, psid$lwage, psid$id, psid$year,
             reliability = 0.65, pilot = "point")
#> Panel Adequacy Report — Measurement Error
#> Design: n=4165, N=595, T=7, d_K=601, rho=0.1443
#> Within reliability lambda_hat = 0.650   ((1-lambda)/lambda = 0.538)
#> Pilot: beta* = 0.7332 -> corrected beta0 = 1.128 (se 0.716)
#> Non-centrality |eta| = 0.848   Threshold (delta=0.05) = 0.652
#> Breakdown threshold = 0.707
#> Implied size of nominal 5% test: 13.6%
#> VERDICT: FLAGGED at delta=0.05
#> Note: POINT PASS, not a certificate (rem-plugin): the corrected pilot at its point estimate is descriptive. Under weak information eta_hat converges to a nondegenerate random multiple (|B|/|beta0|)|eta| of the target -- median close to it, but no concentration. Its sampling variability is a first-order feature of the regime, not a vanishing approximation error. For a size-controlled statement use pilot = "conservative".
#> Note: corrected pilot beta*/lambda_hat is not a consistent point estimate under weak information; reported with its sampling band (Corollary cor-slope)
#> Note: power tax: local power slope attenuated by sqrt(lambda) = 0.81 (Proposition prop-power)
```

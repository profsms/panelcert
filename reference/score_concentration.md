# Realized score concentration

Computes the largest residual-score variance share, its Herfindahl
index, and inverse-Herfindahl effective support size for a two-way
fixed-effect regression. This is a one-outcome warning diagnostic, not
by itself a consistent estimator under unrestricted heteroskedasticity.

## Usage

``` r
score_concentration(y, x, unit, time, controls = NULL)
```

## Arguments

- y, x:

  outcome and regressor vectors.

- unit, time:

  fixed-effect identifiers.

- controls:

  optional numeric nuisance-covariate vector or matrix.

## Value

A list with \`lambda_score\`, \`H_score\`, and \`n_eff_score\`.

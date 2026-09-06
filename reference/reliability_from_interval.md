# Measurement-error SDs from published credible-interval bounds

Computes \`sigma_nu = (codehigh - codelow)/2\`. Use this helper only
when the interval half-width is substantively calibrated as one
measurement-error SD; do not apply it mechanically to
highest-posterior-density bounds.

## Usage

``` r
reliability_from_interval(codelow, codehigh)
```

## Arguments

- codelow, codehigh:

  interval bounds, one per observation

## Value

numeric vector of per-observation measurement-error SDs

## Examples

``` r
reliability_from_interval(c(0.20, 0.31), c(0.24, 0.35))
#> [1] 0.02 0.02
```

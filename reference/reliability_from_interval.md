# Measurement-error SDs from published credible-interval bounds

V-Dem convention: the interval brackets one posterior SD, so \`sigma_nu
= (codehigh - codelow)/2\` (lead application in the accompanying
article).

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

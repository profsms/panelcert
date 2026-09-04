# Design-summary primitive (spec section 2.1)

Pre-outcome design description reused by all four inference/diagnostic
modules and useful standalone: n, N, T, the fixed-effect dimension
\`d_K\` via union-find, \`rho = d_K/n\`, and (if a regressor is
supplied) its within residual variation \`tau_star2 = x' M x\`.

## Usage

``` r
design_summary(
  unit = NULL,
  time = NULL,
  x = NULL,
  fe_levels = NULL,
  controls = NULL
)
```

## Arguments

- unit, time:

  raw identifier vectors

- x:

  optional regressor

- fe_levels:

  optional non-empty list of fixed-effect identifier vectors; when
  supplied, the sparse dummy-matrix rank is used for \`d_K\`

- controls:

  optional numeric nuisance-covariate vector or matrix. These are
  partialled after the fixed effects; \`d_K\` and \`rho\` still describe
  the fixed-effect space alone.

## Value

object of class `DesignSummary`

## Examples

``` r
design_summary(rep(1:20, each = 10), rep(1:10, times = 20))
#> Panel Design Summary
#>   n = 200 obs | N = 20 units | T = 10 periods
#>   d_K = 29 (connected components: 1) | rho = d_K/n = 0.1450
```

# Cost of the estimand-aligned (ETI) analysis of a stepped-wedge design

Variance of the exposure-time-indicator (ETI) estimator of the target
(default TATE), the IT variance of the same design, and both relative to
the IT variance of the balanced complete wedge with the same number of
clusters, periods and covariance.

## Usage

``` r
sw_eti_cost(
  design,
  phi = NULL,
  icc = NULL,
  cac = 1,
  m = NULL,
  estimator = c("GLS", "TWFE"),
  sizes = NULL,
  target = NULL
)
```

## Arguments

- design:

  numeric vector of length `T`: cluster counts (or shares) for start
  periods `2, ..., T`, followed by the number wait-listed. The first
  entry must be positive.

- phi:

  Hussey–Hughes parameter \\\phi = m\rho/(1-\rho)\\ (GLS); `Inf` gives
  the fixed-effects limit.

- icc, cac, m:

  alternative nested-exchangeable specification: within-period ICC,
  cluster autocorrelation and individuals per cluster-period (see
  \[sw_phi_nested()\]). With `sizes`, `m` is not needed.

- estimator:

  `"GLS"` (default) or `"TWFE"` (two-way fixed effects, which does not
  use the covariance).

- sizes:

  optional `G x T` matrix of cluster-period sizes, one row per cluster
  in design order (sequence 2 first, wait-list last). Requires integer
  counts in `design` and, for GLS, `icc` and `cac`.

- target:

  ETI contrast weights (default TATE, `1/L`)

## Value

list with `var_eti`, `var_it`, `var_ref` (balanced complete wedge, IT,
equal sizes at the same covariance; with `sizes`, the mean size),
`eti_cost = var_eti / var_ref` and `it_ratio = var_it / var_ref`.
Variances are for the stated number of clusters on the covariance scale
of the model (GLS: \\V = I + \phi J\\; TWFE: i.i.d. unit errors).

## Examples

``` r
sw_eti_cost(c(6, 6, 6, 6, 0), phi = 2)$eti_cost    # 2.27
#> [1] 2.269942
sw_eti_cost(c(6, 6, 1, 0, 11), phi = 2)$eti_cost   # 1.23
#> [1] 1.227424
```

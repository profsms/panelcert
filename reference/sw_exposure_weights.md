# Exposure weights of the immediate-treatment estimator in a stepped-wedge design

Computes the implicit exposure-time weights \\W_e\\, \\e = 1, \ldots,
T-1\\, of the immediate-treatment (IT) estimator, so that \\E(\hat\beta)
= \sum_e W_e \theta_e\\ under exposure-time heterogeneity. The weights
sum to one.

## Usage

``` r
sw_exposure_weights(
  design,
  phi = NULL,
  icc = NULL,
  cac = 1,
  m = NULL,
  estimator = c("GLS", "TWFE"),
  sizes = NULL
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

## Value

named numeric vector `W1, ..., W(T-1)`

## Examples

``` r
sw_exposure_weights(c(6, 6, 6, 6, 0), phi = 2)          # complete wedge
#>          W1          W2          W3          W4 
#>  0.94285714  0.30000000 -0.07142857 -0.17142857 
sw_exposure_weights(c(4, 4, 4, 4, 8), phi = 2)          # 8 wait-listed
#>          W1          W2          W3          W4 
#>  0.67058824  0.30000000  0.06470588 -0.03529412 
sw_exposure_weights(c(6, 6, 6, 6, 0), estimator = "TWFE")
#>   W1   W2   W3   W4 
#>  1.0  0.3 -0.1 -0.2 
```

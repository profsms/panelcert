# Exhaustive search over stepped-wedge and wait-list allocations

Enumerates every integer allocation of `G` clusters to start periods
`2, ..., T` and a wait-list (at least one cluster starting in period 2),
subject to a rollout cap and a maximum wait-list size, and evaluates the
ETI cost, IT variance ratio and adequacy radius of each (equal
cluster-period sizes). The optimum is over this finite set only.

## Usage

``` r
sw_design_search(
  G,
  T,
  cap = Inf,
  max_wait = G - 1,
  phi = NULL,
  icc = NULL,
  cac = 1,
  m = NULL,
  estimator = c("GLS", "TWFE"),
  objective = c("eti_cost", "radius"),
  max_it_ratio = Inf,
  top = 1,
  alpha = 0.05,
  epsilon = 0.05,
  power = 0.8
)
```

## Arguments

- G:

  number of clusters

- T:

  number of periods

- cap:

  maximum number of clusters starting in any one period

- max_wait:

  maximum number of wait-listed clusters

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

- objective:

  `"eti_cost"` (minimise the ETI cost) or `"radius"` (maximise
  \\\kappa^\*\\ subject to `it_ratio <= max_it_ratio`)

- max_it_ratio:

  IT-variance constraint for `objective = "radius"`

- top:

  number of best designs returned

- alpha, epsilon, power:

  radius calibration (see \[sw_radius()\])

## Value

list with `best` (data frame of the `top` designs), `frontier` (best
objective for each wait-list size bound `0, ..., max_wait`) and
`designs` (all evaluated designs)

## Examples

``` r
# \donttest{
s <- sw_design_search(G = 24, T = 5, cap = 6, max_wait = 12, phi = 2)
s$frontier[s$frontier$max_wait %in% c(0, 4, 6, 8, 12), ]
#>    max_wait start2 start3 start4 start5 wait eti_cost  it_ratio    radius
#> 1         0      6      6      6      6    0 2.269942 1.0000000 0.1332467
#> 5         4      6      6      2      6    4 1.519707 0.9237537 0.1819003
#> 7         6      6      6      0      6    6 1.376921 0.9459459 0.2129412
#> 9         8      6      6      1      3    8 1.290928 0.9563567 0.2491485
#> 13       12      6      6      1      0   11 1.227424 1.0508757 0.3768735
# }
```

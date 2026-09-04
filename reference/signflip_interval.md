# Exact confidence set by test inversion

The contrast scores are affine in \`beta0\`, so the whole randomization
orbit is affine too and the p-value curve over an entire grid costs no
more than a single test.

## Usage

``` r
signflip_interval(
  y,
  x,
  unit = NULL,
  time = NULL,
  alpha = 0.05,
  ngrid = 24001,
  nflips = 99999,
  span = 12,
  method = c("structured", "sparse", "greedy"),
  controls = NULL,
  blocks = NULL
)
```

## Arguments

- y:

  outcome vector.

- x:

  treatment vector.

- unit, time:

  fixed-effect identifiers.

- alpha:

  level.

- ngrid:

  grid resolution.

- nflips:

  number of random sign patterns.

- span:

  half-width of the grid in randomization-scale units.

- method:

  \`"structured"\` (all unit-pair four-cycles, greedy by squared
  loading, then a DFS pass on leftovers; for dense panels), \`"sparse"\`
  (stayer digons then firm-pair four-cycles with the optimal extreme
  pairing; for MOBILITY NETWORKS, where \`"structured"\` is \`O(N^2
  T^2)\` and infeasible), or \`"greedy"\` (digons then DFS-extracted
  cycles, retaining the best of four fixed traversal orders; a
  deterministic lower bound on achievable capture).

- controls:

  optional numeric nuisance-covariate vector or matrix. Exact packing
  currently supports at most one linearly independent control. With one
  control, structured packing uses complete 2-by-3 rectangles and local
  cycle-space projection; other methods combine base cycles in pairs.

- blocks:

  optional design-defined dependence-block ids.

## Value

A list with \`lo\`, \`hi\`, \`contiguous\`, \`grid_truncated\`,
\`beta_tilde\`, \`C\`, \`effective_C\`, \`kappa\`, \`grid\`, and
\`pvalue\`. For the unstudentized statistic the acceptance set is an
interval; \`contiguous\` is a numerical audit. Infinite endpoints are
returned when the grid boundary is reached, rather than presenting a
search edge as a confidence limit.

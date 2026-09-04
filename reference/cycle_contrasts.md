# Nuisance-annihilating contrast system for a two-way design

In a two-way design the annihilating contrasts are exactly the cycle
space of the observation multigraph. Each contrast eliminates both sets
of fixed effects identically – exactly, not asymptotically – which is
what makes the sign-flip test of \[signflip_test()\] exact in finite
samples.

## Usage

``` r
cycle_contrasts(
  x,
  unit,
  time,
  method = c("structured", "sparse", "greedy"),
  controls = NULL
)
```

## Arguments

- x:

  treatment vector.

- unit, time:

  fixed-effect identifiers.

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

## Value

A list with \`rows\`, \`signs\`, \`loadings\`, \`V_n\`, \`kappa\`,
\`max_share\`.

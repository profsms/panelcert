# Cycle capture ratio

\`kappa_C = sum_c (v_c'x)^2 / V_n\`. In Paper A's diffuse comparison it
is the Pitman efficiency only when both power conditions (P1) and (P2)
hold. In the concentrated regime it is a capture diagnostic, and
\`1/sqrt(kappa)\` is only a capture-implied signal/SE ratio. It is
computable from the design before any outcome is examined.

## Usage

``` r
cycle_capture(
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

The capture ratio.

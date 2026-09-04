# Exact sign-flip randomization test

With \`blocks = NULL\`, exactness is stated for independent
observation-level symmetric errors. For clustered errors, pass one block
id per observation: every treatment-loaded support is then checked to be
a union of complete blocks. Conditional on that check, the joint score
law is invariant under independent sign flips; scores need not be
exchangeable or identically distributed. The identity sign pattern is
adjoined.

## Usage

``` r
signflip_test(
  y,
  x,
  unit = NULL,
  time = NULL,
  beta0 = 0,
  nflips = 99999,
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

- beta0:

  null value.

- nflips:

  number of random sign patterns.

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

A list with \`p\`, \`C\`, \`effective_C\`, \`kappa\`, \`beta_tilde\`,
and \`full_enumeration_floor = 2^(1-effective_C)\`. \`min_pvalue\` is
retained as an alias. Global sign reversal duplicates the default
absolute statistic.

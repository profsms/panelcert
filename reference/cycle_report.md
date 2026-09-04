# Paper A adequacy report

Reports the concentration of the identifying variation, the capture the
contrast system achieves, the price of exactness, and the exact
confidence set. The verdict answers: is the Gaussian approximation
underlying conventional inference trustworthy on this design?

## Usage

``` r
cycle_report(
  y,
  x,
  unit,
  time,
  alpha = 0.05,
  delta = 0.05,
  method = c("structured", "sparse", "greedy"),
  nflips = 99999,
  interval = TRUE,
  lambda_max = 0.1,
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

- delta:

  size tolerance (reported, not used in the verdict).

- method:

  \`"structured"\` (all unit-pair four-cycles, greedy by squared
  loading, then a DFS pass on leftovers; for dense panels), \`"sparse"\`
  (stayer digons then firm-pair four-cycles with the optimal extreme
  pairing; for MOBILITY NETWORKS, where \`"structured"\` is \`O(N^2
  T^2)\` and infeasible), or \`"greedy"\` (digons then DFS-extracted
  cycles, retaining the best of four fixed traversal orders; a
  deterministic lower bound on achievable capture).

- nflips:

  number of random sign patterns.

- interval:

  whether to compute the exact confidence set.

- lambda_max:

  concentration warning convention; Paper A's condition is the sequence
  statement \`lambda_n -\> 0\`, not a finite cutoff.

- controls:

  optional numeric nuisance-covariate vector or matrix. Exact packing
  currently supports at most one linearly independent control. With one
  control, structured packing uses complete 2-by-3 rectangles and local
  cycle-space projection; other methods combine base cycles in pairs.

- blocks:

  optional design-defined dependence-block ids.

## Value

An \`AdequacyReport\`.

## Details

With \`blocks = NULL\`, concentration is measured across observations.
With dependence blocks, treatment mass is aggregated as \`sum(xt_i^2)\`
within block and realized score concentration uses squared block scores
\`(sum(xt_i \* u_i))^2\`, as in Paper A's worker–firm application.

\`POINT_PASS\` is a descriptive pass below the finite \`lambda_max\`
heuristic, not a theorem-level certificate; \`FLAGGED\` marks a
concentration warning; and \`INCONCLUSIVE\` means the two-sided
full-enumeration floor \`2^(1-effective_C)\` exceeds \`alpha\`.

Detailed diagnostic caveats remain in \`report\$notes\` but are hidden
from the default display. Use \[show_notes()\] or \`print(report, notes
= TRUE)\` to show them; an \`INCONCLUSIVE\` reason remains visible in
the concise display.

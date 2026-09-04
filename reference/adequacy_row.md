# One-row panel adequacy screen

Composes design and cycle diagnostics into one flat named record. The
base record requires no outcome. Supplying \`y\` appends realized score
concentration and the Paper A verdict.

## Usage

``` r
adequacy_row(
  x,
  unit,
  time,
  y = NULL,
  method = c("structured", "sparse", "greedy"),
  alpha = 0.05,
  controls = NULL
)
```

## Arguments

- x:

  treatment or regressor vector.

- unit, time:

  fixed-effect identifiers.

- y:

  optional outcome vector.

- method:

  cycle packing method: \`"structured"\`, \`"sparse"\`, or \`"greedy"\`.

- alpha:

  nominal level used for the Lei–Bickel feasibility condition and, when
  \`y\` is supplied, the cycle-report verdict.

- controls:

  optional numeric nuisance-covariate vector or matrix.

## Value

A flat named list suitable for conversion to a one-row data frame.

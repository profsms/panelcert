# Validated design-only annihilating contrast system

Construct caller-supplied contrasts for any number of categorical fixed
effects. Supports must be disjoint and every contrast must annihilate
every supplied fixed effect. This covers Paper A's worker–firm–year
application, which is not a two-way cycle packing.

## Usage

``` r
contrast_system(
  x,
  fe_levels,
  rows,
  weights,
  blocks = NULL,
  tol = 1e-10,
  maxit = 10000L
)
```

## Arguments

- x:

  treatment vector.

- fe_levels:

  list of fixed-effect identifier vectors.

- rows:

  list of observation-index vectors, one per support.

- weights:

  list of contrast-weight vectors matching \`rows\`.

- blocks:

  optional dependence-block ids; incompatible supports are rejected.

- tol, maxit:

  numerical controls.

## Value

a \`CycleSystem\` object usable by \[signflip_test()\] and
\[signflip_interval()\].

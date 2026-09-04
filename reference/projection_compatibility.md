# Projection compatibility, evaluated directly

Assumption ass-cluster(iv) of the accompanying article asks that
\\\sum_g \\M a^{(g)} - a^{(g)}\\^2 = o_p(\tau^{\*2})\\, where
\\a^{(g)}\\ is the residualized regressor restricted to cluster \\g\\.
This function computes that ratio as it stands, for a two-way (unit and
time) fixed-effect design.

## Usage

``` r
projection_compatibility(
  xt,
  cluster,
  unit,
  time,
  tau_star2 = sum(xt^2),
  tol = 1e-10,
  maxit = 10000L,
  max_cells = 5e+06
)
```

## Arguments

- xt:

  residualized regressor (the same `M x*` used elsewhere).

- cluster:

  cluster identifier, one per observation.

- unit, time:

  raw fixed-effect identifiers, one per observation.

- tau_star2:

  the within variation `sum(xt^2)`; defaults to `sum(xt^2)`.

- tol, maxit:

  convergence controls for the alternating projections.

- max_cells:

  refuse to allocate more than this many matrix cells.

## Value

A list with \`ratio\` (the quantity above; \`NA\` if the guard tripped),
\`G\`, \`n\` and \`cells\`.

## Details

Unlike the `ratio_ne` shortcut in \[cluster_diagnostics()\], this needs
neither of the balance conditions (N1) and (N2) of Lemma lem-nest(b): it
is the quantity the assumption actually names. Use it when `ratio_ne`
and `max_energy` point in different directions, when the panel is far
from balanced, or whenever the cluster-robust column is load-bearing.

Under Lemma lem-nest(a) the ratio is exactly zero when every
fixed-effect cell sits inside one cluster, so a non-zero value is
entirely the work of the fixed effects that cut across clusters.

Cost is one alternating-projection sweep over an `n * G` matrix. The
`max_cells` guard returns `NA` rather than allocating a matrix larger
than that; raise it deliberately if you want the number anyway.

## Examples

``` r
n <- 200; unit <- rep(1:20, each = 10); time <- rep(1:10, times = 20)
x <- rnorm(n)
xt <- twoway_demean(x, unit, time)
# unit FE nest in unit clusters, the 10 time effects do not
projection_compatibility(xt, unit, unit, time)$ratio
#> [1] 0.05
```

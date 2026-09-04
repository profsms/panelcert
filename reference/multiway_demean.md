# Multiway within transformation by alternating projections

Scalable residualization for applications with more than two categorical
fixed-effect dimensions, including Paper A's worker–firm–year design.

## Usage

``` r
multiway_demean(x, fe_levels, tol = 1e-10, maxit = 10000L)
```

## Arguments

- x:

  numeric vector to residualize.

- fe_levels:

  list of fixed-effect identifier vectors.

- tol, maxit:

  convergence controls.

## Value

the residualized numeric vector.

# Diagonal of the two-way fixed-effect projection

\`(P_K)\_ii\` for each observation, computed exactly from the
(N+T)x(N+T) FE Gram matrix via pseudo-inverse (handles the rank
deficiency \`d_K = N + T - \#components\`).

## Usage

``` r
fe_leverage(unit, time)
```

## Arguments

- unit, time:

  raw identifier vectors

## Value

numeric vector of FE leverages, one per observation

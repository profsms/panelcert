# Worst-case estimand mismatch of the immediate-treatment estimator

\\\Gamma = \sqrt{L}\\\lVert W - v \rVert_2\\: the worst-case bias of the
IT estimator for the target \\v^\top\theta\\ per unit RMS exposure-time
heterogeneity (Proposition 1 of Halkiewicz and Alexeev).

## Usage

``` r
sw_gamma(W = NULL, design = NULL, target = NULL, ...)
```

## Arguments

- W:

  exposure weights (as from \[sw_exposure_weights()\]); alternatively
  supply `design` and the covariance arguments in `...`

- design:

  optional design, passed with `...` to \[sw_exposure_weights()\]

- target:

  target weights \\v\\ summing to one; default TATE (\\v = 1/L\\)

- ...:

  covariance and estimator arguments for \[sw_exposure_weights()\]

## Value

\\\Gamma\\

## Examples

``` r
sw_gamma(design = c(6, 6, 6, 6, 0), phi = 2)
#> [1] 1.747535
```

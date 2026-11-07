# Adequacy radius (heterogeneity tolerance) of a stepped-wedge design

\\\kappa^\* = \eta\_{\alpha,\epsilon} / \\(z\_{1-\alpha/2} +
z\_{power})\Gamma\\\\: the largest RMS exposure-time heterogeneity,
relative to the planned effect, at which the worst-case coverage of the
nominal \\1-\alpha\\ IT interval falls by at most \\\epsilon\\, for a
trial planned with the stated power against a constant effect.

## Usage

``` r
sw_radius(
  W = NULL,
  design = NULL,
  target = NULL,
  alpha = 0.05,
  epsilon = 0.05,
  power = 0.8,
  ...
)
```

## Arguments

- W:

  exposure weights (as from \[sw_exposure_weights()\]); alternatively
  supply `design` and the covariance arguments in `...`

- design:

  optional design, passed with `...` to \[sw_exposure_weights()\]

- target:

  target weights \\v\\ summing to one; default TATE (\\v = 1/L\\)

- alpha, epsilon:

  test level and tolerated coverage loss (defaults 0.05)

- power:

  planned power (default 0.80)

- ...:

  covariance and estimator arguments for \[sw_exposure_weights()\]

## Value

\\\kappa^\*\\ (0.2329 / \\\Gamma\\ at the defaults)

## Examples

``` r
sw_radius(design = c(6, 6, 6, 6, 0), phi = 2)    # 0.133
#> [1] 0.1332467
sw_radius(design = c(4, 4, 4, 4, 8), phi = 2)    # 0.214
#> [1] 0.2143148
```

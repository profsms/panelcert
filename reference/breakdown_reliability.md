# Self-consistent breakdown reliability

The fixed point \\\lambda = \lambda^\dagger(\lambda)\\ when the
corrected pilot \\\beta^\*/\lambda\\ is evaluated at the reliability
being solved for (Definition def-breakdown in the accompanying article).
Closed form \\\lambda^\dagger = t^\*/(t^\* + \eta^\dagger)\\ with \\t^\*
= \|\beta^\*\|\sqrt{\tau^{\*2}}/\sigma\\ — the specification's
conventional t-statistic, so no noise input is needed. Under
cluster-robust standardization pass the variance-inflation factor `psi`
(Remark rem-cluster): \\t^\*\\ is deflated by \\\sqrt{\psi}\\. Uses the
exact-inversion root \\\eta^\dagger\\.

## Usage

``` r
breakdown_reliability(
  beta_star,
  sigma,
  tau_star2,
  alpha = 0.05,
  delta = 0.05,
  psi = 1
)
```

## Arguments

- beta_star:

  attenuated slope from the FE regression

- sigma:

  residual standard deviation

- tau_star2:

  observed within variation of the regressor

- alpha:

  nominal test level

- delta:

  size-distortion tolerance

- psi:

  cluster variance-inflation factor (1 = i.i.d. standardization)

## Value

the breakdown reliability in \[0, 1\]

## Examples

``` r
breakdown_reliability(0.0908759, 0.5095233, 537.2959) # twins: ~0.864
#> [1] 0.8637104
```

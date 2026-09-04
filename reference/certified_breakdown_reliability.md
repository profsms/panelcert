# Certified breakdown reliability

Replaces the reported absolute t-statistic by \\\|t\| +
z\_{1-\gamma\_\beta}\\ in the point-breakdown formula:
\$\$\lambda^\dagger\_{\gamma\_\beta}=
(\|t\|+z\_{1-\gamma\_\beta})/(\|t\|+z\_{1-\gamma\_\beta}+\eta^\dagger).\$\$
Compare the result with a known/consistent within reliability or a lower
confidence bound. If that lower bound has coverage error
\\\gamma\_\lambda\\, the total false-certification bound is
\\\gamma\_\beta+\gamma\_\lambda\\; the second budget belongs to the
bound, not to this threshold formula.

## Usage

``` r
certified_breakdown_reliability(
  beta_star,
  sigma,
  tau_star2,
  alpha = 0.05,
  delta = 0.05,
  gamma_beta = 0.05,
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

- gamma_beta:

  coefficient-uncertainty error budget in (0, 0.5\]

- psi:

  cluster variance-inflation factor (1 = i.i.d. standardization)

## Value

the certified breakdown reliability in \[0, 1\]

## Examples

``` r
certified_breakdown_reliability(0.0908759, 0.5095233, 537.2959)
#> [1] 0.8985669
```

# Measurement-error adequacy diagnostic

Certifies whether naive inference on an already-estimated FE regression
is size-controlled under classical measurement error in the regressor.
The regression is reproduced via Frisch–Waugh from raw data (default
method) or consumed directly from a fitted fixest, plm, or `lm` object —
the model is never re-specified. Supply exactly one noise input:
`sigma_nu` (per-observation or scalar measurement-error SD), `codelow` +
`codehigh` (when the interval half-width is calibrated as one error SD),
or `reliability` (the within reliability \\\hat\lambda\\ directly).

## Usage

``` r
# S3 method for class 'fixest'
eiv_adequacy(object, ...)

# S3 method for class 'plm'
eiv_adequacy(object, ...)

# S3 method for class 'lm'
eiv_adequacy(object, x, unit, time, ...)

eiv_adequacy(object, ...)

# Default S3 method
eiv_adequacy(
  object,
  x,
  unit,
  time,
  sigma_nu = NULL,
  codelow = NULL,
  codehigh = NULL,
  reliability = NULL,
  controls = NULL,
  alpha = 0.05,
  delta = 0.05,
  gamma = 0.05,
  reliability_lower = NULL,
  gamma_lambda = 0,
  pilot = c("conservative", "point", "naive"),
  cluster = c("iid", "crve", "ar1"),
  psi = NULL,
  ...
)
```

## Arguments

- object:

  outcome vector (default method), or a fitted `fixest` / `plm` / `lm`
  model with one regressor and two-way fixed effects

- ...:

  passed between methods

- x:

  observed (mismeasured) regressor; for the `lm` method, the NAME of the
  regressor in the model frame

- unit:

  unit identifiers (any type); for the `lm` method, the name of the unit
  variable in the model frame

- time:

  period identifiers (any type); for the `lm` method, the name of the
  time variable in the model frame

- sigma_nu:

  per-observation (or scalar) measurement-error SD

- codelow:

  lower interval bound (with `codehigh`)

- codehigh:

  upper interval bound (with `codelow`)

- reliability:

  the within reliability \\\hat\lambda\\ in (0, 1\]

- controls:

  optional numeric nuisance-covariate vector or matrix. The target
  regressor, outcome, residual degrees of freedom, and measurement-
  error trace are all partialled with respect to these controls as well
  as the two fixed-effect sets.

- alpha:

  nominal test level

- delta:

  size-distortion tolerance

- gamma:

  coefficient-uncertainty budget \\\gamma\_\beta\\; the name is retained
  for backward compatibility

- reliability_lower:

  optional lower confidence bound for within reliability. If omitted,
  the computed reliability is treated as known/consistent, so
  certification is conditional on that treatment.

- gamma_lambda:

  coverage-error budget for `reliability_lower`

- pilot:

  `"conservative"` (default), `"point"`, or `"naive"`

- cluster:

  standardization of the t-test (cluster-robust extension): `"iid"`
  (default), `"crve"` (by-unit Arellano CR1 variance-inflation), or
  `"ar1"` (parametric within-unit AR(1)); \\\|\eta\|\\ and the breakdown
  are deflated by \\\sqrt{\hat\psi}\\

- psi:

  supply your own variance-inflation factor (overrides `cluster`)

## Value

an object of class `AdequacyReport`

## Details

Correct-by-default honesty machinery: `pilot = "conservative"` (default)
compares a reliability lower bound with the certified breakdown obtained
by replacing \\\|t\|\\ by \\\|t\|+z\_{1-\gamma\_\beta}\\; `"point"` is
the descriptive corrected-pilot verdict; `"naive"` plugs in the
attenuated \\\hat\beta^\*\\ and is anti-conservative — exposed for
comparison only, and labelled as such in the report.

## References

Halkiewicz, S. M. S. Breakdown Reliability for Saturated Fixed-Effect
Inference.

## Examples

``` r
n <- 200; unit <- rep(1:20, each = 10); time <- rep(1:10, times = 20)
x <- rnorm(n) + 0.3 * unit; y <- 0.5 * x + rnorm(n)
eiv_adequacy(y, x, unit, time, reliability = 0.9)
#> Panel Adequacy Report — Measurement Error
#> Design: n=200, N=20, T=10, d_K=29, rho=0.1450
#> Within reliability lambda_hat = 0.900   ((1-lambda)/lambda = 0.111)
#> Pilot: beta* = 0.5327 -> corrected beta0 = 0.5919 (se 0.0847)
#> Certified breakdown = 0.930   reliability lower bound = 0.900   |eta| upper bound = 0.959
#> Non-centrality |eta| = 0.777   Threshold (delta=0.05) = 0.652
#> Breakdown threshold = 0.915
#> Implied size of nominal 5% test: 12.1%
#> VERDICT: FLAGGED at delta=0.05
#> Note: formal certificate (prop-certificate): certified breakdown lambda_dagger_gamma = 0.930 is obtained by replacing |t| with |t| + z_(1-gamma_beta); the comparison uses reliability lower bound ell = 0.900. False certification is at most gamma_beta + gamma_lambda = 0.05 + 0 = 0.05, without requiring independence. Implied size shown is at the point pilot.
#> Note: CONDITIONAL RELIABILITY TREATMENT: no reliability_lower was supplied, so lambda_hat is treated as known/consistent and gamma_lambda = 0. A noisy finite-sample reliability estimate requires a lower confidence bound and its coverage-error budget.
#> Note: corrected pilot beta*/lambda_hat is not a consistent point estimate under weak information; reported with its sampling band (Corollary cor-slope)
#> Note: power tax: local power slope attenuated by sqrt(lambda) = 0.95 (Proposition prop-power)
```

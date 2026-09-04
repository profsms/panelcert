# Module A diagnostic: variance-estimator adequacy under FE saturation

Reproduces the user's FE regression of \`y\` on \`x\` with unit and time
fixed effects via Frisch-Waugh (never re-specified), then reports the
naive / df-corrected / HC0-HC3 variance hierarchy, the leverage
diagnostics, and whether the variance-estimator choice materially
changes inference at tolerance \`delta\`.

## Usage

``` r
# S3 method for class 'fixest'
leverage_report(object, ...)

# S3 method for class 'plm'
leverage_report(object, ...)

# S3 method for class 'lm'
leverage_report(object, x, unit, time, ...)

leverage_report(object, ...)

# Default S3 method
leverage_report(
  object,
  x,
  unit,
  time,
  alpha = 0.05,
  delta = 0.05,
  controls = NULL,
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

  regressor of interest; for the `lm` method, its NAME in the model
  frame

- unit:

  unit identifiers (any type); for the `lm` method, the name of the unit
  variable in the model frame

- time:

  period identifiers (any type); for the `lm` method, the name of the
  time variable in the model frame

- alpha:

  nominal test level

- delta:

  size-distortion tolerance

- controls:

  optional numeric nuisance-covariate vector or matrix.

## Value

an object of class `AdequacyReport`

## Details

Verdict: FLAGGED when the asymptotic HC0/naive over-rejection exceeds
\`alpha + delta\` at this design's saturation \`rho\`, or when
significance at level \`alpha\` flips across the estimators on this
data.

## References

Halkiewicz, S. M. S. Corrected diffuse-regime variance framework for
saturated fixed-effect specifications.

## Examples

``` r
n <- 200; unit <- rep(1:20, each = 10); time <- rep(1:10, times = 20)
x <- rnorm(n); y <- 0.5 * x + rnorm(n)
leverage_report(y, x, unit, time)
#> Panel Adequacy Report — Leverage / Variance (diffuse-regime companion)
#> Design: n=200, N=20, T=10, d_K=29, rho=0.1450
#> Max leverage max_i H_ii = 0.195 | spread hmax/hmin = 1.35
#> Design conditions: lambda_n = 0.0502 (N_eff = 65.8) | max|H_ii - rho| = 0.050
#> Realized score concentration: lambda_score = 0.1619 (N_eff,score = 16.5)
#> SE(beta): df-corrected 0.06731 | HC0 0.06949 | HC2 0.07603 | HC3 0.08319
#> beta_hat = 0.6022   t (HC2) = 7.92
#> Breakdown threshold = 0.296
#> Implied size of nominal 5% test: 7.0%
#> VERDICT: CERTIFIED at delta=0.05
#> Note: HC2 is the recommended default among the HC0--HC3 estimators reported here; it is leverage-adjusted but is not the Kline--Saggio--Solvsten leave-out estimator
#> Note: implied sizes are the conditional-homoskedastic limits of thm:hc (omega_eff^2 = sigma^2); under heteroskedasticity the HCc limits carry the extra factor omega^2/omega_eff^2, omega_eff^2 = (1-rho) omega^2 + rho mu
#> Note: HC3 over-correcting regime (rho = 0.145 > 0.1): HC3 intervals are artificially conservative (implied size 3.4%); HC2 is the preferred member of the HC0--HC3 family reported here
#> Note: sample Hessian V_n = 183.568 estimates (1-rho) n Qbar, not n Qbar (lem:hess(b)); the deflation-corrected primitive is Qhat = V_n/(n(1-rho)) = 1.0735
#> Note: realized score diagnostic (Paper A): lambda_score = 0.1619, N_eff,score = 16.5. This is a one-realization warning statistic, not by itself a consistent population concentration estimate.
```

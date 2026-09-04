# Module B diagnostic from regression summary output

For users who have only fitted-model output (no raw data): supply the
attenuated slope, residual SD, observed within variation, n and d_K,
plus exactly one of \`reliability\` or \`sigma_nu2\` (the mean squared
measurement-error SD).

## Usage

``` r
eiv_adequacy_summary(
  beta_star,
  sigma,
  tau_star2,
  n,
  d_K,
  reliability = NULL,
  sigma_nu2 = NULL,
  N = 0L,
  T = 0L,
  alpha = 0.05,
  delta = 0.05,
  gamma = 0.05,
  reliability_lower = NULL,
  gamma_lambda = 0,
  pilot = c("conservative", "point", "naive"),
  psi = NULL
)
```

## Arguments

- beta_star, sigma, tau_star2:

  regression output: attenuated slope, residual SD, observed within
  variation of the regressor

- n, d_K:

  sample size and fixed-effect dimension

- reliability, sigma_nu2:

  noise input (exactly one): within reliability, or the mean squared
  measurement-error SD

- N, T:

  optional design shape (0 = unknown)

- alpha, delta, gamma, gamma_lambda, reliability_lower, pilot:

  as in \[eiv_adequacy()\]

- psi:

  cluster variance-inflation factor (Remark rem-cluster in the
  accompanying measurement-error article); the `cluster` presets are
  unavailable without raw data

## Value

an object of class `AdequacyReport`

## Examples

``` r
# Repeated-report twins application: both independent-report estimates flag.
eiv_adequacy_summary(0.0908759, 1.0, (0.0219815^-2), n = 147, d_K = 4,
                     reliability = 0.574904, pilot = "point")
#> Panel Adequacy Report — Measurement Error
#> Design: n=147, d_K=4, rho=0.0272 (from summary input)
#> Within reliability lambda_hat = 0.575   ((1-lambda)/lambda = 0.739)
#> Pilot: beta* = 0.09088 -> corrected beta0 = 0.1581 (se 0.0382)
#> Non-centrality |eta| = 3.057   Threshold (delta=0.05) = 0.652
#> Breakdown threshold = 0.864
#> Implied size of nominal 5% test: 86.4%
#> VERDICT: FLAGGED at delta=0.05
#> Note: POINT PASS, not a certificate (rem-plugin): the corrected pilot at its point estimate is descriptive. Under weak information eta_hat converges to a nondegenerate random multiple (|B|/|beta0|)|eta| of the target -- median close to it, but no concentration. Its sampling variability is a first-order feature of the regime, not a vanishing approximation error. For a size-controlled statement use pilot = "conservative".
#> Note: corrected pilot beta*/lambda_hat is not a consistent point estimate under weak information; reported with its sampling band (Corollary cor-slope)
#> Note: far from the threshold (|eta| > 1): the local quadratic approximation is uninformative here; the verdict uses exact inversion (Remark rem-exact-cv)
#> Note: power tax: local power slope attenuated by sqrt(lambda) = 0.76 (Proposition prop-power)
```

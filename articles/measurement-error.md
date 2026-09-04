# Module B: measurement-error adequacy

## The problem

Fixed-effect saturation can leave a noisy regressor with little reliable
within variation. *Breakdown Reliability for Saturated Fixed-Effect
Inference* shows that the resulting size distortion of the naive t-test
is governed by one number, the **within reliability**

`lambda_hat = 1 - mean(sigma_nu^2) * (n - d_K) / tau_star2`,

the share of post-projection variation in the observed regressor that is
signal. The feasible non-centrality is
`|eta| = (|beta0|/sigma) * (1 - lambda) * sqrt(tau_star2)`, and naive
inference is `delta`-adequate iff `|eta| <= eta_dagger(alpha, delta)` —
the **exact-inversion** threshold (about 0.652 at
`alpha = delta = 0.05`). The quadratic closed form
`sqrt(delta/(z*phi(z)))` is its accurate companion; a linear expansion
is spurious (the size is even in `eta`) and must not be used.

## Noise inputs

Three pathways, in decreasing order of directness:

1.  **published measurement-model posteriors** (V-Dem-style):
    `codelow`/`codehigh` interval bounds, converted by
    [`reliability_from_interval()`](https://profsms.github.io/panelcert/reference/reliability_from_interval.md);
2.  **repeated measurements on an identical projected sample**:
    `reliability_from_repeats(first, second)`, using the covariance
    estimator without an equal-error-variance restriction or
    `method = "equal_variance"` when that additional restriction is
    defensible;
3.  **validation-study reliability ratios**:
    `reliability_from_ratio(r, within_sd)`;
4.  **direct** `sigma_nu` or `reliability` values, including a
    sensitivity grid.

Both repeated-report formulas require uncorrelated reporting errors. If
reports share person-specific error, supply a correlated-error
reliability estimate or show a sensitivity range instead of treating
either formula as identified.

## The honesty machinery, demonstrated

The most natural implementation — plugging the attenuated `beta*` from
the regression table into the threshold — is biased toward *certifying
failing specifications*, by exactly the factor `lambda` whose smallness
is the reason to run the diagnostic (the article’s proposition on
pilots). The package defaults to the corrected pilot at its upper
confidence bound; the naive form is opt-in and labelled. On a design
built to fail:

``` r

library(panelcert)
N <- 100; T <- 8
unit <- rep(1:N, each = T); time <- rep(1:T, times = N)
n <- N * T
x_true <- 0.5 * rnorm(n)
x_obs <- x_true + rnorm(n, sd = 0.45)   # within reliability ~ 0.55
y <- 1.0 * x_true + rnorm(n)

r_default <- eiv_adequacy(y, x_obs, unit, time, sigma_nu = 0.45)
r_naive <- eiv_adequacy(y, x_obs, unit, time, sigma_nu = 0.45, pilot = "naive")
c(default = r_default$verdict, naive = r_naive$verdict)
#>   default     naive 
#> "FLAGGED" "FLAGGED"
```

The naive pilot understates `|eta|` by the factor `lambda` and can
certify a specification whose true size is several times nominal — in
the V-Dem application it certifies the legislative-constraints index
whose exact size is 14%. The corrected default flags it.

``` r

r_default
#> Panel Adequacy Report — Measurement Error
#> Design: n=800, N=100, T=8, d_K=107, rho=0.1338
#> Within reliability lambda_hat = 0.508   ((1-lambda)/lambda = 0.968)
#> Pilot: beta* = 0.4988 -> corrected beta0 = 0.9819 (se 0.123)
#> Certified breakdown = 0.936   reliability lower bound = 0.508   |eta| upper bound = 9.310
#> Non-centrality |eta| = 7.717   Threshold (delta=0.05) = 0.652
#> Breakdown threshold = 0.924
#> Implied size of nominal 5% test: 100.0%
#> VERDICT: FLAGGED at delta=0.05
#> Note: formal certificate (prop-certificate): certified breakdown lambda_dagger_gamma = 0.936 is obtained by replacing |t| with |t| + z_(1-gamma_beta); the comparison uses reliability lower bound ell = 0.508. False certification is at most gamma_beta + gamma_lambda = 0.05 + 0 = 0.05, without requiring independence. Implied size shown is at the point pilot.
#> Note: CONDITIONAL RELIABILITY TREATMENT: no reliability_lower was supplied, so lambda_hat is treated as known/consistent and gamma_lambda = 0. A noisy finite-sample reliability estimate requires a lower confidence bound and its coverage-error budget.
#> Note: corrected pilot beta*/lambda_hat is not a consistent point estimate under weak information; reported with its sampling band (Corollary cor-slope)
#> Note: far from the threshold (|eta| > 1): the local quadratic approximation is uninformative here; the verdict uses exact inversion (Remark rem-exact-cv)
#> Note: power tax: local power slope attenuated by sqrt(lambda) = 0.71 (Proposition prop-power)
```

Alongside the verdict, the report carries the **breakdown reliability**
`lambda_dagger` — how reliable the regressor would have to be for naive
inference to pass — and the `sqrt(lambda)` **power tax**, which applies
even to certified specifications.

## The empirical applications, reproduced from bundled data

The V-Dem democracy–growth panel ships with the package as `vdem` (the
measurement-error SDs are published data). The article’s two-pole result
runs directly off it — no download:

``` r

# Aggregate polyarchy: certified (lambda ~ 0.90, implied size ~5.6%)
p <- vdem[complete.cases(vdem$ly, vdem$v2x_polyarchy, vdem$v2x_polyarchy_sd), ]
eiv_adequacy(p$ly, p$v2x_polyarchy, p$iso, p$year, sigma_nu = p$v2x_polyarchy_sd)
#> Panel Adequacy Report — Measurement Error
#> Design: n=8930, N=163, T=59, d_K=221, rho=0.0247
#> Within reliability lambda_hat = 0.898   ((1-lambda)/lambda = 0.113)
#> Pilot: beta* = 0.06096 -> corrected beta0 = 0.06785 (se 0.0324)
#> Certified breakdown = 0.851   reliability lower bound = 0.898   |eta| upper bound = 0.422
#> Non-centrality |eta| = 0.236   Threshold (delta=0.05) = 0.652
#> Breakdown threshold = 0.762
#> Implied size of nominal 5% test: 5.6%
#> VERDICT: FORMALLY CERTIFIED at (alpha, delta, gamma) = (0.05, 0.05, 0.05)
#> Note: formal certificate (prop-certificate): certified breakdown lambda_dagger_gamma = 0.851 is obtained by replacing |t| with |t| + z_(1-gamma_beta); the comparison uses reliability lower bound ell = 0.898. False certification is at most gamma_beta + gamma_lambda = 0.05 + 0 = 0.05, without requiring independence. Implied size shown is at the point pilot.
#> Note: CONDITIONAL RELIABILITY TREATMENT: no reliability_lower was supplied, so lambda_hat is treated as known/consistent and gamma_lambda = 0. A noisy finite-sample reliability estimate requires a lower confidence bound and its coverage-error budget.
#> Note: corrected pilot beta*/lambda_hat is not a consistent point estimate under weak information; reported with its sampling band (Corollary cor-slope)
#> Note: power tax: local power slope attenuated by sqrt(lambda) = 0.95 (Proposition prop-power)

# Legislative constraints: flagged (lambda ~ 0.55, implied size ~14%)
l <- vdem[complete.cases(vdem$ly, vdem$v2xlg_legcon, vdem$v2xlg_legcon_sd), ]
eiv_adequacy(l$ly, l$v2xlg_legcon, l$iso, l$year,
             sigma_nu = l$v2xlg_legcon_sd)$verdict
#> [1] "FLAGGED"
```

The second application ships as `twins`, the public repeated-report
extract from the Ashenfelter–Krueger design. Both schooling reports must
be projected on the identical controls and complete-case sample before
estimating reliability.

``` r

needed <- c("DLHRWAGE", "DEDUC1", "DEDUC2", "DTEN", "DMARRIED", "DUNCOV")
d <- twins[complete.cases(twins[needed]), needed]
fit <- lm(DLHRWAGE ~ DEDUC1 + DTEN + DMARRIED + DUNCOV, data = d)
W <- model.matrix(~ DTEN + DMARRIED + DUNCOV, data = d)
y <- qr.resid(qr(W), d$DLHRWAGE)
x <- qr.resid(qr(W), d$DEDUC1)
z <- qr.resid(qr(W), d$DEDUC2)
lambda_cov <- reliability_from_repeats(x, z)
lambda_equal <- reliability_from_repeats(x, z, method = "equal_variance")

beta_star <- unname(coef(fit)["DEDUC1"])
se_star <- unname(coef(summary(fit))["DEDUC1", "Std. Error"])
tau_star2 <- sum(x^2)
sigma <- se_star * sqrt(tau_star2)
twins_cov <- eiv_adequacy_summary(beta_star, sigma, tau_star2, 147, 4,
                                  reliability = lambda_cov, pilot = "point")
twins_equal <- eiv_adequacy_summary(beta_star, sigma, tau_star2, 147, 4,
                                    reliability = lambda_equal, pilot = "point")
c(lambda_cov = lambda_cov, lambda_equal = lambda_equal,
  breakdown = twins_cov$breakdown,
  size_cov = twins_cov$implied_size, size_equal = twins_equal$implied_size)
#>   lambda_cov lambda_equal    breakdown     size_cov   size_equal 
#>    0.5749037    0.5514170    0.8637104    0.8636693    0.9197285
```

Both repeat models flag: reliability is 0.575 or 0.551 against a
breakdown of 0.864, with implied coverage 13.6% or 8.0%. Rouse’s
correlated-report sensitivity also flags (reliability 0.748 below
breakdown 0.872) but gives much less extreme implied coverage, 67.9%.
The verdict is robust to the reporting- error model; the magnitude is
not. The similarity between the covariance- corrected coefficient and
one directional IV estimate is algebraic, so it is corroboration rather
than external validation.

## Cluster-robust standardization

Applied panels often cluster standard errors by unit. The package
computes the realized Arellano variance-inflation factor with
`cluster = "crve"` and rescales `|eta|` and the breakdown by
`sqrt(psi_hat)`. The sign of `psi_hat - 1` is not fixed: clustering can
soften or worsen the verdict, depending jointly on the regressor and
error dependence. In the V-Dem middle case, legislative constraints is
flagged under i.i.d. standard errors but obtains an uncertified point
pass under the country-clustered ones an applied author would report —

``` r

eiv_adequacy(l$ly, l$v2xlg_legcon, l$iso, l$year,
             sigma_nu = l$v2xlg_legcon_sd,
             pilot = "point", cluster = "crve")$verdict
#> [1] "POINT_PASS"
```

— while the judicial-constraints flag survives clustering (`psi_hat` ~
25 still leaves an implied size of 67%). Pass `psi =` to supply your own
factor (e.g. in
[`eiv_adequacy_summary()`](https://profsms.github.io/panelcert/reference/eiv_adequacy_summary.md),
where raw data is unavailable).

Every application entry in the article is regression-locked this way;
`data(package = "panelcert")` lists the bundled datasets.

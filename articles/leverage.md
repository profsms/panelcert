# Diffuse companion: variance estimators under FE saturation

## The problem

Under the repaired diffuse-regime framework, conventional variance
estimators split in opposite directions as fixed effects saturate
(`rho = d_K/n` grows):

- **HC0** is biased downward by the factor `(1 - rho)`: t-tests
  over-reject, with asymptotic size `2(1 - pnorm(z * sqrt(1 - rho)))` —
  about 9% at `rho = 0.25` and 17% at `rho = 0.5` for a nominal 5% test;
- **HC3** over-corrects by exactly the inverse factor `1/(1 - rho)`:
  size drops below 1% at `rho = 0.5`, and confidence intervals are
  artificially wide;
- **HC1 and HC2** remove the saturation factor in the conditionally
  homoskedastic, uniform-leverage limit. Under heteroskedasticity their
  limits retain the paper’s `omega^2/omega_eff^2` factor.

HC2 is the MacKinnon–White leverage adjustment. It is not the
Kline–Saggio–Solvsten leave-out estimator, which this package does not
implement. None of these asymptotic statements licenses the module when
identifying variation is concentrated or leverage balance fails.

[`leverage_report()`](https://profsms.github.io/panelcert/reference/leverage_report.md)
reproduces your FE regression via Frisch–Waugh, computes the full
hierarchy of standard errors with the exact full-regression leverage
`H_ii = (P_K)_ii + Xt_i^2 / (X'MX)`, and reports whether the estimator
choice materially changes inference on your design.

## A saturated design

``` r

library(panelcert)
N <- 12; T <- 5
unit <- rep(1:N, each = T); time <- rep(1:T, times = N)
n <- N * T
x <- rnorm(n) + 0.3 * unit
y <- 0.8 * x + rnorm(n)
leverage_report(y, x, unit, time)
#> Panel Adequacy Report — Leverage / Variance (diffuse-regime companion)
#> Design: n=60, N=12, T=5, d_K=16, rho=0.2667
#> Max leverage max_i H_ii = 0.355 | spread hmax/hmin = 1.33
#> Design conditions: lambda_n = 0.0885 (N_eff = 24.9) | max|H_ii - rho| = 0.089
#> Realized score concentration: lambda_score = 0.1482 (N_eff,score = 14.3)
#> SE(beta): df-corrected 0.1369 | HC0 0.1174 | HC2 0.1412 | HC3 0.1698
#> beta_hat = 0.9077   t (HC2) = 6.43
#> Breakdown threshold = 0.296
#> Implied size of nominal 5% test: 9.3%
#> VERDICT: CERTIFIED at delta=0.05
#> Note: HC2 is the recommended default among the HC0--HC3 estimators reported here; it is leverage-adjusted but is not the Kline--Saggio--Solvsten leave-out estimator
#> Note: implied sizes are the conditional-homoskedastic limits of thm:hc (omega_eff^2 = sigma^2); under heteroskedasticity the HCc limits carry the extra factor omega^2/omega_eff^2, omega_eff^2 = (1-rho) omega^2 + rho mu
#> Note: HC3 over-correcting regime (rho = 0.267 > 0.1): HC3 intervals are artificially conservative (implied size 2.2%); HC2 is the preferred member of the HC0--HC3 family reported here
#> Note: sample Hessian V_n = 64.24 estimates (1-rho) n Qbar, not n Qbar (lem:hess(b)); the deflation-corrected primitive is Qhat = V_n/(n(1-rho)) = 1.46
#> Note: realized score diagnostic (Paper A): lambda_score = 0.1482, N_eff,score = 14.3. This is a one-realization warning statistic, not by itself a consistent population concentration estimate.
```

At `rho = (N + T - 1)/n` near 0.27, the HC0 and HC3 standard errors
bracket HC2 by visibly asymmetric factors, and the report’s implied size
quantifies what an HC0 user would actually suffer. The
`breakdown threshold` line is the saturation `rho_dagger` (about 0.296
at `alpha = delta = 0.05`) at which HC0-based inference exits the
tolerance.

## Reading the flags

Two notes deserve attention when they appear:

- the **HC3 over-correction flag** fires at the module’s `rho > 0.1`
  reporting convention; HC2 is the preferred member of the HC0–HC3
  family shown here;
- the **leverage non-uniformity flag** (`hmax/hmin` large) marks designs
  — small country-year cells, “anchor” units — where the
  uniform-leverage limit is not licensed. Use a separately implemented
  leave-out method rather than treating any result from this module as a
  certificate.

Fitted-model ingestion works the same way:

``` r

d <- data.frame(y = y, x = x, unit = unit, time = time)
m <- fixest::feols(y ~ x | unit + time, data = d)
leverage_report(m)$verdict
#> [1] "CERTIFIED"
```

## Diffuse-companion application, reproduced from bundled data

The firm panel also used by current Paper A ships as `fscore`: Piotroski
F-Scores and one-year-ahead returns for Warsaw, Budapest and Prague
listings, 2010–2024 (217 firm-years, 19 firms). The preferred
specification — log return on the F-Score, firm and year fixed effects —
runs directly off it:

``` r

library(panelcert)
r <- leverage_report(log1p(fscore$ret), fscore$fscore, fscore$uid, fscore$year)
c(beta = round(r$statistic$beta, 5),
  HC0_HC2 = round(r$statistic$se_hc0 / r$statistic$se_hc2, 3),
  HC3_HC2 = round(r$statistic$se_hc3 / r$statistic$se_hc2, 3),
  rho = round(r$design$rho, 3))
#>     beta  HC0_HC2  HC3_HC2      rho 
#> -0.00147  0.91300  1.09600  0.15200
```

The empirical HC0/HC2 and HC3/HC2 ratios track the theoretical
`sqrt(1 - rho)` and `1/sqrt(1 - rho)` almost exactly, and the F-Score
coefficient is indistinguishable from zero. The most favourable
specification is Poland alone, whose `|t|` is still only 0.61:

``` r

pl <- fscore$country == "Poland"
e <- leverage_report(log1p(fscore$ret[pl]), fscore$fscore[pl],
                     fscore$uid[pl], fscore$year[pl])
round(abs(e$statistic$beta / e$statistic$se_hc2), 2)
#> [1] 0.61
```

The country-year-fixed-effect specification is the leverage cautionary
tale of the flag above: each firm sits in exactly one country, so the
firm + country-year design splits into three disconnected country blocks
and the realized FE dimension is `d_K = 61`, which
[`leverage_report()`](https://profsms.github.io/panelcert/reference/leverage_report.md)
computes from the actual design via union-find rather than a nominal
`19 + 45 - 1` count. Pass `paste(fscore$country, fscore$year)` as the
time id to reproduce it. `data(package = "panelcert")` lists the bundled
datasets.

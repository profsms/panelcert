# TWFE-heterogeneity adequacy for staggered DiD

## The design statistic, before any outcome exists

Under heterogeneous treatment effects, the TWFE coefficient is a
weighted average with weights `w` that need not be convex. The
accompanying article shows the *size distortion* of the TWFE t-test is
governed by

`Gamma = sqrt(N1) * ||w - u||`,

the scaled distance of the weights from uniform — a **pure design
statistic**, computable from the adoption pattern alone. `Gamma = 0`
exactly for clean block designs; it grows with staggering and negative
weights. Crucially, `Gamma` is invariant to replicating the panel while
the absolute weighting bias shrinks — so a design can look safe to
estimand-level diagnostics while its t-test stays arbitrarily distorted.

``` r

library(panelcert)
N <- 48; T <- 12
unit <- rep(1:N, each = T); time <- rep(1:T, times = N)

# a block design: one adoption date, never-treated reservoir
ft_block <- ifelse(unit <= 16, 5, NA)
twfe_design(unit, time, ft_block)$statistic$Gamma
#> [1] 0

# three cohorts, no never-treated reservoir
ft_stag <- ifelse(unit <= 16, 3, ifelse(unit <= 32, 7, 11))
twfe_design(unit, time, ft_stag)
#> Panel Adequacy Report — TWFE Heterogeneity
#> Design: n=576, N=48, T=12, d_K=59, rho=0.1024
#> Design statistic Gamma = 1.541   negative-weight share = 11.1%
#>   restricted ladder: Gamma_gt = 1.541 | Gamma_c+e = 1.541 | Gamma_evt = 1.541 | Gamma_coh = 0.652
#> Breakdown threshold = 0.423
#> VERDICT: INCONCLUSIVE
#> Note: pre-outcome: the saturated group-time class requires c/sigma <= 0.423 (= eta-dagger/Gamma_gt); supply the outcome (twfe_adequacy) to calibrate c/sigma and obtain a boundary-robust certificate
```

Run this before committing to a design: the report’s breakdown value is
the largest heterogeneity-to-noise ratio `c/sigma` the design tolerates.

## The inference layer

With an outcome,
[`twfe_adequacy()`](https://profsms.github.io/panelcert/reference/twfe_adequacy.md)
estimates heterogeneity in a nested ladder: cohort, event time, additive
cohort plus event time, and the saturated post-treatment group-time
class. In a balanced panel the last class spans every treated-cell
profile relevant to TWFE weighting, so `Gamma_gt = Gamma`. The function
maps the estimated profile to two distinct quantities:

- the **worst-case size envelope** maximizes the mean shift over that
  calibrated class. It says what the design can certify uniformly; it is
  an upper bound, not a point prediction of the test’s realized
  rejection probability;
- the **directional plug-in** inserts the estimated signed group-time
  profile. It is orientation-sensitive corroboration, not a uniform
  certificate.

The default scale is direct CR1. The function forms the unit cluster
scores, sets `psi_direct` to their CR1 variance relative to the
homoskedastic scale, and reports the equivalent `q_hat` and `se_cr1`.
Parametric AR(1), i.i.d., and a positive user-supplied `psi` remain
explicit sensitivity routes. The sign of the cluster adjustment is never
assumed.

The heterogeneity pilot is covariance-aware. The cluster-score
covariance of the group-time effects supplies the projected
quadratic-form noise trace, which is removed before the descriptive
point radius is formed. Because that quadratic pilot is nonregular at
zero heterogeneity, it does not determine the formal decision. The
regular lower test applies the reverse triangle inequality to an HC2
cluster-multiplier radius. The upper certificate instead inverts the
noncentral chi-square law of the projected Wald statistic. It is
available only when the score covariance spans the full prespecified
heterogeneity class; a rank failure produces `INCONCLUSIVE`, rather than
an artificial upper bound. HC3 is reported as a leverage sensitivity
check.

Choose the decision class with `heterogeneity_class = "group_time"`,
`"additive"`, `"cohort"`, or `"event"` before inspecting the outcome.
The reported lower and upper bounds are separate one-sided statements,
since only one is used for any verdict. If several classes or outcomes
are screened as a family, divide `gamma` by the number of planned tests.
The optional `q_band` argument reports the more general full
confidence-ball construction with a relative denominator band; that
construction is less powerful and does not change the main verdict.

``` r

delta_g <- c(-0.6, 0, 0.6)                       # heterogeneous cohort effects
eff <- delta_g[ifelse(unit <= 16, 1, ifelse(unit <= 32, 2, 3))]
D <- as.numeric(!is.na(ft_stag) & time >= ft_stag)
y <- 0.1 * unit / N + D * eff + rnorm(N * T, sd = 0.5)
twfe_adequacy(y, unit, time, ft_stag)
#> Panel Adequacy Report — TWFE Heterogeneity
#> Design: n=576, N=48, T=12, d_K=59, rho=0.1024
#> Design statistic Gamma = 1.541   negative-weight share = 11.1%
#>   restricted ladder: Gamma_gt = 1.541 | Gamma_c+e = 1.541 | Gamma_evt = 1.541 | Gamma_coh = 0.652
#> Cluster normalization (direct; psi_hat = 0.751): Gamma_gt,CR = 1.779 | Gamma_c+e,CR = 1.779
#> TWFE beta_hat = 0.07581   sigma = 0.512
#> Trace-debiased point pilots c_S/sigma: gt = 0 | c+e = 0 | evt = 0.368 | coh = 0
#> Point worst-case envelopes: group-time = 5.0% | c+e = 5.0% | cohort 5.0% | event 10.0%
#>   descriptive point-envelope bootstrap (B=999): median 80.7%, 95% [5.0, 100.0]; psi in [0.67, 0.84]
#> Directional plug-in: eta = NA (alignment NA), size NA%
#>   wild directional size: median NA%, 95% [NA, NA]
#> Sign-reversal RMS threshold = 0.04919
#> Worst-case |eta| envelope = 0.000   Threshold (delta=0.05) = 0.652
#> Breakdown threshold = 0.367
#> Worst-case size envelope for nominal 5% test: 5.0%
#> VERDICT: INCONCLUSIVE
#> Note: one-sided procedures unavailable: the group-time estimator does not cover every treated cell
#> Note: point pilot: covariance trace computed from the fitted Rademacher score covariance; 999 multiplier draws calibrate the 95.0% HC2 lower norm radius; the upper procedure uses projected-Wald noncentrality inversion
#> Note: normalization = direct: psi_hat = 0.751; direct CR1 psi = 0.751; AR(1) psi = 0.931 (rho = -0.065)
#> Note: directional plug-in unavailable because the group-time profile does not cover every treated cell
```

The internal group-time estimator uses either not-yet-treated units (the
default, including never treated) or never-treated units only. Select
the published comparison group with `controls = "not_yet"` or
`controls = "never"`. This is a diagnostic implementation, not a
replacement for a full heterogeneity-robust DiD estimator.

## The applications, reproduced from bundled data

The package includes the exact four-column analysis extract for the
Callaway–Sant’Anna minimum-wage application, the balanced-subset Brazil
property-tax reanalysis, and the castle-doctrine and divorce panels. The
minimum-wage data contain 15,988 county-years from 2,284 counties over
seven years; Brazil contains 34,080 municipality-years from 2,840
municipalities. The current diagnostics are:

- minimum wage: an equally weighted treated-cell ATT of about -5.2%, a
  31.6% saturated-class point envelope, and a regular lower bound above
  the threshold, so uniform certification is withheld;
- Brazil property tax: 17.3% of treated-cell weights are negative and
  even the smallest class-specific lower bound exceeds the threshold, so
  all four prespecified classes withhold certification;
- castle doctrine: the prespecified additive homicide class has an HC2
  covariance-aware upper bound near 0.236 and certifies. The saturated
  class is inconclusive because its projected covariance has rank 18
  rather than 19;
- no-fault divorce: a roughly 70% saturated-class point envelope. Its
  lower bound does not clear the threshold, while the saturated upper
  procedure is unavailable because the projected covariance has rank 49
  rather than 167. The fixed-`T` approximation is also strained
  (`G = 49`, `T = 27`).

The last result illustrates why the point envelope and its uncertainty
should both be reported. The package never relabels a worst-case
envelope as realized size.

``` r

data(minimum_wage)
mw <- twfe_adequacy(minimum_wage$y, minimum_wage$uid,
                    minimum_wage$tid, minimum_wage$ft,
                    controls = "never", bootstrap = 999)
c(target_att = mw$statistic$att_target,
  envelope = mw$statistic$size_gt,
  K_lower = mw$statistic$K_lower_gt,
  directional = mw$statistic$size_directional,
  verdict = mw$verdict)

data(brazil)
brazil_report <- twfe_adequacy(brazil$y, brazil$uid, brazil$tid, brazil$ft,
                                bootstrap = 999)
brazil_report$statistic[c("verdict_coh", "verdict_evt",
                          "verdict_cmb", "verdict_gt")]

# Castle certifies for its prespecified additive class; the saturated class is
# retained as a transparent rank-inconclusive sensitivity.
castle_report <- twfe_adequacy(castle$y, castle$uid, castle$tid, castle$ft,
                               heterogeneity_class = "additive")
castle_saturated <- twfe_adequacy(castle$y, castle$uid, castle$tid, castle$ft)
divorce_report <- twfe_adequacy(divorce$y, divorce$uid, divorce$tid, divorce$ft)
c(castle_additive = castle_report$verdict,
  castle_saturated = castle_saturated$verdict,
  divorce_saturated = divorce_report$verdict)
```

`data(package = "panelcert")` lists the bundled datasets. The full
numerical audit is also executable through the package tests and the
sibling `reproduction/reproduce_with_package.jl` harness.

# panelcert: certifying naive panel inference

## What the package does

`panelcert` diagnoses whether conventional fixed-effect inference is
reliable and implements Paper A’s exact nuisance-annihilating contrast
alternative. The diagnostic modules consume a fitted model
([`fixest::feols`](https://lrberge.github.io/fixest/reference/feols.html),
[`plm::plm`](https://rdrr.io/pkg/plm/man/plm.html), `lm`) or raw panel;
the contrast engine works from the design and outcome. Reports use four
plain-language verdicts:

- **CERTIFIED at tolerance delta** — the worst-case size of the
  nominal-alpha test stays below `alpha + delta`;
- **FLAGGED (implied size X%)** — it does not, and the report quantifies
  the failure and points to the appropriate remedy;
- **POINT PASS** — a plug-in or finite-threshold comparison passes, but
  is not a formal confidence-band certificate;
- **INCONCLUSIVE** — a pre-outcome design statistic that needs a pilot
  input.

The current source map is:

| Source | Problem | Headline output |
|----|----|----|
| Paper A | concentrated identifying variation | `lambda_n`, `kappa`, exact sign-flip interval |
| Paper B | classical measurement error | within reliability `lambda_hat` |
| Paper C | TWFE heterogeneity in staggered DiD | design statistic `Gamma` |
| Diffuse companion | variance-estimator choice under FE saturation | HC0–HC3 hierarchy, leverage |

## Bundled data — the paper applications run offline

The small reference panels ship with the package. After
[`library(panelcert)`](https://profsms.github.io/panelcert) the datasets
are available by name:

``` r

library(panelcert)
data(package = "panelcert")$results[, "Item"]   # includes brazil, castle, divorce, minimum_wage
#> [1] "brazil"       "castle"       "divorce"      "fscore"       "grunfeld"    
#> [6] "minimum_wage" "psid"         "twins"        "vdem"
twfe_adequacy(castle$y, castle$uid, castle$tid, castle$ft,
              heterogeneity_class = "additive")$verdict
#> [1] "CERTIFIED"
```

Paper A’s F-score workflow, Papers B and C, and the diffuse companion
use these objects. The public KSS worker–firm testing extract is
intentionally not vendored; Paper A’s replication code downloads it from
the source repository.

## The shared design summary

Every module rests on the same design primitives: the fixed-effect
dimension `d_K` (union-find on the unit–time graph), the saturation
`rho = d_K/n`, and the two-way within transformation.

``` r

library(panelcert)
unit <- rep(1:40, each = 12)
time <- rep(1:12, times = 40)
design_summary(unit, time)
#> Panel Design Summary
#>   n = 480 obs | N = 40 units | T = 12 periods
#>   d_K = 51 (connected components: 1) | rho = d_K/n = 0.1062
```

## A sixty-second tour

Paper A exact inference in a balanced worker-time panel:

``` r

set.seed(3)
n <- length(unit)
x <- rnorm(n)
y <- 0.4 * x + rnorm(n)
cycle_report(y, x, unit, time, nflips = 199, interval = FALSE)
#> Panel Adequacy Report — Concentrated Identifying Variation
#> Design: n=480, N=40, T=12, d_K=51, rho=0.1062
#> Concentration: lambda_n = 0.0287 (N_eff = 153.1)
#> Realized score concentration: lambda_score = 0.0812 (N_eff,score = 52.9)
#> Capture kappa_C = 0.9079 over 115 supports (cycle-space dim 429) | capture-implied SE ratio 1.050x | max share 0.041
#> VERDICT: POINT PASS at delta=0.05 (descriptive — not a certificate)
#> Diagnostic notes hidden (6); call show_notes(report) to display them.
```

Module B on a panel with a noisy regressor (measurement-error SD known,
e.g. from a validation study or a published measurement model):

``` r

x_true <- rnorm(n) + 0.2 * unit          # true regressor
x_obs <- x_true + rnorm(n, sd = 0.3)     # observed with error
y <- 0.6 * x_true + rnorm(n)
eiv_adequacy(y, x_obs, unit, time, sigma_nu = 0.3)
#> Panel Adequacy Report — Measurement Error
#> Design: n=480, N=40, T=12, d_K=51, rho=0.1062
#> Within reliability lambda_hat = 0.916   ((1-lambda)/lambda = 0.091)
#> Pilot: beta* = 0.5414 -> corrected beta0 = 0.5908 (se 0.0516)
#> Certified breakdown = 0.953   reliability lower bound = 0.916   |eta| upper bound = 1.197
#> Non-centrality |eta| = 1.046   Threshold (delta=0.05) = 0.652
#> Breakdown threshold = 0.946
#> Implied size of nominal 5% test: 18.2%
#> VERDICT: FLAGGED at delta=0.05
#> Note: formal certificate (prop-certificate): certified breakdown lambda_dagger_gamma = 0.953 is obtained by replacing |t| with |t| + z_(1-gamma_beta); the comparison uses reliability lower bound ell = 0.916. False certification is at most gamma_beta + gamma_lambda = 0.05 + 0 = 0.05, without requiring independence. Implied size shown is at the point pilot.
#> Note: CONDITIONAL RELIABILITY TREATMENT: no reliability_lower was supplied, so lambda_hat is treated as known/consistent and gamma_lambda = 0. A noisy finite-sample reliability estimate requires a lower confidence bound and its coverage-error budget.
#> Note: corrected pilot beta*/lambda_hat is not a consistent point estimate under weak information; reported with its sampling band (Corollary cor-slope)
#> Note: far from the threshold (|eta| > 1): the local quadratic approximation is uninformative here; the verdict uses exact inversion (Remark rem-exact-cv)
#> Note: power tax: local power slope attenuated by sqrt(lambda) = 0.96 (Proposition prop-power)
```

Module C, before any outcome exists — vetting a staggered adoption
design:

``` r

ft <- ifelse(unit <= 15, 4, ifelse(unit <= 30, 8, NA))
twfe_design(unit, time, ft)
#> Panel Adequacy Report — TWFE Heterogeneity
#> Design: n=480, N=40, T=12, d_K=51, rho=0.1062
#> Design statistic Gamma = 1.035   negative-weight share = 35.7%
#>   restricted ladder: Gamma_gt = 1.035 | Gamma_c+e = 0.960 | Gamma_evt = 0.956 | Gamma_coh = 0.488
#> Breakdown threshold = 0.630
#> VERDICT: INCONCLUSIVE
#> Note: pre-outcome: the saturated group-time class requires c/sigma <= 0.630 (= eta-dagger/Gamma_gt); supply the outcome (twfe_adequacy) to calibrate c/sigma and obtain a boundary-robust certificate
```

The verdict logic is *correct by default*: Module B evaluates the
certificate at the upper confidence bound of the attenuation-corrected
pilot (the naive attenuated pilot is available only as an explicitly
labelled opt-in); Module C uses the shrinkage-corrected cohort
dispersion and computes the cluster variance-inflation factor from the
realized design rather than assuming its direction. These defaults guard
against specific errors the source papers document; the module vignettes
explain each one.

## Where the numbers come from

All thresholds derive from one limit experiment: the studentized
statistic is asymptotically `N(eta, 1)`, the two-sided size is even in
`eta` with leading distortion `z * phi(z) * eta^2` (quadratic — a linear
expansion is spurious and badly anti-conservative), and the operational
threshold `eta_dagger(alpha, delta)` is the exact-inversion root of
`size(eta) = alpha + delta`, about 0.652 at `alpha = delta = 0.05`. The
empirical applications in the source papers (V-Dem democracy–growth,
castle-doctrine, minimum-wage and no-fault-divorce DiD) are reproduced
to reference tolerance by this package’s test suite, in lockstep with
the Julia reference implementation.

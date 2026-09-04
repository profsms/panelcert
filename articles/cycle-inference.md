# Concentrated design: diagnostics to exact inference

## Bundled application

Paper A’s return-panel application is bundled as `fscore`. The
specification uses firm and country-year fixed effects, so every number
below is reproducible without a download.

``` r

library(panelcert)
data(fscore)
cy <- interaction(fscore$country, fscore$year, drop = TRUE)
y <- log1p(fscore$ret)
x <- fscore$fscore
```

## Diagnose concentration

The design summary requires no outcome. `lambda_n` is the largest
observation’s share of within-treatment variation and `n_eff` is its
inverse-Herfindahl effective support.

``` r

d <- design_summary(fscore$uid, cy, x = x)
d
#> Panel Design Summary
#>   n = 217 obs | N = 19 units | T = 45 periods
#>   d_K = 61 (connected components: 3) | rho = d_K/n = 0.2811
#>   within variation tau*^2 = 461.694 | lambda_n = 0.0359037 | N_eff = 83.134
#>   NOTE: design is disconnected (3 components); within comparisons exist only inside each component.
```

The one-call screen adds the greedy lower bound and the designed
packing. It can be used for a prevalence exercise before outcomes are
collected.

``` r

screen <- adequacy_row(x, fscore$uid, cy, method = "structured")
screen[c("Vn", "lambda_n", "n_eff", "kappa_greedy",
         "kappa_designed", "se_price")]
#> $Vn
#> [1] 461.694
#> 
#> $lambda_n
#> [1] 0.03590372
#> 
#> $n_eff
#> [1] 83.1337
#> 
#> $kappa_greedy
#> [1] 0.3741656
#> 
#> $kappa_designed
#> [1] 0.8273879
#> 
#> $se_price
#> [1] 1.099374
```

## Capture and granularity

The structured packing uses disjoint four-cycles. `kappa` is the
observable share of within-treatment variation captured by its
annihilating contrasts; it is a Pitman efficiency only under Paper A’s
two additional power conditions.

``` r

cs <- cycle_contrasts(x, fscore$uid, cy, method = "structured")
c(C = length(cs$rows), kappa = cs$kappa,
  effective_C = length(panelcert:::.active_supports(cs)))
#>           C       kappa effective_C 
#>  43.0000000   0.8273879  40.0000000
```

Every contrast satisfies `v'D = 0`. Consequently `v'x_tilde = v'x`
exactly: packing is ranked on raw `x`, while iterative residualization
is used only for the denominator `V_n`. Both package implementations
regression-test this identity and the selected supports.

## Exact interval

The report adds realized score concentration and the
pass/flag/inconclusive verdict. The interval is then obtained by
sign-flip inversion. The reduced Monte Carlo and grid sizes below keep
vignette builds fast; replication runs should increase them.

``` r

report <- cycle_report(y, x, fscore$uid, cy, method = "structured",
                       nflips = 1999, interval = FALSE)
report
#> Panel Adequacy Report — Concentrated Identifying Variation
#> Design: n=217, N=19, T=45, d_K=61, rho=0.2811
#> Concentration: lambda_n = 0.0359 (N_eff = 83.1)
#> Realized score concentration: lambda_score = 0.1334 (N_eff,score = 26.6)
#> Capture kappa_C = 0.8274 over 43 supports (40 treatment-loaded) (cycle-space dim 156) | capture-implied SE ratio 1.099x | max share 0.065
#> VERDICT: POINT PASS at delta=0.05 (descriptive — not a certificate)
#> Diagnostic notes hidden (6); call show_notes(report) to display them.

ci <- signflip_interval(y, cs, alpha = 0.05, nflips = 1999,
                        ngrid = 2401)
unlist(ci[c("lo", "hi", "beta_tilde", "effective_C", "kappa")])
#>           lo           hi   beta_tilde  effective_C        kappa 
#> -0.040159313  0.026051189 -0.005736055 40.000000000  0.827387880
```

Exactness requires joint sign invariance at the dependence-block level
and supports that contain whole blocks. For clustered errors, pass
`blocks` to
[`cycle_contrasts()`](https://profsms.github.io/panelcert/reference/cycle_contrasts.md),
[`cycle_report()`](https://profsms.github.io/panelcert/reference/cycle_report.md),
and
[`signflip_interval()`](https://profsms.github.io/panelcert/reference/signflip_interval.md);
the package rejects an incompatible packing instead of silently
asserting exactness.

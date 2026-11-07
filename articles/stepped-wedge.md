# Concurrent controls and heterogeneity tolerance in stepped-wedge trials

This vignette recomputes the main design numbers of Halkiewicz and
Alexeev, *Concurrent controls and heterogeneity tolerance in
stepped-wedge trials* (arXiv:2610.35853).

A design is written as cluster counts for start periods `2, ..., T`,
followed by the number of wait-listed clusters (never treated during
follow-up). With 24 clusters and five periods, the balanced complete
wedge is `c(6, 6, 6, 6, 0)`.

``` r

library(panelcert)
designs <- list(complete = c(6, 6, 6, 6, 0), wait4 = c(5, 5, 5, 5, 4),
                wait8 = c(4, 4, 4, 4, 8), wait12 = c(3, 3, 3, 3, 12),
                capped_opt = c(6, 6, 1, 0, 11))
```

## Exposure weights and the adequacy radius

The conventional immediate-treatment (IT) estimator averages the
exposure-specific effects with weights `W`. The adequacy radius `kappa*`
is the RMS exposure-time heterogeneity, relative to the planned effect,
that the IT analysis tolerates before worst-case coverage of its 95%
interval falls to 90% (planned power 80%).

``` r

round(sw_exposure_weights(designs$complete, phi = 2), 3)
#>     W1     W2     W3     W4 
#>  0.943  0.300 -0.071 -0.171
sapply(designs, function(d) round(sw_radius(design = d, phi = 2), 3))
#>   complete      wait4      wait8     wait12 capped_opt 
#>      0.133      0.174      0.214      0.254      0.377
```

The complete wedge tolerates heterogeneity of only about 13% of the
planned effect, and wait-listing raises the radius.

## Cost of the estimand-aligned analysis

The exposure-time-indicator (ETI) analysis targets the time-averaged
treatment effect directly. Its cost is its variance relative to the IT
variance of the complete wedge.

``` r

t(sapply(designs, function(d) {
  ec <- sw_eti_cost(d, phi = 2)
  round(c(eti_cost = ec$eti_cost, it_ratio = ec$it_ratio), 2)
}))
#>            eti_cost it_ratio
#> complete       2.27     1.00
#> wait4          1.62     0.92
#> wait8          1.44     0.93
#> wait12         1.50     1.04
#> capped_opt     1.23     1.05
```

For a fixed rollout shape, the ETI variance first falls and then rises
with the wait-list share (Theorem 3 of the paper).

``` r

q <- seq(0, 0.8, by = 0.05)
cost <- sapply(q, function(x) sw_eti_cost(c(rep((1 - x) / 4, 4), x), phi = 2)$eti_cost)
round(cbind(q, cost), 2)[c(1, 5, 8, 9, 11, 17), ]
#>         q cost
#> [1,] 0.00 2.27
#> [2,] 0.20 1.56
#> [3,] 0.35 1.44
#> [4,] 0.40 1.44
#> [5,] 0.50 1.50
#> [6,] 0.80 2.67
```

## Rollout-capped design search

At most six clusters may start in any period. The frontier gives the
smallest ETI cost for each bound on the number wait-listed.

``` r

s <- sw_design_search(G = 24, T = 5, cap = 6, max_wait = 12, phi = 2)
s$frontier[s$frontier$max_wait %in% c(0, 4, 6, 8, 12), ]
#>    max_wait start2 start3 start4 start5 wait eti_cost  it_ratio    radius
#> 1         0      6      6      6      6    0 2.269942 1.0000000 0.1332467
#> 5         4      6      6      2      6    4 1.519707 0.9237537 0.1819003
#> 7         6      6      6      0      6    6 1.376921 0.9459459 0.2129412
#> 9         8      6      6      1      3    8 1.290928 0.9563567 0.2491485
#> 13       12      6      6      1      0   11 1.227424 1.0508757 0.3768735
```

## Nested-exchangeable correlation and unequal cluster sizes

Under nested-exchangeable correlation, every result holds with the
effective parameter `phi'`. With an ICC of 0.05, cluster autocorrelation
0.5 and 20 individuals per cluster-period:

``` r

sw_phi_nested(icc = 0.05, cac = 0.5, m = 20)
#> [1] 0.3448276
sw_radius(design = designs$complete, icc = 0.05, cac = 0.5, m = 20)
#> [1] 0.1653514
```

The bundled weekend allied-health trial (`haines2017`; CC BY 4.0) has
unequal ward-period sizes. Its IT exposure weights at the actual cell
sizes, under the fitted intracluster correlation 0.0709, give a radius
of 0.124.

``` r

n <- with(haines2017, tapply(rep(1, nrow(haines2017)), list(ward, period), sum))
starts <- tapply(haines2017$start, haines2017$ward, `[`, 1)
ord <- order(starts, as.numeric(names(starts)))
W <- sw_exposure_weights(c(2, 2, 2, 2, 2, 2, 0), icc = 0.0709, cac = 1,
                         sizes = n[ord, ])
round(W, 3)
#>     W1     W2     W3     W4     W5     W6 
#>  0.744  0.395  0.123 -0.038 -0.106 -0.118
round(sw_radius(W), 4)
#> [1] 0.1244
```

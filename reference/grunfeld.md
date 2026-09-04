# Canonical 11-firm Grunfeld investment panel

The complete corrected 1935–1954 panel used for Paper A's concentrated-
variation showcase. The canonical regression is investment on capital
and market value with firm and year fixed effects. Several 10- and
five-firm variants circulate; Kleiber and Zeileis (2010) document their
omissions and transcription errors.

## Usage

``` r
grunfeld
```

## Format

A data frame with 220 firm-year rows and 5 variables:

- invest:

  gross investment in 1947 dollars

- value:

  market value at year end in 1947 dollars

- capital:

  stock of plant and equipment in 1947 dollars

- firm:

  firm name

- year:

  calendar year, 1935–1954

## Source

Public-domain statsmodels Grunfeld data at commit
`57169f2cfc7089141513ed6103f79fa56ab213ac`, accessed 2026-08-02;
reconstructed from Grunfeld (1958) by Kleiber and Zeileis (2010).

## Examples

``` r
applicable(grunfeld$capital, grunfeld$firm, grunfeld$year,
           controls = grunfeld$value)
#> $ok
#> [1] TRUE
#> 
#> $reason
#> [1] "not a binary treatment; the binary granularity floor does not apply."
#> 
# \donttest{
cycle_report(grunfeld$invest, grunfeld$capital, grunfeld$firm,
             grunfeld$year, controls = grunfeld$value, interval = FALSE)
#> Panel Adequacy Report — Concentrated Identifying Variation
#> Design: n=220, N=11, T=20, d_K=30, rho=0.1364
#> Concentration: lambda_n = 0.2065 (N_eff = 16.6)
#> Realized score concentration: lambda_score = 0.7388 (N_eff,score = 1.8)
#> Capture kappa_C = 0.6270 over 32 supports (cycle-space dim 189) | capture-implied SE ratio 1.263x | max share 0.352
#> VERDICT: FLAGGED at delta=0.05
#> Diagnostic notes hidden (7); call show_notes(report) to display them.
# }
```

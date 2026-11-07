# Weekend allied-health stepped-wedge trial (Haines et al. 2017, Trial 1)

Analysis extract of Trial 1 of Haines et al. (2017): withdrawal of
weekend allied-health services from 12 acute wards in two hospitals, a
complete stepped wedge with 7 periods and 2 wards starting in each of
periods 2–7. One row per patient episode. Used by Halkiewicz and Alexeev
(arXiv:2610.35853) as the application with visible exposure-time
heterogeneity; the exposure weights of the IT estimator at the trial's
actual cell sizes give an adequacy radius of 0.124.

## Usage

``` r
haines2017
```

## Format

A data frame with 14,834 rows and 6 variables:

- ward:

  ward identifier (cluster)

- site:

  hospital

- period:

  study period 1–7

- start:

  first period in which the ward was treated (sequence)

- trt:

  1 if the weekend service had been withdrawn (intervention)

- acute_los:

  acute length of stay in days (analysed on the log scale)

## Source

Haines TP, Bowles K-A, Mitchell D, et al. (2017). Impact of
disinvestment from weekend allied health services across acute medical
and surgical wards: 2 stepped-wedge cluster randomised controlled
trials. *PLOS Medicine* 14(10): e1002412,
[doi:10.1371/journal.pmed.1002412](https://doi.org/10.1371/journal.pmed.1002412)
; S2 Data,
[doi:10.1371/journal.pmed.1002412.s002](https://doi.org/10.1371/journal.pmed.1002412.s002)
. Licence: CC BY 4.0 (not covered by the package's MIT licence). Only
the six variables above are retained; see `data-raw/make_sw_trials.py`
and `inst/DATA_SOURCES.md`.

## Examples

``` r
n <- with(haines2017, tapply(rep(1, nrow(haines2017)), list(ward, period), sum))
starts <- tapply(haines2017$start, haines2017$ward, `[`, 1)
ord <- order(starts, as.numeric(names(starts)))
sw_radius(design = c(2, 2, 2, 2, 2, 2, 0), icc = 0.0709, cac = 1,
          sizes = n[ord, ])
#> [1] 0.1243615
```

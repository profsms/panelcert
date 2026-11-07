# 99DOTS digital adherence stepped-wedge trial, Uganda (Cattamanchi et al. 2021)

Patient-level public data of the 99DOTS stepped-wedge trial: 18 health
facilities in six randomisation blocks, eight trial months, treatment
success as the primary outcome. Used by Halkiewicz and Alexeev
(arXiv:2610.35853) as the application with no detectable heterogeneity
but low power to rule out heterogeneity that would damage IT inference.

## Usage

``` r
dots99
```

## Format

A data frame with 2,224 rows and 13 variables:

- block:

  intervention group (randomisation block 1–6)

- facility:

  health facility 1–18 (cluster)

- facility_type:

  facility type

- month:

  trial month 0–7

- period:

  study period: Control, Buffer, or Intervention

- sex, age, hiv, disease_class, retreatment:

  patient characteristics

- treated_successfully, persisting, not_ltfu:

  outcomes (Yes/No)

## Source

Cattamanchi A, Crowder R, Kityamuwesi A, et al. (2021). Digital
adherence technology for tuberculosis treatment supervision: a
stepped-wedge cluster-randomized trial in Uganda. *PLOS Medicine* 18(5):
e1003628,
[doi:10.1371/journal.pmed.1003628](https://doi.org/10.1371/journal.pmed.1003628)
; S1 Data,
[doi:10.1371/journal.pmed.1003628.s002](https://doi.org/10.1371/journal.pmed.1003628.s002)
. Licence: CC BY 4.0 (not covered by the package's MIT licence). Column
names are snake_case; commas inside text labels are replaced by
semicolons. See `inst/DATA_SOURCES.md`.

## Examples

``` r
table(dots99$block, dots99$period)
#>    
#>     Buffer period Control Intervention
#>   1            47      30          192
#>   2            33      56          164
#>   3            64     150          231
#>   4            36     148          116
#>   5            69     303          136
#>   6            62     335           52
```

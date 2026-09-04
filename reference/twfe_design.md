# Pre-outcome TWFE design-statistic ladder

The design statistic `Gamma = sqrt(N1)||w - u||` and its restricted
variants `Gamma_coh`, `Gamma_evt`, `Gamma_c+e`, and `Gamma_gt` (Prop.
prop-restricted), the negative-weight share, and the breakdown ratio
from the adoption pattern alone. In a balanced panel,
`Gamma_gt = Gamma`. Always-treated units are dropped.

## Usage

``` r
twfe_design(unit, time, first_treat, alpha = 0.05, delta = 0.05)
```

## Arguments

- unit, time:

  raw identifiers (time numeric)

- first_treat:

  adoption time per observation (constant within unit; NA = never)

- alpha, delta:

  level and size tolerance

## Value

an `AdequacyReport`

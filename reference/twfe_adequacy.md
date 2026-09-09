# TWFE-heterogeneity adequacy (covariance-aware, wild-bootstrap)

Restricted design-statistic ladder, direct CR1 rescaling
`Gamma_{S,CR} = Gamma_S/sqrt(psi)`, covariance-aware pilots `c_S/sigma`,
saturated group-time and restricted point envelopes, and the signed
directional plug-in. A fixed-design cluster-multiplier radius gives the
regular one-sided lower bound. The upper certificate inverts the
noncentral chi-square law of the projected Wald statistic and is issued
only when the score covariance has full rank in the prespecified class.
Point-envelope bootstrap percentiles remain descriptive. Always-treated
units are dropped.

## Usage

``` r
# S3 method for class 'fixest'
twfe_adequacy(object, first_treat, ...)

# S3 method for class 'plm'
twfe_adequacy(object, first_treat, ...)

twfe_adequacy(object, ...)

# Default S3 method
twfe_adequacy(
  object,
  unit,
  time,
  first_treat,
  alpha = 0.05,
  delta = 0.05,
  cluster = c("direct", "ar1", "iid"),
  psi = NULL,
  controls = c("not_yet", "never"),
  heterogeneity_class = c("group_time", "additive", "cohort", "event"),
  bootstrap = 999L,
  seed = 20260715L,
  gamma = 0.05,
  q_band = NULL,
  ...
)
```

## Arguments

- object:

  outcome vector

- ...:

  passed between methods

- unit, time, first_treat:

  as in \[twfe_design()\]

- alpha, delta:

  level and size tolerance

- cluster:

  normalization: `"direct"` (default, the reported CR1 score scale),
  `"ar1"`, or `"iid"`

- psi:

  optional positive user-supplied variance-inflation factor; when
  supplied it overrides `cluster`

- controls:

  comparison group for group-time effects: not-yet-treated (including
  never-treated) or never-treated only

- heterogeneity_class:

  prespecified class used for the report verdict: saturated group-time,
  additive cohort-plus-event-time, cohort, or event-time

- bootstrap:

  number of wild-cluster draws (\>0 enables the covariance correction,
  descriptive point summaries, and one-sided procedures)

- seed:

  RNG seed for the wild bootstrap

- gamma:

  error probability for each reported one-sided procedure

- q_band:

  optional relative half-width for the denominator band in the more
  general full confidence-ball construction. It does not affect the main
  decision-specific verdict.

## Value

an `AdequacyReport`

# TWFE-heterogeneity adequacy (covariance-aware, wild-bootstrap)

Restricted design-statistic ladder, direct CR1 rescaling
`Gamma_{S,CR} = Gamma_S/sqrt(psi)`, covariance-aware pilots `c_S/sigma`,
saturated group-time and restricted point envelopes, and the signed
directional plug-in. A fixed-design cluster-multiplier radius yields
boundary-robust lower and upper bounds on the population envelope;
point-envelope bootstrap percentiles remain descriptive. Always-treated
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
  bootstrap = 999L,
  seed = 20260715L,
  gamma = 0.05,
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

- bootstrap:

  number of wild-cluster draws (\>0 enables the covariance correction,
  descriptive point summaries, and projected-norm bounds)

- seed:

  RNG seed for the wild bootstrap

- gamma:

  error probability for the projected-norm confidence set

## Value

an `AdequacyReport`

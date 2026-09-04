# TWFE design-statistic ladder

Convenience accessor returning the unrestricted, cohort, event-time,
combined additive, and saturated group-time design statistics from the
adoption pattern.

## Usage

``` r
twfe_gammas(unit, time, first_treat)
```

## Arguments

- unit, time:

  raw identifiers (time numeric)

- first_treat:

  adoption time per observation (constant within unit; NA = never)

## Value

A list with \`unr\`, \`coh\`, \`evt\`, \`cmb\`, \`gt\`, \`neg_share\`,
\`N1\`, and \`n_w\`.

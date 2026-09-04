# Pre-flight applicability check for Paper A

For a binary treatment with \`n1\` treated and \`n0\` untreated
observations, no more than \`min(n1, n0)\` disjoint contrasts can carry
a nonzero loading. A two-sided level-\`alpha\` test therefore requires
\`min(n1, n0) \>= 1 + log2(1/alpha)\`. Other designs pass this
structural check; packing quality is assessed separately by
\[cycle_report()\].

## Usage

``` r
applicable(x, unit, time, alpha = 0.05, controls = NULL)
```

## Arguments

- x:

  treatment or regressor vector.

- unit, time:

  fixed-effect identifiers.

- alpha:

  nominal test level.

- controls:

  optional numeric nuisance-covariate vector or matrix.

## Value

A named list with logical \`ok\` and a plain-language \`reason\`.

# Checkable conditions behind the cluster-robust layer

Assumption ass-cluster and Lemma lem-nest of the accompanying article.
All are computable from the design alone, before any outcome is
examined.

## Usage

``` r
cluster_diagnostics(xt, cluster, fe_levels, tau_star2)
```

## Arguments

- xt:

  residualized regressor.

- cluster:

  cluster identifier, one per observation.

- fe_levels:

  named list of fixed-effect id vectors, one per FE dimension.

- tau_star2:

  the within variation \`x\*'M x\*\`.

## Value

A list with \`G\`, \`max_size\`, \`min_size\`, \`max_energy\`, \`d_ne\`,
\`ratio_ne\` and \`nested\`.

## Details

Condition (iv) of ass-cluster – projection compatibility – is the one
with no counterpart in the i.i.d. theory, and it fails silently when the
fixed effects cut across clusters in a high-dimensional way.

Lemma lem-nest(b) bounds it by \`varpi_n \* max_g A_g\`, which is
unconditional. The familiar reduction to \`d_ne / G -\> 0\` is \*not\*:
it holds only under two balance conditions, which this function reports
rather than assumes.

- (N1) spectral: \`lambda_min^+(Lambda_n) ~ n / d_ne\`. Exact in a
  balanced two-way panel, where \`Lambda_n = G (I_T - 11'/T)\`, but
  equal cell counts alone do not control the smallest non-zero
  eigenvalue of the residualized Gram matrix. Not computed here.

- (N2) energy: \`max_g A_g = O(tau\*2 / G)\`, reported as `max_energy`.
  If the residualized signal concentrates in \`sqrt(G)\` clusters then
  \`max_energy ~ G^-1/2\` and the requirement becomes \`d_ne / sqrt(G)
  -\> 0\` instead.

So read `ratio_ne` and `max_energy` together: a small `ratio_ne` carries
no warrant on its own. When the two disagree, or when (N1) is in doubt,
use \[projection_compatibility()\], which evaluates \`sum_g \|\|M
a^(g) - a^(g)\|\|^2 / tau\*2\` directly and needs neither condition.

In the V-Dem application country effects nest in country clusters while
the 59 year effects do not, against \`G = 163\`, so \`ratio_ne ~ 0.36\`
and the direct diagnostic is load-bearing. The repeated-report twins
application samples independent pairs and therefore uses the baseline
i.i.d. theorem; it does not invoke this cluster condition or a \`psi\`
rescaling.

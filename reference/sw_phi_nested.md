# Effective Hussey–Hughes parameter under nested-exchangeable correlation

For cross-sectional sampling with `m` individuals in every
cluster-period, within-period intracluster correlation `icc` and cluster
autocorrelation `cac`, the covariance of a cluster's period means is
proportional to \\I + \phi' J\\ with \$\$\phi' = m\\r\\\rho /
\\m(1-r)\rho + 1 - \rho\\,\$\$ so every Hussey–Hughes result holds with
\\\phi\\ replaced by \\\phi'\\ (Corollary 1 of Halkiewicz and Alexeev).
At `cac = 1`, \\\phi' = m\rho/(1-\rho)\\.

## Usage

``` r
sw_phi_nested(icc, cac, m)
```

## Arguments

- icc:

  within-period intracluster correlation \\\rho \in \[0, 1)\\

- cac:

  cluster autocorrelation \\r \in \[0, 1\]\\

- m:

  individuals per cluster-period

## Value

the effective parameter \\\phi'\\ (vectorised)

## Examples

``` r
sw_phi_nested(icc = 0.05, cac = 0.5, m = 20)
#> [1] 0.3448276
sw_phi_nested(icc = 0.05, cac = 1, m = 20)   # = 20 * 0.05 / 0.95
#> [1] 1.052632
```

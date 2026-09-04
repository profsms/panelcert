# Finite-n non-centrality mapping

Paper B eq-eta-n under the local drift \`sigma_nu^2 = c^2/n\`: \`eta_n =
-beta0 c^2 sqrt(1-rho) / (sigma sqrt(tau^2 + c^2))\`. The
\`sqrt(1-rho)\` factor is the common-scale shrinkage implied by \`X'M_K
X/(nQ_K) -\>p 1-rho\`; it was absent from the pre-correction formula.

## Usage

``` r
eta_finite_n(beta0, c2, sigma, tau2, rho)
```

## Arguments

- beta0, c2, sigma, tau2, rho:

  model primitives.

## Value

The non-centrality \`eta_n\`.

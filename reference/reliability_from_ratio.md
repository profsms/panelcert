# Measurement-error SD implied by an external reliability ratio

With the documented default \`scale = "observed"\`, \`within_sd\` is the
SD of the observed regressor and \`sigma_nu = sqrt(1-r) \* within_sd\`.
Set \`scale = "signal"\` only when \`within_sd\` is the latent-signal
SD; then \`sigma_nu = sqrt((1-r)/r) \* within_sd\`.

## Usage

``` r
reliability_from_ratio(r, within_sd, scale = c("observed", "signal"))
```

## Arguments

- r:

  reliability ratio in (0, 1\]

- within_sd:

  within-SD on the scale selected by \`scale\`

- scale:

  \`"observed"\` (default) or \`"signal"\`

## Value

scalar measurement-error SD

## Examples

``` r
reliability_from_ratio(0.75, within_sd = 0.15)
#> [1] 0.075
```

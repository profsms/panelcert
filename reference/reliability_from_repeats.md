# Reliability from two measurements of the same regressor

Apply the identical projection and complete-case sample to both
measurements before calling this helper. The default covariance
estimator is \`cov(first_measure, second_measure) / var(first_measure)\`
and identifies the first measurement's reliability when the two
classical reporting errors are uncorrelated; it does not require equal
error variances. The \`"equal_variance"\` estimator is \`1 -
var(first_measure - second_measure) / (2 \* var(first_measure))\` and
additionally imposes equal reporting-error variances.

## Usage

``` r
reliability_from_repeats(
  first_measure,
  second_measure,
  method = c("covariance", "equal_variance")
)
```

## Arguments

- first_measure, second_measure:

  numeric vectors containing two projected measurements on the same
  complete-case sample.

- method:

  \`"covariance"\` (default) or \`"equal_variance"\`.

## Value

A scalar estimate of the first measurement's reliability.

## Details

Correlation across reporting errors invalidates both formulas and must
be modeled or examined as a sensitivity. The estimate is returned
without truncation so assumption or sampling failures remain visible;
\[eiv_adequacy()\] requires a reliability in \`(0, 1\]\`.

## Examples

``` r
x <- c(-2, -1, 1, 2)
z <- c(-1.8, -1.2, 0.9, 2.1)
reliability_from_repeats(x, z)
#> [1] 0.99
reliability_from_repeats(x, z, method = "equal_variance")
#> [1] 0.995
```

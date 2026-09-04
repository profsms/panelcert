# Ashenfelter–Krueger repeated-report twins extract

Public teaching extract from the Ashenfelter–Krueger identical-twins
design. Each twin reports both siblings' schooling, producing two
measurements of the within-pair schooling difference. The article's
controlled calculation uses the 147 complete rows on \`DLHRWAGE\`,
\`DEDUC1\`, \`DEDUC2\`, \`DTEN\`, \`DMARRIED\`, and \`DUNCOV\`, applying
the identical control projection and sample to both schooling reports.

## Usage

``` r
twins
```

## Format

A data frame with 183 twin-pair rows and 16 variables:

- DLHRWAGE:

  within-pair difference in log hourly wages

- DEDUC1:

  schooling difference from each twin's self-report

- AGE:

  age of the twin pair

- AGESQ:

  squared age

- HRWAGEH:

  hourly wage of the H-labelled twin

- WHITEH:

  white indicator for the H-labelled twin

- MALEH:

  male indicator for the H-labelled twin

- EDUCH:

  years of schooling of the H-labelled twin

- HRWAGEL:

  hourly wage of the L-labelled twin

- WHITEL:

  white indicator for the L-labelled twin

- MALEL:

  male indicator for the L-labelled twin

- EDUCL:

  years of schooling of the L-labelled twin

- DEDUC2:

  schooling difference from the co-twin reports

- DTEN:

  within-pair tenure difference

- DMARRIED:

  within-pair married-status difference

- DUNCOV:

  within-pair union-coverage difference

## Source

\`RbyExample::twins\` version 0.0.100, attributed to Ashenfelter and
Krueger (1994), \*American Economic Review\* 84(5), 1157–1173. Upstream
package license: GPL (\>= 2).

## Examples

``` r
keep <- stats::complete.cases(twins[c("DLHRWAGE", "DEDUC1", "DEDUC2",
                                      "DTEN", "DMARRIED", "DUNCOV")])
W <- stats::model.matrix(~ DTEN + DMARRIED + DUNCOV, data = twins[keep, ])
x <- qr.resid(qr(W), twins$DEDUC1[keep])
z <- qr.resid(qr(W), twins$DEDUC2[keep])
reliability_from_repeats(x, z)
#> [1] 0.5749037
```

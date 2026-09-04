# Piotroski F-Score / Visegrad firm panel (Paper A application)

Hand-collected firm–year panel of Piotroski F-Scores and one-year-ahead
returns for firms on the Warsaw, Budapest and Prague exchanges,
2010–2024, consolidated from the three exchange production files. Paper
A uses the firm and country–year specification for a complete
exact-inference workflow. The diffuse-regime variance companion also
uses seven outcome/subsample/FE variants to compare HC0–HC3 behavior.

## Usage

``` r
fscore
```

## Format

A data frame with 217 firm-year rows and 6 variables:

- uid:

  firm ticker (unit id)

- year:

  fiscal year (time id)

- country:

  Poland, Hungary or Czech Republic

- status:

  `"active"` or `"delisted"`

- fscore:

  Piotroski F-Score, 0–9 integer composite (the regressor)

- ret:

  one-year-ahead return (outcome)

## Source

Hand-collected from Warsaw (WSE), Budapest (BSE) and Prague (PSE)
exchange filings; the panel of Paper A's empirical application.

## Examples

``` r
# Paper A specification: log return on F-Score, firm + country-year effects
cy <- interaction(fscore$country, fscore$year, drop = TRUE)
cycle_capture(fscore$fscore, fscore$uid, cy)
#> [1] 0.8273879
```

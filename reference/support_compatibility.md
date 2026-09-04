# Check dependence-block compatibility of contrast supports

Paper A Assumption ass:sym(iii): each treatment-loaded support must
contain all or none of every dependence block, and a block may not be
split across supports. Zero-loading supports are irrelevant to the test
statistic.

## Usage

``` r
support_compatibility(cs, blocks)
```

## Arguments

- cs:

  a \`CycleSystem\`.

- blocks:

  one block id per observation.

## Value

list with \`compatible\`, \`incompatible_blocks\`, \`nblocks\`, and
\`effective_C\`.

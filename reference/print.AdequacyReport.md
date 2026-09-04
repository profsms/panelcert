# Print an adequacy report

Renders the plain-language verdict block (spec section 2.2): design
line, pathology-specific statistics, non-centrality vs threshold,
implied size, and the VERDICT. Cycle-inference caveats remain in
\`x\$notes\` and are shown only when requested.

## Usage

``` r
# S3 method for class 'AdequacyReport'
print(x, notes = x$pathology != "cycle_inference", ...)
```

## Arguments

- x:

  an `AdequacyReport`

- notes:

  logical; whether to print detailed diagnostic notes. Defaults to
  \`FALSE\` for cycle-inference reports and \`TRUE\` for other reports.

- ...:

  unused

## Value

`x`, invisibly

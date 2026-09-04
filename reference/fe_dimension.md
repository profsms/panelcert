# Fixed-effect dimension of a two-way design

Union-find on the bipartite unit-time graph: \`d_K = N + T - \#connected
components\`.

## Usage

``` r
fe_dimension(uid, tid, N, T)
```

## Arguments

- uid, tid:

  integer codes 1..N / 1..T (one per observation)

- N, T:

  number of units / periods

## Value

list with \`d_K\` and \`ncomponents\`

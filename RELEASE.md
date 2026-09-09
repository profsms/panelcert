# Release checklist

1.  Confirm R-CMD-check, coverage, pkgdown, and cross-language parity
    are green on `main`.
2.  Confirm `R CMD check --as-cran` reports zero errors and warnings.
    Before the first CRAN acceptance, the sole permitted note is CRAN’s
    administrative `New submission`; subsequent releases must report
    zero notes.
3.  Add `https://github.com/profsms/panelcert` to the `packages.json`
    file in the `profsms` r-universe configuration repository.
4.  Tag and push `v0.8.0`; create a GitHub release from `NEWS.md`.
5.  Verify the r-universe binary and source builds, then make its
    install command primary in the README.
6.  Submit to CRAN after the public r-universe build and
    reverse-dependency check are clean.

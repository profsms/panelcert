# panelcert: Inference-Adequacy Diagnostics for Panel-Data Designs

Implements exact nuisance-annihilating contrast inference for
fixed-effect designs with concentrated identifying variation, and
diagnoses when conventional inference is unreliable because of leverage,
classical measurement error, or staggered-DiD treatment-effect
heterogeneity. Reports design and realized-score concentration, capture,
exact sign-flip tests and confidence sets, HC0-HC3 behavior,
exact-normal measurement-error thresholds, and TWFE design exposure,
direct cluster-score normalization, saturated group-time envelopes,
regular multiplier lower tests, covariance-aware projected-Wald upper
certificates, class-specific verdicts, and directional diagnostics.

## Details

Inference and diagnostics: \[cycle_report()\] and \[contrast_system()\]
(Paper A, exact contrast inference), \[eiv_adequacy()\] (Paper B,
measurement error), \[twfe_design()\] / \[twfe_adequacy()\]
(staggered-DiD heterogeneity), and \[leverage_report()\] (the
diffuse-regime variance companion). Shared primitives include
\[design_summary()\], \[twoway_demean()\], \[multiway_demean()\], and
\[fe_leverage()\]. Reports use a plain-language verdict; see
\`vignette("panelcert")\` and \`vignette("cycle-inference")\`.

## See also

Useful links:

- <https://github.com/profsms/panelcert>

- Report bugs at <https://github.com/profsms/panelcert/issues>

## Author

**Maintainer**: Stanislaw M. S. Halkiewicz <stashal@o2.pl>
([ORCID](https://orcid.org/0009-0000-7344-7522))

Authors:

- Stanislaw M. S. Halkiewicz <stashal@o2.pl>
  ([ORCID](https://orcid.org/0009-0000-7344-7522))

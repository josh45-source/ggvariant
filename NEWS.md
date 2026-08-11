# ggvariant (development version)

* `plot_variant_spectrum()`'s `context` and `genome` arguments, and
  `read_vcf()`'s `info_fields` argument, now abort with an informative
  error instead of being silently ignored. None of the three were actually
  implemented, despite being documented (#1, #2).
* Removed stray scratch files (`Mutation 1.pdf`, `Rplot*.png`/`.pdf`, etc.)
  that were shipping inside the package tarball and triggering an
  `R CMD check` WARNING and NOTE.
* Added `R-CMD-check` and `test-coverage` GitHub Actions workflows, and
  fixed the Codecov badge, which pointed at a nonexistent `main` branch.

# ggvariant 0.1.0

* Initial CRAN release.
* `read_vcf()` and `coerce_variants()` read VCF files or plain data frames
  into a tidy `gvf` object.
* `plot_lollipop()`, `plot_consequence_summary()`, and
  `plot_variant_spectrum()` provide `ggplot2`-native variant visualizations.
* `gv_palette()` and `theme_ggvariant()` provide built-in colour palettes
  and a shared plot theme.

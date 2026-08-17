# Changelog

## ggvariant 0.2.0

### New features

- Added
  [`plot_oncoprint()`](https://josh45-source.github.io/ggvariant/reference/plot_oncoprint.md)
  /
  [`plot_waterfall()`](https://josh45-source.github.io/ggvariant/reference/plot_oncoprint.md):
  a gene-by-sample mutation matrix (the same visualisation under the
  `ComplexHeatmap`/cBioPortal and `GenVisR`/`maftools` community names,
  respectively). Genes are ranked by descending number of altered
  samples; samples are ordered by the memo-sort/cascade algorithm for
  the characteristic staircase pattern. Multi-hit cells (a gene mutated
  more than once in one sample) render as a distinct `"Multi_Hit"`
  category. Supports an optional clinical annotation track from a
  sample-metadata data frame. Suggested by Nour-al-dain Marzouka (#N).
- Added
  [`plot_tmb()`](https://josh45-source.github.io/ggvariant/reference/plot_tmb.md):
  a per-sample tumour mutational burden bar chart, raw counts or
  normalised to mutations/Mb.
- [`read_vcf()`](https://josh45-source.github.io/ggvariant/reference/read_vcf.md)’s
  `gvf` objects now have real
  [`print.gvf()`](https://josh45-source.github.io/ggvariant/reference/gvf-methods.md)
  and
  [`summary.gvf()`](https://josh45-source.github.io/ggvariant/reference/gvf-methods.md)
  methods, showing a compact header and consequence/ sample/chromosome
  breakdowns instead of a raw `data.frame` dump.

### Bug fixes

- **[`read_vcf()`](https://josh45-source.github.io/ggvariant/reference/read_vcf.md)
  no longer attaches a variant to a homozygous-reference (`"0/0"`)
  sample.** `.pivot_samples()`’s presence check only excluded literal
  missing-genotype codes, so every sample was pivoted in as carrying
  every variant regardless of its actual genotype on any multi-sample
  VCF. This affected every per-sample dimension of
  [`read_vcf()`](https://josh45-source.github.io/ggvariant/reference/read_vcf.md)’s
  output: `plot_lollipop(color_by = "sample")`,
  [`plot_consequence_summary()`](https://josh45-source.github.io/ggvariant/reference/plot_consequence_summary.md)’s
  per-sample bars, and both new sample-aware plots above. Predates this
  release entirely (present since 0.1.0)
  ([\#3](https://github.com/josh45-source/ggvariant/issues/3)).
- [`plot_variant_spectrum()`](https://josh45-source.github.io/ggvariant/reference/plot_variant_spectrum.md)’s
  `context`/`genome` arguments and
  [`read_vcf()`](https://josh45-source.github.io/ggvariant/reference/read_vcf.md)’s
  `info_fields` argument now abort with an informative error instead of
  being silently ignored. None of the three were actually implemented,
  despite being documented
  ([\#1](https://github.com/josh45-source/ggvariant/issues/1),
  [\#2](https://github.com/josh45-source/ggvariant/issues/2)).
- A VCF missing its `#CHROM` header line, or containing a data line with
  the wrong number of tab-separated fields, now aborts with a clear
  message instead of failing deep inside on an opaque base-R error.
- A nonexistent file path now aborts with a clear message instead of
  [`normalizePath()`](https://rdrr.io/r/base/normalizePath.html)’s
  base-R error.

### Performance

- [`read_vcf()`](https://josh45-source.github.io/ggvariant/reference/read_vcf.md)’s
  `.parse_ann_csq()` (ANN/CSQ INFO field parsing) is now vectorised over
  the whole INFO column instead of looping row by row. ~8.75x faster in
  isolation on a synthetic 100,000-record VCF (8.05s -\> 0.92s); full
  [`read_vcf()`](https://josh45-source.github.io/ggvariant/reference/read_vcf.md)
  on the same file drops from 11.62s to 3.68s.

### Documentation

- Added a pkgdown site. The introductory vignette was renamed to
  `ggvariant.Rmd` to activate the “Get started” navbar convention, and
  three new articles (Gallery, interactive plots, customising with
  ggplot2) live in `vignettes/articles/` rather than the shipped
  vignettes, per R Packages’ guidance for graphics-heavy content.
- Removed stray scratch files (`Mutation 1.pdf`, `Rplot*.png`/`.pdf`,
  etc.) that were shipping inside the package tarball and triggering an
  `R CMD check` WARNING and NOTE.
- Added `R-CMD-check`, `test-coverage`, and `pkgdown` GitHub Actions
  workflows, and fixed the Codecov badge, which pointed at a nonexistent
  `main` branch.
- Test coverage increased from 68.12% to over 92%, including the first
  tests of malformed/edge-case input and the first `vdiffr` visual
  regression snapshots.

## ggvariant 0.1.0

CRAN release: 2026-02-27

- Initial CRAN release.
- [`read_vcf()`](https://josh45-source.github.io/ggvariant/reference/read_vcf.md)
  and
  [`coerce_variants()`](https://josh45-source.github.io/ggvariant/reference/coerce_variants.md)
  read VCF files or plain data frames into a tidy `gvf` object.
- [`plot_lollipop()`](https://josh45-source.github.io/ggvariant/reference/plot_lollipop.md),
  [`plot_consequence_summary()`](https://josh45-source.github.io/ggvariant/reference/plot_consequence_summary.md),
  and
  [`plot_variant_spectrum()`](https://josh45-source.github.io/ggvariant/reference/plot_variant_spectrum.md)
  provide `ggplot2`-native variant visualizations.
- [`gv_palette()`](https://josh45-source.github.io/ggvariant/reference/gv_palette.md)
  and
  [`theme_ggvariant()`](https://josh45-source.github.io/ggvariant/reference/theme_ggvariant.md)
  provide built-in colour palettes and a shared plot theme.

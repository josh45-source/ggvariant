# Package index

## Package overview

- [`ggvariant-package`](https://josh45-source.github.io/ggvariant/reference/ggvariant-package.md)
  [`ggvariant`](https://josh45-source.github.io/ggvariant/reference/ggvariant-package.md)
  : ggvariant: Tidy, ggplot2-Native Visualization for Genomic Variants

## Reading and coercing variants

Get variant data into the tidy `gvf` format every `ggvariant` plot
function accepts.

- [`read_vcf()`](https://josh45-source.github.io/ggvariant/reference/read_vcf.md)
  : Read a VCF file into a tidy variant data frame

- [`coerce_variants()`](https://josh45-source.github.io/ggvariant/reference/coerce_variants.md)
  :

  Coerce a plain data frame to a `gvf` object

- [`print(`*`<gvf>`*`)`](https://josh45-source.github.io/ggvariant/reference/gvf-methods.md)
  [`summary(`*`<gvf>`*`)`](https://josh45-source.github.io/ggvariant/reference/gvf-methods.md)
  [`print(`*`<summary.gvf>`*`)`](https://josh45-source.github.io/ggvariant/reference/gvf-methods.md)
  :

  Print and summarise `gvf` objects

## Plots

ggplot2-native plotting functions, one per visualisation type.

- [`plot_lollipop()`](https://josh45-source.github.io/ggvariant/reference/plot_lollipop.md)
  : Lollipop plot of variants along a gene
- [`plot_consequence_summary()`](https://josh45-source.github.io/ggvariant/reference/plot_consequence_summary.md)
  : Consequence summary bar chart
- [`plot_variant_spectrum()`](https://josh45-source.github.io/ggvariant/reference/plot_variant_spectrum.md)
  : Mutational spectrum (SBS) bar chart
- [`plot_oncoprint()`](https://josh45-source.github.io/ggvariant/reference/plot_oncoprint.md)
  [`plot_waterfall()`](https://josh45-source.github.io/ggvariant/reference/plot_oncoprint.md)
  : Oncoprint / waterfall plot of variants across samples
- [`plot_tmb()`](https://josh45-source.github.io/ggvariant/reference/plot_tmb.md)
  : Tumour mutational burden bar chart

## Palettes and themes

Built-in colour palettes and the shared plot theme, exported for reuse
when customising a plot’s `ggplot2` layers directly.

- [`gv_palette()`](https://josh45-source.github.io/ggvariant/reference/gv_palette.md)
  : ggvariant colour palettes
- [`theme_ggvariant()`](https://josh45-source.github.io/ggvariant/reference/theme_ggvariant.md)
  : ggvariant ggplot2 theme

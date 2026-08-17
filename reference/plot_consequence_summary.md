# Consequence summary bar chart

Summarises variant consequences (e.g. missense, frameshift, synonymous)
across one or more samples, producing a stacked or grouped bar chart.

## Usage

``` r
plot_consequence_summary(
  variants,
  samples = NULL,
  group_by = c("consequence", "gene"),
  top_n = 10L,
  position = c("stack", "fill", "dodge"),
  palette = NULL,
  flip = FALSE,
  interactive = FALSE
)
```

## Arguments

- variants:

  A `gvf` object from
  [`read_vcf()`](https://josh45-source.github.io/ggvariant/reference/read_vcf.md)
  or
  [`coerce_variants()`](https://josh45-source.github.io/ggvariant/reference/coerce_variants.md),
  or any `data.frame` with columns `pos`, `consequence`, and optionally
  `gene` and `sample`.

- samples:

  Character vector of sample names to include. `NULL` (default) uses all
  samples. Ignored if there is no `sample` column.

- group_by:

  `"consequence"` (default) stacks bars by consequence per sample;
  `"gene"` stacks by gene per consequence.

- top_n:

  Integer. For `group_by = "gene"`, show only the top N genes by total
  variant count. Default `10`.

- position:

  `"stack"` (default) or `"fill"` (proportional) or `"dodge"`.

- palette:

  Named character vector of colours for each consequence/sample
  category. `NULL` uses the built-in `ggvariant` palette.

- flip:

  Logical. If `TRUE`, flips coordinates for horizontal bars. Default
  `FALSE`.

- interactive:

  Logical. If `TRUE`, returns a `plotly` interactive plot (requires the
  `plotly` package).

## Value

A `ggplot` object.

## See also

[`plot_lollipop()`](https://josh45-source.github.io/ggvariant/reference/plot_lollipop.md),
[`plot_variant_spectrum()`](https://josh45-source.github.io/ggvariant/reference/plot_variant_spectrum.md),
[`gv_palette()`](https://josh45-source.github.io/ggvariant/reference/gv_palette.md)

Other ggvariant plots:
[`plot_lollipop()`](https://josh45-source.github.io/ggvariant/reference/plot_lollipop.md),
[`plot_oncoprint()`](https://josh45-source.github.io/ggvariant/reference/plot_oncoprint.md),
[`plot_tmb()`](https://josh45-source.github.io/ggvariant/reference/plot_tmb.md),
[`plot_variant_spectrum()`](https://josh45-source.github.io/ggvariant/reference/plot_variant_spectrum.md)

## Examples

``` r
vcf_file <- system.file("extdata", "example.vcf", package = "ggvariant")
variants <- read_vcf(vcf_file)
#> ℹ Reading VCF: example.vcf
#> ✔ Reading VCF: example.vcf [10ms]
#> 
#> Loaded 19 variant records across 7 chromosomes.

# Consequence counts per sample
plot_consequence_summary(variants)


# Proportional bars
plot_consequence_summary(variants, position = "fill")


# Top 10 genes coloured by consequence
plot_consequence_summary(variants, group_by = "gene", top_n = 10)

```

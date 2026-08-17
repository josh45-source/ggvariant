# Tumour mutational burden bar chart

Plots per-sample tumour mutational burden (TMB): the number of mutations
carried by each sample, optionally normalised to mutations per megabase.

## Usage

``` r
plot_tmb(
  variants,
  samples = NULL,
  mb_size = NULL,
  sample_order = NULL,
  palette = NULL,
  interactive = FALSE
)
```

## Arguments

- variants:

  A `gvf` object from
  [`read_vcf()`](https://josh45-source.github.io/ggvariant/reference/read_vcf.md)
  or
  [`coerce_variants()`](https://josh45-source.github.io/ggvariant/reference/coerce_variants.md),
  or any `data.frame` with `sample` and `pos` columns.

- samples:

  Character vector of sample names to include. `NULL` (default) uses
  every sample present in `variants`.

- mb_size:

  Numeric. Size, in megabases, of the sequenced region used to normalise
  mutation counts into mutations/Mb. `NULL` (default) plots raw mutation
  counts instead.

- sample_order:

  Character vector giving an explicit left-to-right sample order, e.g.
  to align with
  [`plot_oncoprint()`](https://josh45-source.github.io/ggvariant/reference/plot_oncoprint.md)'s
  memo-sorted column order. `NULL` (default) orders samples by
  descending TMB.

- palette:

  Single hex colour string for the bars. `NULL` uses a built-in default.

- interactive:

  Logical. Returns a `plotly` object if `TRUE`.

## Value

A `ggplot` object (or a `plotly` object when `interactive = TRUE`).

## Details

The per-sample mutation count is computed by an internal helper shared
with any future
[`plot_oncoprint()`](https://josh45-source.github.io/ggvariant/reference/plot_oncoprint.md)
TMB marginal, so the two always agree. Pass `mb_size` (the size, in
megabases, of the sequenced region) to convert raw counts into
mutations/Mb, the conventional TMB unit; leave it `NULL` to plot raw
counts.

## References

Chalmers ZR, Connelly CF, Fabrizio D, et al. (2017). Analysis of 100,000
human cancer genomes reveals the landscape of tumor mutational burden.
*Genome Medicine*, 9(1), 34.
[doi:10.1186/s13073-017-0424-2](https://doi.org/10.1186/s13073-017-0424-2)

## See also

[`plot_oncoprint()`](https://josh45-source.github.io/ggvariant/reference/plot_oncoprint.md),
[`gv_palette()`](https://josh45-source.github.io/ggvariant/reference/gv_palette.md)

Other ggvariant plots:
[`plot_consequence_summary()`](https://josh45-source.github.io/ggvariant/reference/plot_consequence_summary.md),
[`plot_lollipop()`](https://josh45-source.github.io/ggvariant/reference/plot_lollipop.md),
[`plot_oncoprint()`](https://josh45-source.github.io/ggvariant/reference/plot_oncoprint.md),
[`plot_variant_spectrum()`](https://josh45-source.github.io/ggvariant/reference/plot_variant_spectrum.md)

## Examples

``` r
vcf_file <- system.file("extdata", "example.vcf", package = "ggvariant")
variants <- read_vcf(vcf_file)
#> ℹ Reading VCF: example.vcf
#> ✔ Reading VCF: example.vcf [13ms]
#> 
#> Loaded 19 variant records across 7 chromosomes.

# Raw mutation counts, samples ordered by descending TMB
plot_tmb(variants)


# Normalised to mutations/Mb for a 38 Mb exome
plot_tmb(variants, mb_size = 38)


# Aligned to a specific sample order, e.g. from plot_oncoprint()
plot_tmb(variants, sample_order = c("TUMOR_S1", "TUMOR_S2"))

```

# Mutational spectrum (SBS) bar chart

Plots the single-base substitution (SBS) spectrum — the relative
frequency of each of the 6 substitution classes (C\>A, C\>G, C\>T, T\>A,
T\>C, T\>G) — optionally broken down by trinucleotide context.

## Usage

``` r
plot_variant_spectrum(
  variants,
  sample = NULL,
  context = FALSE,
  genome = NULL,
  facet_by_sample = FALSE,
  palette = NULL,
  normalize = TRUE,
  interactive = FALSE
)
```

## Arguments

- variants:

  A `gvf` object or compatible `data.frame` containing SNVs. Indels are
  automatically excluded.

- sample:

  Character. Sample name to filter on. `NULL` uses all variants pooled
  (or facets by sample if `facet_by_sample = TRUE`).

- context:

  Logical. 96-trinucleotide context bars are not yet implemented;
  passing `TRUE` aborts with an error. Default `FALSE`, which produces
  the 6-class SBS spectrum. See
  <https://github.com/josh45-source/ggvariant/issues/1>.

- genome:

  Not yet implemented; passing a non-`NULL` value aborts with an error.
  Reserved for future `BSgenome`-based trinucleotide context extraction.
  See <https://github.com/josh45-source/ggvariant/issues/1>.

- facet_by_sample:

  Logical. If `TRUE`, facets the plot by sample. Default `FALSE`.

- palette:

  Named character vector with names matching substitution classes
  (`"C>A"`, `"C>G"`, etc.). `NULL` uses COSMIC-style colours.

- normalize:

  Logical. If `TRUE` (default), shows relative proportions. If `FALSE`,
  shows raw counts.

- interactive:

  Logical. If `TRUE`, returns a `plotly` interactive plot (requires the
  `plotly` package).

## Value

A `ggplot` object.

## References

Alexandrov LB, Kim J, Haradhvala NJ, et al.; PCAWG Consortium (2020).
The repertoire of mutational signatures in human cancer. *Nature*,
578(7793), 94-101.
[doi:10.1038/s41586-020-1943-3](https://doi.org/10.1038/s41586-020-1943-3)

Blokzijl F, Janssen R, van Boxtel R, Cuppen E (2018).
MutationalPatterns: comprehensive genome-wide analysis of mutational
processes. *Genome Medicine*, 10(1), 33.
[doi:10.1186/s13073-018-0539-0](https://doi.org/10.1186/s13073-018-0539-0)

## See also

[`plot_lollipop()`](https://josh45-source.github.io/ggvariant/reference/plot_lollipop.md),
[`plot_consequence_summary()`](https://josh45-source.github.io/ggvariant/reference/plot_consequence_summary.md),
[`gv_palette()`](https://josh45-source.github.io/ggvariant/reference/gv_palette.md)

Other ggvariant plots:
[`plot_consequence_summary()`](https://josh45-source.github.io/ggvariant/reference/plot_consequence_summary.md),
[`plot_lollipop()`](https://josh45-source.github.io/ggvariant/reference/plot_lollipop.md),
[`plot_oncoprint()`](https://josh45-source.github.io/ggvariant/reference/plot_oncoprint.md),
[`plot_tmb()`](https://josh45-source.github.io/ggvariant/reference/plot_tmb.md)

## Examples

``` r
vcf_file <- system.file("extdata", "example.vcf", package = "ggvariant")
variants <- read_vcf(vcf_file)
#> ℹ Reading VCF: example.vcf
#> ✔ Reading VCF: example.vcf [10ms]
#> 
#> Loaded 19 variant records across 7 chromosomes.

# Basic 6-class SBS spectrum
plot_variant_spectrum(variants)
#> Excluded 2 non-SNV records from spectrum plot.


# Faceted by sample
plot_variant_spectrum(variants, facet_by_sample = TRUE)
#> Excluded 2 non-SNV records from spectrum plot.

```

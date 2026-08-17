# Read a VCF file into a tidy variant data frame

Parses a standard VCF (v4.x) file and returns a tidy `data.frame` (a
`gvf` object) that all `ggvariant` plotting functions accept. For users
who already have variant data in a plain `data.frame` or `tibble`, see
[`coerce_variants()`](https://josh45-source.github.io/ggvariant/reference/coerce_variants.md).

## Usage

``` r
read_vcf(path, samples = NULL, pass_only = TRUE, info_fields = NULL)
```

## Arguments

- path:

  Path to a `.vcf` or `.vcf.gz` file.

- samples:

  Character vector of sample names to retain. `NULL` (default) keeps all
  samples.

- pass_only:

  Logical. If `TRUE` (default), only variants with `FILTER` equal to
  `"PASS"` or `"."` are retained.

- info_fields:

  Not yet implemented; passing a non-`NULL` value aborts with an error.
  Reserved for future INFO field expansion. See
  <https://github.com/josh45-source/ggvariant/issues/2>.

## Value

A `gvf` (genomic variant frame) — a `data.frame` with columns:

- chrom:

  Chromosome (character)

- pos:

  Position (integer)

- ref:

  Reference allele

- alt:

  Alternate allele (multi-allelic sites are split into rows)

- qual:

  QUAL score (numeric)

- filter:

  FILTER field

- sample:

  Sample name (NA for single-sample VCFs without GT field)

- consequence:

  Variant consequence if ANN/CSQ INFO field is present

- gene:

  Gene symbol if ANN/CSQ INFO field is present

## References

Danecek P, Auton A, Abecasis G, et al.; 1000 Genomes Project Analysis
Group (2011). The variant call format and VCFtools. *Bioinformatics*,
27(15), 2156-2158.
[doi:10.1093/bioinformatics/btr330](https://doi.org/10.1093/bioinformatics/btr330)

Cingolani P, Platts A, Wang LL, et al. (2012). A program for annotating
and predicting the effects of single nucleotide polymorphisms, SnpEff:
SNPs in the genome of *Drosophila melanogaster* strain w1118; iso-2;
iso-3. *Fly*, 6(2), 80-92.
[doi:10.4161/fly.19695](https://doi.org/10.4161/fly.19695)

McLaren W, Gil L, Hunt SE, et al. (2016). The Ensembl Variant Effect
Predictor. *Genome Biology*, 17(1), 122.
[doi:10.1186/s13059-016-0974-4](https://doi.org/10.1186/s13059-016-0974-4)

## See also

[`coerce_variants()`](https://josh45-source.github.io/ggvariant/reference/coerce_variants.md),
[`plot_lollipop()`](https://josh45-source.github.io/ggvariant/reference/plot_lollipop.md),
[`plot_consequence_summary()`](https://josh45-source.github.io/ggvariant/reference/plot_consequence_summary.md)

Other ggvariant input:
[`coerce_variants()`](https://josh45-source.github.io/ggvariant/reference/coerce_variants.md)

## Examples

``` r
vcf_file <- system.file("extdata", "example.vcf", package = "ggvariant")
variants <- read_vcf(vcf_file)
#> ℹ Reading VCF: example.vcf
#> ✔ Reading VCF: example.vcf [11ms]
#> 
#> Loaded 19 variant records across 7 chromosomes.
head(variants)
#> <gvf: 6 variants, 1 sample, 3 chromosomes, 4 genes>
#>   chrom      pos ref alt qual filter        consequence  gene   sample
#> 1 chr17  7577120   C   T  250   PASS   missense_variant  TP53 TUMOR_S1
#> 3 chr17  7578210   G   T   95   PASS        stop_gained  TP53 TUMOR_S1
#> 4 chr17  7579472   A   G  310   PASS synonymous_variant  TP53 TUMOR_S1
#> 6 chr13 32914437   G   A  175   PASS frameshift_variant BRCA2 TUMOR_S1
#> 7 chr17 41244000  AT   A  310   PASS frameshift_variant BRCA1 TUMOR_S1
#> 9  chr7 55242465   C   A   90   PASS synonymous_variant  EGFR TUMOR_S1
```

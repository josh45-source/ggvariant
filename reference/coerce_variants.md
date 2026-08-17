# Coerce a plain data frame to a `gvf` object

If you already have variant data in a `data.frame` (e.g. exported from
Excel, a database, or another tool), use this function to prepare it for
use with `ggvariant` plotting functions.

## Usage

``` r
coerce_variants(
  x,
  chrom = "chrom",
  pos = "pos",
  ref = "ref",
  alt = "alt",
  consequence = "consequence",
  gene = "gene",
  sample = "sample"
)
```

## Arguments

- x:

  A `data.frame` or `tibble`.

- chrom:

  Column name containing chromosome (default `"chrom"`).

- pos:

  Column name containing position (default `"pos"`).

- ref:

  Column name containing reference allele (default `"ref"`).

- alt:

  Column name containing alternate allele (default `"alt"`).

- consequence:

  Column name containing variant consequence annotation, e.g.
  `"Missense_Mutation"`. If `NULL`, consequence is inferred from REF/ALT
  lengths.

- gene:

  Column name containing gene symbol (default `"gene"`).

- sample:

  Column name containing sample identifier (default `"sample"`).

## Value

A `gvf` object.

## See also

[`read_vcf()`](https://josh45-source.github.io/ggvariant/reference/read_vcf.md)

Other ggvariant input:
[`read_vcf()`](https://josh45-source.github.io/ggvariant/reference/read_vcf.md)

## Examples

``` r
df <- data.frame(
  chromosome = c("chr1", "chr1", "chr7"),
  position   = c(100200, 100350, 55249071),
  ref_allele = c("A", "G", "C"),
  alt_allele = c("T", "A", "T"),
  variant_class = c("missense_variant", "synonymous_variant", "missense_variant"),
  hugo_symbol = c("GENE1", "GENE1", "EGFR"),
  tumor_sample = c("S1", "S2", "S2")
)

variants <- coerce_variants(df,
  chrom       = "chromosome",
  pos         = "position",
  ref         = "ref_allele",
  alt         = "alt_allele",
  consequence = "variant_class",
  gene        = "hugo_symbol",
  sample      = "tumor_sample"
)
```

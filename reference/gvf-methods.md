# Print and summarise `gvf` objects

`print.gvf()` shows a compact header (variant / sample / chromosome /
gene counts) followed by a truncated preview, rather than dumping the
full underlying `data.frame`. `summary.gvf()` returns consequence,
per-sample, and per-chromosome breakdowns as a `summary.gvf` object,
printed via its own method.

## Usage

``` r
# S3 method for class 'gvf'
print(x, ..., n = 6L)

# S3 method for class 'gvf'
summary(object, ...)

# S3 method for class 'summary.gvf'
print(x, ...)
```

## Arguments

- x, object:

  A `gvf` object.

- ...:

  Passed to further methods; currently unused.

- n:

  Integer. Number of rows to preview in `print.gvf()`. Default `6`.

## Value

`print.gvf()` returns `x`, invisibly. `summary.gvf()` returns a
`summary.gvf` object.

## Examples

``` r
vcf_file <- system.file("extdata", "example.vcf", package = "ggvariant")
variants <- read_vcf(vcf_file)
#> ℹ Reading VCF: example.vcf
#> ✔ Reading VCF: example.vcf [15ms]
#> 
#> Loaded 19 variant records across 7 chromosomes.
variants
#> <gvf: 19 variants, 2 samples, 7 chromosomes, 8 genes>
#>   chrom      pos ref alt qual filter        consequence  gene   sample
#> 1 chr17  7577120   C   T  250   PASS   missense_variant  TP53 TUMOR_S1
#> 3 chr17  7578210   G   T   95   PASS        stop_gained  TP53 TUMOR_S1
#> 4 chr17  7579472   A   G  310   PASS synonymous_variant  TP53 TUMOR_S1
#> 6 chr13 32914437   G   A  175   PASS frameshift_variant BRCA2 TUMOR_S1
#> 7 chr17 41244000  AT   A  310   PASS frameshift_variant BRCA1 TUMOR_S1
#> 9  chr7 55242465   C   A   90   PASS synonymous_variant  EGFR TUMOR_S1
#> # ... 13 more rows
summary(variants)
#> gvf summary: 19 variants
#> 
#> Consequence breakdown:
#> 
#>    missense_variant  frameshift_variant         stop_gained splice_site_variant 
#>                   9                   3                   3                   2 
#>  synonymous_variant 
#>                   2 
#> 
#> Per-sample counts:
#> 
#> TUMOR_S1 TUMOR_S2 
#>       10        9 
#> 
#> Chromosome distribution:
#> 
#> chr17 chr12 chr13  chr3  chr7  chr9 chr10 
#>     7     3     2     2     2     2     1 
```

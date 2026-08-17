## R CMD check results

0 errors | 0 warnings | 0 notes

## Summary of changes

ggvariant 0.2.0 is a minor release. Highlights:

  * New `plot_oncoprint()` / `plot_waterfall()` (gene-by-sample mutation
    matrix) and `plot_tmb()` (tumour mutational burden) plots.
  * New `print.gvf()` / `summary.gvf()` methods.
  * Fixed a correctness bug (present since the 0.1.0 release on CRAN):
    `read_vcf()`'s internal sample-pivoting incorrectly attached a variant
    to every sample regardless of its actual genotype, rather than only to
    samples with a non-reference allele. This affected every per-sample
    plot on a multi-sample VCF with mixed genotypes.
  * `plot_variant_spectrum()`'s `context`/`genome` arguments and
    `read_vcf()`'s `info_fields` argument, previously documented but
    silently ignored, now abort with an informative error rather than
    doing nothing.
  * Performance: VCF INFO-field parsing is now vectorised
    (~8.75x faster in isolation on a 100k-record synthetic VCF).

Full details in NEWS.md.

## Downstream dependencies

There are no downstream dependencies for this package.

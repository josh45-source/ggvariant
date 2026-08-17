# Related work

This article surveys other R and Bioconductor tools that address parts
of the same problem space as `ggvariant`, and states where `ggvariant`
sits relative to each. None of the alternatives here are deficient; each
was built for a different point in the workflow, and two of them predate
`ggvariant` by roughly a decade.

## Reading VCF into R

vcfR (Knaus & Grünwald 2017) and VariantAnnotation (Obenchain et
al. 2014) both read VCF files into R and both predate `ggvariant` by
about ten years. vcfR represents a VCF as an S4 `vcfR` object with base
R plotting methods for quality-control diagnostics such as coverage,
depth, and quality distributions; it is not built around `ggplot2` and
its output is not a data frame. VariantAnnotation is part of the core
Bioconductor variant-calling infrastructure: it parses VCF into
`VCF`/`CollapsedVCF` objects built on `GRanges`, designed to
interoperate with the rest of Bioconductor’s genomic-ranges ecosystem
(annotation lookups, overlap operations, reference genome access) rather
than to produce plots directly.

[`read_vcf()`](https://josh45-source.github.io/ggvariant/reference/read_vcf.md)
covers a narrower slice of the same job: it returns a plain tidy data
frame (a `gvf` object — a `data.frame` with an added class attribute,
nothing more), trading the richer Bioconductor object model for a format
any `dplyr`/`tidyr`/base R workflow already knows how to manipulate, and
for direct compatibility with `ggvariant`’s own plotting functions.

## Cancer genomics summary plots

maftools (Mayakonda et al. 2018), GenVisR (Skidmore et al. 2016), and
ComplexHeatmap (Gu et al. 2016) all provide oncoprint/waterfall-style
summary plots, tracing back to the visualisation introduced by
cBioPortal (Gao et al. 2013; Cerami et al. 2012). maftools works from
MAF-format input and is tightly coupled to that cancer-genomics-specific
data model. ComplexHeatmap’s `oncoPrint()` is general-purpose heatmap
infrastructure that happens to support the oncoprint layout, at the cost
of a larger, more general API surface aimed at heatmaps in general, not
variant data specifically.

[`plot_oncoprint()`](https://josh45-source.github.io/ggvariant/reference/plot_oncoprint.md)
(aliased
[`plot_waterfall()`](https://josh45-source.github.io/ggvariant/reference/plot_oncoprint.md))
implements the same memo-sort gene/sample ordering these tools use, but
takes the same tidy `gvf` input as every other `ggvariant` function and
returns a plain `ggplot` object, so it composes with `ggplot2` layers
the way the rest of the package does.

## Mutational signature analysis

MutationalPatterns (Blokzijl et al. 2018) computes and visualises
mutational signatures against the full 96-class trinucleotide context
(Alexandrov et al. 2020), including signature refitting against
reference catalogues, and depends on a `BSgenome` reference-genome
package to extract that context.
[`plot_variant_spectrum()`](https://josh45-source.github.io/ggvariant/reference/plot_variant_spectrum.md)
does one narrower thing: it plots the unweighted 6-class
single-base-substitution spectrum with no genome dependency. 96-class
context support is planned but not yet implemented (`context`/`genome`
currently abort rather than doing anything); until then,
[`plot_variant_spectrum()`](https://josh45-source.github.io/ggvariant/reference/plot_variant_spectrum.md)
is not a substitute for MutationalPatterns’ signature-analysis
capability, only for the simpler 6-class summary.

## Where ggvariant sits

`ggvariant`’s distinguishing feature is not that it does more than these
tools — several of them do substantially more. It is a small API
surface: read a VCF file or a plain data frame, get back a `gvf` tidy
data frame, call one function per plot type, get back an ordinary
`ggplot` object. There is no Bioconductor dependency and no bespoke S4
class to learn. That trade-off is real: less capability, in exchange for
less to learn and full `ggplot2` composability.

## References

- Knaus BJ, Grünwald NJ (2017). vcfR: a package to manipulate and
  visualize variant call format data in R. *Molecular Ecology
  Resources*, 17(1), 44-53.
  [doi:10.1111/1755-0998.12549](https://doi.org/10.1111/1755-0998.12549)
- Obenchain V, Lawrence M, Carey V, Gogarten S, Shannon P, Morgan M
  (2014). VariantAnnotation: a Bioconductor package for exploration and
  annotation of genetic variants. *Bioinformatics*, 30(14), 2076-2078.
  [doi:10.1093/bioinformatics/btu168](https://doi.org/10.1093/bioinformatics/btu168)
- Mayakonda A, Lin DC, Assenov Y, Plass C, Koeffler HP (2018). Maftools:
  efficient and comprehensive analysis of somatic variants in cancer.
  *Genome Research*, 28(11), 1747-1756.
  [doi:10.1101/gr.239244.118](https://doi.org/10.1101/gr.239244.118)
- Skidmore ZL, Wagner AH, Lesurf R, Campbell KM, Kunisaki J, Griffith
  OL, Griffith M (2016). GenVisR: Genomic Visualizations in R.
  *Bioinformatics*, 32(19), 3012-3014.
  [doi:10.1093/bioinformatics/btw325](https://doi.org/10.1093/bioinformatics/btw325)
- Gu Z, Eils R, Schlesner M (2016). Complex heatmaps reveal patterns and
  correlations in multidimensional genomic data. *Bioinformatics*,
  32(18), 2847-2849.
  [doi:10.1093/bioinformatics/btw313](https://doi.org/10.1093/bioinformatics/btw313)
- Gao J, Aksoy BA, Dogrusoz U, et al. (2013). Integrative analysis of
  complex cancer genomics and clinical profiles using the cBioPortal.
  *Science Signaling*, 6(269), pl1.
  [doi:10.1126/scisignal.2004088](https://doi.org/10.1126/scisignal.2004088)
- Cerami E, Gao J, Dogrusoz U, et al. (2012). The cBio Cancer Genomics
  Portal: an open platform for exploring multidimensional cancer
  genomics data. *Cancer Discovery*, 2(5), 401-404.
  [doi:10.1158/2159-8290.CD-12-0095](https://doi.org/10.1158/2159-8290.CD-12-0095)
- Blokzijl F, Janssen R, van Boxtel R, Cuppen E (2018).
  MutationalPatterns: comprehensive genome-wide analysis of mutational
  processes. *Genome Medicine*, 10(1), 33.
  [doi:10.1186/s13073-018-0539-0](https://doi.org/10.1186/s13073-018-0539-0)
- Alexandrov LB, Kim J, Haradhvala NJ, et al.; PCAWG Consortium (2020).
  The repertoire of mutational signatures in human cancer. *Nature*,
  578(7793), 94-101.
  [doi:10.1038/s41586-020-1943-3](https://doi.org/10.1038/s41586-020-1943-3)

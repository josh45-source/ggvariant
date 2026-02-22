# ggvariant <img src="man/figures/logo.png" align="right" height="120" alt=""/>

<!-- badges: start -->
[![R-CMD-check](https://github.com/josh45-source/ggvariant/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/josh45-source/ggvariant/actions)
[![Codecov](https://app.codecov.io/gh/josh45-source/ggvariant/branch/main/graph/badge.svg)](https://app.codecov.io/gh/josh45-source/ggvariant)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

> **Publication-ready genomic variant plots in a few lines of R.**

`ggvariant` fills a gap in the Bioconductor/CRAN ecosystem: a simple,
ggplot2-native package that takes you directly from a **VCF file or data
frame** to beautiful, customisable variant visualisations — without wrestling
with complex APIs or writing 50-line wrangling scripts.

---

## The gap this fills

| Need | Existing options | The problem |
|---|---|---|
| Lollipop plots | Gviz, karyoploteR | Very steep learning curve |
| Consequence summaries | maftools | Tightly coupled to MAF / cancer genomics |
| Mutational spectra | MutationalPatterns | Heavyweight; requires BSgenome |
| General VCF → ggplot | — | **Nothing simple exists** |

`ggvariant` gives both wet-lab biologists and experienced bioinformaticians
the same clean, ggplot2-idiomatic entry point.

---

## Installation

```r
# Install from GitHub (once released)
# install.packages("remotes")
remotes::install_github("josh45-source/ggvariant")
```

---

## Quick start

```r
library(ggvariant)

# 1. Load a VCF file
variants <- read_vcf("my_variants.vcf")

# 2. If you have a data frame instead (e.g. from Excel)
variants <- coerce_variants(my_df,
  chrom = "Chr", pos = "Position",
  ref   = "Ref", alt = "Alt",
  gene  = "Gene", sample = "SampleID"
)
```

---

## Core plots

### Lollipop plot — variants along a gene

```r
plot_lollipop(variants, gene = "TP53")
```

Add protein domain annotations:

```r
tp53_domains <- data.frame(
  name  = c("Transactivation", "DNA-binding", "Tetramerization"),
  start = c(1,   102, 323),
  end   = c(67,  292, 356)
)

plot_lollipop(variants, gene = "TP53", domains = tp53_domains)
```

Colour by sample instead of consequence:

```r
plot_lollipop(variants, gene = "TP53", color_by = "sample")
```

---

### Consequence summary

```r
# Stacked bar by sample
plot_consequence_summary(variants)

# Proportional
plot_consequence_summary(variants, position = "fill")

# Top 10 mutated genes
plot_consequence_summary(variants, group_by = "gene", top_n = 10)
```

---

### Mutational spectrum

```r
# 6-class SBS spectrum
plot_variant_spectrum(variants)

# Faceted by sample
plot_variant_spectrum(variants, facet_by_sample = TRUE)

# Raw counts, not proportions
plot_variant_spectrum(variants, normalize = FALSE)
```

---

## Interactive plots

All plot functions accept `interactive = TRUE` to return a `plotly` object
for sharing with collaborators who don't use R:

```r
plot_lollipop(variants, gene = "BRCA1", interactive = TRUE)
```

---

## Customisation

Because every function returns a standard `ggplot` object, you can layer
on any `ggplot2` or extension code:

```r
library(ggplot2)

plot_lollipop(variants, gene = "KRAS") +
  scale_colour_brewer(palette = "Set2") +
  theme(legend.position = "bottom") +
  labs(subtitle = "KRAS mutations in cohort X")
```

Access palettes directly:

```r
gv_palette("consequence")   # named hex vector
gv_palette("spectrum")      # COSMIC SBS colours
```

---

## Design philosophy

- **Minimal code** — one function call per plot type
- **Two entry points** — VCF files *and* plain data frames
- **ggplot2-native** — every plot is a `ggplot` object; extend freely
- **Opinionated defaults** — looks good out of the box; no mandatory config
- **Progressive disclosure** — simple API for novices, full control for experts

---

## Package structure

```
ggvariant/
├── R/
│   ├── ggvariant-package.R       # Package documentation
│   ├── read_vcf.R                # read_vcf() and coerce_variants()
│   ├── plot_lollipop.R           # plot_lollipop()
│   ├── plot_functions.R          # plot_consequence_summary(), plot_variant_spectrum()
│   └── utils.R                   # Theme, palettes, shared helpers
├── tests/
│   └── testthat/
│       └── test-core.R           # Unit tests
├── inst/
│   └── extdata/
│       └── example.vcf           # Bundled example VCF
├── DESCRIPTION
└── NAMESPACE
```

---

## Roadmap

- [ ] `plot_oncoprint()` — sample × gene mutation matrix
- [ ] `plot_copy_number()` — CNV segment visualisation
- [ ] `plot_rainfall()` — kataegis / mutation density along genome
- [ ] `plot_tmb()` — tumour mutation burden comparison across cohorts
- [ ] BSgenome integration for automatic trinucleotide context extraction
- [ ] Shiny module for non-coding users

---

## Contributing

Pull requests are welcome. Please open an issue first to discuss proposed
changes. All contributions should include tests.

## License

MIT

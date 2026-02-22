#' Consequence summary bar chart
#'
#' Summarises variant consequences (e.g. missense, frameshift, synonymous)
#' across one or more samples, producing a stacked or grouped bar chart.
#'
#' @param variants A `gvf` object or compatible `data.frame`.
#' @param samples Character vector of sample names to include. `NULL` (default)
#'   uses all samples. Ignored if there is no `sample` column.
#' @param group_by `"consequence"` (default) stacks bars by consequence per
#'   sample; `"gene"` stacks by gene per consequence.
#' @param top_n Integer. For `group_by = "gene"`, show only the top N genes by
#'   total variant count. Default `10`.
#' @param position `"stack"` (default) or `"fill"` (proportional) or
#'   `"dodge"`.
#' @param palette Named character vector of colours. `NULL` uses built-in.
#' @param flip Logical. If `TRUE`, flips coordinates for horizontal bars.
#'   Default `FALSE`.
#' @param interactive Logical. Returns a `plotly` object if `TRUE`.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' vcf_file <- system.file("extdata", "example.vcf", package = "ggvariant")
#' variants <- read_vcf(vcf_file)
#'
#' # Consequence counts per sample
#' plot_consequence_summary(variants)
#'
#' # Proportional bars
#' plot_consequence_summary(variants, position = "fill")
#'
#' # Top 10 genes coloured by consequence
#' plot_consequence_summary(variants, group_by = "gene", top_n = 10)
#'
#' @export
plot_consequence_summary <- function(variants,
                                     samples    = NULL,
                                     group_by   = c("consequence", "gene"),
                                     top_n      = 10L,
                                     position   = c("stack", "fill", "dodge"),
                                     palette    = NULL,
                                     flip       = FALSE,
                                     interactive = FALSE) {

  group_by <- match.arg(group_by)
  position <- match.arg(position)
  variants <- .prepare_variants(variants)

  if (!is.null(samples) && "sample" %in% colnames(variants))
    variants <- variants[variants$sample %in% samples, ]

  pal <- palette %||% .consequence_palette()

  if (group_by == "consequence") {
    # X = sample (or "All" if no sample column), fill = consequence
    if (!"sample" %in% colnames(variants) || all(is.na(variants$sample))) {
      variants$sample <- "All variants"
    }
    counts <- aggregate(
      pos ~ sample + consequence,
      data    = variants,
      FUN     = length
    )
    colnames(counts)[3] <- "n"
    counts$consequence  <- .standardise_consequence(counts$consequence)

    # Order samples by total count
    sample_order <- names(sort(tapply(counts$n, counts$sample, sum),
                                decreasing = TRUE))
    counts$sample <- factor(counts$sample, levels = sample_order)

    p <- ggplot2::ggplot(counts,
           ggplot2::aes(
             x    = .data$sample,
             y    = .data$n,
             fill = .data$consequence
           )) +
      ggplot2::geom_col(position = position, width = 0.7) +
      ggplot2::scale_fill_manual(values = pal, name = "Consequence",
                                 na.value = "grey70") +
      ggplot2::labs(
        title = "Variant consequence summary",
        x     = "Sample",
        y     = if (position == "fill") "Proportion" else "Count"
      )

  } else {
    # X = gene, fill = consequence — show top N genes
    if (!"gene" %in% colnames(variants) || all(is.na(variants$gene))) {
      cli::cli_abort("No {.val gene} column found; cannot use {.code group_by = 'gene'}.")
    }
    counts <- aggregate(
      pos ~ gene + consequence,
      data = variants[!is.na(variants$gene), ],
      FUN  = length
    )
    colnames(counts)[3] <- "n"
    counts$consequence  <- .standardise_consequence(counts$consequence)

    gene_totals  <- sort(tapply(counts$n, counts$gene, sum), decreasing = TRUE)
    top_genes    <- names(gene_totals)[seq_len(min(top_n, length(gene_totals)))]
    counts       <- counts[counts$gene %in% top_genes, ]
    counts$gene  <- factor(counts$gene, levels = rev(top_genes))

    p <- ggplot2::ggplot(counts,
           ggplot2::aes(
             x    = .data$gene,
             y    = .data$n,
             fill = .data$consequence
           )) +
      ggplot2::geom_col(position = position, width = 0.7) +
      ggplot2::scale_fill_manual(values = pal, name = "Consequence",
                                 na.value = "grey70") +
      ggplot2::labs(
        title = paste0("Top ", top_n, " mutated genes"),
        x     = "Gene",
        y     = if (position == "fill") "Proportion" else "Count"
      ) +
      ggplot2::coord_flip()
    flip <- FALSE  # already flipped
  }

  p <- p + .ggvariant_theme()
  if (flip) p <- p + ggplot2::coord_flip()

  if (interactive) {
    .require_plotly()
    return(plotly::ggplotly(p))
  }
  p
}


#' Mutational spectrum (SBS) bar chart
#'
#' Plots the single-base substitution (SBS) spectrum — the relative frequency
#' of each of the 6 substitution classes (C>A, C>G, C>T, T>A, T>C, T>G) —
#' optionally broken down by trinucleotide context.
#'
#' @param variants A `gvf` object or compatible `data.frame` containing SNVs.
#'   Indels are automatically excluded.
#' @param sample Character. Sample name to filter on. `NULL` uses all variants
#'   pooled (or facets by sample if `facet_by_sample = TRUE`).
#' @param context Logical. If `TRUE`, shows 96-trinucleotide context bars
#'   (requires a `context` column or a reference genome via `genome`).
#'   Default `FALSE`.
#' @param genome A `BSgenome` object or genome abbreviation string (e.g.
#'   `"hg38"`) used to extract trinucleotide context when `context = TRUE`
#'   and no `context` column is present. Requires the `BSgenome` and
#'   `Biostrings` packages.
#' @param facet_by_sample Logical. If `TRUE`, facets the plot by sample.
#'   Default `FALSE`.
#' @param palette Named character vector with names matching substitution
#'   classes (`"C>A"`, `"C>G"`, etc.). `NULL` uses COSMIC-style colours.
#' @param normalize Logical. If `TRUE` (default), shows relative proportions.
#'   If `FALSE`, shows raw counts.
#' @param interactive Logical. Returns a `plotly` object if `TRUE`.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' vcf_file <- system.file("extdata", "example.vcf", package = "ggvariant")
#' variants <- read_vcf(vcf_file)
#'
#' # Basic 6-class SBS spectrum
#' plot_variant_spectrum(variants)
#'
#' # Faceted by sample
#' plot_variant_spectrum(variants, facet_by_sample = TRUE)
#'
#' @export
plot_variant_spectrum <- function(variants,
                                  sample          = NULL,
                                  context         = FALSE,
                                  genome          = NULL,
                                  facet_by_sample = FALSE,
                                  palette         = NULL,
                                  normalize       = TRUE,
                                  interactive     = FALSE) {

  variants <- .prepare_variants(variants)

  # Keep SNVs only
  snvs <- variants[nchar(variants$ref) == 1 & nchar(variants$alt) == 1 &
                   !is.na(variants$ref) & !is.na(variants$alt), ]
  n_dropped <- nrow(variants) - nrow(snvs)
  if (n_dropped > 0)
    cli::cli_inform("Excluded {n_dropped} non-SNV record{?s} from spectrum plot.")
  if (nrow(snvs) == 0)
    cli::cli_abort("No SNVs remaining after filtering.")

  if (!is.null(sample) && "sample" %in% colnames(snvs))
    snvs <- snvs[snvs$sample == sample, ]

  # Classify substitution class (pyrimidine-normalised)
  snvs$subst <- .classify_substitution(snvs$ref, snvs$alt)

  pal <- palette %||% .spectrum_palette()

  # Aggregate
  group_vars <- if (facet_by_sample && "sample" %in% colnames(snvs))
    c("sample", "subst") else "subst"

  counts <- aggregate(
    pos ~ .,
    data = snvs[, c(group_vars, "pos"), drop = FALSE],
    FUN  = length
  )
  colnames(counts)[ncol(counts)] <- "n"

  if (normalize) {
    if ("sample" %in% group_vars) {
      counts$freq <- ave(counts$n, counts$sample,
                         FUN = function(x) x / sum(x))
    } else {
      counts$freq <- counts$n / sum(counts$n)
    }
    y_var <- "freq"
    y_lab <- "Proportion"
  } else {
    y_var <- "n"
    y_lab <- "Count"
  }

  # Ensure all 6 classes are present
  all_subst <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
  counts$subst <- factor(counts$subst, levels = all_subst)

  p <- ggplot2::ggplot(counts,
         ggplot2::aes(
           x    = .data$subst,
           y    = .data[[y_var]],
           fill = .data$subst
         )) +
    ggplot2::geom_col(width = 0.7, show.legend = FALSE) +
    ggplot2::scale_fill_manual(values = pal, drop = FALSE) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.05)),
      labels = if (normalize) scales::percent_format(accuracy = 1) else scales::comma
    ) +
    ggplot2::labs(
      title = "Mutational spectrum",
      x     = "Substitution class",
      y     = y_lab
    ) +
    .ggvariant_theme()

  if (facet_by_sample && "sample" %in% colnames(snvs)) {
    p <- p + ggplot2::facet_wrap(~ sample, scales = "free_y")
  }

  if (interactive) {
    .require_plotly()
    return(plotly::ggplotly(p))
  }
  p
}

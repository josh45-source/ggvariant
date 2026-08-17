#' Oncoprint / waterfall plot of variants across samples
#'
#' Draws a gene-by-sample mutation matrix ("oncoprint" in the
#' `ComplexHeatmap`/cBioPortal vocabulary, "waterfall plot" in the
#' `GenVisR`/`maftools` vocabulary — the same visualisation under two
#' names). `plot_waterfall()` is an alias for `plot_oncoprint()`; both
#' produce identical output.
#'
#' @details
#' **Gene order.** Genes are ranked by the number of distinct samples
#' carrying at least one mutation in that gene (ties broken by total
#' mutation count, then alphabetically), and the `top_n` most frequently
#' altered genes are shown, most-altered at the top.
#'
#' **Sample order.** Samples are ordered using the memo-sort ("cascade")
#' algorithm: the displayed genes, already ranked most-to-least altered,
#' are treated as bits of a binary number (the most-altered gene the most
#' significant bit). Each sample's mutation pattern across the displayed
#' genes is scored as that binary number, and samples are sorted by
#' descending score (ties broken alphabetically by sample name). This
#' greedily groups samples that share mutations in the top genes together,
#' producing the characteristic left-to-right staircase pattern. Only the
#' displayed genes contribute to the score; genes excluded by `top_n`/
#' `genes` have no effect on sample order.
#'
#' **Multi-hit cells.** A gene mutated more than once in the same sample
#' cannot be represented by a single consequence colour, so it is shown as
#' a distinct `"Multi_Hit"` category instead of either consequence.
#'
#' @param variants A `gvf` object from [read_vcf()] or [coerce_variants()],
#'   or any `data.frame` with `gene`, `sample`, and `consequence` columns.
#' @param top_n Integer. Show the `top_n` most frequently altered genes.
#'   Mutually exclusive with `genes`; if both are `NULL`, defaults to `10`.
#' @param genes Character vector of specific genes to show, in place of
#'   `top_n`. Mutually exclusive with `top_n`. Genes not present in
#'   `variants` are shown as fully unaltered rows.
#' @param samples Character vector of sample names to include as columns.
#'   `NULL` (default) uses every sample present in `variants`.
#' @param annotation Optional `data.frame` of sample-level metadata drawn as
#'   annotation tracks below the mutation matrix. Must contain a `sample`
#'   column matching sample identifiers in `variants`, plus one or more
#'   additional columns to display, one per track. `NULL` (default) omits
#'   annotation tracks.
#' @param palette Named character vector of colours keyed by consequence.
#'   `NULL` uses the built-in `gv_palette("consequence")`. If it does not
#'   already contain a `"Multi_Hit"` entry, one is added automatically.
#' @param interactive Logical. Returns a `plotly` object if `TRUE`.
#'
#' @return A `ggplot` object (or a `plotly` object when `interactive = TRUE`).
#'
#' @examples
#' vcf_file <- system.file("extdata", "example.vcf", package = "ggvariant")
#' variants <- read_vcf(vcf_file)
#'
#' # Top 5 most-altered genes
#' plot_oncoprint(variants, top_n = 5)
#'
#' # A specific gene panel instead of top_n
#' plot_oncoprint(variants, genes = c("TP53", "BRCA1", "BRCA2"))
#'
#' # With a clinical annotation track
#' clinical <- data.frame(
#'   sample = c("TUMOR_S1", "TUMOR_S2"),
#'   stage  = c("III", "IV")
#' )
#' plot_oncoprint(variants, top_n = 5, annotation = clinical)
#'
#' @references
#' Gao J, Aksoy BA, Dogrusoz U, et al. (2013). Integrative analysis of
#' complex cancer genomics and clinical profiles using the cBioPortal.
#' *Science Signaling*, 6(269), pl1. \doi{10.1126/scisignal.2004088}
#'
#' Skidmore ZL, Wagner AH, Lesurf R, et al. (2016). GenVisR: Genomic
#' Visualizations in R. *Bioinformatics*, 32(19), 3012-3014.
#' \doi{10.1093/bioinformatics/btw325}
#'
#' Gu Z, Eils R, Schlesner M (2016). Complex heatmaps reveal patterns and
#' correlations in multidimensional genomic data. *Bioinformatics*,
#' 32(18), 2847-2849. \doi{10.1093/bioinformatics/btw313}
#'
#' Mayakonda A, Lin DC, Assenov Y, Plass C, Koeffler HP (2018). Maftools:
#' efficient and comprehensive analysis of somatic variants in cancer.
#' *Genome Research*, 28(11), 1747-1756. \doi{10.1101/gr.239244.118}
#'
#' @family ggvariant plots
#' @seealso [plot_tmb()], [plot_consequence_summary()], [gv_palette()]
#' @export
plot_oncoprint <- function(variants,
                           top_n       = NULL,
                           genes       = NULL,
                           samples     = NULL,
                           annotation  = NULL,
                           palette     = NULL,
                           interactive = FALSE) {

  if (!is.null(top_n) && !is.null(genes))
    cli::cli_abort("Specify only one of {.arg top_n} or {.arg genes}, not both.")
  if (is.null(top_n) && is.null(genes)) top_n <- 10L

  variants <- .prepare_variants(variants)
  if (!all(c("gene", "sample") %in% colnames(variants)))
    cli::cli_abort(
      "{.arg variants} must have {.val gene} and {.val sample} columns."
    )

  variants <- variants[!is.na(variants$gene) & !is.na(variants$sample), ,
                        drop = FALSE]

  all_samples <- sort(unique(variants$sample))
  sample_universe <- if (!is.null(samples)) intersect(samples, all_samples) else
    all_samples
  if (length(sample_universe) == 0)
    cli::cli_abort("No samples remaining after filtering.")

  variants <- variants[variants$sample %in% sample_universe, , drop = FALSE]

  candidate_genes <- if (!is.null(genes)) unique(genes) else
    unique(variants$gene)
  if (length(candidate_genes) == 0)
    cli::cli_abort("No genes to plot.")

  gene_rank  <- .rank_genes(variants, candidate_genes)
  gene_order <- if (!is.null(genes)) gene_rank else
    utils::head(gene_rank, top_n)

  sample_order <- .memo_sort_samples(
    variants[variants$gene %in% gene_order, , drop = FALSE],
    gene_order, sample_universe
  )

  cells <- .oncoprint_cells(variants, gene_order, sample_order)

  pal <- palette %||% .consequence_palette()
  if (!"Multi_Hit" %in% names(pal)) pal <- c(pal, Multi_Hit = "grey15")

  n_genes  <- length(gene_order)
  y_breaks <- rev(seq_len(n_genes))
  y_labels <- gene_order

  p <- ggplot2::ggplot(
         cells,
         ggplot2::aes(x = .data$x, y = .data$y, fill = .data$label)
       ) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.5) +
    ggplot2::scale_fill_manual(
      values = pal, name = "Consequence",
      na.value = "grey92", na.translate = FALSE
    )

  ann_tracks <- character(0)
  if (!is.null(annotation)) {
    ann <- .oncoprint_annotation(annotation, sample_order)
    ann_tracks <- ann$tracks
    p <- p +
      ggplot2::geom_tile(
        data = ann$cells,
        ggplot2::aes(x = .data$x, y = .data$y, fill = I(.data$hex)),
        inherit.aes = FALSE, colour = "white", linewidth = 0.5
      ) +
      ggplot2::geom_text(
        data = ann$cells,
        ggplot2::aes(x = .data$x, y = .data$y, label = .data$value),
        inherit.aes = FALSE, size = 2.2, na.rm = TRUE
      )
  }

  n_tracks     <- length(ann_tracks)
  y_breaks_all <- c(y_breaks, if (n_tracks) -seq_len(n_tracks))
  y_labels_all <- c(y_labels, ann_tracks)

  p <- p +
    ggplot2::scale_x_discrete(expand = c(0, 0)) +
    ggplot2::scale_y_continuous(
      breaks = y_breaks_all, labels = y_labels_all, expand = c(0, 0)
    ) +
    ggplot2::labs(title = "Oncoprint", x = NULL, y = NULL) +
    .ggvariant_theme() +
    ggplot2::theme(
      panel.grid  = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5)
    )

  if (interactive) {
    .require_plotly()
    return(plotly::ggplotly(p))
  }
  p
}

#' @rdname plot_oncoprint
#' @export
plot_waterfall <- plot_oncoprint


# ── Internal helpers ──────────────────────────────────────────────────────────

.rank_genes <- function(variants, candidate_genes) {
  sub <- variants[variants$gene %in% candidate_genes, , drop = FALSE]
  n_samples <- tapply(sub$sample, sub$gene, function(s) length(unique(s)))
  n_muts    <- tapply(sub$gene,   sub$gene, length)

  stats_df <- data.frame(
    gene      = candidate_genes,
    n_samples = as.integer(n_samples[candidate_genes]),
    n_muts    = as.integer(n_muts[candidate_genes]),
    stringsAsFactors = FALSE
  )
  stats_df$n_samples[is.na(stats_df$n_samples)] <- 0L
  stats_df$n_muts[is.na(stats_df$n_muts)]       <- 0L

  ord <- order(-stats_df$n_samples, -stats_df$n_muts, stats_df$gene)
  stats_df$gene[ord]
}

.memo_sort_samples <- function(variants, gene_order, sample_universe) {
  mat <- matrix(FALSE, nrow = length(gene_order), ncol = length(sample_universe),
                dimnames = list(gene_order, sample_universe))
  if (nrow(variants) > 0) {
    ridx <- match(variants$gene, gene_order)
    cidx <- match(variants$sample, sample_universe)
    keep <- !is.na(ridx) & !is.na(cidx)
    mat[cbind(ridx[keep], cidx[keep])] <- TRUE
  }
  weights <- 2^(length(gene_order) - seq_along(gene_order))
  score   <- as.numeric(weights %*% (mat * 1))
  sample_universe[order(-score, sample_universe)]
}

.oncoprint_cells <- function(variants, gene_order, sample_order) {
  sub  <- variants[variants$gene %in% gene_order &
                    variants$sample %in% sample_order, , drop = FALSE]
  full <- expand.grid(gene = gene_order, sample = sample_order,
                       stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)

  if (nrow(sub) > 0) {
    key   <- paste(sub$gene, sub$sample, sep = "\r")
    n_tab <- table(key)
    first_idx <- match(names(n_tab), key)
    hits <- data.frame(
      key         = names(n_tab),
      n           = as.integer(n_tab),
      consequence = sub$consequence[first_idx],
      stringsAsFactors = FALSE
    )
    hits$label <- ifelse(hits$n > 1L, "Multi_Hit",
                          .standardise_consequence(hits$consequence))
    full$key   <- paste(full$gene, full$sample, sep = "\r")
    full$label <- hits$label[match(full$key, hits$key)]
    full$key   <- NULL
  } else {
    full$label <- NA_character_
  }

  n_genes <- length(gene_order)
  gene_y  <- stats::setNames(rev(seq_len(n_genes)), gene_order)
  full$y  <- gene_y[full$gene]
  full$x  <- factor(full$sample, levels = sample_order)
  full
}

.oncoprint_annotation <- function(annotation, sample_order) {
  if (!is.data.frame(annotation))
    cli::cli_abort("{.arg annotation} must be a data frame.")
  if (!"sample" %in% colnames(annotation))
    cli::cli_abort("{.arg annotation} must have a {.val sample} column.")
  tracks <- setdiff(colnames(annotation), "sample")
  if (length(tracks) == 0)
    cli::cli_abort(
      "{.arg annotation} must have at least one column besides {.val sample}."
    )

  rows <- lapply(seq_along(tracks), function(i) {
    col  <- tracks[i]
    vals <- as.character(annotation[[col]])
    uniq <- sort(unique(stats::na.omit(vals)))
    cols <- stats::setNames(.domain_palette(length(uniq)), uniq)
    data.frame(
      sample = annotation$sample,
      track  = col,
      value  = vals,
      hex    = unname(cols[vals]),
      y      = -i,
      stringsAsFactors = FALSE
    )
  })
  cells <- do.call(rbind, rows)
  cells <- cells[cells$sample %in% sample_order, , drop = FALSE]
  cells$x <- factor(cells$sample, levels = sample_order)
  cells$hex[is.na(cells$hex)] <- "grey95"
  list(cells = cells, tracks = tracks)
}

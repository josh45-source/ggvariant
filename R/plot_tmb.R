#' Tumour mutational burden bar chart
#'
#' Plots per-sample tumour mutational burden (TMB): the number of mutations
#' carried by each sample, optionally normalised to mutations per megabase.
#'
#' @details
#' The per-sample mutation count is computed by an internal helper shared
#' with any future [plot_oncoprint()] TMB marginal, so the two always agree.
#' Pass `mb_size` (the size, in megabases, of the sequenced region) to
#' convert raw counts into mutations/Mb, the conventional TMB unit; leave it
#' `NULL` to plot raw counts.
#'
#' @param variants A `gvf` object from [read_vcf()] or [coerce_variants()],
#'   or any `data.frame` with `sample` and `pos` columns.
#' @param samples Character vector of sample names to include. `NULL`
#'   (default) uses every sample present in `variants`.
#' @param mb_size Numeric. Size, in megabases, of the sequenced region used
#'   to normalise mutation counts into mutations/Mb. `NULL` (default) plots
#'   raw mutation counts instead.
#' @param sample_order Character vector giving an explicit left-to-right
#'   sample order, e.g. to align with [plot_oncoprint()]'s memo-sorted
#'   column order. `NULL` (default) orders samples by descending TMB.
#' @param palette Single hex colour string for the bars. `NULL` uses a
#'   built-in default.
#' @param interactive Logical. Returns a `plotly` object if `TRUE`.
#'
#' @return A `ggplot` object (or a `plotly` object when `interactive = TRUE`).
#'
#' @examples
#' vcf_file <- system.file("extdata", "example.vcf", package = "ggvariant")
#' variants <- read_vcf(vcf_file)
#'
#' # Raw mutation counts, samples ordered by descending TMB
#' plot_tmb(variants)
#'
#' # Normalised to mutations/Mb for a 38 Mb exome
#' plot_tmb(variants, mb_size = 38)
#'
#' # Aligned to a specific sample order, e.g. from plot_oncoprint()
#' plot_tmb(variants, sample_order = c("TUMOR_S1", "TUMOR_S2"))
#'
#' @references
#' Chalmers ZR, Connelly CF, Fabrizio D, et al. (2017). Analysis of
#' 100,000 human cancer genomes reveals the landscape of tumor
#' mutational burden. *Genome Medicine*, 9(1), 34.
#' \doi{10.1186/s13073-017-0424-2}
#'
#' @family ggvariant plots
#' @seealso [plot_oncoprint()], [gv_palette()]
#' @export
plot_tmb <- function(variants,
                     samples      = NULL,
                     mb_size      = NULL,
                     sample_order = NULL,
                     palette      = NULL,
                     interactive  = FALSE) {

  tmb <- .compute_tmb(variants, samples = samples, mb_size = mb_size)

  if (!is.null(sample_order)) {
    missing <- setdiff(sample_order, tmb$sample)
    if (length(missing))
      cli::cli_abort(
        "{.arg sample_order} contains sample{?s} not present after \\
         filtering: {.val {missing}}."
      )
    tmb <- tmb[match(sample_order, tmb$sample), , drop = FALSE]
    tmb$sample <- factor(tmb$sample, levels = sample_order)
  } else {
    tmb <- tmb[order(-tmb$tmb, tmb$sample), , drop = FALSE]
    tmb$sample <- factor(tmb$sample, levels = tmb$sample)
  }

  bar_colour <- palette %||% "#4E79A7"
  y_lab <- if (!is.null(mb_size)) "Mutations / Mb" else "Mutation count"

  p <- ggplot2::ggplot(tmb, ggplot2::aes(x = .data$sample, y = .data$tmb)) +
    ggplot2::geom_col(fill = bar_colour, width = 0.7) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.05)),
      labels = scales::comma
    ) +
    ggplot2::labs(title = "Tumour mutational burden", x = "Sample", y = y_lab) +
    .ggvariant_theme() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5)
    )

  if (interactive) {
    .require_plotly()
    return(plotly::ggplotly(p))
  }
  p
}


# ── Internal helpers ──────────────────────────────────────────────────────────

.compute_tmb <- function(variants, samples = NULL, mb_size = NULL) {
  variants <- .prepare_variants(variants)
  if (!"sample" %in% colnames(variants))
    cli::cli_abort("{.arg variants} must have a {.val sample} column.")

  variants <- variants[!is.na(variants$sample), , drop = FALSE]
  if (!is.null(samples))
    variants <- variants[variants$sample %in% samples, , drop = FALSE]
  if (nrow(variants) == 0)
    cli::cli_abort("No mutations remaining after filtering.")

  counts <- as.data.frame(table(sample = variants$sample),
                           stringsAsFactors = FALSE)
  colnames(counts) <- c("sample", "n")
  counts$sample <- as.character(counts$sample)
  counts <- counts[counts$n > 0, , drop = FALSE]

  if (!is.null(mb_size)) {
    if (!is.numeric(mb_size) || length(mb_size) != 1L || mb_size <= 0)
      cli::cli_abort("{.arg mb_size} must be a single positive number.")
    counts$tmb <- counts$n / mb_size
  } else {
    counts$tmb <- counts$n
  }
  counts
}

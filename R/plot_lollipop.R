#' Lollipop plot of variants along a gene
#'
#' Draws a lollipop (stem-and-dot) diagram showing variant positions along
#' a gene, coloured by consequence. Optionally overlays protein domain
#' annotations when domain boundaries are supplied.
#'
#' @param variants A `gvf` object from [read_vcf()] or [coerce_variants()],
#'   or any `data.frame` with columns `pos`, `consequence`, and optionally
#'   `gene` and `sample`.
#' @param gene Character. Gene to filter on. If `NULL` and `variants` contains
#'   a `gene` column, the most-mutated gene is chosen automatically.
#' @param domains A `data.frame` with columns `name`, `start`, `end` (amino
#'   acid positions) for domain annotation. `NULL` (default) omits domains.
#' @param color_by Column name to use for dot colour. Default `"consequence"`.
#'   Set to `"sample"` to colour by sample instead.
#' @param palette Named character vector of colours for each consequence/sample
#'   category. `NULL` uses the built-in `ggvariant` palette.
#' @param protein_length Integer. Total length of the protein in amino acids,
#'   used to scale the x-axis. If `NULL`, inferred from `max(pos)`.
#' @param stack_dots Logical. If `TRUE` (default), dots at the same position
#'   are stacked vertically (beeswarm-style) rather than overlapping.
#' @param title Character. Plot title. Defaults to the gene name.
#' @param interactive Logical. If `TRUE`, returns a `plotly` interactive plot
#'   (requires the `plotly` package).
#'
#' @return A `ggplot` object (or a `plotly` object when `interactive = TRUE`).
#'
#' @examples
#' vcf_file <- system.file("extdata", "example.vcf", package = "ggvariant")
#' variants <- read_vcf(vcf_file)
#'
#' # Basic lollipop for the most-mutated gene
#' plot_lollipop(variants)
#'
#' # Specific gene
#' plot_lollipop(variants, gene = "TP53")
#'
#' # With domain annotation
#' tp53_domains <- data.frame(
#'   name  = c("Transactivation", "DNA-binding", "Tetramerization"),
#'   start = c(1, 102, 323),
#'   end   = c(67, 292, 356)
#' )
#' plot_lollipop(variants, gene = "TP53", domains = tp53_domains)
#'
#' @export
plot_lollipop <- function(variants,
                          gene           = NULL,
                          domains        = NULL,
                          color_by       = "consequence",
                          palette        = NULL,
                          protein_length = NULL,
                          stack_dots     = TRUE,
                          title          = NULL,
                          interactive    = FALSE) {

  variants <- .prepare_variants(variants)

  # Gene selection
  if (is.null(gene)) {
    if ("gene" %in% colnames(variants) && any(!is.na(variants$gene))) {
      gene <- names(sort(table(variants$gene), decreasing = TRUE))[1]
      cli::cli_inform("No gene specified; using most-mutated gene: {.val {gene}}")
    }
  }
  if (!is.null(gene) && "gene" %in% colnames(variants)) {
    variants <- variants[variants$gene == gene & !is.na(variants$gene), ]
  }
  if (nrow(variants) == 0)
    cli::cli_abort("No variants remaining after filtering for gene {.val {gene}}.")

  # Resolve colour column
  if (!color_by %in% colnames(variants))
    cli::cli_abort("Column {.val {color_by}} not found in variant data.")

  pal <- palette %||% .consequence_palette()

  # Stack / jitter y positions
  if (stack_dots) {
    variants <- .stack_positions(variants)
  } else {
    variants$y_pos <- 1L
  }

  x_max <- protein_length %||% max(variants$pos, na.rm = TRUE)
  plot_title <- title %||% if (!is.null(gene)) paste(gene, "variants") else "Variant lollipop"

  p <- ggplot2::ggplot(variants,
         ggplot2::aes(x = .data$pos, y = .data$y_pos,
                      colour = .data[[color_by]],
                      text   = .lollipop_tooltip(variants))) +
    # Stems
    ggplot2::geom_segment(
      ggplot2::aes(xend = .data$pos, yend = 0),
      linewidth = 0.4, colour = "grey70"
    ) +
    # Dots
    ggplot2::geom_point(size = 3.5, alpha = 0.9) +
    # Protein backbone
    ggplot2::annotate("rect",
      xmin = 0, xmax = x_max, ymin = -0.15, ymax = 0.15,
      fill = "#CCCCCC", colour = NA
    ) +
    ggplot2::scale_colour_manual(
      values = pal,
      name   = .label_case(color_by),
      na.value = "grey50"
    ) +
    ggplot2::scale_x_continuous(
      limits = c(0, x_max * 1.02),
      expand = ggplot2::expansion(0),
      labels = scales::comma
    ) +
    ggplot2::scale_y_continuous(breaks = NULL) +
    ggplot2::labs(
      title = plot_title,
      x     = "Amino acid position",
      y     = NULL
    ) +
    .ggvariant_theme()

  # Domain annotation layer
  if (!is.null(domains)) {
    .check_domain_df(domains)
    domain_colors <- .domain_palette(nrow(domains))
    p <- p +
      ggplot2::geom_rect(
        data = domains,
        ggplot2::aes(
          xmin = .data$start, xmax = .data$end,
          ymin = -0.4, ymax = 0.4,
          fill = .data$name
        ),
        inherit.aes = FALSE, alpha = 0.7, colour = NA
      ) +
      ggplot2::scale_fill_manual(
        values = domain_colors,
        name   = "Domain"
      ) +
      ggplot2::geom_text(
        data = domains,
        ggplot2::aes(
          x     = (.data$start + .data$end) / 2,
          y     = 0,
          label = .data$name
        ),
        inherit.aes = FALSE,
        size = 2.8, fontface = "bold", colour = "white"
      )
  }

  if (interactive) {
    .require_plotly()
    return(plotly::ggplotly(p, tooltip = "text"))
  }
  p
}

# ── Helpers ───────────────────────────────────────────────────────────────────

.stack_positions <- function(df) {
  df <- df[order(df$pos), ]
  df$y_pos <- stats::ave(df$pos, df$pos, FUN = function(x) seq_along(x))
  df$y_pos <- as.integer(df$y_pos)
  df
}

.check_domain_df <- function(d) {
  required <- c("name", "start", "end")
  missing  <- setdiff(required, colnames(d))
  if (length(missing))
    cli::cli_abort(
      "{.arg domains} must have columns {.val {required}}. \\
       Missing: {.val {missing}}."
    )
}

.lollipop_tooltip <- function(df) {
  paste0(
    "<b>", df$gene %||% "", "</b><br>",
    "Pos: ", df$pos, "<br>",
    "Change: ", df$ref, ">", df$alt, "<br>",
    "Consequence: ", df$consequence,
    if ("sample" %in% colnames(df)) paste0("<br>Sample: ", df$sample) else ""
  )
}

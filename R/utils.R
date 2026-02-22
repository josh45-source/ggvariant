#' ggvariant colour palettes
#'
#' Access the built-in colour palettes used by `ggvariant` plot functions.
#'
#' @param type One of `"consequence"` (default), `"spectrum"`, or `"domain"`.
#' @param n Integer. For `"domain"`, the number of colours to generate.
#'
#' @return A named character vector of hex colour codes.
#'
#' @examples
#' gv_palette("consequence")
#' gv_palette("spectrum")
#'
#' @export
gv_palette <- function(type = c("consequence", "spectrum", "domain"), n = 8L) {
  type <- match.arg(type)
  switch(type,
    consequence = .consequence_palette(),
    spectrum    = .spectrum_palette(),
    domain      = .domain_palette(n)
  )
}


#' ggvariant ggplot2 theme
#'
#' A clean, publication-ready theme based on `theme_minimal`. Applied
#' automatically by all `ggvariant` plot functions; export it to customise
#' further.
#'
#' @param base_size Base font size in pt. Default `12`.
#' @param base_family Base font family. Default `""` (system sans-serif).
#'
#' @return A `ggplot2` theme object.
#'
#' @examples
#' library(ggplot2)
#' ggplot(mtcars, aes(mpg, wt)) + geom_point() + theme_ggvariant()
#'
#' @export
theme_ggvariant <- function(base_size = 12, base_family = "") {
  ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      plot.title        = ggplot2::element_text(face = "bold", size = base_size * 1.15,
                                                 margin = ggplot2::margin(b = 8)),
      axis.title        = ggplot2::element_text(size = base_size * 0.9, colour = "grey30"),
      axis.text         = ggplot2::element_text(size = base_size * 0.85),
      legend.title      = ggplot2::element_text(face = "bold", size = base_size * 0.9),
      legend.text       = ggplot2::element_text(size = base_size * 0.85),
      panel.grid.major  = ggplot2::element_line(colour = "grey92", linewidth = 0.4),
      panel.grid.minor  = ggplot2::element_blank(),
      panel.border      = ggplot2::element_blank(),
      strip.text        = ggplot2::element_text(face = "bold"),
      plot.margin       = ggplot2::margin(12, 12, 12, 12)
    )
}


# ── Internal colour helpers ───────────────────────────────────────────────────

.consequence_palette <- function() {
  # Inspired by Ensembl VEP consequence colours
  c(
    "missense_variant"         = "#FD8D3C",
    "missense"                 = "#FD8D3C",
    "Missense_Mutation"        = "#FD8D3C",
    "synonymous_variant"       = "#74C476",
    "synonymous"               = "#74C476",
    "Silent"                   = "#74C476",
    "frameshift_variant"       = "#E31A1C",
    "frameshift_deletion"      = "#E31A1C",
    "frameshift_insertion"     = "#FC4E2A",
    "stop_gained"              = "#800026",
    "Nonsense_Mutation"        = "#800026",
    "stop_lost"                = "#BD0026",
    "splice_site_variant"      = "#9E9AC8",
    "Splice_Site"              = "#9E9AC8",
    "inframe_insertion"        = "#6BAED6",
    "inframe_deletion"         = "#2171B5",
    "In_Frame_Del"             = "#2171B5",
    "In_Frame_Ins"             = "#6BAED6",
    "5_prime_UTR_variant"      = "#BCBDDC",
    "3_prime_UTR_variant"      = "#DADAEB",
    "intron_variant"           = "#D9D9D9",
    "SNV"                      = "#BDBDBD",
    "deletion"                 = "#41AB5D",
    "insertion"                = "#A1D99B",
    "MNV"                      = "#F768A1"
  )
}

.spectrum_palette <- function() {
  # COSMIC SBS colour scheme
  c(
    "C>A" = "#03BDEE",
    "C>G" = "#010101",
    "C>T" = "#E52926",
    "T>A" = "#CAC9C9",
    "T>C" = "#A0CE63",
    "T>G" = "#ECC7C5"
  )
}

.domain_palette <- function(n) {
  base <- c(
    "#4E79A7", "#F28E2B", "#E15759", "#76B7B2",
    "#59A14F", "#EDC948", "#B07AA1", "#FF9DA7",
    "#9C755F", "#BAB0AC"
  )
  if (n <= length(base)) base[seq_len(n)]
  else grDevices::colorRampPalette(base)(n)
}


# ── Internal substitution helpers ─────────────────────────────────────────────

.classify_substitution <- function(ref, alt) {
  # Pyrimidine-normalise: if ref is purine (A/G), take complement
  purines <- ref %in% c("A", "G")
  ref_norm <- ref
  alt_norm <- alt
  ref_norm[purines] <- .complement(ref[purines])
  alt_norm[purines] <- .complement(alt[purines])
  paste0(ref_norm, ">", alt_norm)
}

.complement <- function(bases) {
  comp <- c(A = "T", T = "A", C = "G", G = "C")
  comp[bases]
}

.standardise_consequence <- function(x) {
  # Map common aliases to a unified label for display
  lut <- c(
    "Missense_Mutation"  = "missense_variant",
    "Silent"             = "synonymous_variant",
    "Nonsense_Mutation"  = "stop_gained",
    "Splice_Site"        = "splice_site_variant",
    "Frame_Shift_Del"    = "frameshift_deletion",
    "Frame_Shift_Ins"    = "frameshift_insertion",
    "In_Frame_Del"       = "inframe_deletion",
    "In_Frame_Ins"       = "inframe_insertion"
  )
  # Preserve unmatched values
  x <- ifelse(x %in% names(lut), lut[x], x)
  x
}


# ── Shared internal utilities ─────────────────────────────────────────────────

# Null-coalescing operator
`%||%` <- function(a, b) if (!is.null(a)) a else b

.prepare_variants <- function(x) {
  if (inherits(x, "gvf")) return(x)
  if (is.data.frame(x)) {
    required <- c("pos", "ref", "alt")
    missing  <- setdiff(required, colnames(x))
    if (length(missing))
      cli::cli_abort(
        "Data frame must contain columns {.val {required}}. Missing: {.val {missing}}."
      )
    class(x) <- c("gvf", "data.frame")
    return(x)
  }
  cli::cli_abort(
    "{.arg variants} must be a {.cls gvf} or {.cls data.frame}, not {.cls {class(x)}}."
  )
}

.ggvariant_theme <- function() theme_ggvariant()

.label_case <- function(x) {
  gsub("_", " ", paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x))))
}

.require_plotly <- function() {
  if (!requireNamespace("plotly", quietly = TRUE))
    cli::cli_abort(
      "Package {.pkg plotly} is required for interactive plots. \\
       Install it with {.code install.packages('plotly')}."
    )
}

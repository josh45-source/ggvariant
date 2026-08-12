#' Print and summarise `gvf` objects
#'
#' @description
#' `print.gvf()` shows a compact header (variant / sample / chromosome /
#' gene counts) followed by a truncated preview, rather than dumping the
#' full underlying `data.frame`. `summary.gvf()` returns consequence,
#' per-sample, and per-chromosome breakdowns as a `summary.gvf` object,
#' printed via its own method.
#'
#' @param x,object A `gvf` object.
#' @param n Integer. Number of rows to preview in `print.gvf()`. Default `6`.
#' @param ... Passed to further methods; currently unused.
#'
#' @return `print.gvf()` returns `x`, invisibly. `summary.gvf()` returns a
#'   `summary.gvf` object.
#'
#' @examples
#' vcf_file <- system.file("extdata", "example.vcf", package = "ggvariant")
#' variants <- read_vcf(vcf_file)
#' variants
#' summary(variants)
#'
#' @name gvf-methods
#' @export
print.gvf <- function(x, ..., n = 6L) {
  n_var    <- nrow(x)
  n_sample <- if ("sample" %in% colnames(x))
    length(unique(stats::na.omit(x$sample))) else 0L
  n_chrom  <- if ("chrom" %in% colnames(x))
    length(unique(stats::na.omit(x$chrom))) else 0L
  n_gene   <- if ("gene" %in% colnames(x))
    length(unique(stats::na.omit(x$gene))) else 0L

  cat(sprintf(
    "<gvf: %d variant%s, %d sample%s, %d chromosome%s, %d gene%s>\n",
    n_var,    if (n_var    == 1L) "" else "s",
    n_sample, if (n_sample == 1L) "" else "s",
    n_chrom,  if (n_chrom  == 1L) "" else "s",
    n_gene,   if (n_gene   == 1L) "" else "s"
  ))
  print.data.frame(utils::head(as.data.frame(x), n))
  if (n_var > n)
    cat(sprintf(
      "# ... %d more row%s\n", n_var - n, if (n_var - n == 1L) "" else "s"
    ))
  invisible(x)
}

#' @rdname gvf-methods
#' @export
summary.gvf <- function(object, ...) {
  out <- list(
    n_variants  = nrow(object),
    consequence = if ("consequence" %in% colnames(object))
      sort(table(object$consequence), decreasing = TRUE) else table(character()),
    per_sample  = if ("sample" %in% colnames(object))
      sort(table(object$sample), decreasing = TRUE) else table(character()),
    chromosome  = if ("chrom" %in% colnames(object))
      sort(table(object$chrom), decreasing = TRUE) else table(character())
  )
  class(out) <- "summary.gvf"
  out
}

#' @rdname gvf-methods
#' @export
print.summary.gvf <- function(x, ...) {
  cat(sprintf(
    "gvf summary: %d variant%s\n\n",
    x$n_variants, if (x$n_variants == 1L) "" else "s"
  ))
  cat("Consequence breakdown:\n")
  print(x$consequence)
  cat("\nPer-sample counts:\n")
  print(x$per_sample)
  cat("\nChromosome distribution:\n")
  print(x$chromosome)
  invisible(x)
}


# ── Internal helpers ──────────────────────────────────────────────────────────

validate_gvf <- function(x) {
  required <- c("chrom", "pos", "ref", "alt", "consequence", "gene", "sample")
  missing  <- setdiff(required, colnames(x))
  if (length(missing))
    cli::cli_abort(
      "Malformed {.cls gvf} object: missing column{?s} {.val {missing}}."
    )
  if (!is.numeric(x$pos))
    cli::cli_abort("Malformed {.cls gvf} object: {.field pos} must be numeric.")
  invisible(x)
}

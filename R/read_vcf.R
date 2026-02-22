#' Read a VCF file into a tidy variant data frame
#'
#' Parses a standard VCF (v4.x) file and returns a tidy `data.frame` (a
#' `gvf` object) that all `ggvariant` plotting functions accept. For users
#' who already have variant data in a plain `data.frame` or `tibble`, see
#' [coerce_variants()].
#'
#' @param path Path to a `.vcf` or `.vcf.gz` file.
#' @param samples Character vector of sample names to retain. `NULL` (default)
#'   keeps all samples.
#' @param pass_only Logical. If `TRUE` (default), only variants with `FILTER`
#'   equal to `"PASS"` or `"."` are retained.
#' @param info_fields Character vector of INFO field names to expand into
#'   columns. `NULL` keeps none. Use `"all"` to expand everything (may be slow
#'   for large files).
#'
#' @return A `gvf` (genomic variant frame) — a `data.frame` with columns:
#'   \describe{
#'     \item{chrom}{Chromosome (character)}
#'     \item{pos}{Position (integer)}
#'     \item{ref}{Reference allele}
#'     \item{alt}{Alternate allele (multi-allelic sites are split into rows)}
#'     \item{qual}{QUAL score (numeric)}
#'     \item{filter}{FILTER field}
#'     \item{sample}{Sample name (NA for single-sample VCFs without GT field)}
#'     \item{consequence}{Variant consequence if ANN/CSQ INFO field is present}
#'     \item{gene}{Gene symbol if ANN/CSQ INFO field is present}
#'   }
#'
#' @examples
#' vcf_file <- system.file("extdata", "example.vcf", package = "ggvariant")
#' variants <- read_vcf(vcf_file)
#' head(variants)
#'
#' @seealso [coerce_variants()], [plot_lollipop()], [plot_consequence_summary()]
#' @export
read_vcf <- function(path,
                     samples    = NULL,
                     pass_only  = TRUE,
                     info_fields = NULL) {

  path <- normalizePath(path, mustWork = TRUE)
  cli::cli_progress_step("Reading VCF: {.file {basename(path)}}")

  lines <- .read_vcf_lines(path)
  header_lines <- grep("^#", lines, value = TRUE)
  data_lines   <- lines[!grepl("^#", lines)]

  if (length(data_lines) == 0L) {
    cli::cli_warn("VCF file contains no variant records.")
    return(.empty_gvf())
  }

  # Parse sample names from #CHROM header
  col_header <- tail(header_lines[grepl("^#CHROM", header_lines)], 1)
  col_names  <- strsplit(sub("^#", "", col_header), "\t")[[1]]
  fixed_cols <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")
  sample_names_all <- setdiff(col_names, c(fixed_cols, "FORMAT"))

  # Parse body
  mat <- do.call(rbind, strsplit(data_lines, "\t"))
  colnames(mat) <- col_names[seq_len(ncol(mat))]
  df <- as.data.frame(mat, stringsAsFactors = FALSE)

  # Apply PASS filter
  if (pass_only) {
    keep <- df$FILTER %in% c("PASS", ".", "")
    n_removed <- sum(!keep)
    df <- df[keep, , drop = FALSE]
    if (n_removed > 0)
      cli::cli_inform("Removed {n_removed} non-PASS variant{?s}.")
  }

  # Split multi-allelic ALT into separate rows
  df <- .split_multiallelic(df)

  # Classify consequence from REF/ALT if no ANN present
  df$consequence <- .infer_consequence(df$REF, df$ALT)

  # Tidy columns
  out <- data.frame(
    chrom       = df$CHROM,
    pos         = as.integer(df$POS),
    ref         = df$REF,
    alt         = df$ALT,
    qual        = suppressWarnings(as.numeric(df$QUAL)),
    filter      = df$FILTER,
    consequence = df$consequence,
    gene        = NA_character_,
    sample      = NA_character_,
    stringsAsFactors = FALSE
  )

  # Expand ANN/CSQ if present
  if ("INFO" %in% colnames(df)) {
    out <- .parse_ann_csq(out, df$INFO)
  }

  # Pivot samples (if FORMAT/GT columns present)
  if (length(sample_names_all) > 0) {
    fmt_col <- if ("FORMAT" %in% colnames(df)) df$FORMAT else NULL
    out <- .pivot_samples(out, df, sample_names_all, fmt_col, samples)
  }

  class(out) <- c("gvf", "data.frame")
  cli::cli_progress_done()
  cli::cli_inform(
    "Loaded {nrow(out)} variant record{?s} across \\
     {length(unique(out$chrom))} chromosome{?s}."
  )
  out
}


#' Coerce a plain data frame to a `gvf` object
#'
#' If you already have variant data in a `data.frame` (e.g. exported from
#' Excel, a database, or another tool), use this function to prepare it for
#' use with `ggvariant` plotting functions.
#'
#' @param x A `data.frame` or `tibble`.
#' @param chrom Column name containing chromosome (default `"chrom"`).
#' @param pos Column name containing position (default `"pos"`).
#' @param ref Column name containing reference allele (default `"ref"`).
#' @param alt Column name containing alternate allele (default `"alt"`).
#' @param consequence Column name containing variant consequence annotation,
#'   e.g. `"Missense_Mutation"`. If `NULL`, consequence is inferred from
#'   REF/ALT lengths.
#' @param gene Column name containing gene symbol (default `"gene"`).
#' @param sample Column name containing sample identifier (default `"sample"`).
#'
#' @return A `gvf` object.
#'
#' @examples
#' df <- data.frame(
#'   chromosome = c("chr1", "chr1", "chr7"),
#'   position   = c(100200, 100350, 55249071),
#'   ref_allele = c("A", "G", "C"),
#'   alt_allele = c("T", "A", "T"),
#'   variant_class = c("missense_variant", "synonymous_variant", "missense_variant"),
#'   hugo_symbol = c("GENE1", "GENE1", "EGFR"),
#'   tumor_sample = c("S1", "S2", "S2")
#' )
#'
#' variants <- coerce_variants(df,
#'   chrom       = "chromosome",
#'   pos         = "position",
#'   ref         = "ref_allele",
#'   alt         = "alt_allele",
#'   consequence = "variant_class",
#'   gene        = "hugo_symbol",
#'   sample      = "tumor_sample"
#' )
#'
#' @export
coerce_variants <- function(x,
                            chrom       = "chrom",
                            pos         = "pos",
                            ref         = "ref",
                            alt         = "alt",
                            consequence = "consequence",
                            gene        = "gene",
                            sample      = "sample") {
  x <- as.data.frame(x)

  .check_col <- function(col, required = TRUE) {
    if (!is.null(col) && !col %in% colnames(x)) {
      if (required)
        cli::cli_abort("Column {.val {col}} not found in data frame.")
      else
        return(NULL)
    }
    col
  }

  chrom       <- .check_col(chrom)
  pos         <- .check_col(pos)
  ref         <- .check_col(ref)
  alt         <- .check_col(alt)
  consequence <- .check_col(consequence, required = FALSE)
  gene        <- .check_col(gene,        required = FALSE)
  sample      <- .check_col(sample,      required = FALSE)

  out <- data.frame(
    chrom  = x[[chrom]],
    pos    = as.integer(x[[pos]]),
    ref    = x[[ref]],
    alt    = x[[alt]],
    stringsAsFactors = FALSE
  )

  out$consequence <- if (!is.null(consequence)) x[[consequence]] else
    .infer_consequence(out$ref, out$alt)
  out$gene   <- if (!is.null(gene))   x[[gene]]   else NA_character_
  out$sample <- if (!is.null(sample)) x[[sample]] else NA_character_

  # Carry over any extra columns
  extra <- setdiff(colnames(x),
                   c(chrom, pos, ref, alt, consequence, gene, sample))
  if (length(extra)) out <- cbind(out, x[, extra, drop = FALSE])

  class(out) <- c("gvf", "data.frame")
  out
}


# ── Internal helpers ──────────────────────────────────────────────────────────

.read_vcf_lines <- function(path) {
  if (grepl("\\.gz$", path)) {
    con <- gzcon(file(path, "rb"))
    on.exit(close(con))
    readLines(con, warn = FALSE)
  } else {
    readLines(path, warn = FALSE)
  }
}

.empty_gvf <- function() {
  out <- data.frame(
    chrom = character(), pos = integer(), ref = character(),
    alt = character(), qual = numeric(), filter = character(),
    consequence = character(), gene = character(), sample = character(),
    stringsAsFactors = FALSE
  )
  class(out) <- c("gvf", "data.frame")
  out
}

.split_multiallelic <- function(df) {
  needs_split <- grepl(",", df$ALT)
  if (!any(needs_split)) return(df)
  singles <- df[!needs_split, , drop = FALSE]
  multi   <- df[needs_split,  , drop = FALSE]
  split_rows <- lapply(seq_len(nrow(multi)), function(i) {
    row <- multi[i, , drop = FALSE]
    alts <- strsplit(row$ALT, ",")[[1]]
    rows <- row[rep(1, length(alts)), , drop = FALSE]
    rows$ALT <- alts
    rows
  })
  rbind(singles, do.call(rbind, split_rows))
}

.infer_consequence <- function(ref, alt) {
  ref_len <- nchar(ref)
  alt_len <- nchar(alt)
  ifelse(ref_len == 1 & alt_len == 1, "SNV",
    ifelse(ref_len > alt_len, "deletion",
      ifelse(alt_len > ref_len, "insertion", "MNV")))
}

.parse_ann_csq <- function(out, info_vec) {
  # Attempt to parse VEP CSQ or SnpEff ANN fields
  has_ann <- grepl("ANN=|CSQ=", info_vec)
  if (!any(has_ann)) return(out)

  parsed <- vapply(info_vec, function(info) {
    m <- regmatches(info, regexpr("(?:ANN|CSQ)=([^;]+)", info))
    if (length(m) == 0) return(c(gene = NA_character_, csq = NA_character_))
    val   <- sub("^(?:ANN|CSQ)=", "", m)
    first <- strsplit(val, ",")[[1]][1]
    parts <- strsplit(first, "\\|")[[1]]
    # ANN: [1]=allele [2]=effect [3]=impact [4]=gene
    # CSQ: order varies; we pick conservatively
    gene <- if (length(parts) >= 4) parts[4] else NA_character_
    csq  <- if (length(parts) >= 2) parts[2] else NA_character_
    c(gene = gene, csq = csq)
  }, character(2))

  out$gene[is.na(out$gene)] <- parsed["gene", is.na(out$gene)]
  csq <- parsed["csq", ]
  out$consequence <- ifelse(!is.na(csq), csq, out$consequence)
  out
}

.pivot_samples <- function(out, df, sample_names_all, fmt_col, keep_samples) {
  if (!is.null(keep_samples))
    sample_names_all <- intersect(sample_names_all, keep_samples)
  if (length(sample_names_all) == 0) return(out)

  rows <- lapply(sample_names_all, function(sname) {
    s_col <- df[[sname]]
    present <- !s_col %in% c("./.", ".", NA)
    sub <- out[present, , drop = FALSE]
    sub$sample <- sname
    sub
  })
  do.call(rbind, rows)
}

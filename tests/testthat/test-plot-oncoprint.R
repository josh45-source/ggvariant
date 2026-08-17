vcf_path <- system.file("extdata", "example.vcf", package = "ggvariant")

test_that("plot_oncoprint orders genes by descending sample frequency", {
  skip_if(!nzchar(vcf_path))
  vf <- read_vcf(vcf_path)
  ranked <- ggvariant:::.rank_genes(vf, unique(vf$gene))
  # TP53 is mutated in both samples across 4 records -> most-altered gene
  expect_equal(ranked[1], "TP53")
  # frequency must be non-increasing along the ranking
  n_samples <- vapply(ranked, function(g) {
    length(unique(vf$sample[vf$gene == g]))
  }, integer(1))
  expect_true(all(diff(n_samples) <= 0))
})

test_that("plot_oncoprint memo-sorts samples into a deterministic staircase", {
  skip_if(!nzchar(vcf_path))
  vf <- read_vcf(vcf_path)
  genes <- c("TP53", "BRCA1", "BRCA2")
  order1 <- ggvariant:::.memo_sort_samples(
    vf[vf$gene %in% genes, ], genes, c("TUMOR_S1", "TUMOR_S2")
  )
  order2 <- ggvariant:::.memo_sort_samples(
    vf[vf$gene %in% genes, ], genes, c("TUMOR_S2", "TUMOR_S1")
  )
  # sample_universe order must not affect the resulting sort
  expect_equal(order1, order2)
  expect_length(order1, 2L)
})

test_that("plot_oncoprint marks a gene mutated more than once per sample as Multi_Hit", {
  skip_if(!nzchar(vcf_path))
  vf <- read_vcf(vcf_path)
  cells <- ggvariant:::.oncoprint_cells(vf, "TP53", c("TUMOR_S1", "TUMOR_S2"))
  # TP53 has 3 records for TUMOR_S1 and 2 for TUMOR_S2 in the example VCF
  expect_equal(cells$label[cells$sample == "TUMOR_S1"], "Multi_Hit")
  expect_equal(cells$label[cells$sample == "TUMOR_S2"], "Multi_Hit")
})

test_that("plot_oncoprint errors when both top_n and genes are supplied", {
  skip_if(!nzchar(vcf_path))
  vf <- read_vcf(vcf_path)
  expect_error(
    plot_oncoprint(vf, top_n = 3, genes = "TP53"),
    regexp = "only one of"
  )
})

test_that("plot_oncoprint handles a single sample", {
  skip_if(!nzchar(vcf_path))
  vf <- read_vcf(vcf_path)
  p <- plot_oncoprint(vf, top_n = 5, samples = "TUMOR_S1")
  expect_s3_class(p, "ggplot")
  built <- ggplot2::ggplot_build(p)
  expect_length(unique(built$data[[1]]$x), 1L)
})

test_that("plot_oncoprint handles a single gene", {
  skip_if(!nzchar(vcf_path))
  vf <- read_vcf(vcf_path)
  p <- plot_oncoprint(vf, genes = "TP53")
  expect_s3_class(p, "ggplot")
})

test_that("plot_oncoprint errors cleanly when no samples remain after filtering", {
  skip_if(!nzchar(vcf_path))
  vf <- read_vcf(vcf_path)
  expect_error(
    plot_oncoprint(vf, top_n = 3, samples = "NOT_A_SAMPLE"),
    regexp = "No samples remaining"
  )
})

test_that("a gene absent from the data is shown as a fully unaltered row", {
  skip_if(!nzchar(vcf_path))
  vf <- read_vcf(vcf_path)
  cells <- ggvariant:::.oncoprint_cells(vf, "NOTAGENE", c("TUMOR_S1", "TUMOR_S2"))
  expect_true(all(is.na(cells$label)))
})

test_that("plot_waterfall is an alias of plot_oncoprint", {
  skip_if(!nzchar(vcf_path))
  vf <- read_vcf(vcf_path)
  expect_identical(
    plot_oncoprint(vf, top_n = 5)$data,
    plot_waterfall(vf, top_n = 5)$data
  )
})

test_that("plot_oncoprint draws an annotation track from sample metadata", {
  skip_if(!nzchar(vcf_path))
  vf <- read_vcf(vcf_path)
  clinical <- data.frame(sample = c("TUMOR_S1", "TUMOR_S2"), stage = c("III", "IV"))
  p <- plot_oncoprint(vf, top_n = 5, annotation = clinical)
  built <- ggplot2::ggplot_build(p)
  expect_length(built$data, 3L)  # matrix tiles + annotation tiles + annotation text
})

test_that("plot_oncoprint errors when annotation lacks a sample column", {
  skip_if(!nzchar(vcf_path))
  vf <- read_vcf(vcf_path)
  expect_error(
    plot_oncoprint(vf, top_n = 3, annotation = data.frame(x = 1)),
    regexp = "sample"
  )
})

test_that("plot_oncoprint snapshot is stable", {
  skip_if(!nzchar(vcf_path))
  testthat::skip_if_not_installed("vdiffr")
  vf <- read_vcf(vcf_path)
  vdiffr::expect_doppelganger("oncoprint-top5", plot_oncoprint(vf, top_n = 5))
})

test_that("an unmapped consequence term renders as a visible Other cell, not NA", {
  df <- data.frame(
    chrom = "chr1", pos = 1:4, ref = "A", alt = "T",
    gene   = c("TP53", "TP53", "BRCA1", "BRCA1"),
    sample = c("S01", "S02", "S01", "S02"),
    consequence = c("missense_variant", "upstream_gene_variant",
                    "intergenic_variant", "stop_gained"),
    stringsAsFactors = FALSE
  )
  vf <- coerce_variants(df)
  p <- plot_oncoprint(vf, top_n = 2)
  built <- ggplot2::ggplot_build(p)

  # the two unmapped terms are standardised to "Other" ...
  expect_equal(sum(built$plot$data$label == "Other", na.rm = TRUE), 2L)
  # ... and get a real, non-missing fill colour, not the NA used for
  # genuinely unmutated cells
  other_fill <- built$data[[1]]$fill[built$plot$data$label == "Other"]
  expect_false(anyNA(other_fill))
  expect_equal(unique(other_fill), unname(gv_palette("consequence")["Other"]))

  # x-axis sample labels must still render regardless of the unmapped terms
  x_labels <- built$layout$panel_params[[1]]$x$get_labels()
  expect_setequal(x_labels, c("S01", "S02"))
})

test_that("gene row order is a factor tied to the frequency ranking, independent of input row order", {
  skip_if(!nzchar(vcf_path))
  vf <- read_vcf(vcf_path)

  expected_order <- ggvariant:::.rank_genes(vf, unique(vf$gene))
  expected_order <- utils::head(expected_order, 5)

  set.seed(123)
  vf_shuffled <- vf[sample(nrow(vf)), ]

  p_normal   <- plot_oncoprint(vf, top_n = 5)
  p_shuffled <- plot_oncoprint(vf_shuffled, top_n = 5)

  built_normal   <- ggplot2::ggplot_build(p_normal)
  built_shuffled <- ggplot2::ggplot_build(p_shuffled)

  expect_s3_class(built_normal$plot$data$gene, "factor")
  expect_equal(levels(built_normal$plot$data$gene), expected_order)
  expect_equal(levels(built_shuffled$plot$data$gene), expected_order)

  # TP53 (most-altered) must sit at the top (max y); the least-altered
  # displayed gene must sit at the bottom (min y), in both orderings
  d <- built_shuffled$plot$data
  expect_equal(as.character(d$gene[d$y == max(d$y)][1]), expected_order[1])
  expect_equal(as.character(d$gene[d$y == min(d$y)][1]), utils::tail(expected_order, 1))
})

test_that("plot_oncoprint aborts exactly once, cleanly, on genuinely empty input", {
  skip_if(!nzchar(vcf_path))
  vf <- read_vcf(vcf_path)
  expect_error(plot_oncoprint(vf[0, ]), class = "rlang_error")
  err <- tryCatch(plot_oncoprint(vf[0, ]), error = function(e) e)
  expect_null(err$parent)
})

test_that("subsetting the input to a single sample still plots", {
  skip_if(!nzchar(vcf_path))
  vf <- read_vcf(vcf_path)
  p <- plot_oncoprint(vf[vf$sample == "TUMOR_S1", ])
  expect_s3_class(p, "ggplot")
})

test_that("subsetting the input to a single gene still plots", {
  skip_if(!nzchar(vcf_path))
  vf <- read_vcf(vcf_path)
  p <- plot_oncoprint(vf[vf$gene == "TP53", ])
  expect_s3_class(p, "ggplot")
})

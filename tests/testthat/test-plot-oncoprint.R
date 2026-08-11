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

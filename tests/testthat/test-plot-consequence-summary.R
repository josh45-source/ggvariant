vcf_path <- system.file("extdata", "example.vcf", package = "ggvariant")

test_that("plot_consequence_summary group_by = 'gene' shows top_n genes", {
  skip_if(!nzchar(vcf_path))
  vf <- read_vcf(vcf_path)
  p <- plot_consequence_summary(vf, group_by = "gene", top_n = 3)
  expect_s3_class(p, "ggplot")
  expect_length(levels(p$data$gene), 3L)
})

test_that("plot_consequence_summary group_by = 'gene' errors with no gene column", {
  df <- data.frame(pos = 1:3, ref = "A", alt = "T", consequence = "SNV")
  expect_error(
    plot_consequence_summary(df, group_by = "gene"),
    regexp = "No.*gene.*column"
  )
})

test_that("plot_consequence_summary supports position = 'fill'", {
  skip_if(!nzchar(vcf_path))
  vf <- read_vcf(vcf_path)
  p <- plot_consequence_summary(vf, position = "fill")
  expect_equal(p$labels$y, "Proportion")
})

test_that("plot_consequence_summary supports position = 'dodge'", {
  skip_if(!nzchar(vcf_path))
  vf <- read_vcf(vcf_path)
  p <- plot_consequence_summary(vf, position = "dodge")
  expect_s3_class(p, "ggplot")
  built <- ggplot2::ggplot_build(p)
  expect_s3_class(built$plot$layers[[1]]$position, "PositionDodge")
})

test_that("plot_consequence_summary applies the samples filter", {
  skip_if(!nzchar(vcf_path))
  vf <- read_vcf(vcf_path)
  p <- plot_consequence_summary(vf, samples = "TUMOR_S1")
  expect_equal(levels(droplevels(p$data$sample)), "TUMOR_S1")
})

test_that("plot_consequence_summary flip = TRUE flips coordinates", {
  skip_if(!nzchar(vcf_path))
  vf <- read_vcf(vcf_path)
  p <- plot_consequence_summary(vf, flip = TRUE)
  expect_s3_class(p$coordinates, "CoordFlip")
})

test_that("plot_consequence_summary accepts a custom palette", {
  skip_if(!nzchar(vcf_path))
  vf <- read_vcf(vcf_path)
  pal <- c(missense_variant = "#111111", frameshift_variant = "#222222",
           stop_gained = "#333333", synonymous_variant = "#444444",
           splice_site_variant = "#555555")
  p <- plot_consequence_summary(vf, palette = pal)
  built <- ggplot2::ggplot_build(p)
  expect_true(all(built$data[[1]]$fill %in% pal))
})

test_that("plot_consequence_summary snapshot is stable", {
  skip_if(!nzchar(vcf_path))
  testthat::skip_if_not_installed("vdiffr")
  vf <- read_vcf(vcf_path)
  vdiffr::expect_doppelganger(
    "consequence-summary-default", plot_consequence_summary(vf)
  )
})

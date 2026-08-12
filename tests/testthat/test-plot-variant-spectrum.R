vcf_path <- system.file("extdata", "example.vcf", package = "ggvariant")

test_that("plot_variant_spectrum applies the sample filter", {
  skip_if(!nzchar(vcf_path))
  vf <- read_vcf(vcf_path)
  p <- plot_variant_spectrum(vf, sample = "TUMOR_S1")
  expect_s3_class(p, "ggplot")
  expect_lte(sum(p$data$n), sum(vf$sample == "TUMOR_S1"))
})

test_that("plot_variant_spectrum normalize = FALSE plots raw counts", {
  skip_if(!nzchar(vcf_path))
  vf <- read_vcf(vcf_path)
  p <- plot_variant_spectrum(vf, normalize = FALSE)
  expect_equal(p$labels$y, "Count")
  expect_true("n" %in% colnames(p$data))
  expect_false("freq" %in% colnames(p$data))
})

test_that("plot_variant_spectrum errors when no SNVs remain after filtering", {
  df <- data.frame(
    pos = 1:2, ref = c("AT", "GC"), alt = c("A", "G"),
    consequence = "deletion", stringsAsFactors = FALSE
  )
  expect_error(plot_variant_spectrum(df), regexp = "No SNVs remaining")
})

test_that("plot_variant_spectrum reports dropped non-SNV records", {
  skip_if(!nzchar(vcf_path))
  vf <- read_vcf(vcf_path)
  expect_message(plot_variant_spectrum(vf), regexp = "Excluded")
})

test_that("plot_variant_spectrum snapshot is stable", {
  skip_if(!nzchar(vcf_path))
  testthat::skip_if_not_installed("vdiffr")
  vf <- read_vcf(vcf_path)
  vdiffr::expect_doppelganger("variant-spectrum-default", plot_variant_spectrum(vf))
})

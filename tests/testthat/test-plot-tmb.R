vcf_path <- system.file("extdata", "example.vcf", package = "ggvariant")

test_that("plot_tmb orders samples by descending TMB by default", {
  skip_if(!nzchar(vcf_path))
  vf <- read_vcf(vcf_path)
  p <- plot_tmb(vf)
  tmb <- ggvariant:::.compute_tmb(vf)
  expect_equal(levels(p$data$sample), tmb$sample[order(-tmb$tmb, tmb$sample)])
})

test_that("plot_tmb normalises by mb_size", {
  skip_if(!nzchar(vcf_path))
  vf <- read_vcf(vcf_path)
  raw  <- ggvariant:::.compute_tmb(vf)
  norm <- ggvariant:::.compute_tmb(vf, mb_size = 10)
  expect_equal(norm$tmb, raw$n / 10)
})

test_that("plot_tmb rejects a non-positive mb_size", {
  skip_if(!nzchar(vcf_path))
  vf <- read_vcf(vcf_path)
  expect_error(plot_tmb(vf, mb_size = -5), regexp = "positive")
  expect_error(plot_tmb(vf, mb_size = 0), regexp = "positive")
})

test_that("plot_tmb honours an explicit sample_order", {
  skip_if(!nzchar(vcf_path))
  vf <- read_vcf(vcf_path)
  p <- plot_tmb(vf, sample_order = c("TUMOR_S2", "TUMOR_S1"))
  expect_equal(levels(p$data$sample), c("TUMOR_S2", "TUMOR_S1"))
})

test_that("plot_tmb errors when sample_order references an unknown sample", {
  skip_if(!nzchar(vcf_path))
  vf <- read_vcf(vcf_path)
  expect_error(
    plot_tmb(vf, sample_order = c("TUMOR_S1", "NOPE")),
    regexp = "not present"
  )
})

test_that("plot_tmb handles a single sample", {
  skip_if(!nzchar(vcf_path))
  vf <- read_vcf(vcf_path)
  p <- plot_tmb(vf, samples = "TUMOR_S1")
  expect_s3_class(p, "ggplot")
  expect_equal(nrow(p$data), 1L)
})

test_that("plot_tmb errors cleanly when no mutations remain after filtering", {
  skip_if(!nzchar(vcf_path))
  vf <- read_vcf(vcf_path)
  expect_error(
    plot_tmb(vf, samples = "NOT_A_SAMPLE"),
    regexp = "No mutations remaining"
  )
})

test_that(".compute_tmb agrees with a manual per-sample tally", {
  skip_if(!nzchar(vcf_path))
  vf <- read_vcf(vcf_path)
  tmb <- ggvariant:::.compute_tmb(vf)
  manual <- table(vf$sample)
  expect_equal(
    tmb$n[match(names(manual), tmb$sample)],
    as.integer(manual)
  )
})

test_that("plot_tmb snapshot is stable", {
  skip_if(!nzchar(vcf_path))
  testthat::skip_if_not_installed("vdiffr")
  vf <- read_vcf(vcf_path)
  vdiffr::expect_doppelganger("tmb-default", plot_tmb(vf))
})

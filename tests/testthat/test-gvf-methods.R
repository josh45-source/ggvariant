vcf_path <- system.file("extdata", "example.vcf", package = "ggvariant")

test_that("print.gvf shows a compact header and truncated preview", {
  skip_if(!nzchar(vcf_path))
  vf <- read_vcf(vcf_path)
  out <- capture.output(print(vf))
  expected_header <- sprintf(
    "<gvf: %d variants, %d samples, %d chromosomes, %d genes>",
    nrow(vf), length(unique(vf$sample)), length(unique(vf$chrom)),
    length(unique(vf$gene))
  )
  expect_equal(out[1], expected_header)
  expect_true(any(grepl("more row", out)))
})

test_that("print.gvf respects the n argument", {
  skip_if(!nzchar(vcf_path))
  vf <- read_vcf(vcf_path)
  out <- capture.output(print(vf, n = 2))
  expect_true(any(grepl(paste0(nrow(vf) - 2, " more rows"), out, fixed = TRUE)))
})

test_that("print.gvf returns its input invisibly", {
  skip_if(!nzchar(vcf_path))
  vf <- read_vcf(vcf_path)
  ret <- withVisible(print(vf))
  expect_false(ret$visible)
  expect_identical(ret$value, vf)
})

test_that("summary.gvf returns consequence, sample, and chromosome breakdowns", {
  skip_if(!nzchar(vcf_path))
  vf <- read_vcf(vcf_path)
  s <- summary(vf)
  expect_s3_class(s, "summary.gvf")
  expect_equal(s$n_variants, nrow(vf))
  expect_equal(sum(s$consequence), nrow(vf))
  expect_equal(sum(s$per_sample), nrow(vf))
  expect_equal(sum(s$chromosome), nrow(vf))
  expect_equal(
    unname(s$per_sample[["TUMOR_S1"]]), sum(vf$sample == "TUMOR_S1")
  )
})

test_that("print.summary.gvf prints all three breakdown sections", {
  skip_if(!nzchar(vcf_path))
  vf <- read_vcf(vcf_path)
  out <- capture.output(print(summary(vf)))
  expect_true(any(grepl("Consequence breakdown", out)))
  expect_true(any(grepl("Per-sample counts", out)))
  expect_true(any(grepl("Chromosome distribution", out)))
})

test_that("validate_gvf aborts on a gvf object missing a required column", {
  skip_if(!nzchar(vcf_path))
  vf <- read_vcf(vcf_path)
  bad <- vf
  bad$gene <- NULL
  class(bad) <- c("gvf", "data.frame")
  expect_error(ggvariant:::validate_gvf(bad), regexp = "missing column")
})

test_that("validate_gvf aborts on a non-numeric pos column", {
  skip_if(!nzchar(vcf_path))
  vf <- read_vcf(vcf_path)
  bad <- vf
  bad$pos <- as.character(bad$pos)
  class(bad) <- c("gvf", "data.frame")
  expect_error(ggvariant:::validate_gvf(bad), regexp = "pos.*numeric")
})

test_that("read_vcf() output passes validate_gvf()", {
  skip_if(!nzchar(vcf_path))
  vf <- read_vcf(vcf_path)
  expect_silent(ggvariant:::validate_gvf(vf))
})

test_that("coerce_variants() output passes validate_gvf()", {
  df <- data.frame(chrom2 = "chr1", pos2 = 100L, ref2 = "A", alt2 = "T")
  cv <- coerce_variants(df, chrom = "chrom2", pos = "pos2",
                         ref = "ref2", alt = "alt2")
  expect_silent(ggvariant:::validate_gvf(cv))
})

test_that("print.gvf and summary.gvf handle a single-row gvf", {
  df <- data.frame(chrom2 = "chr1", pos2 = 100L, ref2 = "A", alt2 = "T")
  cv <- coerce_variants(df, chrom = "chrom2", pos = "pos2",
                         ref = "ref2", alt = "alt2")
  expect_output(print(cv), "<gvf: 1 variant,", fixed = TRUE)
  s <- summary(cv)
  expect_equal(s$n_variants, 1L)
})

test_that("read_vcf returns a gvf object", {
  vcf_path <- system.file("extdata", "example.vcf", package = "ggvariant")
  skip_if(!nzchar(vcf_path), "example VCF not found")
  vf <- read_vcf(vcf_path)
  expect_s3_class(vf, "gvf")
  expect_true(all(c("chrom", "pos", "ref", "alt", "consequence") %in% colnames(vf)))
  expect_gt(nrow(vf), 0)
})

test_that("read_vcf filters non-PASS variants", {
  vcf_path <- system.file("extdata", "example.vcf", package = "ggvariant")
  skip_if(!nzchar(vcf_path))
  vf <- read_vcf(vcf_path, pass_only = TRUE)
  expect_true(all(vf$filter %in% c("PASS", ".", "")))
})

test_that("coerce_variants works with custom column names", {
  df <- data.frame(
    chrom2 = c("chr1", "chr2"),
    pos2   = c(1000L, 2000L),
    ref2   = c("A", "G"),
    alt2   = c("T", "C"),
    stringsAsFactors = FALSE
  )
  vf <- coerce_variants(df,
    chrom = "chrom2", pos = "pos2", ref = "ref2", alt = "alt2"
  )
  expect_s3_class(vf, "gvf")
  expect_equal(nrow(vf), 2L)
  expect_equal(vf$ref, c("A", "G"))
})

test_that("coerce_variants errors on missing required columns", {
  df <- data.frame(a = 1, b = 2)
  expect_error(coerce_variants(df), regexp = "not found")
})

test_that("plot_lollipop returns a ggplot", {
  vcf_path <- system.file("extdata", "example.vcf", package = "ggvariant")
  skip_if(!nzchar(vcf_path))
  vf <- read_vcf(vcf_path)
  p  <- plot_lollipop(vf, gene = "TP53")
  expect_s3_class(p, "ggplot")
})

test_that("plot_consequence_summary returns a ggplot", {
  vcf_path <- system.file("extdata", "example.vcf", package = "ggvariant")
  skip_if(!nzchar(vcf_path))
  vf <- read_vcf(vcf_path)
  p  <- plot_consequence_summary(vf)
  expect_s3_class(p, "ggplot")
})

test_that("plot_variant_spectrum returns a ggplot", {
  vcf_path <- system.file("extdata", "example.vcf", package = "ggvariant")
  skip_if(!nzchar(vcf_path))
  vf <- read_vcf(vcf_path)
  p  <- plot_variant_spectrum(vf)
  expect_s3_class(p, "ggplot")
})

test_that("gv_palette returns named colour vectors", {
  pal <- gv_palette("consequence")
  expect_type(pal, "character")
  expect_named(pal)
  expect_true(all(grepl("^#[0-9A-Fa-f]{6}$", pal)))
})

test_that(".classify_substitution pyrimidine-normalises purines", {
  # A>G is the reverse complement of C>T
  expect_equal(ggvariant:::.classify_substitution("A", "G"), "T>C")
  expect_equal(ggvariant:::.classify_substitution("C", "T"), "C>T")
})

test_that(".infer_consequence classifies correctly", {
  expect_equal(ggvariant:::.infer_consequence("A", "T"), "SNV")
  expect_equal(ggvariant:::.infer_consequence("AT", "A"), "deletion")
  expect_equal(ggvariant:::.infer_consequence("A", "ATG"), "insertion")
})

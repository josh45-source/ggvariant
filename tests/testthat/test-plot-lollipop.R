vcf_path <- system.file("extdata", "example.vcf", package = "ggvariant")

test_that("plot_lollipop draws domain annotations", {
  skip_if(!nzchar(vcf_path))
  vf <- read_vcf(vcf_path)
  domains <- data.frame(
    name  = c("Transactivation", "DNA-binding"),
    start = c(1, 102),
    end   = c(67, 292)
  )
  p <- plot_lollipop(vf, gene = "TP53", domains = domains)
  expect_s3_class(p, "ggplot")
  built <- ggplot2::ggplot_build(p)
  expect_gt(length(built$data), 2L)  # stems + dots + domain rect + domain text
})

test_that("plot_lollipop errors on a malformed domains data frame", {
  skip_if(!nzchar(vcf_path))
  vf <- read_vcf(vcf_path)
  expect_error(
    plot_lollipop(vf, gene = "TP53", domains = data.frame(x = 1)),
    regexp = "domains"
  )
})

test_that("plot_lollipop colours by sample when color_by = 'sample'", {
  skip_if(!nzchar(vcf_path))
  vf <- read_vcf(vcf_path)
  sample_pal <- c(TUMOR_S1 = "#1B9E77", TUMOR_S2 = "#D95F02")
  p <- plot_lollipop(vf, gene = "TP53", color_by = "sample",
                      palette = sample_pal)
  built <- ggplot2::ggplot_build(p)
  expect_true(all(built$data[[2]]$colour %in% sample_pal))
  colour_scale <- Filter(function(s) "colour" %in% s$aesthetics, p$scales$scales)
  expect_equal(colour_scale[[1]]$name, "Sample")
})

test_that("plot_lollipop errors when color_by references a missing column", {
  skip_if(!nzchar(vcf_path))
  vf <- read_vcf(vcf_path)
  expect_error(
    plot_lollipop(vf, gene = "TP53", color_by = "not_a_column"),
    regexp = "not_a_column"
  )
})

test_that("plot_lollipop errors when filtering removes every variant", {
  skip_if(!nzchar(vcf_path))
  vf <- read_vcf(vcf_path)
  expect_error(
    plot_lollipop(vf, gene = "NOTAGENE"),
    regexp = "No variants remaining"
  )
})

test_that("plot_lollipop with stack_dots = FALSE does not stack y positions", {
  skip_if(!nzchar(vcf_path))
  vf <- read_vcf(vcf_path)
  p <- plot_lollipop(vf, gene = "TP53", stack_dots = FALSE)
  expect_true(all(p$data$y_pos == 1L))
})

test_that("plot_lollipop accepts a custom palette", {
  skip_if(!nzchar(vcf_path))
  vf <- read_vcf(vcf_path)
  custom_pal <- c(missense_variant = "#123456", stop_gained = "#654321",
                   synonymous_variant = "#ABCDEF")
  p <- plot_lollipop(vf, gene = "TP53", palette = custom_pal)
  built <- ggplot2::ggplot_build(p)
  expect_true(all(built$data[[2]]$colour %in% c(custom_pal, "grey50")))
})

test_that("plot_lollipop accepts a custom title", {
  skip_if(!nzchar(vcf_path))
  vf <- read_vcf(vcf_path)
  p <- plot_lollipop(vf, gene = "TP53", title = "My custom title")
  expect_equal(p$labels$title, "My custom title")
})

test_that("plot_lollipop labels the x-axis as genomic position by default", {
  skip_if(!nzchar(vcf_path))
  vf <- read_vcf(vcf_path)
  p <- plot_lollipop(vf, gene = "TP53")
  # pos is a raw genomic coordinate unless the caller has already rescaled
  # it themselves; labelling it "Amino acid position" would be factually
  # wrong in the default case, so the axis must read "Genomic position".
  expect_equal(p$labels$x, "Genomic position")
})

test_that("plot_lollipop labels the x-axis as amino acid position when protein_length is supplied", {
  skip_if(!nzchar(vcf_path))
  vf <- read_vcf(vcf_path)
  # An explicit protein_length is the caller's signal that pos has already
  # been rescaled into protein coordinates -- only then is "Amino acid
  # position" accurate.
  p <- plot_lollipop(vf, gene = "TP53", protein_length = 393)
  expect_equal(p$labels$x, "Amino acid position")
})

test_that("plot_lollipop snapshot is stable", {
  skip_if(!nzchar(vcf_path))
  testthat::skip_if_not_installed("vdiffr")
  vf <- read_vcf(vcf_path)
  vdiffr::expect_doppelganger("lollipop-tp53", plot_lollipop(vf, gene = "TP53"))
})

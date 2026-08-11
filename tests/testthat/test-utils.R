test_that("gv_palette('spectrum') returns the 6 COSMIC SBS classes", {
  pal <- gv_palette("spectrum")
  expect_type(pal, "character")
  expect_named(pal)
  expect_setequal(names(pal), c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"))
  expect_true(all(grepl("^#[0-9A-Fa-f]{6}$", pal)))
})

test_that("gv_palette('domain') generates n colours", {
  pal5  <- gv_palette("domain", n = 5)
  pal15 <- gv_palette("domain", n = 15)
  expect_length(pal5, 5L)
  expect_length(pal15, 15L)   # beyond the 10 base colours -> ramp-generated
  expect_true(all(grepl("^#[0-9A-Fa-f]{6}$", pal5)))
  expect_true(all(grepl("^#[0-9A-Fa-f]{6}$", pal15)))
})

test_that("theme_ggvariant returns a ggplot2 theme", {
  th <- theme_ggvariant()
  expect_s3_class(th, "theme")
  expect_s3_class(th, "gg")
})

test_that("theme_ggvariant honours base_size and base_family", {
  th <- theme_ggvariant(base_size = 16, base_family = "serif")
  expect_equal(th$text$size, 16)
  expect_equal(th$text$family, "serif")
})

test_that("theme_ggvariant composes onto a ggplot", {
  library(ggplot2)
  p <- ggplot(mtcars, aes(mpg, wt)) + geom_point() + theme_ggvariant()
  expect_s3_class(p, "ggplot")
  built <- ggplot2::ggplot_build(p)
  expect_s3_class(built, "ggplot_built")
})

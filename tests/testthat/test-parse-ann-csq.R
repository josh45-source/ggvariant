test_that(".parse_ann_csq extracts gene and consequence from ANN", {
  out <- data.frame(gene = NA_character_, consequence = "SNV",
                     stringsAsFactors = FALSE)
  res <- ggvariant:::.parse_ann_csq(
    out, "ANN=T|missense_variant|MODERATE|TP53"
  )
  expect_equal(res$gene, "TP53")
  expect_equal(res$consequence, "missense_variant")
})

test_that(".parse_ann_csq extracts gene and consequence from CSQ", {
  out <- data.frame(gene = NA_character_, consequence = "SNV",
                     stringsAsFactors = FALSE)
  res <- ggvariant:::.parse_ann_csq(
    out, "CSQ=A|frameshift_variant|HIGH|BRCA2"
  )
  expect_equal(res$gene, "BRCA2")
  expect_equal(res$consequence, "frameshift_variant")
})

test_that(".parse_ann_csq takes only the first comma-separated annotation", {
  out <- data.frame(gene = NA_character_, consequence = "SNV",
                     stringsAsFactors = FALSE)
  res <- ggvariant:::.parse_ann_csq(
    out,
    "ANN=T|missense_variant|MODERATE|TP53,A|synonymous_variant|LOW|TP53"
  )
  expect_equal(res$gene, "TP53")
  expect_equal(res$consequence, "missense_variant")
})

test_that(".parse_ann_csq leaves rows with no ANN/CSQ field unchanged", {
  out <- data.frame(gene = NA_character_, consequence = "SNV",
                     stringsAsFactors = FALSE)
  res <- ggvariant:::.parse_ann_csq(out, "DP=45;AF=0.5")
  expect_true(is.na(res$gene))
  expect_equal(res$consequence, "SNV")
})

test_that(".parse_ann_csq handles a mixed vector, including NA INFO", {
  out <- data.frame(
    gene = rep(NA_character_, 4),
    consequence = rep("SNV", 4),
    stringsAsFactors = FALSE
  )
  info <- c(
    "ANN=T|missense_variant|MODERATE|TP53",
    "DP=10",
    NA_character_,
    "ANN=A|stop_gained|HIGH|BRCA1"
  )
  res <- ggvariant:::.parse_ann_csq(out, info)
  expect_equal(res$gene, c("TP53", NA, NA, "BRCA1"))
  expect_equal(res$consequence, c("missense_variant", "SNV", "SNV", "stop_gained"))
})

test_that("read_vcf() gene/consequence extraction matches the bundled example", {
  vcf_path <- system.file("extdata", "example.vcf", package = "ggvariant")
  skip_if(!nzchar(vcf_path))
  vf <- read_vcf(vcf_path)
  expect_equal(sum(vf$gene == "TP53"), 8L)
  expect_equal(sum(vf$consequence == "missense_variant"), 16L)
})

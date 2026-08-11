test_that("read_vcf errors clearly on a nonexistent path", {
  expect_error(read_vcf("does_not_exist_at_all.vcf"), regexp = "File not found")
})

test_that("read_vcf errors clearly when the #CHROM header is missing", {
  tmp <- tempfile(fileext = ".vcf")
  writeLines(c(
    "##fileformat=VCFv4.2",
    "chr1\t100\t.\tA\tT\t50\tPASS\t."
  ), tmp)
  on.exit(unlink(tmp))
  expect_error(read_vcf(tmp), regexp = "no.*#CHROM.*header")
})

test_that("read_vcf errors clearly on a data line with the wrong field count", {
  tmp <- tempfile(fileext = ".vcf")
  writeLines(c(
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "chr1\t100\t.\tA\tT\t50\tPASS\t.",
    "chr1\t200\t.\tA\tT\t50\tPASS"
  ), tmp)
  on.exit(unlink(tmp))
  expect_error(read_vcf(tmp), regexp = "wrong number of tab-separated fields")
})

test_that("read_vcf returns an empty gvf for a header-only VCF", {
  tmp <- tempfile(fileext = ".vcf")
  writeLines(c(
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
  ), tmp)
  on.exit(unlink(tmp))
  expect_warning(vf <- read_vcf(tmp), regexp = "no variant records")
  expect_s3_class(vf, "gvf")
  expect_equal(nrow(vf), 0L)
})

test_that("read_vcf reads a gzipped VCF", {
  tmp <- tempfile(fileext = ".vcf.gz")
  con <- gzfile(tmp, "w")
  writeLines(c(
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "chr1\t100\t.\tA\tT\t50\tPASS\tANN=T|missense_variant|MODERATE|TP53"
  ), con)
  close(con)
  on.exit(unlink(tmp))
  vf <- read_vcf(tmp)
  expect_equal(nrow(vf), 1L)
  expect_equal(vf$gene, "TP53")
})

test_that("read_vcf handles a VCF with no ANN/CSQ INFO field", {
  tmp <- tempfile(fileext = ".vcf")
  writeLines(c(
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "chr1\t100\t.\tA\tT\t50\tPASS\tDP=45"
  ), tmp)
  on.exit(unlink(tmp))
  vf <- read_vcf(tmp)
  expect_equal(nrow(vf), 1L)
  expect_true(is.na(vf$gene))
  expect_equal(vf$consequence, "SNV")
})

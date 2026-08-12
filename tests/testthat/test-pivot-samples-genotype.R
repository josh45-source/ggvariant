test_that(".gt_has_alt distinguishes reference from alternate genotypes", {
  expect_false(ggvariant:::.gt_has_alt("0/0"))
  expect_true(ggvariant:::.gt_has_alt("0/1"))
  expect_true(ggvariant:::.gt_has_alt("1/1"))
  expect_true(ggvariant:::.gt_has_alt("1/2"))
})

test_that(".gt_has_alt handles phased genotypes", {
  expect_true(ggvariant:::.gt_has_alt("0|1"))
  expect_false(ggvariant:::.gt_has_alt("0|0"))
})

test_that(".gt_has_alt treats missing genotypes as absent", {
  expect_false(ggvariant:::.gt_has_alt("./."))
  expect_false(ggvariant:::.gt_has_alt("."))
  expect_false(ggvariant:::.gt_has_alt(".|0"))
})

test_that(".gt_has_alt handles a partially-called genotype with an alt allele", {
  expect_true(ggvariant:::.gt_has_alt("1|."))
})

test_that(".gt_has_alt handles haploid calls", {
  expect_true(ggvariant:::.gt_has_alt("1"))
  expect_false(ggvariant:::.gt_has_alt("0"))
})

test_that(".gt_has_alt reads only the GT subfield of a multi-field FORMAT column", {
  expect_true(ggvariant:::.gt_has_alt("0/1:10,5:15:40"))
  expect_false(ggvariant:::.gt_has_alt("0/0:12,0:12:36"))
})

test_that("read_vcf() does not attach a variant to a homozygous-reference sample", {
  tmp <- tempfile(fileext = ".vcf")
  writeLines(c(
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2",
    "chr1\t100\t.\tA\tT\t50\tPASS\tANN=T|missense_variant|MODERATE|GENE1\tGT\t0/0\t0/1"
  ), tmp)
  on.exit(unlink(tmp))
  vf <- read_vcf(tmp)
  expect_equal(nrow(vf), 1L)
  expect_equal(vf$sample, "S2")
})

test_that("read_vcf() excludes samples with a missing genotype", {
  tmp <- tempfile(fileext = ".vcf")
  writeLines(c(
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2",
    "chr1\t100\t.\tA\tT\t50\tPASS\t.\tGT\t./.\t0/1"
  ), tmp)
  on.exit(unlink(tmp))
  vf <- read_vcf(tmp)
  expect_equal(nrow(vf), 1L)
  expect_equal(vf$sample, "S2")
})

test_that("read_vcf() correctly splits a multi-field FORMAT column", {
  tmp <- tempfile(fileext = ".vcf")
  writeLines(c(
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2",
    "chr1\t100\t.\tA\tT\t50\tPASS\t.\tGT:DP:GQ\t0/0:20:99\t0/1:18:80"
  ), tmp)
  on.exit(unlink(tmp))
  vf <- read_vcf(tmp)
  expect_equal(nrow(vf), 1L)
  expect_equal(vf$sample, "S2")
})

test_that("read_vcf() on the bundled example VCF assigns EGFR to the correct sample only", {
  vcf_path <- system.file("extdata", "example.vcf", package = "ggvariant")
  skip_if(!nzchar(vcf_path))
  vf <- read_vcf(vcf_path)
  egfr <- vf[vf$gene == "EGFR", ]
  expect_equal(nrow(egfr), 2L)
  expect_setequal(egfr$sample, c("TUMOR_S1", "TUMOR_S2"))
  # each EGFR record belongs to exactly one sample, not both
  expect_equal(
    egfr$sample[egfr$consequence == "missense_variant"], "TUMOR_S2"
  )
  expect_equal(
    egfr$sample[egfr$consequence == "synonymous_variant"], "TUMOR_S1"
  )
})

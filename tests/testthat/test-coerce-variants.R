test_that("coerce_variants infers consequence when consequence column is absent", {
  df <- data.frame(
    chrom2 = c("chr1", "chr1"),
    pos2   = c(1000L, 2000L),
    ref2   = c("A", "AT"),
    alt2   = c("T", "A"),
    stringsAsFactors = FALSE
  )
  vf <- coerce_variants(df, chrom = "chrom2", pos = "pos2",
                         ref = "ref2", alt = "alt2")
  expect_equal(vf$consequence, c("SNV", "deletion"))
})

test_that("coerce_variants defaults gene to NA when gene column is absent", {
  df <- data.frame(chrom2 = "chr1", pos2 = 100L, ref2 = "A", alt2 = "T")
  vf <- coerce_variants(df, chrom = "chrom2", pos = "pos2",
                         ref = "ref2", alt = "alt2")
  expect_true(is.na(vf$gene))
})

test_that("coerce_variants defaults sample to NA when sample column is absent", {
  df <- data.frame(chrom2 = "chr1", pos2 = 100L, ref2 = "A", alt2 = "T")
  vf <- coerce_variants(df, chrom = "chrom2", pos = "pos2",
                         ref = "ref2", alt = "alt2")
  expect_true(is.na(vf$sample))
})

test_that("coerce_variants uses an explicit consequence column when present", {
  df <- data.frame(
    chrom2 = "chr1", pos2 = 100L, ref2 = "A", alt2 = "T",
    class2 = "Missense_Mutation", stringsAsFactors = FALSE
  )
  vf <- coerce_variants(df, chrom = "chrom2", pos = "pos2", ref = "ref2",
                         alt = "alt2", consequence = "class2")
  expect_equal(vf$consequence, "Missense_Mutation")
})

test_that("coerce_variants carries over columns not mapped to a known role", {
  df <- data.frame(
    chrom2 = "chr1", pos2 = 100L, ref2 = "A", alt2 = "T",
    depth  = 45L, vaf = 0.3, stringsAsFactors = FALSE
  )
  vf <- coerce_variants(df, chrom = "chrom2", pos = "pos2",
                         ref = "ref2", alt = "alt2")
  expect_true(all(c("depth", "vaf") %in% colnames(vf)))
  expect_equal(vf$depth, 45L)
  expect_equal(vf$vaf, 0.3)
})

test_that("coerce_variants does not carry over columns already mapped to a role", {
  df <- data.frame(
    chrom2 = "chr1", pos2 = 100L, ref2 = "A", alt2 = "T",
    gene2  = "TP53", stringsAsFactors = FALSE
  )
  vf <- coerce_variants(df, chrom = "chrom2", pos = "pos2", ref = "ref2",
                         alt = "alt2", gene = "gene2")
  expect_equal(sum(colnames(vf) == "gene"), 1L)
  expect_false("gene2" %in% colnames(vf))
})

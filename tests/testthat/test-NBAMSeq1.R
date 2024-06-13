test_that("NBAMSeq1 returns same object as NBAMSeq works", {
  set.seed(1234)
  sim <- makeExample(n = 4, m = 10)
  gsd <- suppressWarnings(NBAMSeq(sim))
  gsd_r <- results(gsd, name = "pheno")
  set.seed(1234)
  gsd1 <- suppressWarnings(NBAMSeq1(sim, AIC = TRUE))
  gsd1_r <- results1(gsd1, name = "pheno")

  expect_identical(gsd_r, gsd1_r)
})

test_that("bs = 're' works", {
  set.seed(1234)
  sim1 <- makeExample(n = 4, m = 10)
  sim1$id <- factor(seq_along(sim1$pheno))
  sim1@design <- ~ s(pheno) + s(id, bs = "re")
  expect_error(NBAMSeq(sim1), regexp = "the model matrix is not full rank")
  expect_no_error(suppressWarnings(NBAMSeq1(sim1)))
})

test_that("bam works", {
  set.seed(1234)
  sim1 <- makeExample(n = 4, m = 10)
  sim1$id <- factor(seq_along(sim1$pheno))
  sim1@design <- ~ s(pheno) + s(id, bs = "re")
  expect_no_error(suppressWarnings(NBAMSeq1(sim1, use_bam = TRUE)))
})

test_that("parallel", {
  skip(message = "Skipping because of duration")
  set.seed(1234)
  sim <- makeExample(n = 4, m = 10)
  param <- BiocParallel::SnowParam(workers = 2)
  expect_no_error(suppressWarnings(NBAMSeq1(sim, BPPARAM = param)))
})

test_that("failed to fit genes are removed", {
  set.seed(1234)
  sim1 <- makeExample(n = 4, m = 10)
  sim1_assay <- assay(sim1)
  sim1_assay[1, ] <- 0 # Expect gene 1 to be removed
  sim1 <- suppressWarnings(
    NBAMSeqDataSet(countData = sim1_assay, colData = colData(sim1), design = ~ s(pheno))
  )

  fit <- suppressWarnings(NBAMSeq1(sim1))
  # One gene removed = 3
  expect_equal(nrow(fit), 3)
})

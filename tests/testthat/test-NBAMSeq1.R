test_that("bam works", {
  sim <- makeExample(n = 3, m = 10)
  set.seed(1234)
  gsd <- suppressWarnings(NBAMSeq(sim))

  set.seed(1234)
  gsd1 <- suppressWarnings(NBAMSeq1(sim))
  
  expect_identical(gsd, gsd1)
})

test_that("bs = 're' works", {
  sim1 <- makeExample(n = 3, m = 10)
  sim1$id <- factor(seq_along(sim$pheno))
  sim1@design <- ~s(pheno) + s(id, bs = "re")
  expect_error(NBAMSeq(sim1), regexp = "the model matrix is not full rank")
  expect_no_error(suppressWarnings(NBAMSeq1(sim1)))
})

test_that("parallel", {
  skip(message = "Skipping because of duration")
  sim <- makeExample(n = 3, m = 10)
  set.seed(1234)
  param <- BiocParallel::SnowParam(workers = 2)
  expect_no_error(suppressWarnings(NBAMSeq1(sim, BPPARAM = param)))  
})

test_that("failed to fit genes are removed", {
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

test_that("bam fitlin works", {
  skip(message = "Skipping because of duration")
  sim <- makeExample(n = 3, m = 10)
  set.seed(1234)
  gsd_lin <- suppressWarnings(NBAMSeq(sim, fitlin = TRUE))
  
  set.seed(1234)
  gsd1_lin <- suppressWarnings(NBAMSeq1(sim, fitlin = TRUE))
  
  expect_identical(gsd_lin, gsd1_lin)
})

# sim1 <- makeExample(n = 4, m = 10)
# sim1$id <- factor(seq_along(sim$pheno))
# sim1$f <- factor(rbinom(n = seq_along(sim$pheno), size = 1, prob= 0.5))
# sim1@design <- ~s(pheno) + s(id, bs = "re") + f + s(pheno, by = f)
# 
# library(microbenchmark)
# NBAMSeq1(sim1)
# microbenchmark(
#   NBAMSeq1(sim1),
#   NBAMSeq1(sim1, bam = TRUE), 
#   times = 3
# )
# 
# debug(NBAMSeq1)
# fit <- NBAMSeq1(sim1, bam = TRUE)
# fit |> rowData()
# fit |> results(name = "pheno")

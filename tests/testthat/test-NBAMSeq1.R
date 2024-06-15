test_that("NBAMSeq1 returns same object as NBAMSeq works", {
  set.seed(1234)
  sim <- makeExample(n = 4, m = 10)
  gsd <- suppressWarnings(NBAMSeq(sim))
  gsd_r <- wrap_results(gsd, name = "pheno")
  set.seed(1234)
  gsd1 <- suppressWarnings(NBAMSeq1(sim, alpha_bound = 10, use_bam = F))
  gsd1_r <- topTable1(gsd1, coef = "s(pheno)")[row.names(gsd_r), c("edf", "P.Value")]

  expect_identical(gsd_r, gsd1_r)
})

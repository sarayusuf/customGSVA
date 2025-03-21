library(testthat)
library(customGSVA)

test_that("Custom GSVA matches original GSVA", {
  expr <- matrix(rnorm(1000), nrow = 100, ncol = 10)
  rownames(expr) <- paste0("Gene", 1:100)
  gset.idx.list <- list(Pathway1 = c("Gene1", "Gene2", "Gene3"))
  
  es_custom_gsva <- gsva(expr, gset.idx.list, method = "gsva")
  param_gsva <- GSVA::gsvaParam(expr, gset.idx.list, kcdf = "Gaussian")
  es_original_gsva <- GSVA::gsva(param_gsva)
  
  expect_equal(es_custom_gsva, es_original_gsva)
})
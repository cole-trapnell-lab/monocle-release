library(monocle)
library(HSMMSingleCell)
context("relative2abs works properly")

data(HSMM_expr_matrix)
data(HSMM_gene_annotation)
data(HSMM_sample_sheet)

pd <- new("AnnotatedDataFrame", data = HSMM_sample_sheet)
fd <- new("AnnotatedDataFrame", data = HSMM_gene_annotation)
HSMM <- newCellDataSet(as.matrix(HSMM_expr_matrix),
                       phenoData = pd,
                       featureData = fd,
                       lowerDetectionLimit = 0.1,
                       expressionFamily = tobit(Lower = 0.1))

#write test code for this: 
test_that("relative2abs works with valid input and return_all = FALSE", {
  expect_error(relative2abs(HSMM, method = "num_genes", verbose = TRUE), NA)
})

test_that("relative2abs works with valid input and return_all = TRUE", {
  expect_error(relative2abs(HSMM, method = "num_genes", verbose = TRUE, return_all = TRUE), NA)
})

test_that("throws error if infinite parameter is present", {
  expect_error(relative2abs(HSMM, method = "num_genes", volume = Inf), "Your input parameters should not contain either null or infinite numbers")
})

test_that("throws error if concentration detection limit is too low", {
  expect_error(relative2abs(HSMM, method = "num_genes", detection_threshold = 0), "concentration detection limit should be between 0.01430512 and 7500")
})

test_that("throws error if concentration detection limit is too high", {
  expect_error(relative2abs(HSMM, method = "num_genes", detection_threshold = 7501), "concentration detection limit should be between 0.01430512 and 7500")
})

test_that("throws error if ERCC_controls is null and ERCC_annotation is not", {
  expect_error(relative2abs(HSMM, method = "num_genes", ERCC_controls = NULL, ERCC_annotation = 1), "If you want to transform the data to copy number with your spikein data, please provide both of ERCC_controls and ERCC_annotation data frame...")
})

test_that("throws error if ERCC_annotation is null and ERCC_controls is not", {
  expect_error(relative2abs(HSMM, method = "num_genes", ERCC_controls = 1, ERCC_annotation = NULL), "If you want to transform the data to copy number with your spikein data, please provide both of ERCC_controls and ERCC_annotation data frame...")
})

## TO DO:
## ADD TEST WHERE ERCC_ANNOTATION AND ERCC_CONTROLS ARE NOT NULL


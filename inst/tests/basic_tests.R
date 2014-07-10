load_small_test <- function(){
  sample_sheet_small <- read.delim("../../data/sample_sheet_small.txt", row.names=1)
  fpkm_matrix_small <- read.delim("../../data/fpkm_matrix_small.txt")
  pd <- new("AnnotatedDataFrame", data = sample_sheet_small)
  small_set <- new("CellDataSet", exprs = as.matrix(fpkm_matrix_small), phenoData = pd)
  small_set
}

test_that("CellDataSet objects can be built from scratch", {
  small_set <- load_small_test()
  expect_false(is.null(small_set))
})

test_that("Complete workflow for small data set is OK", {
  small_set <- load_small_test()
  expect_false(is.null(small_set))
  
  small_set <- detectGenes(small_set, 0.1)
  expect_false(is.null(fData(small_set)$num_detected))
  
  small_set <- setOrderingFilter(small_set, row.names(fData(small_set)))
  expect_false(is.null(fData(small_set)$use_for_ordering))
  
})


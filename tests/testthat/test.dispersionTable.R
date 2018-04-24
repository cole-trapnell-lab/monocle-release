library(monocle)
library(HSMMSingleCell)
context("dispersionTable validate input")

#write test code for this: 
test_that("dispersionTable validates input", {
  data(HSMM_expr_matrix)
  data(HSMM_gene_annotation)
  data(HSMM_sample_sheet)
  
  pd <- new("AnnotatedDataFrame", data = HSMM_sample_sheet)
  fd <- new("AnnotatedDataFrame", data = HSMM_gene_annotation)
  HSMM <- newCellDataSet(as.matrix(HSMM_expr_matrix),
                         phenoData = pd, 
                         featureData = fd,
                         lowerDetectionLimit=1,
                         expressionFamily=tobit())  
  expect_error(dispersionTable(HSMM))

})
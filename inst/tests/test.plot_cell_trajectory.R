library(monocle)
library(HSMMSingleCell)
context("plot_cell_trajectory")

test_that("plot_cell_trajectory() is able to scale markers linearly and logarithmically"){
  data(HSMM_expr_matrix)
  data(HSMM_gene_annotation)
  data(HSMM_sample_sheet)
  pd <- new("AnnotatedDataFrame", data = HSMM_sample_sheet)
  fd <- new("AnnotatedDataFrame", data = HSMM_gene_annotation)
  HSMM <- newCellDataSet(as.matrix(HSMM_expr_matrix),   
                         phenoData = pd, 
                         featureData = fd)
  
  rpc_matrix <- relative2abs(HSMM)
  HSMM <- newCellDataSet(as(as.matrix(rpc_matrix), "sparseMatrix"),
                         phenoData = pd, 
                         featureData = fd,
                         lowerDetectionLimit=1,
                         expressionFamily=negbinomial.size())
  
  HSMM <- estimateSizeFactors(HSMM)
  HSMM <- estimateDispersions(HSMM)
  HSMM <- clusterCells(HSMM, num_clusters=2)
  
  expect_error(plot_cell_trajectory(HSMM, markers = "MYH3"), NA)
  expect_error(plot_cell_trajectory(HSMM, markers = "MYH3", markers_linear = TRUE), NA)
}
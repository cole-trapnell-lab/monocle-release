library(monocle)
library(HSMMSingleCell)
context("plotCellTrajectory")

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

test_that("plot_cell_trajectory() throws an error if reduceDimensionality hasn't been called yet", {
  expect_error(plot_cell_trajectory(HSMM))
})

test_that("plot_cell_trajectory() throws an error if the CellDataSet's reduction method is unrecognizable", {
  HSMM <- estimateSizeFactors(HSMM)
  HSMM <- estimateDispersions(HSMM)
  HSMM <- clusterCells(HSMM, num_clusters=2)
  HSMM@dim_reduce_type <- "asdasd"
  expect_error(plot_cell_trajectory(HSMM), "Error: unrecognized dimensionality reduction method.")
})

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

#test_that("plot_cell_trajectory() is able to scale markers linearly and logarithmically", {
#  HSMM <- estimateSizeFactors(HSMM)
#  HSMM <- estimateDispersions(HSMM)
#  HSMM <- clusterCells(HSMM, num_clusters=2)
#  expect_error(plot_cell_trajectory(HSMM, markers = "MYH3"), NA)
 # expect_error(plot_cell_trajectory(HSMM, markers = "MYH3", markers_linear = TRUE), NA)
#})

#test_that("plot_cell_trajectory() is able to run without the presence of markers", {
#  HSMM <- estimateSizeFactors(HSMM)
#  HSMM <- estimateDispersions(HSMM)
#  HSMM <- clusterCells(HSMM, num_clusters=2)
#  expect_error(plot_cell_trajectory(HSMM), NA)
#})

#test_that("plot_cell_trajectory() is able to show cell names", {
#  HSMM <- estimateSizeFactors(HSMM)
#  HSMM <- estimateDispersions(HSMM)
#  HSMM <- clusterCells(HSMM, num_clusters=2)
#  expect_error(plot_cell_trajectory(HSMM, show_cell_names = TRUE), NA)
#})

#test_that("plot_cell_trajectory() is able to show cell state numbers", {
#  HSMM <- estimateSizeFactors(HSMM)
#  HSMM <- estimateDispersions(HSMM)
#  HSMM <- clusterCells(HSMM, num_clusters=2)
#  expect_error(plot_cell_trajectory(HSMM, show_state_number = TRUE), NA)
#})












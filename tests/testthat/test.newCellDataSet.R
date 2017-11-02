library(monocle)
library(HSMMSingleCell)
data(HSMM_expr_matrix)
data(HSMM_gene_annotation)
data(HSMM_sample_sheet)
pd <- new("AnnotatedDataFrame", data = HSMM_sample_sheet)
fd <- new("AnnotatedDataFrame", data = HSMM_gene_annotation)

context("newCellDataSet")

  test_that("newCellDataSet is able to create a cell data set with phenoData and featureData inputted from user", { 
      expect_error(newCellDataSet(as.matrix(HSMM_expr_matrix), phenoData = pd, featureData = fd), NA)
  })
  
  test_that("newCellDataSet is able to create a cell data set with only phenoData inputted from user", { 
      expect_warning(newCellDataSet(as.matrix(HSMM_expr_matrix), phenoData = pd), "Warning: featureData must contain a column verbatim named 'gene_short_name' for certain functions")
  })
  
  test_that("newCellDataSet is able to create a cell data set with only featureData inputted from user", { 
    expect_error(newCellDataSet(as.matrix(HSMM_expr_matrix), featureData = fd), NA)
  })
  
  test_that("newCellDataSet throws warning if user inputs a featureData table that doesn't contain a column named 'gene_short_name'", {
    HSMM_gene_annotation$gene_short_name <- NULL
    fd <- new("AnnotatedDataFrame", data = HSMM_gene_annotation)
    expect_warning(newCellDataSet(as.matrix(HSMM_expr_matrix), featureData = fd), "Warning: featureData must contain a column verbatim named 'gene_short_name' for certain functions")
  })
  
  test_that("newCellData throws error if cellData is not a matrix or a sparse matrix", {
    expect_error(newCellDataSet(pd, phenoData = pd, featureData = fd))
  })
  
  
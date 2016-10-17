library(HSMMSingleCell)
library(monocle)
context("detectGenes")

test_that("detectGenes() reports non-NULL count vectors", {
  small_set <- load_HSMM_markers()
  expect_false(is.null(small_set))
  
  small_set <- detectGenes(small_set)
  
  expect_false(is.null(fData(small_set)$num_cells_expressed))
  expect_false(is.null(pData(small_set)$num_genes_expressed))
  
  small_set <- setOrderingFilter(small_set, row.names(fData(small_set)))
  expect_false(is.null(fData(small_set)$use_for_ordering))
})


test_that("detectGenes() reports expected frequencies of the classic muscle markers", {
  small_set <- load_HSMM_markers()
  expect_false(is.null(small_set))
  
  small_set <- detectGenes(small_set)
  
  expect_equal(fData(small_set)$num_cells_expressed, c(211, 117, 131, 238, 148, 109, 126, 241,  82, 105,  51, 122,  7, 117,  82, 271, 211, 117, 181, 72, 271))
  expect_equal(head(pData(small_set)$num_genes_expressed, n=10), c(13, 14, 15,  9, 12, 16, 14,  9,  9, 12))
  expect_equal(sum(pData(small_set)$num_genes_expressed), 3010)
  #expect_false(is.null(pData(small_set)$num_genes_expressed))
})



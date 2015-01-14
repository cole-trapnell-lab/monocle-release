library(monocle)
context("CellDataSet-setup")

test_that("CellDataSet objects can be built from scratch", {
  small_set <- load_HSMM()
  expect_false(is.null(small_set))
})

test_that("CellDataSet objects can be subsetted", {
  small_set <- load_HSMM_markers()
  expect_false(is.null(small_set))
  expect_equivalent(dim(small_set), c(21, 271))
})




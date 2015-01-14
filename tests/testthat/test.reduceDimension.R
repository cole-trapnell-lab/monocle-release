library(monocle)
context("reduceDimension")

test_that("reduceDimension() correctly distills expression to the plane", {
  set.seed(123)
  small_set <- load_HSMM_markers()
  
  small_set <- setOrderingFilter(small_set, row.names(fData(small_set)))
  
  small_set <- suppressMessages(reduceDimension(small_set, use_irlba=FALSE))
  
  expect_false(is.null(reducedDimA(small_set)))
  expect_equivalent(dim(reducedDimA(small_set)), c(21, 2))
  
  expect_false(is.null(reducedDimW(small_set)))
  expect_equivalent(dim(reducedDimW(small_set)), c(2, 2))
  
  expect_false(is.null(reducedDimS(small_set)))
  expect_equivalent(dim(reducedDimS(small_set)), c(2, 271))
  
  expect_false(is.null(reducedDimK(small_set)))
  expect_equivalent(dim(reducedDimK(small_set)), c(21, 2))
  
  })

test_that("reduceDimension() properly validates its input",{
  
})
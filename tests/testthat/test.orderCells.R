library(monocle)
context("orderCells")

test_that("orderCells() correctly orders cells by State and Pseudotime", {
  set.seed(123)
  small_set <- load_HSMM_markers()
  
  small_set <- setOrderingFilter(small_set, row.names(fData(small_set)))
  
  small_set <- suppressMessages(reduceDimension(small_set, use_irlba=FALSE))
  
  small_set <- suppressMessages(orderCells(small_set))
  
  expect_false(is.null(pData(small_set)$Pseudotime))
  expect_false(is.null(pData(small_set)$State))
  
  expect_equivalent(levels(pData(small_set)$State), "1")
  
  expect_less_than(abs(max(pData(small_set)$Pseudotime) -  124.8243), 1e-4)
})

test_that("orderCells() properly validates its input",{
  
  small_set <- load_HSMM_markers()
  
  small_set <- setOrderingFilter(small_set, row.names(fData(small_set)))
  
  results <- evaluate_promise(orderCells(small_set))
  
  expect_equals(results$error, "Error: reducedDimS is NULL. Did you call reducedDimension()?")
  
  small_set <- setOrderingFilter(small_set, row.names(fData(small_set)))
  small_set <- suppressMessages(reduceDimension(small_set, use_irlba=FALSE))
  
  results <- evaluate_promise(orderCells(small_set, num_paths=1000))
  
  expect_equals(results$error, "Error: num_paths must be fewer than the number of cells")
  
  results <- evaluate_promise(orderCells(small_set, num_paths=0))
  
  expect_equals(results$error, "Error: num_paths must be > 0")
  
  results <- evaluate_promise(orderCells(small_set, root_cell ="ARGHYBLARGH"))
  
  expect_equals(results$error, "Error: root_cell must be the name of a cell in the CellDataSet")
})
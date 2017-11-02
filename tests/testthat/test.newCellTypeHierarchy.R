context("newCellTypeHierarchy")

test_that("newCellTypeHierarchy creates a newCellTypeHierarchy", {
  cth <- newCellTypeHierarchy()
  expect_equal(typeof(cth), "S4")
  expect_equal(typeof(cth@classificationTree), "list")
})
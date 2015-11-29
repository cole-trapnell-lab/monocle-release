library(monocle)
context("relative2abs")

test_that("relative2abs() recovers the absolute transcript counts from the relative abundance", {

	HSMM <- load_HSMM()
	rpc_matrix <- relative2abs(HSMM)

	expect_equal(round(as.numeric(rpc_matrix[1, 1:5])), 
             c(7, 0, 2, 0, 1))
}
)
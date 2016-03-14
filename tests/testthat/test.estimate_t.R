library(monocle)
context("estimate_t")

test_that("estimate_t() reports the estimated mode of the relative abundance", {

	HSMM <- load_HSMM()

	#test the t_estimate with non-default values
	t_estimate <- as.numeric(estimate_t(exprs(HSMM), relative_expr_thresh = 0.2))

	expect_equal(length(t_estimate), as.vector(ncol(HSMM)))
	expect_equal(round(as.numeric(t_estimate[1:5])), 
             c(1, 11, 26, 21, 52))
})
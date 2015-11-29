library(monocle)
context("estimate_t")

test_that("estimate_t() reports the estimated mode of the relative abundance", {

	HSMM <- load_HSMM()
	t_estimate <- as.numeric(estimate_t(exprs(HSMM)))

	expect_equal(round(as.numeric(t_estimate[1:5])), 
             c(1, 11, 25, 21, 51))
})
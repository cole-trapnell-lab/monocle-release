library(monocle)
context("estimate_t")

test_that("estimate_t() reports the estimated mode of the relative abundance"), {

	HSMM <- load_HSMM()
	rpc_matrix <- relative2abs(HSMM)
}
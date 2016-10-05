library(monocle)
context("scale_pseudotime")

#write test code for this: 
test_that("scale_pseudotime() scales pseudotime for each lineages up to 100", {
	#test the scaling on the HSMM data: 
	HSMM <- load_HSMM()
    HSMM <- reduceDimension(HSMM, use_vst = F, pseudo_expr = 1)
	HSMM <- orderCells(HSMM, num_paths=2, scale_pseudotime = T)
	 
	expect_equal(range(pData(HSMM)$Pseudotime), c(0, 100))

	#test the scaling on the lung data: 
	lung <- load_lung()
	lung <- orderCells(lung, num_paths=2, scale_pseudotime = T)
	expect_equal(range(pData(lung)$Pseudotime), c(0, 100))

})
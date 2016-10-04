library(monocle)
context("BEAM")

<<<<<<< 7afe2e5a61fe398816056758e84b77f3d44ceb50
test_that("BEAM() reports valid branchTest, bifurcationTmePoint detection results for markers in lung dataset", {
=======
test_that("BEAM() reports valid branchTest, bifurcationTimePoint detection results for markers in lung dataset", {
>>>>>>> First commit

	set.seed(123)

	lung <- load_lung()

	BEAM_res <- BEAM(lung)

	expect_false(is.null(lung))

	expect_equal(colnames(BEAM_res)[1:4], c("status", "family", "pval", "qval"))

	# test the branchTest
	expect_gt(sum(BEAM_res$qval < 0.01), 100)

	# test the bifurcation time point detection
	expect_equal(BEAM_res['ENSMUSG00000000031.9', 'Bifurcation_time_point'], 27)
	
	# test on the non-default parameters: 
	BEAM_res2 <- BEAM(lung, relative_expr = F, stretch = F, weighted = F, 
		pseudocount = 1, lineage_labels = c('AT1', 'AT2'), q_thrsld = 0.01, 
		cores = detectCores(), verbose = T)

	expect_equal(colnames(BEAM_res2)[1:4], c("status", "family", "pval", "qval"))

	# test the branchTest
	expect_gt(sum(BEAM_res2$qval < 0.01), 100)
<<<<<<< 7afe2e5a61fe398816056758e84b77f3d44ceb50

=======
	
>>>>>>> First commit
})



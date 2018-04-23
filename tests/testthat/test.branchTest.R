library(monocle)
context("branchTest")

test_that("branchTest() report valid branch test results for markers in lung dataset", {

	set.seed(123)

	lung <- load_lung()

	#test for the default behavior
	branchTest_res <- branchTest(lung)
	#test output from the branchTest
	expect_equal(colnames(branchTest_res)[1:4], c("status", "family", "pval", "qval"))
	expect_gt(sum(branchTest_res$qval < 0.01), 100)

	#test with different parameters
	branchTest_res2 <- branchTest(lung, fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3)*Lineage",
	                             reducedModelFormulaStr = "~sm.ns(Pseudotime, df = 3)", 
	                             lineage_states =  c(2, 3), 
	                             relative_expr = F,
	                             stretch = F,
	                             pseudocount = 1,
	                             cores = detectCores(), 
	                             weighted = F, 
	                             lineage_labels = c('AT2', 'AT1'))

	expect_false(is.null(lung))

	#test output from the branchTest
	expect_equal(colnames(branchTest_res2)[1:4], c("status", "family", "pval", "qval"))
	expect_gt(sum(branchTest_res2$qval < 0.01), 100)
	
})
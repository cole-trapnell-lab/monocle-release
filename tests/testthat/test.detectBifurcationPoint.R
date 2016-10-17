library(monocle)
context("detectBifurcationPoint")

test_that("detectBifurcationPoint() reports the bifurcation time point detected for the marker genes in lung dataset", {

	set.seed(123)

	lung <- load_lung()

	#make a newCellDataSet object with the smoothed data? 
	ILRs_res <- calILRs(cds = lung, 
				  trajectory_states = c(2, 3), 
				  lineage_labels = NULL, 
				  stretch = T, 
				  cores = 1, 
				  trend_formula = "~sm.ns(Pseudotime, df = 3)*Lineage",
				  ILRs_limit = 3, 
				  relative_expr = T, 
				  weighted = T, 
				  pseudocount = 0, 
				  return_all = T)

	BifurcationTimePoint_res <- detectBifurcationPoint(str_log_df = ILRs_res$str_norm_div_df,
	  lineage_states = c(2, 3), 
	  stretch = T, 
	  cores = 1, 
	  trend_formula = "~sm.ns(Pseudotime, df = 3)*Lineage", 
	  relative_expr = T, 
	  weighted = T, 
	  pseudocount = 0)

	expect_false(is.null(lung))
	expect_equal(as.vector(BifurcationTimePoint_res[1:4]), c(27, -45, 8, -36))

})
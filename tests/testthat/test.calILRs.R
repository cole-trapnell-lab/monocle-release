library(monocle)
context("calILRs")

test_that("calILRs() calculates the Instant Log Ratio between two lineages", {
	lung <- load_lung()

	ILRs_res <- calILRs(cds = lung, 
  			  trajectory_states = c(2, 3), 
  			  lineage_labels = c('AT1', 'AT2'), 
  			  stretch = T, 
  			  cores = detectCores(), 
  			  trend_formula = "~sm.ns(Pseudotime, df = 3)*Lineage",
  			  ILRs_limit = 3, 
  			  relative_expr = T, 
  			  weighted = T, 
  			  pseudocount = 0, 
  			  return_all = T)
	expect_equal(names(ILRs_res), c("str_logfc_df", "norm_str_logfc_df", "str_norm_div_df", "str_raw_div_df", "str_branchA_expression_curve_matrix", "str_branchB_expression_curve_matrix"))
	expect_equal(dim(ILRs_res$str_logfc_df), c(as.vector(nrow(lung)), 100))
	expect_equal(dim(ILRs_res$norm_str_logfc_df), c(as.vector(nrow(lung)), 100))
	expect_equal(dim(ILRs_res$str_norm_div_df), c(as.vector(nrow(lung)), 100))
	expect_equal(dim(ILRs_res$str_raw_div_df), c(as.vector(nrow(lung)), 100))
	expect_equal(dim(ILRs_res$str_branchA_expression_curve_matrix), c(as.vector(nrow(lung)), 100))
	expect_equal(dim(ILRs_res$str_branchB_expression_curve_matrix), c(as.vector(nrow(lung)), 100))
	
	#test on the variable value: not all NA
	
	ILRs_res <- calILRs(cds = lung, 
	                    trajectory_states = c(2, 3), 
	                    lineage_labels = c('AT1', 'AT2'), 
	                    stretch = F, 
	                    cores = detectCores(), 
	                    trend_formula = "~sm.ns(Pseudotime, df = 3)*Lineage",
	                    ILRs_limit = 1, 
	                    relative_expr = F, 
	                    weighted = F, 
	                    pseudocount = 1, 
	                    return_all = F)
	expect_equal(dim(ILRs_res),  c(as.vector(nrow(lung)), 100))
	
	expect_false(is.null(lung))
})
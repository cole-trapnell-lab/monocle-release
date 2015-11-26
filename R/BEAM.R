#' Peform the beam analysis test
#'
#' Perform BEAM analysis
#'
#' @param cds a CellDataSet object upon which to perform this operation
#' @param fullModelFormulaStr a formula string specifying the full model in differential expression tests (i.e. likelihood ratio tests) for each gene/feature.
#' @param reducedModelFormulaStr a formula string specifying the reduced model in differential expression tests (i.e. likelihood ratio tests) for each gene/feature.
#' @param lineage_states ids for the immediate branch lineage which obtained from lineage construction based on MST
#' @param relative_expr a logic flag to determine whether or not the relative gene expression should be used
#' @param stretch  a logic flag to determine whether or not each lineage should be stretched
#' @param pseudocount pseudo count added before fitting the spline curves 
#' @param weighted  A logic flag to determine whether or not we should use the navie logLikelihood weight scheme for the duplicated progenitor cells
#' @param cores the number of cores to be used while testing each gene for differential expression
#' @param lineage_labels the name for each lineage, for example, AT1 or AT2  
#' @return a data frame containing the p values and q-values from the likelihood ratio tests on the parallel arrays of models.
#' @export
#'
BEAM <- function(cds, fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3)*Lineage", 
					reducedModelFormulaStr = "~sm.ns(Pseudotime, df = 3)", 
					lineage_states = c(2, 3), 
					relative_expr = TRUE, 
					stretch = TRUE, 
					pseudocount = 0, 
					weighted = T, 
					q_thrsld = 0.05, 
					lineage_labels = NULL, 
					cores = 1, 
					...) {

	branchTest_res <- branchTest(cds, fullModelFormulaStr = fullModelFormulaStr,
	                       reducedModelFormulaStr = reducedModelFormulaStr, 
	                       lineage_states = lineage_states, 
	                       relative_expr = relative_expr,
	                       stretch = stretch,
	                       pseudocount = pseudocount,
	                       cores = cores, 
	                       weighted = weighted, 
	                       lineage_labels = lineage_labels, ...)
	cmbn_df <- branchTest_res[, 1:4] 

# 	if(calABC) {
# 		ABCs_res <- calABCs(cds, trajectory_type = "Lineage", 
# 						  trajectory_states = lineage_states,
# 						  branchTest = FALSE, 
# 						  relative_expr = relative_expr, 
# 						  trend_formula = fullModelFormulaStr,
# 						  stretch = stretch, 
# 						  pseudocount = pseudocount, 
# 						  cores = cores, 
# 						  weighted = weighted, 
# 						  lineage_labels = lineage_labels,
# 						  ...)
# 
# 		cmbn_df <- cbind(cmbn_df, ABCs_res[, 1])
# 	}

	#make a newCellDataSet object with the smoothed data? 
		ILRs_res <- calILRs(cds = cds, 
					  trajectory_type = "Lineage", 
					  trajectory_states = lineage_states, 
					  lineage_labels = lineage_labels, 
					  stretch = stretch, 
					  cores = cores, 
					  trend_formula = fullModelFormulaStr,
					  ILRs_limit = 3, 
					  relative_expr = relative_expr, 
					  weighted = weighted, 
					  pseudocount = pseudocount, 
					  ...)

		cmbn_df <- cbind(cmbn_df, ILRs_res[, 1])

		BifurcationTimePoint_res <- detectBifurcationPoint(str_log_df = ILRs_res$str_norm_div_df,
		  lineage_states = lineage_states, 
		  stretch = stretch, 
		  cores = cores, 
		  trend_formula = fullModelFormulaStr, 
		  relative_expr = relative_expr, 
		  weighted = weighted, 
		  pseudocount = pseudocount, 
			...)

		cmbn_df <- cbind(cmbn_df, BifurcationTimePoint_res[, 1])

	# if(draw_branched_kinetics) {
	# 	plot_genes_branched_pseudotime(cds, 
	# 		lineage_states = lineage_states, 
	# 		lineage_labels = lineage_labels,
	# 		method = "fitting", 
	# 		stretch = TRUE, 
	# 		min_expr = NULL, 
	# 		cell_size = 0.75,
	# 		nrow = NULL, 
	# 		ncol = 1, 
	# 		panel_order = NULL, 
	# 		color_by = "State",
	# 		cell_color_by = "State",
	# 		trajectory_color_by = "State", 
	# 		trend_formula = fullModelFormulaStr, 
	# 		reducedModelFormulaStr = reducedModelFormulaStr, 
	# 		label_by_short_name = TRUE,
	# 		weighted = TRUE, 
	# 		add_ABC = FALSE, 
	# 		add_pval = FALSE,
	# 		normalize = TRUE,
	# 		bifurcation_time = NULL, 
	# 		#gene_pairs = NULL,
	# 	...)
	# }

	# if(draw_branched_heatmap) {
	# 	plot_genes_branched_heatmap(cds_subset, 
	# 	  num_clusters = 6,
	# 	  ABC_df = NULL, 
	# 	  branchTest_df = NULL, 
	# 	  lineage_labels = lineage_labels, 
	# 	  stretch = T, 
	# 	  scaling = T,
	# 	  norm_method = c("vstExprs", "log"), 
	# 	  use_fitting_curves = T, 
	# 	  dist_method = NULL, 
	# 	  hclust_method = "ward", 
	# 	  heatmap_height = 3, 
	# 	  heatmap_width = 4,
	# 	  ABC_lowest_thrsd = 0, 
	# 	  ABC_highest_thrsd = 2,
	# 	  qval_lowest_thrsd = 1, 
	# 	  qval_highest_thrsd = 5,
	# 	  hmcols = NULL, 
	# 	  Cell_type_color = c('#979797', '#F05662', '#7990C8'), 
	# 	  trend_formula = '~sm.ns(Pseudotime, df=3) * Lineage',
	# 	  pseudo_cnt = 0, 
	# 	  add_annotation_row = NULL,
	# 	  add_annotation_col = NULL,
	# 	  show_rownames = F, 
	# 	  cores = cores,
	# 	  use_gene_short_name = F,
	# 	  file_name = 'branched_heatmap.pdf')
	# }


	fd <- fData(cds)

	#combined dataframe: 
	
	fData(cds) <- cbind(cmbn_df, fd)


	return(cds)
}

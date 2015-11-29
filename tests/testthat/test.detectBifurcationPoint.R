library(monocle)
context("detectBifurcationPoint")

test_that("detectBifurcationPoint() reports the bifurcation time point detected for the marker genes in lung dataset"), {

	set.seed(123)

	baseLoc <- system.file(package="monocle")
	#baseLoc <- './inst'
	extPath <- file.path(baseLoc, "extdata")
	load(file.path(extPath, "lung_phenotype_data.RData"))
	load(file.path(extPath, "lung_exprs_data.RData"))
	load(file.path(extPath, "lung_feature_data.RData"))
	lung_exprs_data <- lung_exprs_data[,row.names(lung_phenotype_data)]

	pd <- new("AnnotatedDataFrame", data = lung_phenotype_data)
	fd <- new("AnnotatedDataFrame", data = lung_feature_data)

	# Now, make a new CellDataSet using the RNA counts
	lung <- newCellDataSet(as.matrix(lung_exprs_data), 
	                       phenoData = pd, 
	                       featureData = fd,
	                       lowerDetectionLimit=1,
	                       expressionFamily=negbinomial())

	lung <- estimateSizeFactors(lung)
	lung <- estimateDispersions(lung)

	pData(lung)$Total_mRNAs <- colSums(exprs(lung))
	lung <- detectGenes(lung, min_expr = 1)
	expressed_genes <- row.names(subset(fData(lung), num_cells_expressed >= 5))
	ordering_genes <- expressed_genes
	lung <- setOrderingFilter(lung, ordering_genes)
	lung <- reduceDimension(lung, use_vst = F, pseudo_expr = 1)
	lung <- orderCells(lung, num_paths=2)

	#make a newCellDataSet object with the smoothed data? 
	ILRs_res <- calILRs(cds = lung, 
				  trajectory_type = "Lineage", 
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

}
library(monocle)
context("BEAM")

test_that("BEAM() reports valid branch test results for markers in lung dataset", {

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
	plot_spanning_tree(lung, color_by="Time")

	BEAM_res <- BEAM(lung)

	expect_false(is.null(lung))

	expect_equal(colnames(BEAM_res)[1:4], c("status", "family", "pval", "qval"))

	# test the branchTest
	expect_equal(sum(BEAM_res$qval < 0.01), 131)

	# test the bifurcation time point detection
	expect_equal(BEAM_res['ENSMUSG00000000031.9', 'Bifurcation_time_point'], 27)
	
	BEAM_res <- BEAM(lung, relative_expr = F, stretch = F, weighted = F, 
		pseudocount = 1, lineage_labels = c('AT1', 'AT2'), q_thrsld = 0.01, 
		cores = detectCores(), verbose = T)


})



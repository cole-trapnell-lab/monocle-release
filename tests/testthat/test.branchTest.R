library(monocle)
context("branchTest")

test_that("branchTest() report valid branch test results for markers in lung dataset", {

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

	branchTest_res <- branchTest(lung)

	expect_false(is.null(lung))
	expect_equal(colnames(branchTest_res)[1:4], c("status", "family", "pval", "qval"))

	expect_equal(sum(branchTest_res$qval < 0.01), 131)
})
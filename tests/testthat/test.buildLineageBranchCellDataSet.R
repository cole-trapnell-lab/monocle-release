library(monocle)
context("buildLineageBranchCellDataSet")

test_that("buildLineageBranchCellDataSet() creates a new CellDataSet object with duplicated progenitor cells from lung dataset", {
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

    lung_new <- buildLineageBranchCellDataSet(cds = lung, lineage_states = c(2, 3),
                                            lineage_labels = NULL, method = 'fitting', stretch = T,
                                            weighted = T)
	
	# test the dimension of new built cds object
	expect_equal(as.vector(dim(lung_new)), c(218, 248))

	# test the stretched pseudotime 
	expect_equal(range(pData(lung_new)$Pseudotime), c(0, 100))

	# test the weighted of duplicated cells 
	expect_equal(unique(pData(lung_new[, ])$weight), c(0.5, 1))
	
}
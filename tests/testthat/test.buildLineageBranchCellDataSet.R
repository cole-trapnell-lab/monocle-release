library(monocle)
context("buildLineageBranchCellDataSet")

test_that("buildLineageBranchCellDataSet() creates a new CellDataSet object with duplicated progenitor cells from lung dataset", {
	set.seed(123)

	lung <- load_lung()

    lung_new <- buildLineageBranchCellDataSet(cds = lung, lineage_states = c(2, 3),
                                            lineage_labels = c('AT1', 'AT2'), method = 'fitting', stretch = T,
                                            weighted = T)
	# test the dimension of new built cds object
	expect_equal(as.vector(dim(lung_new)), c(218, 185 + sum(pData(lung)$State == 1)))

	# test the stretched pseudotime 
	expect_equal(range(pData(lung_new)$Pseudotime), c(0, 100))

	# test the weighted of duplicated cells 
	expect_equal(unique(sort(pData(lung_new[, ])$weight)), c(0.5, 1))
	
	#use un-default parameters
	lung_new2 <- buildLineageBranchCellDataSet(cds = lung, lineage_states = c(2, 3),
                                        lineage_labels = c('AT1', 'AT2'), method = 'fitting', stretch = F,
                                        weighted = F)
	# test the dimension of new built cds object
	expect_equal(as.vector(dim(lung_new2)), c(218, 185 + sum(pData(lung)$State == 1)))

	# test the stretched pseudotime 
	expect_less_than(range(pData(lung_new2)$Pseudotime)[2], 100)

	# test the weighted of duplicated cells 
	expect_equal(unique(sort(pData(lung_new2[, ])$weight)), 1)
	

})
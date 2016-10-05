library(monocle)
context("calABCs")

test_that("calABCs() calculates the Area Between Curves for two lineages", {
	set.seed(123)

  lung <- load_lung()
  ABCs_res <- calABCs(lung[1:8, ],
                 trajectory_states = c(2, 3),
                 branchTest = FALSE, 
                 relative_expr = T, 
                 trend_formula = "~sm.ns(Pseudotime, df = 3)*Lineage",
                 stretch = F, 
                 pseudocount = 0, 
                 cores = 1, 
                 weighted = F, 
                 verbose = F,
                 min_expr = 0.5, 
                 integer_expression = FALSE, 
                 num = 500, 
                 lineage_labels = NULL,
                 ABC_method = c('integral'))

    expect_false(is.null(lung))
    expect_equal(colnames(ABCs_res)[1], "ABCs")       
})
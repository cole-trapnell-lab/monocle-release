library(monocle)
context("differentialExpression")

test_that("differentialGeneTest() reports valid test results for markers", {
  skip("Skipping this for now to save time while running tests")
  
  set.seed(123)
  
  small_set <- load_HSMM_markers()
  expect_false(is.null(small_set))
  
  diff_test_res <- differentialGeneTest(small_set, 
                                        fullModelFormulaStr="expression~Media")
  expect_false(is.null(small_set))
  expect_equal(nrow(diff_test_res), length(get_classic_muscle_markers()))
  expect_equal(colnames(diff_test_res), c("status", "pval", "qval"))
  
  # Attach the HUGO symbols and other featureData for these genes
  diff_test_res <- merge(fData(small_set), diff_test_res, by="row.names")
  
  sig_genes <- subset(diff_test_res, qval < 0.01)
  expect_equal(sig_genes$gene_short_name,  c("TNNT1",
                                             "MYH3",
                                             "CCND1",
                                             "MYF5",
                                             "TNNC1",
                                             "TNNT2",
                                             "MYOG",
                                             "CDK2",
                                             "MYH2",
                                             "ID1",
                                             "CCNB1",
                                             "PDGFRA",
                                             "TPM1",
                                             "CCNB2",
                                             "ANPEP",
                                             "CDK1",
                                             "TPM2"))
})

test_that("differentialGeneTest() properly validates its input",{
  
})
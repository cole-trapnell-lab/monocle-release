library(monocle)
library(HSMMSingleCell)
context("detectGenes functions properly")

data(HSMM_expr_matrix)
data(HSMM_gene_annotation)
data(HSMM_sample_sheet)

pd <- new("AnnotatedDataFrame", data = HSMM_sample_sheet)
fd <- new("AnnotatedDataFrame", data = HSMM_gene_annotation)
HSMM <- newCellDataSet(as.matrix(HSMM_expr_matrix),
                       phenoData = pd,
                       featureData = fd,
                       lowerDetectionLimit = 0.1,
                       expressionFamily = tobit(Lower = 0.1))

HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)
HSMM <- detectGenes(HSMM, min_expr = 0.1)
HSMM <- HSMM[,pData(HSMM)$Total_mRNAs < 1e6]
upper_bound <- 10^(mean(log10(pData(HSMM)$Total_mRNAs)) + 2*sd(log10(pData(HSMM)$Total_mRNAs)))
lower_bound <- 10^(mean(log10(pData(HSMM)$Total_mRNAs)) - 2*sd(log10(pData(HSMM)$Total_mRNAs)))
HSMM <- HSMM[,pData(HSMM)$Total_mRNAs > lower_bound & pData(HSMM)$Total_mRNAs < upper_bound]
cth <- newCellTypeHierarchy()

MYF5_id <- row.names(subset(fData(HSMM), gene_short_name == "MYF5"))

#write test code for this: 

test_that("test addCellType works properly", {
  expect_error(cth <- addCellType(cth, "Myoblast", classify_func = function(x) {x[MYF5_id,] >= 1}), NA)
})

test_that("test addCellType throws error when same cell name is used", {
  cth <- addCellType(cth, "Myoblast1", classify_func = function(x) {x[MYF5_id,] >= 1})
  expect_error(cth <- addCellType(cth, "Myoblast1", classify_func = function(x) {x[MYF5_id,] >= 1}), NA)
})



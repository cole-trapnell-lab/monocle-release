library(monocle)
library(HSMMSingleCell)
context("classifyCells functions properly")

data(HSMM_expr_matrix)
data(HSMM_gene_annotation)
data(HSMM_sample_sheet)

pd <- new("AnnotatedDataFrame", data = HSMM_sample_sheet)
fd <- new("AnnotatedDataFrame", data = HSMM_gene_annotation)
HSMM <- newCellDataSet(as.matrix(HSMM_expr_matrix),   
                       phenoData = pd, 
                       featureData = fd,
                       lowerDetectionLimit=0.1,
                       expressionFamily=tobit(Lower=0.1))


rpc_matrix <- relative2abs(HSMM, method = "num_genes")


HSMM <- newCellDataSet(as(as.matrix(rpc_matrix), "sparseMatrix"),
                       phenoData = pd, 
                       featureData = fd,
                       lowerDetectionLimit=0.5,
                       expressionFamily=negbinomial.size())

HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)
HSMM <- detectGenes(HSMM, min_expr = 0.1)
HSMM <- HSMM[,pData(HSMM)$Total_mRNAs < 1e6]
upper_bound <- 10^(mean(log10(pData(HSMM)$Total_mRNAs)) + 2*sd(log10(pData(HSMM)$Total_mRNAs)))
lower_bound <- 10^(mean(log10(pData(HSMM)$Total_mRNAs)) - 2*sd(log10(pData(HSMM)$Total_mRNAs)))
HSMM <- HSMM[,pData(HSMM)$Total_mRNAs > lower_bound & pData(HSMM)$Total_mRNAs < upper_bound]
cth <- newCellTypeHierarchy()

MYF5_id <- row.names(subset(fData(HSMM), gene_short_name == "MYF5"))
ANPEP_id <- row.names(subset(fData(HSMM), gene_short_name == "ANPEP"))

cth <- newCellTypeHierarchy()
cth <- addCellType(cth, "Myoblast", classify_func=function(x) {x[MYF5_id,] >= 1})
cth <- addCellType(cth, "Fibroblast", classify_func=function(x)
{x[MYF5_id,] < 1 & x[ANPEP_id,] > 1})

#write test code for this: 

test_that("test classifyCells throws error if frequency_thresh & enrichment_thresh are both NULL and character string is passed to method", {
  expect_error(expect_error(classifyCells(HSMM, cth, frequency_thresh = NULL, enrichment_thresh = NULL, test = c("t", "e", "s", "t")), "Error: to use classifyCells in grouped mode, you must also set frequency_thresh"))
})

test_that("test classifyCells works properly 1", {
  expect_error(classifyCells(HSMM, cth, 0.1), NA)
})

test_that("test classifyCells works when enrichment_thresh is passed and frequency_thresh is null", {
  expect_error(classifyCells(HSMM, cth, enrichment_thresh = 0.1), NA)
})

test_that("test classifyCells works when frequency_thresh is passed and enrichment_thresh is null", {
  expect_error(classifyCells(HSMM, cth, frequency_thresh = 0.1))
})





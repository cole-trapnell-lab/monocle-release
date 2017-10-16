library(monocle)
library(HSMMSingleCell)
context("addCellType functions properly")

data(HSMM_expr_matrix)
data(HSMM_gene_annotation)
data(HSMM_sample_sheet)

pd <- new("AnnotatedDataFrame", data = HSMM_sample_sheet)
fd <- new("AnnotatedDataFrame", data = HSMM_gene_annotation)

# First create a CellDataSet from the relative expression levels
HSMM <- newCellDataSet(as.matrix(HSMM_expr_matrix),   
                       phenoData = pd, 
                       featureData = fd,
                       lowerDetectionLimit=0.1,
                       expressionFamily=tobit(Lower=0.1))

# Next, use it to estimate RNA counts
rpc_matrix <- relative2abs(HSMM, method = "num_genes")

# Now, make a new CellDataSet using the RNA counts
HSMM <- newCellDataSet(as(as.matrix(rpc_matrix), "sparseMatrix"),
                       phenoData = pd, 
                       featureData = fd,
                       lowerDetectionLimit=0.5,
                       expressionFamily=negbinomial.size())

HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)
HSMM <- detectGenes(HSMM, min_expr = 0.1)
HSMM <- HSMM[,pData(HSMM)$Total_mRNAs < 1e6]
cth <- newCellTypeHierarchy()

MYF5_id <- row.names(subset(fData(HSMM), gene_short_name == "MYF5"))

#write test code for this: 

test_that("test addCellType works properly", {
  expect_error(cth <- addCellType(cth, "Myoblast", classify_func = function(x) {x[MYF5_id,] >= 1}), NA)
})

test_that("test addCellType throws error when same cell name is used", {
  cth <- addCellType(cth, "Myoblast1", classify_func = function(x) {x[MYF5_id,] >= 1})
  expect_error(cth <- addCellType(cth, "Myoblast1", classify_func = function(x) {x[MYF5_id,] >= 1}), "cell type Myoblast1 already exists.")
})



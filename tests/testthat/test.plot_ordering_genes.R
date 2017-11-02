library(monocle)
library(HSMMSingleCell)
context("plot_ordering_genes functions properly")

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

cth <- newCellTypeHierarchy()

MYF5_id <- row.names(subset(fData(HSMM), gene_short_name == "MYF5"))
ANPEP_id <- row.names(subset(fData(HSMM), gene_short_name == "ANPEP"))

cth <- newCellTypeHierarchy()
cth <- addCellType(cth, "Myoblast", classify_func=function(x) {x[MYF5_id,] >= 1})
cth <- addCellType(cth, "Fibroblast", classify_func=function(x)
{x[MYF5_id,] < 1 & x[ANPEP_id,] > 1})

HSMM <- classifyCells(HSMM, cth, 0.1)

disp_table <- dispersionTable(HSMM)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
HSMM <- setOrderingFilter(HSMM, unsup_clustering_genes$gene_id)

test_that("plot_ordering_genes functions properly under normal conditions", 
          expect_error(plot_ordering_genes(HSMM), NA))

test89 <- class(plot_ordering_genes(HSMM))[2]

test_that("plot_ordering_genes returns ggplot", expect_equal(test89, "ggplot"))

test_that("plot_ordering_genes fails if CellDataSet is not passed through", 
          expect_error(plot_ordering_genes(7), "Error input object is not of type 'CellDataSet'"))

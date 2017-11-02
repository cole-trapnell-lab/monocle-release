library(monocle)
library(HSMMSingleCell)
context("reduceDimension is functioning properly")

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
print(head(fData(HSMM)))
expressed_genes <- row.names(subset(fData(HSMM), num_cells_expressed >= 10))

pData(HSMM)$Total_mRNAs <- Matrix::colSums(exprs(HSMM))


HSMM <- HSMM[,pData(HSMM)$Total_mRNAs < 1e6]

HSMM <- detectGenes(HSMM, min_expr = 0.1)

L <- log(exprs(HSMM[expressed_genes,]))

melted_dens_df <- melt(Matrix::t(scale(Matrix::t(L))))

MYF5_id <- row.names(subset(fData(HSMM), gene_short_name == "MYF5"))
ANPEP_id <- row.names(subset(fData(HSMM), gene_short_name == "ANPEP"))

# HSMM <- classifyCells(HSMM, cth, 0.1)

disp_table <- dispersionTable(HSMM)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
HSMM <- setOrderingFilter(HSMM, unsup_clustering_genes$gene_id)

test_that("reduceDimension works properly", expect_error(HSMM <- reduceDimension(HSMM, max_components=2, num_dim = 6, 
                                                                                 reduction_method = 'tSNE', verbose = T), NA))
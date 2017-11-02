library(monocle)
library(HSMMSingleCell)
context("plot_cell_clusters is functioning properly")

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
expressed_genes <- row.names(subset(fData(HSMM), num_cells_expressed >= 10))
pData(HSMM)$Total_mRNAs <- Matrix::colSums(exprs(HSMM))

HSMM <- detectGenes(HSMM, min_expr = 0.1)

L <- log(exprs(HSMM[expressed_genes,]))

melted_dens_df <- melt(Matrix::t(scale(Matrix::t(L))))

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

plot_pc_variance_explained(HSMM, return_all = F) # norm_method = 'log',

HSMM <- reduceDimension(HSMM, max_components=2, num_dim = 6, 
                        reduction_method = 'tSNE', verbose = T) 
HSMM <- clusterCells(HSMM,
                     num_clusters=2)

test_that("plot_cell_clusters functions in vignette", expect_error(plot_cell_clusters(HSMM, 1, 2, color = "CellType"), NA))
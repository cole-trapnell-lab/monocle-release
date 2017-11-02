library(monocle)
library(HSMMSingleCell)
context("differentialGeneTest functions properly")

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


HSMM <- HSMM[,pData(HSMM)$Total_mRNAs < 1e6]
upper_bound <- 10^(mean(log10(pData(HSMM)$Total_mRNAs)) + 2*sd(log10(pData(HSMM)$Total_mRNAs)))
lower_bound <- 10^(mean(log10(pData(HSMM)$Total_mRNAs)) - 2*sd(log10(pData(HSMM)$Total_mRNAs)))
qplot(Total_mRNAs, data=pData(HSMM), color=Hours, geom="density") + 
  geom_vline(xintercept=lower_bound) + 
  geom_vline(xintercept=upper_bound)

HSMM <- HSMM[,pData(HSMM)$Total_mRNAs > lower_bound & 
               pData(HSMM)$Total_mRNAs < upper_bound]								  
HSMM <- detectGenes(HSMM, min_expr = 0.1)

L <- log(exprs(HSMM[expressed_genes,]))

# Standardize each gene, so that they are all on the same scale,
# Then melt the data with plyr so we can plot it easily"
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

HSMM <- reduceDimension(HSMM, max_components=2, num_dim = 6, 
                        reduction_method = 'tSNE', verbose = T) 
HSMM <- clusterCells(HSMM,
                     num_clusters=2)

marker_diff <- markerDiffTable(HSMM[expressed_genes,], 
                               cth, 
                               residualModelFormulaStr="~Media + num_genes_expressed",
                               cores=detectCores())

set.seed(0)

candidate_clustering_genes <- row.names(subset(marker_diff, qval < 0.01))
marker_spec <- calculateMarkerSpecificity(HSMM[candidate_clustering_genes,], cth)
head(selectTopMarkers(marker_spec, 3))

semisup_clustering_genes <- unique(selectTopMarkers(marker_spec, 500)$gene_id)
HSMM <- setOrderingFilter(HSMM, semisup_clustering_genes)

HSMM <- reduceDimension(HSMM, max_components=2, num_dim = 3, norm_method = 'log', reduction_method = 'tSNE', 
                        residualModelFormulaStr="~Media + num_genes_expressed", verbose = T) 
HSMM <- clusterCells(HSMM, num_clusters=2) 

HSMM_myo <- HSMM[,pData(HSMM)$CellType == "Myoblast"]	
HSMM_myo <- estimateDispersions(HSMM_myo)

disp_table <- dispersionTable(HSMM_myo)
ordering_genes <- subset(disp_table, 
                         mean_expression >= 0.5 & 
                           dispersion_empirical >= 1 * dispersion_fit)$gene_id

test_that("differentialGeneTest functions under standard conditions (conditions set in vignette)", 
          expect_error(differentialGeneTest(HSMM_myo[expressed_genes,], fullModelFormulaStr = "~Media"), NA))

test_that("Error is thrown if cds is not of type 'CellDataSet'", 
          expect_error(differentialGeneTest(8, fullModelFormulaStr = "~Media"), NA))

holdPseudotime <- pData(HSMM_myo[expressed_genes,])$Pseudotime
pData(HSMM_myo)$Pseudotime <- Inf

test_that("Error is thrown if Inf is located in pData", 
          expect_error(differentialGeneTest(HSMM_myo, fullModelFormulaStr = "~Media"), "Error: Inf, NaN, or NA values were located in pData of cds in columns mentioned in model terms"))

pData(HSMM_myo)$Pseudotime <- NaN

test_that("Error is thrown if NaN is located in pData", expect_error(differentialGeneTest(HSMM_myo, fullModelFormulaStr = "~Media"), "Error: Inf, NaN, or NA values were located in pData of cds in columns mentioned in model terms"))

pData(HSMM_myo)$Pseudotime <- NA

test_that("Error is thrown if NA is located in pData", expect_error(differentialGeneTest(HSMM_myo, fullModelFormulaStr = "~Media"), "Error: Inf, NaN, or NA values were located in pData of cds in columns mentioned in model terms"))

pData(HSMM_myo)$Pseudotime <- holdPseudotime
pData(HSMM_myo)$Size_Factor <- NULL

test_that("Error is thrown if i) relative_expr is set to true, ii) expression family is 'negbinomial' or 'negbinomial.size', iii) size factors of cds is null", 
  expect_error(differentialGeneTest(HSMM_myo, fullModelFormulaStr = "~Media"), "Error: to call this function with relative_expr==TRUE, you must first call estimateSizeFactors() on the CellDataSet."))

pData(HSMM_myo)$Size_Factor <- NA
test_that("Error is thrown if i) relative_expr is set to true, ii) expression family is 'negbinomial' or 'negbinomial.size', iii) sum of size factors of cds is NA",
  expect_error(differentialGeneTest(HSMM_myo, fullModelFormulaStr = "~Media"), "Error: to call this function with relative_expr==TRUE, you must first call estimateSizeFactors() on the CellDataSet."))

test_that("differentialGeneTest functions properly when cores is larger than 1", expect_error(differentialGeneTest(HSMM_myo[expressed_genes,], fullModelFormulaStr = "~Media", cores = 2), NA))

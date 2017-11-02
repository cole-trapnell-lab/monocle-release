library(monocle)
library(HSMMSingleCell)
context("markerDiffTable is functioning properly")

pd <- new("AnnotatedDataFrame", data = HSMM_sample_sheet)
fd <- new("AnnotatedDataFrame", data = HSMM_gene_annotation)
HSMM <- newCellDataSet(as.matrix(HSMM_expr_matrix), phenoData = pd, featureData = fd)

HSMM <- newCellDataSet(as(umi_matrix, "sparseMatrix"),
                       phenoData = pd,
                       featureData = fd,
                       lowerDetectionLimit = 0.5,
                       expressionFamily = negbinomial.size())

cellranger_pipestance_path <- "/path/to/your/pipeline/output/directory"
gbm <- load_cellranger_matrix(cellranger_pipestance_path)
gbm_cds <- newCellDataSet(exprs(gbm),
                          phenoData = new("AnnotatedDataFrame", data = pData(gbm)),
                          phenoData = new("AnnotatedDataFrame", data = fData(gbm)),
                          lowerDetectionLimit = 0.5,
                          expressionFamily = negbinomial.size())

HSMM <- detectGenes(HSMM, min_expr = 0.1)
expressed_genes <- row.names(subset(fData(HSMM), num_cells_expressed >= 10))
valid_cells <- row.names(subset(pData(HSMM),
                                Cells.in.Well == 1 &
                                  Control == FALSE &
                                  Clump == FALSE &
                                  Debris == FALSE &
                                  Mapped.Fragments > 1000000))
HSMM <- HSMM[,valid_cells]

pData(HSMM)$Total_mRNAs <- Matrix::colSums(exprs(HSMM))

HSMM <- HSMM[,pData(HSMM)$Total_mRNAs < 1e6]

upper_bound <- 10^(mean(log10(pData(HSMM)$Total_mRNAs)) +
                     2*sd(log10(pData(HSMM)$Total_mRNAs)))
lower_bound <- 10^(mean(log10(pData(HSMM)$Total_mRNAs)) -
                     2*sd(log10(pData(HSMM)$Total_mRNAs)))

HSMM <- HSMM[,pData(HSMM)$Total_mRNAs > lower_bound &
               pData(HSMM)$Total_mRNAs < upper_bound]
HSMM <- detectGenes(HSMM, min_expr = 0.1)

# Log-transform each value in the expression matrix.
L <- log(exprs(HSMM[expressed_genes,]))

# Standardize each gene, so that they are all on the same scale,
# Then melt the data with plyr so we can plot it easily
melted_dens_df <- melt(Matrix::t(scale(Matrix::t(L))))

# Plot the distribution of the standardized gene expression values.
qplot(value, geom = "density", data = melted_dens_df) +
  stat_function(fun = dnorm, size = 0.5, color = 'red') +
  xlab("Standardized log(FPKM)") +
  ylab("Density")

MYF5_id <- row.names(subset(fData(HSMM), gene_short_name == "MYF5"))
ANPEP_id <- row.names(subset(fData(HSMM), gene_short_name == "ANPEP"))

cth <- newCellTypeHierarchy()
cth <- addCellType(cth, "Myoblast", classify_func = function(x) { x[MYF5_id,] >= 1 })
cth <- addCellType(cth, "Fibroblast", classify_func = function(x)
{ x[MYF5_id,] < 1 & x[ANPEP_id,] > 1 })

HSMM <- classifyCells(HSMM, cth, 0.1)

test_that("markerDiffTable works properly in vignette", 
          expect_error(markerDiffTable(HSMM[expressed_genes,],
                                       cth,
                                       residualModelFormulaStr = "~Media + num_genes_expressed",
                                       cores = 1), NA))

test_that("markerDiffTable works when 'balanced set to true'", 
          expect_error(markerDiffTable(HSMM[expressed_genes,],
                                       cth,
                                       balanced = TRUE,
                                       residualModelFormulaStr = "~Media + num_genes_expressed",
                                       cores = 1), NA))

test_that("markerDiffTable works when 'verbose set to true'", 
          expect_error(markerDiffTable(HSMM[expressed_genes,],
                                       cth,
                                       verbose = TRUE,
                                       residualModelFormulaStr = "~Media + num_genes_expressed",
                                       cores = 1), NA))

test_that("markerDiffTable throws error if cds is not of type 'CellDataSet'", 
          expect_error(markerDiffTable(cth,
                                       cth,
                                       balanced = TRUE,
                                       residualModelFormulaStr = "~Media + num_genes_expressed",
                                       cores = 1), "Error cds is not of type 'CellDataSet'"))

test_that("markerDiffTable throws error if cds is not of type 'CellDataSet'", 
          expect_error(markerDiffTable(HSMM[expressed_genes,],
                                       HSMM[expressed_genes],
                                       balanced = TRUE,
                                       residualModelFormulaStr = "~Media + num_genes_expressed",
                                       cores = 1), "Error cth is not of type 'CellTypeHierarchy'"))

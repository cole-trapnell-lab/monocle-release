library(HSMMSingleCell)
library(monocle)
library(DDRTree)
context("fitModels")

test_that("fitModels() properly validates its input",{
  load("HSMM_gene_annotation.rda")
  load("HSMM_expr_matrix.rda")
  load("HSMM_sample_sheet.rda")
  pd <- new("AnnotatedDataFrame", data = HSMM_sample_sheet)
  fd <- new("AnnotatedDataFrame", data = HSMM_gene_annotation)
  HSMM <- newCellDataSet(as.matrix(HSMM_expr_matrix), phenoData = pd, featureData = fd)
  rpc_matrix <- relative2abs(HSMM, cores = detectCores())
  HSMM <- newCellDataSet(as(as.matrix(rpc_matrix), "sparseMatrix"),
                         phenoData = pd,
                         featureData = fd,
                         lowerDetectionLimit = 1,
                         expressionFamily = negbinomial.size())
  HSMM <- estimateSizeFactors(HSMM)
  HSMM <- estimateDispersions(HSMM)
  HSMM <- detectGenes(HSMM, min_expr = 0.1)
  expressed_genes <- row.names(subset(fData(HSMM), num_cells_expressed >= 10))
  pData(HSMM)$Total_mRNAs <- Matrix::colSums(exprs(HSMM))
  HSMM <- HSMM[,row.names(subset(pData(HSMM), Total_mRNAs >= 10000 & Total_mRNAs <= 400000))]
  HSMM <- detectGenes(HSMM, min_expr = 0.1)
  L <- log(exprs(HSMM[expressed_genes,])) 
  MYF5_id <- row.names(subset(fData(HSMM), gene_short_name == "MYF5"))
  ANPEP_id <- row.names(subset(fData(HSMM), gene_short_name == "ANPEP"))
  cth <- newCellTypeHierarchy()
  cth <- addCellType(cth, "Myoblast", classify_func = function(x) {x[MYF5_id,] >= 1})
  cth <- addCellType(cth, "Fibroblast", classify_func = function(x) {x[MYF5_id,] < 1 & x[ANPEP_id,] > 1})
  HSMM <- classifyCells(HSMM, cth, 0.1)
  
  disp_table <- dispersionTable(HSMM)
  unsup_clustering_genes <- subset(disp_table, mean_expression >= 2 & dispersion_empirical >= 1 * dispersion_fit)
  HSMM <- setOrderingFilter(HSMM, unsup_clustering_genes$gene_id)
  
  
  HSMM <- clusterCells(HSMM, num_clusters = 2)
  
  HSMM <- clusterCells(HSMM, residualModelFormulaStr = "~Media + num_genes_expressed", num_clusters = 2)
  
  marker_diff <- markerDiffTable(HSMM[expressed_genes,],
                                 cth,
                                 residualModelFormulaStr = "~Media",
                                 cores = detectCores())
  candidate_clustering_genes <- row.names(subset(marker_diff, qval < 0.05))
  marker_spec <- calculateMarkerSpecificity(HSMM[candidate_clustering_genes,], cth)
  
  semisup_clustering_genes <- unique(selectTopMarkers(marker_spec, 10)$gene_id)
  HSMM <- setOrderingFilter(HSMM, semisup_clustering_genes)
  
  HSMM <- clusterCells(HSMM,
                       num_clusters = 2,
                       clustering_genes = semisup_clustering_genes,
                       residualModelFormulaStr = "~Media + num_genes_expressed")
  
  HSMM <- clusterCells(HSMM, 
                       num_clusters = 2,
                       cell_type_hierarchy = cth,
                       frequency_thresh = 0.1,
                       clustering_genes = row.names(subset(marker_diff, qval < 0.05)),
                       residualModelFormulaStr = "~Media + num_genes_expressed")
  
  HSMM_myo <- HSMM[,pData(HSMM)$CellType == "Myoblast"]
  HSMM_myo <- estimateDispersions(HSMM_myo)
  
  disp_table <- dispersionTable(HSMM_myo)
  ordering_genes <- subset(disp_table,
                           mean_expression >= 1 &
                             dispersion_empirical >= 2 * dispersion_fit)$gene_id
  HSMM_myo <- setOrderingFilter(HSMM_myo, ordering_genes)
  HSMM_myo <- reduceDimension(HSMM_myo, max_components = 2)
  HSMM_myo <- orderCells(HSMM_myo, reverse = FALSE)
  
  HSMM_filtered <- HSMM_myo[expressed_genes,]
  my_genes <- row.names(subset(fData(HSMM_filtered),
                               gene_short_name %in% c("CDK1", "MEF2C", "MYH3")))
  cds_subset <- HSMM_filtered[my_genes,]
  
  CCNB2_id <- row.names(subset(fData(HSMM), gene_short_name == "CCNB2"))
  MYH3_id <- row.names(subset(fData(HSMM), gene_short_name == "MYH3"))
  
  cth <- newCellTypeHierarchy()
  cth <- addCellType(cth, "Cycling myoblast", classify_func = function(x) {x[CCNB2_id,] >= 1})
  cth <- addCellType(cth, "Myotube", classify_func = function(x) {x[MYH3_id,] >= 1})
  
  marker_diff <- markerDiffTable(HSMM_myo[expressed_genes,], cth, cores = detectCores())
  semisup_clustering_genes <- row.names(subset(marker_diff, qval < 0.01))
  length(semisup_clustering_genes)
  
  HSMM_myo <- setOrderingFilter(HSMM_myo, semisup_clustering_genes)
  
  HSMM_myo <- reduceDimensions(HSMM_myo, max_components = 2)
  HSMM_myo <- orderCells(HSMM_myo, reverse = FALSE)
  
  HSMM_filtered <- HSMM_myo[expressed_genes,]
  my_genes <- row.names(subset(fData(HSMM_filtered), 
                               gene_short_name %in% c("CDK1", "MEF2C", "MYOG")))
  cds_subset <- HSMM_filtered[my_genes,]
  
  marker_genes <- row.names(subset(fData(HSMM_myo),
                                   gene_short_name %in% c("MEF2C", "MEF2D", "MYF5",
                                                          "ANPEP", "PDGFRA", "MYOG",
                                                          "TPM1", "TPM2", "MYH2",
                                                          "MYH3", "NCAM1", "TNNT1",
                                                          "TNNT2", "TNNC1", "CDK1",
                                                          "CDK2", "CCNB1", "CCNB2",
                                                          "CCND1", "CCNA1", "ID1")))
  diff_test_res <- differentialGeneTest(HSMM_myo[marker_genes,],
                                        fullModelFormulaStr = "~Media")
  sig_genes <- subset(diff_test_res, qval < 0.1)
  
  MYOG_ID1 <- HSMM_myo[row.names(subset(fData(HSMM_myo),
                                        gene_short_name %in% c("MYOG", "CCNB2"))),]
  
  to_be_tested <- row.names(subset(fData(HSMM),
                                   gene_short_name %in% c("UBC", "NCAM1", "ANPEP")))
  cds_subset <- HSMM[to_be_tested,]
  diff_test_res <- differentialGeneTest(cds_subset, fullModelFormulaStr = "~CellType")
  
  full_model_fits <- fitModel(cds_subset, modelFormulaStr = "~CellType")
  reduced_model_fits <- fitModel(cds_subset, modelFormulaStr = "~1")
  diff_test_res <- compareModels(full_model_fits, reduced_model_fits)
  expect_equal(colnames(diff_test_res), c("status", "family", "pval", "qval"))
  expect_equal(levels(diff_test_res$status),"OK")
  expect_equal(levels(diff_test_res$family), "negbinomial.size")
  expect_equal(diff_test_res$pval, c(1.956459e-32, 8.9011122e-01, 1.490169e-45))
  expect_equal(diff_test_res$qval, c(2.934689e-32, 8.9011122e-01, 4.470506e-45))
})


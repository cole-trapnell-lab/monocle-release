library(HSMMSingleCell)
library(monocle)
context("clusterGenes")

test_that("clusterGenes() properly validates its input",{
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
  
  marker_diff <- markerDiffTable(HSMM_myo[expressed_genes,], cth, cores = detectCores())
  semisup_clustering_genes <- row.names(subset(marker_diff, qval < 0.01))
  
  HSMM_myo <- setOrderingFilter(HSMM_myo, semisup_clustering_genes)
  
  HSMM_myo <- reduceDimensions(HSMM_myo, max_components = 2)
  HSMM_myo <- orderCells(HSMM_myo, reverse = FALSE)
  
  HSMM_filtered <- HSMM_myo[expressed_genes,]
  
  full_model_fits <- fitModel(HSMM[sample(nrow(fData(HSMM_filtered)), 100),], modelFormulaStr="~sm.ns(Pseudotime)")
  expression_curve_matrix <- responseMatrix(full_model_fits)
  clusters <- clusterGenes(expression_curve_matrix, k=4)
   
  expect_equal(names(clusters), c("medoids", "id.med", "clustering", "objective", "isolation", "clusinfo", "silinfo", "diss", "call", "exprs"))
  for(i in 1:length(clusters$id.med)) {
   expect_gte(clusters$id.med[i], 1)
  }
  expect_equal(substring(clusters$medoids, 0, 8), c("ENSG0000", "ENSG0000", "ENSG0000", "ENSG0000")) 
})
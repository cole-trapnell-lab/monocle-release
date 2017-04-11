library(HSMMSingleCell)
library(monocle)
context("clusterCells")

test_that("clusterCells() properly validates its input",{
  lung <- load_lung()
  lung <- clusterCells(lung, num_clusters = 2)
  lung_pData <- pData(lung)
  
  expect_equal(colnames(lung_pData), c("file", "total_mass", "internal_scale", "external_scale", "median_transcript_frags", "BioSample", "age", "genotype", "Sample.Name", "SRA.Sample", "MBases", "MBytes", 
                                       "SRA.Study", "BioProject", "source_name", "strain", "tissue", "Assay.Type", "Center.Name", "Platform", "Consent", "Time", "Size_Factor", "Total_mRNAs", "endogenous_RNA", 
                                       "Pseudotime", "State", "Parent", "num_genes_expressed", "Cluster"))
  for(i in 1:length(lung_pData$Pseudotime)) {
    expect_gte(lung_pData$Pseudotime[i], 0)
    expect_lte(lung_pData$Pseudotime[i], 100)
  }
  expect_equal(levels(lung_pData$Cluster), c("1", "2"))
  expect_equal(levels(lung_pData$State), c("1", "2", "3"))
  for(i in 1:length(lung_pData$Size_Factor)) {
    expect_gte(lung_pData$Size_Factor[i], 0.2796534)
    expect_lte(lung_pData$Size_Factor[i], 5)
  }
  for(i in 1:length(lung_pData$num_genes_expressed)) {
    expect_gte(lung_pData$num_genes_expressed[i], 10)
    expect_lte(lung_pData$num_genes_expressed[i], 196)
  }
  for(i in 1:length(lung_pData$Total_mRNAs)) {
  expect_gte(lung_pData$Total_mRNAs[i], 50.42859)
  expect_lte(lung_pData$Total_mRNAs[i], 16184.71)
  }
  
  pd <- new("AnnotatedDataFrame", data = HSMM_sample_sheet)
  fd <- new("AnnotatedDataFrame", data = HSMM_gene_annotation)
  HSMM <- newCellDataSet(as.matrix(HSMM_expr_matrix), phenoData = pd, featureData = fd)
  HSMM <- estimateSizeFactors(HSMM)
  HSMM <- estimateDispersions(HSMM)
  HSMM <- detectGenes(HSMM, min_expr = 0.1)
  cth <- newCellTypeHierarchy()
  cth <- addCellType(cth, "Myoblast", classify_func = function(x) {x[MYF5_id,] >= 1})
  cth <- addCellType(cth, "Fibroblast", classify_func = function(x)
    {x[MYF5_id,] < 1 & x[ANPEP_id,] > 1})
  HSMM <- classifyCells(HSMM, cth, 0.1)
  disp_table <- dispersionTable(HSMM)
  unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)
  HSMM <- setOrderingFilter(HSMM, unsup_clustering_genes$gene_id)
  
  #"DDRTree", "ICA", 'tSNE', "SimplePPT", 'L1-graph', 'SGL-tree'
  
  HSMM_DDRTree <- reduceDimension(HSMM, max_components = 2, num_dim = 5, reduction_method = 'DDRTree', verbose = T)
  HSMM_ICA <- reduceDimension(HSMM, max_components = 2, num_dim = 5, reduction_method = 'ICA', verbose = T)
  HSMM_tSNE <- reduceDimension(HSMM, max_components = 2, num_dim = 5, reduction_method = 'tSNE', verbose = T)
  HSMM_SimplePPT <- reduceDimension(HSMM, max_components = 2, num_dim = 5, reduction_method = 'SimplePPT', verbose = T)
  HSMM_L1-graph <- reduceDimension(HSMM, max_components = 2, num_dim = 5, reduction_method = 'L1-graph', verbose = T)
  HSMM_SGL-tree <- reduceDimension(HSMM, max_components = 2, num_dim = 5, reduction_method = 'SGL-tree', verbose = T)
  
  HSMM_DDRTree <- clusterCells(HSMM_DDRTree, num_clusters = 2, clustering_genes = unsup_clustering_genes$gene_id)
  HSMM_ICA <- clusterCells(HSMM_ICA, num_clusters = 2, clustering_genes = unsup_clustering_genes$gene_id)
  HSMM_tSNE <- clusterCells(HSMM_tSNE, num_clusters = 2, clustering_genes = unsup_clustering_genes$gene_id)
  HSMM_SimplePPT <- clusterCells(HSMM_SimplePPT, num_clusters = 2, clustering_genes = unsup_clustering_genes$gene_id)
  HSMM_L1-graph <- clusterCells(HSMM_L1-graph, num_clusters = 2, clustering_genes = unsup_clustering_genes$gene_id)
  
})
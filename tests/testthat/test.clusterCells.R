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
})
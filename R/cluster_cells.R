#' Cluster cells into a specified number of groups.
#' 
#' Unsupervised clustering of cells is a common step in many single-cell expression
#' workflows. In an experiment containing a mixture of cell types, each cluster might
#' correspond to a different cell type. This method takes a CellDataSet as input
#' along with a requested number of clusters, clusters them with an unsupervised 
#' algorithm, and then returns the CellDataSet with the cluster assignments stored in
#' the pData table.
#' 
#' @param cds the CellDataSet upon which to perform this operation
#' @param num_clusters number of desired cell clusters
#' @param num_reduced_dims number of dimensions to reduce the dataset to
#' @param residualModelFormulaStr A model formula specifying the effects to subtract from the data before clustering.
#' @param ddrtree_gamma gamma parameter for DDRTree
#' @param verbose Verbose parameter for DDRTree
#' @return an updated CellDataSet object, in which phenoData contains values for Cluster for each cell
#' @export
clusterCells <- function(cds, 
                         num_clusters=2, 
                         num_reduced_dims=10, 
                         residualModelFormulaStr=NULL,
                         ddrtree_gamma=100,
			                   verbose = F) {
  
  # disp_table <- dispersionTable(cds)
  # ordering_genes <- row.names(subset(disp_table, dispersion_empirical >= 2 * dispersion_fit))
  # cds <- setOrderingFilter(cds, ordering_genes)
  
  cds <- reduceDimension(cds, 
                         max_components=num_reduced_dims, 
                         residualModelFormulaStr=residualModelFormulaStr,
                         use_vst=T, 
                         method = "DDRTree",
			                   pseudo_expr=0, 
                         verbose=verbose,
                         maxIter=20,
                         param.gamma=ddrtree_gamma,
                         ncenter=num_clusters)
  pData(cds)$Cluster <- as.factor(cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex)

  cds

}

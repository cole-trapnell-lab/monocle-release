#' Cluster cells to nearest DDRTree center.
#' @param cds the CellDataSet upon which to perform this operation
#' @param max_components number of dimensions to reduce the dataset to
#' @param ncenter number of desired cell clusters
#' @return an updated CellDataSet object, in which phenoData contains values for Cluster for each cell
#' @export
clusterCells <- function(cds, max_comp=18, ncent=15, ...) {

  cds <- reduceDimension(cds, 
                         max_components=max_comp, 
                         ncenter=ncent,
			 ...)


  pData(cds)$Cluster <- as.factor(cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex)

  cds

}

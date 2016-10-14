#' Clusters genes by pseudotime trend.
#'
#' This function takes a matrix of expression values and performs k-means 
#' clustering on the genes. 
#'
#' @param expr_matrix A matrix of expression values to cluster together. Rows 
#' are genes, columns are cells.
#' @param k How many clusters to create
#' @param method The distance function to use during clustering
#' @param ... Extra parameters to pass to pam() during clustering
#' @return a pam cluster object
#' @importFrom cluster pam
#' @importFrom stats as.dist cor
#' @export
#' @examples
#' \dontrun{
#' full_model_fits <- fitModel(HSMM[sample(nrow(fData(HSMM_filtered)), 100),],  
#'    modelFormulaStr="~sm.ns(Pseudotime)")
#' expression_curve_matrix <- responseMatrix(full_model_fits)
#' clusters <- clusterGenes(expression_curve_matrix, k=4)
#' plot_clusters(HSMM_filtered[ordering_genes,], clusters)
#' }
clusterGenes<-function(expr_matrix, k, method=function(x){as.dist((1 - cor(Matrix::t(x)))/2)}, ...){
  expr_matrix <- expr_matrix[rowSums(is.na(expr_matrix)) == 0,]
  expr_matrix <- expr_matrix[is.nan(rowSums(expr_matrix)) == FALSE,] 
  expr_matrix[is.na(expr_matrix)] <- 0
  n<-method(expr_matrix)
  clusters<-cluster::pam(n,k, ...)
  class(clusters)<-"list"
  clusters$exprs<-expr_matrix
  clusters
}

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
#' @param cell_type_hierarchy the CellTypeHierarchy that divides the cells from cds into different types of cells
#' @param max_components number of dimensions to project the data into via \code{\link{reduceDimension}()}
#' @param frequency_thresh When a CellTypeHierarchy is provided, cluster cells will impute cell types in clusters that are composed of at least this much of exactly one cell type.
#' @param clustering_genes a vector of genes used to differentiate between the cell types in the CellTypeHierarchy
#' @param residualModelFormulaStr A model formula specifying the effects to subtract from the data before clustering.
#' @param param.gamma gamma parameter for DDRTree
#' @param verbose Verbose parameter for DDRTree
#' @param ... Additional arguments passed to \code{\link{reduceDimension}()}
#' @importFrom Biobase fData pData pData<-
#' @return an updated CellDataSet object, in which phenoData contains values for Cluster for each cell
#' @export
clusterCells <- function(cds, 
                         num_clusters, 
                         cell_type_hierarchy=NULL,
                         frequency_thresh=0.10,
                         clustering_genes=NULL,
                         max_components=10, 
                         residualModelFormulaStr=NULL,
                         param.gamma=100,
                         verbose = F, 
                         ...) {
  
  # disp_table <- dispersionTable(cds)
  # ordering_genes <- row.names(subset(disp_table, dispersion_empirical >= 2 * dispersion_fit))
  # cds <- setOrderingFilter(cds, ordering_genes)
  use_for_ordering <- NA
  if (is.null(fData(cds)$use_for_ordering) == FALSE)
    old_ordering_genes <- row.names(subset(fData(cds), use_for_ordering)) 
  else
    old_ordering_genes <- NULL
  
  if (is.null(clustering_genes) == FALSE) 
    cds <- setOrderingFilter(cds, clustering_genes)
  
  cds <- reduceDimension(cds, 
                         max_components=max_components, 
                         residualModelFormulaStr=residualModelFormulaStr,
                         reduction_method = "DDRTree",
                         verbose=verbose,
                         param.gamma=param.gamma,
                         ncenter=num_clusters, 
                         ...)
  pData(cds)$Cluster <- as.factor(cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex)
  
  if (is.null(old_ordering_genes) == FALSE)
    cds <- setOrderingFilter(cds, old_ordering_genes)
  
  if (is.null(cell_type_hierarchy) == FALSE)
    cds <- classifyCells(cds, cell_type_hierarchy, frequency_thresh, "Cluster")
  
  return(cds)
}

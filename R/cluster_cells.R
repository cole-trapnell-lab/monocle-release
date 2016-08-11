#' Cluster cells into a specified number of groups based on .
#' 
#' Unsupervised clustering of cells is a common step in many single-cell expression
#' workflows. In an experiment containing a mixture of cell types, each cluster might
#' correspond to a different cell type. This method takes a CellDataSet as input
#' along with a requested number of clusters, clusters them with an unsupervised 
#' algorithm, and then returns the CellDataSet with the cluster assignments stored in
#' the pData table. When number of clusters is set to NULL (num_clusters = NULL), 
#' the decision plot as introduced in the above citation will be plotted and the 
#' users are required to click on the decision plot to select the rho and delta 
#' to determine the number of clusters to cluster.  
#' 
#' 
#' @param cds the CellDataSet upon which to perform this operation
#' @param skip_rho_sigma A logic flag to determine whether or not you want to skip the calculation of rho / sigma 
#' @param variance_explained Variance explained by the PCA components, used to select the component numbers for tSNE 
#' @param inspect_rho_sigma A logical flag to determine whether or not you want to interactively select the rho and sigma for assigning up clusters
#' @param num_clusters Number of clusters you wanted. If you specify this argument, the rho and sigma used to define the threshold for density cluster
#' will automatically calculated based on the top num_cluster product of rho and sigma. 
#' @param frequency_thresh When a CellTypeHierarchy is provided, cluster cells will impute cell types in clusters that are composed of at least this much of exactly one cell type.
#' @param verbose Verbose parameter for DDRTree
#' @param ... Additional arguments passed to \code{\link{densityClust}()}
#' @return an updated CellDataSet object, in which phenoData contains values for Cluster for each cell
#' @import densityClust
#' @references Rodriguez, A., & Laio, A. (2014). Clustering by fast search and find of density peaks. Science, 344(6191), 1492-1496. doi:10.1126/science.1242072
#' @export

clusterCells_Density_Peak <- function(cds, 
                                      skip_rho_sigma = F, 
                                      variance_explained = 0.8, 
                                      inspect_rho_sigma = F, 
                                      rho_val = NULL, 
                                      delta_val = NULL, 
                                      cell_type_hierarchy=NULL,
                                      frequency_thresh=0.10,
                                      clustering_genes=NULL,
                                      verbose = F, 
                                      ...) {
  tsne_data <- reducedDimA(cds)
  if(ncol(tsne_data) != ncol(cds))
    stop("reduced dimension space doesn't match the dimension of the CellDataSet object")

  dataDist <- dist(t(tsne_data)) #calculate distances between cells

  if(skip_rho_sigma 
    & !is.null(cds@auxClusteringData[["tSNE"]]$densityPeak) 
    & !is.null(pData(cds)$Cluster)
    & !is.null(pData(cds)$peaks)
    & !is.null(pData(cds)$halo)
    & !is.null(pData(cds)$delta)
    & !is.null(pData(cds)$rho)) {
    dataClust <- cds@auxClusteringData[["tSNE"]]$densityPeak
    dataClust$rho <- pData(cds)$rho
    dataClust$delta <- pData(cds)$delta
    dataClust$distance <- dataDist
    dataClust$peaks <- pData(cds)$peaks
    dataClust$clusters <- pData(cds)$clusters
    dataClust$halo <- pData(cds)$halo

    # res <- list(rho=rho, delta=delta, distance=distance, dc=dc, threshold=c(rho=NA, delta=NA), peaks=NA, clusters=NA, halo=NA)
    dataClust <- dataClust[c('rho', 'delta', 'distance', 'dc', 'threshold', 'peaks', 'clusters', 'halo')]
    class(dataClust) <- 'densityCluster'
  }
  else{

    #finally use density peak to determine the number of clusters
    if (verbose) 
        message("Run densityPeak algorithm to automatically cluster cells based on distance of cells on tSNE components...")

    dataClust <- densityClust::densityClust(dataDist, ...) #gaussian = F
  }
  #automatically find the rho / sigma based on the number of cells you want: 
  if(!is.null(rho_val) & !is.null(delta_val)){
    if(verbose)
      message('Use the user provided rho and delta for assigning the density peaks and clusters')
  }
  else 
  {
    if(verbose)
      message('Use 0.9 of the delta and 0.95 of the rho as the cutoff for assigning density peaks and clusters')

    rho_val <- quantile(dataClust$rho, probs = 0.95)
    delta_val <- quantile(dataClust$delta, probs = 0.90)
  }

  #automatically pick up the rho and delta values: 
  if(inspect_rho_sigma == F)
    dataClust <- densityClust::findClusters(dataClust, rho = rho_val, delta = delta_val)
  else {
   if (verbose) 
      message("Please click on the decision plot to select rho and delta for density peak clustering...")
   
    dataClust <- densityClust::findClusters(dataClust)
  }
  
  #find the number of clusters: 
  #cluster_num <- length(unique(dataClust$clusters))
  
  pData(cds)$Cluster <- dataClust$clusters
  pData(cds)$peaks <- F
  pData(cds)$peaks[dataClust$peaks] <- T
  pData(cds)$halo <- dataClust$halo
  pData(cds)$delta <- dataClust$delta
  pData(cds)$rho <- dataClust$rho
  
  if (is.null(cell_type_hierarchy) == FALSE)
    cds <- classifyCells(cds, cell_type_hierarchy, frequency_thresh, "Cluster")
  
  cds@auxClusteringData[["tSNE"]]$densityPeak <- dataClust[c("dc", "threshold")] #, "peaks"
  
  return(cds)

}
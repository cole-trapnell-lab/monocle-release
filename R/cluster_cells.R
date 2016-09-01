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
#' @param number_clusters Number of clusters. The algorithm use 0.5 of the rho as the threshold of rho and the delta 
#' corresponding to the number_clusters sample with the highest delta as the density peaks and for assigning clusters
#' @param inspect_rho_sigma A logical flag to determine whether or not you want to interactively select the rho and sigma for assigning up clusters
#' @param rho_threshold The threshold of local density (rho) used to select the density peaks 
#' @param delta_threshold The threshold of local distance (delta) used to select the density peaks 
#' @param peaks A numeric vector indicates the index of density peaks used for clustering. This vector should be retrieved from the decision plot with caution. No checking involved.  
#' will automatically calculated based on the top num_cluster product of rho and sigma. 
#' @param gaussian A logic flag passed to densityClust function in desnityClust package to determine whether or not Gaussian kernel will be used for calculating the local density
#' @param frequency_thresh When a CellTypeHierarchy is provided, cluster cells will impute cell types in clusters that are composed of at least this much of exactly one cell type.
#' @param verbose Verbose parameter for DDRTree
#' @param ... Additional arguments passed to \code{\link{densityClust}()}
#' @return an updated CellDataSet object, in which phenoData contains values for Cluster for each cell
#' @import densityClust
#' @references Rodriguez, A., & Laio, A. (2014). Clustering by fast search and find of density peaks. Science, 344(6191), 1492-1496. doi:10.1126/science.1242072
#' @export

clusterCells_Density_Peak <- function(cds, 
                                      skip_rho_sigma = F, 
                                      number_clusters = NULL, 
                                      inspect_rho_sigma = F, 
                                      rho_threshold = NULL, 
                                      delta_threshold = NULL, 
                                      peaks = NULL,
                                      gaussian = T, 
                                      cell_type_hierarchy=NULL,
                                      frequency_thresh=NULL,
                                      enrichment_thresh=NULL,
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

    dataClust <- densityClust::densityClust(dataDist, gaussian = gaussian, verbose = verbose) #gaussian = F
  }
  #automatically find the rho / sigma based on the number of cells you want: 
  if(!is.null(rho_threshold) & !is.null(delta_threshold)){
    if(verbose)
      message('Use the user provided rho and delta for assigning the density peaks and clusters')
  }
  else 
  {
    if(is.null(number_clusters)) {
      if(verbose)
        message('Use 0.95 of the delta and 0.95 of the rho as the cutoff for assigning density peaks and clusters')

      rho_threshold <- quantile(dataClust$rho, probs = 0.95)
      delta_threshold <- quantile(dataClust$delta, probs = 0.95)
    }
    else {
        if(verbose)
          message(paste('Use 0.5 of the rho as the cutoff and first', number_clusters , 'samples with highest delta as the density peaks and for assigning clusters'))
      delta_rho_df <- data.frame(delta = dataClust$delta, rho = dataClust$rho)
      rho_valid_threshold <- quantile(dataClust$rho, probs = 0.5)
      delta_rho_df <- subset(delta_rho_df, rho > rho_valid_threshold) 
      threshold_ind <- order(dataClust$delta)[number_clusters]
      delta_threshold <- delta_rho_df$delta[threshold_ind] 
      rho_threshold <- delta_rho_df$rho[threshold_ind]
    }
  }

  #automatically pick up the rho and delta values: 
  if(inspect_rho_sigma == F)
    dataClust <- densityClust::findClusters(dataClust, rho = rho_threshold, delta = delta_threshold, peaks = peaks, verbose = verbose)
  else {
   if (verbose) 
      message("Please click on the decision plot to select rho and delta for density peak clustering...")
   
    dataClust <- densityClust::findClusters(dataClust)
  }
  
  #find the number of clusters: 
  #cluster_num <- length(unique(dataClust$clusters))
  
  pData(cds)$Cluster <- factor(dataClust$clusters)
  pData(cds)$peaks <- F
  pData(cds)$peaks[dataClust$peaks] <- T
  pData(cds)$halo <- dataClust$halo
  pData(cds)$delta <- dataClust$delta
  pData(cds)$rho <- dataClust$rho
  
  if (is.null(cell_type_hierarchy) == FALSE)
    cds <- classifyCells(cds, cell_type_hierarchy, frequency_thresh, enrichment_thresh, "Cluster")
  
  cds@auxClusteringData[["tSNE"]]$densityPeak <- dataClust[c("dc", "threshold")] #, "peaks"
  
  return(cds)

}
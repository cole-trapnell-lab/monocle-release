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
  clusters$exprs <- expr_matrix
  clusters
}

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
#' @param num_clusters Number of clusters. The algorithm use 0.5 of the rho as the threshold of rho and the delta 
#' corresponding to the number_clusters sample with the highest delta as the density peaks and for assigning clusters
#' @param inspect_rho_sigma A logical flag to determine whether or not you want to interactively select the rho and sigma for assigning up clusters
#' @param rho_threshold The threshold of local density (rho) used to select the density peaks 
#' @param delta_threshold The threshold of local distance (delta) used to select the density peaks 
#' @param peaks A numeric vector indicates the index of density peaks used for clustering. This vector should be retrieved from the decision plot with caution. No checking involved.  
#' will automatically calculated based on the top num_cluster product of rho and sigma. 
#' @param gaussian A logic flag passed to densityClust function in desnityClust package to determine whether or not Gaussian kernel will be used for calculating the local density
#' @param frequency_thresh When a CellTypeHierarchy is provided, cluster cells will impute cell types in clusters that are composed of at least this much of exactly one cell type.
#' @param verbose Verbose parameter for DDRTree
#' @param cell_type_hierarchy A data structure used for organizing functions that can be used for organizing cells 
#' @param enrichment_thresh fraction to be multipled by each cell type percentage. Only used if frequency_thresh is NULL, both cannot be NULL
#' @param clustering_genes a vector of feature ids (from the CellDataSet's featureData) used for ordering cells
#' @param method method for clustering cells. By default, we use density peak clustering algorithm for clustering. 
#' The other method is based on DDRTree. 
#' @param ... Additional arguments passed to \code{\link{densityClust}()}
#' @return an updated CellDataSet object, in which phenoData contains values for Cluster for each cell
#' @importFrom densityClust densityClust findClusters
#' @references Rodriguez, A., & Laio, A. (2014). Clustering by fast search and find of density peaks. Science, 344(6191), 1492-1496. doi:10.1126/science.1242072
#' @export

clusterCells <- function(cds, 
                         skip_rho_sigma = F, 
                         num_clusters = NULL, 
                         inspect_rho_sigma = F, 
                         rho_threshold = NULL, 
                         delta_threshold = NULL, 
                         peaks = NULL,
                         gaussian = T, 
                         cell_type_hierarchy=NULL,
                         frequency_thresh=NULL,
                         enrichment_thresh=NULL,
                         clustering_genes=NULL,
                         method = c('densityPeak', 'DDRTree'),
                         verbose = F, 
                         ...) {
  method <- match.arg(method)
  
  if(method == 'DDRTree') { # this option will be removed in future
    ################### DDRTREE START ###################
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
                           max_components=2, #
                           residualModelFormulaStr=NULL,
                           reduction_method = "DDRTree",
                           verbose=verbose,
                           param.gamma=100,
                           ncenter=num_clusters, 
                           ...)
    pData(cds)$Cluster <- as.factor(cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex)
    
    if (is.null(old_ordering_genes) == FALSE)
      cds <- setOrderingFilter(cds, old_ordering_genes)
    
    if (is.null(cell_type_hierarchy) == FALSE)
      cds <- classifyCells(cds, cell_type_hierarchy, frequency_thresh, enrichment_thresh, "Cluster")
    
    return(cds)
   ################### DDRTREE END ###################
  } else if(method == 'densityPeak'){
    ################### DENSITYPEAK START ###################
    set.seed(2017)
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
      
    } else {
      #finally use density peak to determine the number of clusters
      if (verbose) {
        message("Run densityPeak algorithm to automatically cluster cells based on distance of cells on tSNE components...")
      }
      dataClust <- densityClust::densityClust(dataDist, gaussian = gaussian) #gaussian = F
    }
    #automatically find the rho / sigma based on the number of cells you want: 
    if(!is.null(rho_threshold) & !is.null(delta_threshold)){
      if(verbose) {
        message('Use the user provided rho and delta for assigning the density peaks and clusters')
      }
    } else {
      if(is.null(num_clusters)) {
        if(verbose) {
          message('Use 0.95 of the delta and 0.95 of the rho as the cutoff for assigning density peaks and clusters')
        }
        
        rho_threshold <- quantile(dataClust$rho, probs = 0.95)
        delta_threshold <- quantile(dataClust$delta, probs = 0.95)
      } else {
        if(verbose) {
          message(paste('Use 0.5 of the rho as the cutoff and first', num_clusters , 'samples with highest delta as the density peaks and for assigning clusters'))
        }

        delta_rho_df <- data.frame(delta = dataClust$delta, rho = dataClust$rho)
        rho_valid_threshold <- quantile(dataClust$rho, probs = 0.5)
        delta_rho_df <- subset(delta_rho_df, rho > rho_valid_threshold) 
        threshold_ind <- order(delta_rho_df$delta, decreasing = T)[num_clusters + 1]
        candidate_peaks <- subset(delta_rho_df, delta >= delta_rho_df$delta[threshold_ind])
        #delta_threshold <- min(candidate_peaks$delta) - .Machine$double.eps
        delta_threshold <- delta_rho_df$delta[threshold_ind] #- 10*.Machine$double.eps 
        rho_threshold <- min(candidate_peaks$rho) - 10*.Machine$double.eps
        #head(subset(delta_rho_df, delta > delta_threshold & rho > rho_threshold))
        #delta_threshold <- delta_rho_df$delta[threshold_ind] - .Machine$double.eps 
        #rho_threshold <- delta_rho_df$rho[threshold_ind] - .Machine$double.eps
      }
    }
    
    #automatically pick up the rho and delta values: 

    if(inspect_rho_sigma == F) {
      dataClust <- densityClust::findClusters(dataClust, rho = rho_threshold, delta = delta_threshold, peaks=peaks)
    } else {
      if (verbose) 
        message("Please click on the decision plot to select rho and delta for density peak clustering...")
      
      dataClust <- densityClust::findClusters(dataClust)
    }
    
    #find the number of clusters: 
    #cluster_num <- length(unique(dataClust$clusters))
    
    pData(cds)$Cluster <- as.character(dataClust$clusters)
    pData(cds)$peaks <- F
    pData(cds)$peaks[dataClust$peaks] <- T
    pData(cds)$halo <- dataClust$halo
    pData(cds)$delta <- dataClust$delta
    pData(cds)$rho <- dataClust$rho
    
    if (is.null(cell_type_hierarchy) == FALSE) {
      cds <- classifyCells(cds, cell_type_hierarchy, frequency_thresh, enrichment_thresh, "Cluster")
    }
    
    cds@auxClusteringData[["tSNE"]]$densityPeak <- dataClust[c("dc", "threshold")] #, "peaks"
    
    return(cds)
    
  ################### DENSITYPEAK END ###################
  } else {
    stop('Cluster method ', method, ' is not implemented')
  }
}

#' #' Select features for constructing the developmental trajectory 
#' #' 
#' #' @param cds the CellDataSet upon which to perform this operation
#' #' @param num_cells_expressed A logic flag to determine whether or not you want to skip the calculation of rho / sigma 
#' #' @param num_dim Number of clusters. The algorithm use 0.5 of the rho as the threshold of rho and the delta 
#' #' corresponding to the number_clusters sample with the highest delta as the density peaks and for assigning clusters
#' #' @param rho_threshold The threshold of local density (rho) used to select the density peaks 
#' #' @param delta_threshold The threshold of local distance (delta) used to select the density peaks 
#' #' @param qval_threshold A logic flag passed to densityClust function in desnityClust package to determine whether or not Gaussian kernel will be used for calculating the local density
#' #' @param verbose Verbose parameter for DDRTree
#' #' @return an list contains a updated CellDataSet object after clustering and tree construction as well as a vector including the selected top significant differentially expressed genes across clusters of cells  
#' #' @references Rodriguez, A., & Laio, A. (2014). Clustering by fast search and find of density peaks. Science, 344(6191), 1492-1496. doi:10.1126/science.1242072
#' #' @export
#' dpFeature <- function(cds, num_cells_expressed = 5, num_dim = 5, rho_threshold = NULL, delta_threshold = NULL, qval_threshold = 0.01, verbose = F){
#'   #1. determine how many pca dimension you want:
#'   cds <- detectGenes(cds)
#'   fData(cds)$use_for_ordering[fData(cds)$num_cells_expressed > num_cells_expressed] <- T
#'   
#'   if(is.null(num_dim)){
#'     lung_pc_variance <- plot_pc_variance_explained(cds, return_all = T)
#'     ratio_to_first_diff <- diff(lung_pc_variance$variance_explained) / diff(lung_pc_variance$variance_explained)[1]
#'     num_dim <- which(ratio_to_first_diff < 0.1) + 1
#'   }
#'   
#'   #2. run reduceDimension with tSNE as the reduction_method
#'   # absolute_cds <- setOrderingFilter(absolute_cds, quake_id)
#'   cds <- reduceDimension(cds, return_all = F, max_components=2, norm_method = 'log', reduction_method = 'tSNE', num_dim = num_dim,  verbose = verbose)
#'   
#'   #3. initial run of clusterCells_Density_Peak
#'   cds <- clusterCells_Density_Peak(cds, rho_threshold = rho_threshold, delta_threshold = delta_threshold, verbose = verbose)
#'   
#'   #perform DEG test across clusters: 
#'   cds@expressionFamily <- negbinomial.size()
#'   pData(cds)$Cluster <- factor(pData(cds)$Cluster)
#'   clustering_DEG_genes <- differentialGeneTest(cds, fullModelFormulaStr = '~Cluster', cores = detectCores() - 2)
#'   clustering_DEG_genes_subset <- lung_clustering_DEG_genes[fData(cds)$num_cells_expressed > num_cells_expressed, ]
#'   
#'   #use all DEG gene from the clusters
#'   clustering_DEG_genes_subset <- clustering_DEG_genes_subset[order(clustering_DEG_genes_subset$qval), ]
#'   ordering_genes <- row.names(subset(clustering_DEG_genes, qval < qval_threshold))
#'   
#'   cds <- setOrderingFilter(cds, ordering_genes = lung_ordering_genes)
#'   cds <- reduceDimension(cds, norm_method = 'log', verbose = T)
#'   plot_cell_trajectory(cds, color_by = 'Time')
#'   
#'   return(list(new_cds = cds, ordering_genes = ordering_genes))
#' }
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
#' @param variance_explained Variance explained by the PCA components, used to select the component numbers for tSNE 
#' @param inspect_rho_sigma A logical flag to determine whether or not you want to interactively select the rho and sigma for assigning up clusters
#' @param frequency_thresh When a CellTypeHierarchy is provided, cluster cells will impute cell types in clusters that are composed of at least this much of exactly one cell type.
#' @param residualModelFormulaStr A model formula specifying the effects to subtract from the data before clustering.
#' @param param.gamma gamma parameter for DDRTree
#' @param verbose Verbose parameter for DDRTree
#' @param ... Additional arguments passed to \code{\link{reduceDimension}()}
#' @return an updated CellDataSet object, in which phenoData contains values for Cluster for each cell
#' @import densityClust
#' @references Rodriguez, A., & Laio, A. (2014). Clustering by fast search and find of density peaks. Science, 344(6191), 1492-1496. doi:10.1126/science.1242072
#' @export
clusterCells_Density_Peak <- function(cds, 
                                      variance_explained = 0.8, 
                                      inspect_rho_sigma = F, 
                                      num_clusters = NULL,
                                      cell_type_hierarchy=NULL,
                                      frequency_thresh=0.10,
                                      clustering_genes=NULL,
                                      max_components=2, 
                                      residualModelFormulaStr=NULL,
                                      param.gamma=100,
                                      verbose = F, 
                                      ...) {
  
  # disp_table <- dispersionTable(cds)
  # ordering_genes <- row.names(subset(disp_table, dispersion_empirical >= 2 * dispersion_fit))
  # cds <- setOrderingFilter(cds, ordering_genes)
  if (is.null(fData(cds)$use_for_ordering) == FALSE)
    old_ordering_genes <- row.names(subset(fData(cds), use_for_ordering)) 
  else
    old_ordering_genes <- NULL
  
  if (is.null(clustering_genes) == FALSE) 
    cds <- setOrderingFilter(cds, clustering_genes)
  
  FM <- normalize_expr_data(cds, ...)
  
  #FM <- FM[unlist(sparseApply(FM, 1, sd, convert_to_dense=TRUE)) > 0, ]
  xm <- Matrix::rowMeans(FM)
  xsd <- sqrt(Matrix::rowMeans((FM - xm)^2))
  FM <- FM[xsd > 0,]
   
  FM <- as.matrix(Matrix::t(scale(Matrix::t(FM))))
  
  if (nrow(FM) == 0) {
    stop("Error: all rows have standard deviation zero")
  }
  
  #first perform PCA 
  if (verbose) 
      message("Remove noise by PCA ...")

  pca_res <- prcomp(t(FM), center = T, scale = T)
  std_dev <- pca_res$sdev 
  pr_var <- std_dev^2
  prop_varex <- pr_var/sum(pr_var)
  num_dim <- min(which(cumsum(prop_varex) > variance_explained))
  topDim_pca <- pca_res$x[, 1:num_dim]
  
  #perform the model formula transformation right before tSNE: 
  if (is.null(residualModelFormulaStr) == FALSE) {
    if (verbose) 
      message("Removing batch effects")
    X.model_mat <- sparse.model.matrix(as.formula(residualModelFormulaStr), 
                                       data = pData(cds), drop.unused.levels = TRUE)
    
    fit <- limma::lmFit(topDim_pca, X.model_mat, ...)
    beta <- fit$coefficients[, -1, drop = FALSE]
    beta[is.na(beta)] <- 0
    topDim_pca <- as.matrix(FM) - beta %*% t(X.model_mat[, -1])
  }else{
    X.model_mat <- NULL
  }

  #then run tSNE 
  if (verbose) 
      message("Reduce dimension by tSNE ...")

  tsne_res <- Rtsne::Rtsne(as.matrix(topDim_pca), dims = max_components, pca = F,...)
  
  tsne_data <- tsne_res$Y[, 1:max_components]
  row.names(tsne_data) <- colnames(tsne_data)

  #finally use density peak to determine the number of clusters
  if (verbose) 
      message("Run densityPeak algorithm to automatically cluster cells based on distance of cells on tSNE components...")

  dataDist <- dist(tsne_data)
  dataClust <- densityClust::densityClust(dataDist, gaussian = F)
  
  #automatically find the rho / sigma based on the number of cells you want: 
  if(!is.null(num_clusters)){
    gamma <- dataClust$rho * dataClust$delta
    ind <- order(gamma)[num_clusters]
    rho_val <- dataClust$rho[ind]
    delta_val <- dataClust$delta[ind]
  }
  else 
  {
    rho_val <- quantile(dataClust$rho, probs = 0.90)
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
  # cluster_num <- length(unique(dataClust$clusters))
  
  pData(cds)$Cluster <- as.factor(dataClust$clusters)
  pData(cds)$peaks <- as.factor(dataClust$peaks)
  pData(cds)$halo <- as.factor(dataClust$halo)
  pData(cds)$delta <- as.factor(dataClust$delta)
  pData(cds)$rho <- as.factor(dataClust$rho)

  if (is.null(old_ordering_genes) == FALSE)
    cds <- setOrderingFilter(cds, old_ordering_genes)
  
  if (is.null(cell_type_hierarchy) == FALSE)
    cds <- classifyCells(cds, cell_type_hierarchy, frequency_thresh, "Cluster")
  
  reducedDimA(cds) <- t(tsne_data) #this may move to the auxClusteringData environment
  pData(cds)$State <- as.factor(dataClust$clusters)

  #set the important information from densityClust to certain part of the cds object: 
  cds@auxClusteringData[["tSNE"]]$pca_components_used <- num_dim
  cds@auxClusteringData[["tSNE"]]$reduced_dimension <- t(tsne_data) 
  cds@auxClusteringData[["tSNE"]]$densityPeak <- dataClust[c("dc", "threshold")] #, "peaks"
  
  return(cds)
}

#add the new plot for plotting the clusters: 
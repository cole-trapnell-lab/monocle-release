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
#' algorithm (by default, density peak clustering), and then returns the CellDataSet with the 
#' cluster assignments stored in the pData table. When number of clusters is set 
#' to NULL (num_clusters = NULL), the decision plot as introduced in the reference 
#' will be plotted and the users are required to check the decision plot to select 
#' the rho and delta to determine the number of clusters to cluster. When the dataset 
#' is big, for example > 50 k, we recommend the user to use the Louvain clustering 
#' algorithm which is inspired from phenograph paper. Note Louvain doesn't support the 
#' num_cluster argument but the k (number of k-nearest neighbors) is relevant to the final 
#' clustering number. The deafult implementation of Louvain clustering (when res is set to be NULL)
#' is based on the Rphenograph package but updated based on our requirement (for example,
#' changed the jaccard_coeff function as well as adding louvain_iter argument, etc.). We also 
#' support setting resolutioin parameter when performing louvain clustering using the louvain package 
#' from python. With different res values, the users can obtain different granularity of the data.  
#' 
#' @param cds the CellDataSet upon which to perform this operation
#' @param use_pca Whether or not to cluster cells based on top PCA component. Default to be FALSE. 
#' @param skip_rho_sigma A logic flag to determine whether or not you want to skip the calculation of rho / sigma 
#' @param num_clusters Number of clusters. The algorithm use 0.5 of the rho as the threshold of rho and the delta 
#' corresponding to the number_clusters sample with the highest delta as the density peaks and for assigning clusters
#' @param inspect_rho_sigma A logical flag to determine whether or not you want to interactively select the rho and sigma for assigning up clusters
#' @param rho_threshold The threshold of local density (rho) used to select the density peaks 
#' @param delta_threshold The threshold of local distance (delta) used to select the density peaks 
#' @param peaks A numeric vector indicates the index of density peaks used for clustering. This vector should be retrieved from the decision plot with caution. No checking involved.  
#' will automatically calculated based on the top num_cluster product of rho and sigma. 
#' @param gaussian A logic flag passed to densityClust function in desnityClust package to determine whether or not Gaussian kernel will be used for calculating the local density
#' @param clustering_genes a vector of feature ids (from the CellDataSet's featureData) used for ordering cells
#' @param k number of kNN used in creating the k nearest neighbor graph for Louvain clustering. The number of kNN is related to the resolution of the clustering result, bigger number of kNN gives low resolution and vice versa. Default to be 20
#' @param louvain_iter Integer number of iterations used for Louvain clustering. The clustering result gives the largest modularity score will be used as the final clustering result.  Default to be 1. Note that if louvain_iter is large than 1, the `seed` argument will be ignored.  
#' @param weight A logic argument to determine whether or not we will use Jaccard coefficent for two nearest neighbors (based on the overlapping of their kNN) as the weight used for Louvain clustering. Default to be FALSE.
#' @param res Resolution parameter for the louvain clustering. Values between 0 and 1e-2 are good, bigger values give you more clusters. Default is set to be `seq(0, 1e-4, length.out = 5)`. 
#' @param method method for clustering cells. Three methods are available, including densityPeak, louvian and DDRTree. By default, we use density peak clustering algorithm for clustering. For big datasets (like data with 50 k cells or so), we recommend using the louvain clustering algorithm. 
#' @param random_seed  the seed used by the random number generator in louvain-igraph package. This argument will be ignored if louvain_iter is larger than 1.    
#' @param verbose Verbose A logic flag to determine whether or not we should print the running details. 
#' @param cores number of cores computer should use to execute function
#' @param ... Additional arguments passed to \code{\link{densityClust}()}
#' @return an updated CellDataSet object, in which phenoData contains values for Cluster for each cell
#' @importFrom densityClust densityClust findClusters
#' @importFrom igraph graph.data.frame cluster_louvain modularity membership
#' @import ggplot2
#' @importFrom RANN nn2
#' @references Rodriguez, A., & Laio, A. (2014). Clustering by fast search and find of density peaks. Science, 344(6191), 1492-1496. doi:10.1126/science.1242072
#' @references Vincent D. Blondel, Jean-Loup Guillaume, Renaud Lambiotte, Etienne Lefebvre: Fast unfolding of communities in large networks. J. Stat. Mech. (2008) P10008
#' @references Jacob H. Levine and et.al. Data-Driven Phenotypic Dissection of AML Reveals Progenitor-like Cells that Correlate with Prognosis. Cell, 2015. 
#' 
#' @useDynLib monocle
#' 
#' @export

clusterCells <- function(cds, 
                         use_pca = FALSE,
                         skip_rho_sigma = F, 
                         num_clusters = NULL, 
                         inspect_rho_sigma = F, 
                         rho_threshold = NULL, 
                         delta_threshold = NULL, 
                         peaks = NULL,
                         gaussian = T, 
                         clustering_genes=NULL,
                         k = 20, 
                         louvain_iter = 1, 
                         weight = FALSE,
                         res = NULL, 
                         method = c('densityPeak', 'louvain'),
                         random_seed = 0L, 
                         verbose = F, 
                         cores=1,
                         ...) {
  method <- match.arg(method)
  
  if(ncol(cds) > 500000) {
    if(method %in% c('densityPeak', "DDRTree")) {
      warning('Number of cells in your data is larger than 50 k, clusterCells with densityPeak or DDRTree may crash. Please try to use the newly added Louvain clustering algorithm!')
    }
  }
  
  if(method == 'densityPeak'){ 
    ################### DENSITYPEAK START ###################
    set.seed(2017)
    if(use_pca) {
      tsne_data <- t(cds@normalized_data_projection)
    } else {
      tsne_data <- reducedDimA(cds)
    }
    if(nrow(tsne_data) == 0) {
      message('ReduceDimension is not applied to this dataset. We are using the normalized reduced space obtained from preprocessCDS to cluster cells...')
      tsne_data <- t(cds@normalized_data_projection)
    }
    
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
      dataClust <- dataClust[c('rho', 'delta', 'distance', 'dc', 'threshold', 'peaks', 'clusters', 'halo', 'nearest_higher_density_neighbor')]
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
          message(paste('Select top ', num_clusters , 'samples with highest delta as the density peaks and for assigning clusters'))
        }
        
        delta_rho_df <- data.frame("delta" = dataClust$delta, "rho" = dataClust$rho)
        rho_threshold <- 0 
        delta_threshold <- sort(delta_rho_df$delta, decreasing = T)[num_clusters] - .Machine$double.eps
        # rho_valid_threshold <- 0 #quantile(dataClust$rho, probs = 0.05)
        # delta_rho_df <- subset(delta_rho_df, rho > rho_valid_threshold) 
        # threshold_ind <- order(delta_rho_df$delta, decreasing = T)[num_clusters + 1]
        # candidate_peaks <- subset(delta_rho_df, delta >= delta_rho_df$delta[threshold_ind])
        #delta_threshold <- min(candidate_peaks$delta) - .Machine$double.eps
        # delta_threshold <- delta_rho_df$delta[threshold_ind] - 10*.Machine$double.eps 
        # rho_threshold <- min(candidate_peaks$rho) - 10*.Machine$double.eps
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
    
    pData(cds)$Cluster <- factor(dataClust$clusters)
    pData(cds)$peaks <- F
    pData(cds)$peaks[dataClust$peaks] <- T
    pData(cds)$halo <- dataClust$halo
    pData(cds)$delta <- dataClust$delta
    pData(cds)$rho <- dataClust$rho
    pData(cds)$nearest_higher_density_neighbor <- dataClust$nearest_higher_density_neighbor
    
    cds@auxClusteringData[["tSNE"]]$densityPeak <- dataClust[c("dc", "threshold")] #, "peaks"
    
    return(cds)
    
    ################### DENSITYPEAK END ###################
  }  else if(method == 'louvain'){
    if(use_pca) {
      data <- cds@normalized_data_projection
    } else {
      data <- t(reducedDimA(cds))
    }
    if(nrow(data) == 0) {
      message('ReduceDimension is not applied to this dataset. We are using the normalized reduced space obtained from preprocessCDS to cluster cells...')
      data <- cds@normalized_data_projection
      louvain_res <- louvain_clustering(data = data, pd = pData(cds), k = k, weight = weight, louvain_iter = louvain_iter, resolution = res, random_seed = random_seed, verbose = verbose, ...)
    } else {
      if(!('louvain_res' %in% names(cds@auxOrderingData[[cds@dim_reduce_type]]))) {
        louvain_res <- louvain_clustering(data = data, pd = pData(cds), k = k, weight = weight, louvain_iter = louvain_iter, resolution = res, random_seed = random_seed, verbose = verbose, ...)
      } else {
        louvain_res <- cds@auxOrderingData[[cds@dim_reduce_type]]$louvain_res     
      }
    }
    cluster_graph_res <- compute_louvain_connected_components(louvain_res$g, louvain_res$optim_res, verbose = verbose)
    louvain_component <-  components(cluster_graph_res$cluster_g)$membership[louvain_res$optim_res$membership]
    names(louvain_component) <- V(louvain_res$g)$name
    louvain_component <- as.factor(louvain_component)
    pData(cds)$louvain_component <- louvain_component
    
    pData(cds)$Cluster <- factor(igraph::membership(louvain_res$optim_res)) 
    
    cds@auxClusteringData[["louvian"]] <- list(louvain_res = louvain_res)
    
    return(cds)
  }
  else {
    stop('Cluster method ', method, ' is not implemented')
  }
}

#' function to run louvain clustering algorithm 
#'
#' @param data low dimensional space used to perform graph clustering 
#' @param pd the dataframe of the phenotype from the cell dataset (pData(cds))
#' @param k number of nearest neighbors used for Louvain clustering 
#' @param weight whether or not to calculate the weight for each edge in the kNN graph 
#' @param louvain_iter the number of iteraction for louvain clustering 
#' @param resolution resolution of clustering result, specifiying the granularity of clusters.
#' Default to not use resolution and the standard igraph louvain clustering algorithm will be used. 
#' @param random_seed  the seed used by the random number generator in louvain-igraph package  
#' @param verbose Whether to emit verbose output during dimensionality reduction
#' @param ... extra arguments used to run louvain_R
#' @import igraph
#' @return a list with four elements (g (igraph object for the kNN graph), coord (coordinates of the graph with 
#' layout_component, if the number of cells is less than 3000), edge_links (the data frame to plot the edges of 
#' the igraph, if the number of cells is less than 3000) and optim_res (the louvain clustering result)). 
#' 
louvain_clustering <- function(data, pd, k = 20, weight = F, louvain_iter = 1, resolution = NULL, random_seed = 0L, verbose = F, ...) {
  extra_arguments <- list(...)
  cell_names <- row.names(pd)
  if(!identical(cell_names, row.names(pd)))
    stop("phenotype and row name from the data doesn't match")
  
  if (is.data.frame(data))
    data <- as.matrix(data)
  if (!is.matrix(data))
    stop("Wrong input data, should be a data frame of matrix!")
  if (k < 1) {
    stop("k must be a positive integer!")
  } else if (k > nrow(data) - 2) {
    k <- nrow(data) - 2
    warning("RANN counts the point itself, k must be smaller than\nthe total number of points - 1 (all other points) - 1 (itself)!")
  }
  if (verbose) {
    message("Run kNN based graph clustering starts:", "\n", "  -Input data of ",
            nrow(data), " rows and ", ncol(data), " columns",
            "\n", "  -k is set to ", k)
  }
  if (verbose) {
    cat("  Finding nearest neighbors...")
  }
  t1 <- system.time(tmp <- RANN::nn2(data, data, k +
                                       1, searchtype = "standard"))
  neighborMatrix <- tmp[[1]][, -1]
  distMatrix <- tmp[[2]][, -1]
  if (verbose) {
    cat("DONE ~", t1[3], "s\n", " Compute jaccard coefficient between nearest-neighbor sets ...")
  }
  t2 <- system.time(links <- monocle:::jaccard_coeff(neighborMatrix,
                                                     weight))
  if (verbose) {
    cat("DONE ~", t2[3], "s\n", " Build undirected graph from the weighted links ...")
  }
  links <- links[links[, 1] > 0, ]
  relations <- as.data.frame(links)
  colnames(relations) <- c("from", "to", "weight")
  
  relations$from <- cell_names[relations$from]
  relations$to <- cell_names[relations$to]
  t3 <- system.time(g <- igraph::graph.data.frame(relations, directed = FALSE))
  if (verbose) {
    cat("DONE ~", t3[3], "s\n", " Run louvain clustering on the graph ...\n")
  }
  t_start <- Sys.time()
  Qp <- -1
  optim_res <- NULL
  best_max_resolution <- 'No resolution'
  
  if(louvain_iter >= 2) {
    random_seed <- NULL
  }
    
  for (iter in 1:louvain_iter) {
    if(verbose) {
      cat("Running louvain iteration ", iter, "...\n")
    }
    if(!is.null(resolution)) {
      for(i in 1:length(resolution)) {
        cur_resolution <- resolution[i]
        louvain_args <- c(list(X = igraph::get.adjacency(g), res = as.numeric(cur_resolution), random_seed = random_seed, verbose = verbose),
                          extra_arguments[names(extra_arguments) %in% 
                                            c("python_home", "partition_method", "initial_membership", "weights", "node_sizes", 'return_all')])
        Q <- do.call(louvain_R, louvain_args)  
        Qt <- max(Q$modularity)
        if(verbose) {
          message('Current iteration is ', iter, '; current resolution is ', cur_resolution, '; Modularity is ', Qt, '; Number of clusters are ', max(Q$membership))
        }
        if (Qt > Qp) {
          optim_res <- Q
          Qp <- Qt
          best_max_resolution <- cur_resolution
        }
      }
    } else {
      Q <- igraph::cluster_louvain(g)
    }
    if (is.null(optim_res)) {
      Qp <- max(Q$modularity)
      optim_res <- Q
    }
    else {
      Qt <- max(Q$modularity)
      if (Qt > Qp) {
        optim_res <- Q
        Qp <- Qt
      }
    }
  }
  if(verbose) 
    message('Maximal modularity is ', Qp, '; corresponding resolution is ', best_max_resolution)
  t_end <- Sys.time()
  if (verbose) {
    message("\nRun kNN based graph clustering DONE, totally takes ", t_end -
              t_start, " s.")
    cat("  -Number of clusters:", length(unique(igraph::membership(optim_res))), "\n")
  }
  
  if(igraph::vcount(g) < 3000) {
    # coord <- igraph::layout_components(g) 
    # coord <- as.data.frame(coord)
    # colnames(coord) <- c('x', 'y')
    # row.names(coord) <- row.names(pd)
    # coord <- cbind(coord, pd)
    # 
    # edge_links <- cbind(coord[relations$from, 1:2], coord[relations$to, 1:2])
    # edge_links <- as.data.frame(edge_links)
    # colnames(edge_links) <- c('x_start', 'x_end', 'y_start', 'y_end')
    # edge_links$weight <- relations[, 3]
    
    coord <- NULL
    edge_links <- NULL
  } else {
    coord <- NULL
    edge_links <- NULL
  }
  
  V(g)$names <- as.character(V(g))
  return(list(g = g, relations = relations, distMatrix = distMatrix, coord = coord, edge_links = edge_links, optim_res = optim_res))
}

#' @importFrom igraph as_adjacency_matrix
#' @importFrom stats pnorm
compute_louvain_connected_components <- function(g, optim_res, qval_thresh=0.05, verbose = FALSE){
  cell_membership <- as.factor(igraph::membership(optim_res))
  membership_matrix = sparse.model.matrix( ~ cell_membership + 0)
  num_links = t(membership_matrix) %*% as_adjacency_matrix(g) %*% membership_matrix
  diag(num_links) = 0
  louvain_modules = levels(cell_membership)
  
  cluster_mat <- matrix(0, nrow = length(louvain_modules), ncol = length(louvain_modules)) # a matrix storing the overlapping clusters between louvain clusters which is based on the spanning tree
  enrichment_mat <- matrix(0, nrow = length(louvain_modules), ncol = length(louvain_modules)) # a matrix storing the overlapping clusters between louvain clusters which is based on the spanning tree
  
  overlapping_threshold <- 1e-5
  
  edges_per_module = rowSums(num_links)
  total_edges = sum(num_links)
  
  theta <- (as.matrix(edges_per_module) / total_edges) %*% t(edges_per_module / total_edges)
  var_null_num_links <- theta * (1 - theta) / total_edges
  num_links_ij <- num_links / total_edges - theta
  # use mcmapply below (https://stackoverflow.com/questions/7395397/how-to-apply-function-over-each-matrix-elements-indices)
  # tmp <- data.frame(mrow=c(row(num_links)),   # straightens out the arguments
  #          mcol=c(col(num_links)), 
  #          m.f.res= mcmapply(function(r, c) pnorm(num_links_ij[r, c], 0, sqrt(var_null_num_links[r, c]), lower.tail = FALSE), row(num_links), col(num_links), mc.cores = cores  ) )
  # cluster_mat <- as.matrix(dcast(tmp, mrow ~ mcol)[, -1])
  cluster_mat <- pnorm_over_mat(as.matrix(num_links_ij), var_null_num_links) # much faster c++ version 
  
  enrichment_mat <- num_links_ij
  num_links <- num_links_ij / total_edges
  
  cluster_mat = matrix(p.adjust(cluster_mat), nrow=length(louvain_modules), ncol=length(louvain_modules))

  sig_links <- as.matrix(num_links)
  sig_links[cluster_mat > qval_thresh] = 0
  diag(sig_links) = 0
  
  cluster_g <- igraph::graph_from_adjacency_matrix(sig_links, weighted = T, mode = 'undirected')
  louvain_modules <- igraph::cluster_louvain(cluster_g)
  
  # return also the layout coordinates and the edges link for the graph of clusters

  coord <- igraph::layout_components(cluster_g) 
  coord <- as.data.frame(coord)
  colnames(coord) <- c('x', 'y')
  row.names(coord) <- 1:nrow(coord)
  coord$Cluster <- 1:nrow(coord)
  coord$louvain_cluster <- as.character(igraph::membership(louvain_modules))
  
  edge_links <- NULL
  if(length(E(cluster_g)) > 0) { # run this only when there is edges
    edge <- get.data.frame(cluster_g)
    edge <- as.data.frame(edge)
    colnames(edge) <- c('start', 'end', 'weight')
    edge_links <- cbind(coord[edge$start, 1:2], coord[edge$end, 1:2])
    edge_links <- as.data.frame(edge_links)
    colnames(edge_links) <- c('x_start', 'x_end', 'y_start', 'y_end')
    edge_links$weight <- edge[, 3]
  }
  
  list(cluster_g = cluster_g, cluster_optim_res = optim_res, num_links = num_links, cluster_mat = cluster_mat, enrichment_mat = enrichment_mat, cluster_coord = coord, edge_links = edge_links)
}

#' function to run build asymmetric kNN graph 
#'
#' @param data low dimensional space used to perform graph clustering 
#' @param k number of nearest neighbors used for Louvain clustering 
#' @param return_graph whether or not to return the kNN graph instead of the asymmetric adjacency matrix 
#' @return Either a sparse asymmetric adjacent matrix or a corresponding directed weighted kNN graph 
#' 
build_asym_kNN_graph <- function(data, k = 20, dist_type = c('raw', 'euclidean', 'cosine'), return_graph = F) {
  # build an asymmetric kNN graph -- replace with the louvain_clustering one 
  nbrs <- RANN::nn2(data, k = k + 1)
  N <- nrow(data)
  distances <- nbrs$nn.dists[, -1]
  indices =  nbrs$nn.idx[, -1]
  
  rows <- rep(0, N * k)
  cols <- rep(0, N * k)
  dists <- rep(0, N * k)
  location <- 1
  
  for(i in 1:N) {
    inds <- location:(location + k - 1)
    rows[inds] <- i
    cols[inds] <- indices[i, ]
    if(dist_type == 'cosine') {
      if(N < 3000) {
        tmp <- lapply(1:k, function(j) {
          x <- data[i, ]
          y <- data[indices[i, j], ] 
          cosine(x, y) 
        })

        dists[inds] <- unlist(tmp)
      } else {
        message('Cosine distance only supports running on 3K cells at most! Using euclidean distance instead!')
        dists[inds] <- distances[i, ]
      }
    } else if(dist_type == 'euclidean') {
      dists[inds] <- distances[i, ]
    } else if(dist_type == 'raw') {
      dists[inds] <- 1
    }

    location <- location + k
  }
  adj_mat <- sparseMatrix(rows, cols, x = dists^2, dims = c(N, N))
  dimnames(adj_mat) <- list(row.names(data), row.names(data))

  if(return_graph) {
    g <- igraph::graph_from_adjacency_matrix(adj_mat, mode = 'direct', weighted = T)
    return(g)
  } else {
    adj_mat
  }
}

#' function to calculate cosine distance 
cosine <- function (x , y ) {
  cp <- t( x ) %*% y
  normx <- sqrt (t( x ) %*% x )
  normy <- sqrt (t( y ) %*% y )
  cp / normx / normy
}

# cosine <- function ( x ) {
#   x <- t(x)
#   cp <- crossprod ( x )
#   rtdg <- sqrt ( diag ( cp ) )
#   cos <- cp / tcrossprod ( rtdg )
#   return (cos )
# }

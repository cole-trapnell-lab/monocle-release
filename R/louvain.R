#' Cluster cells based on louvain community detection algorithm. 
#' 
#' @description This function is a wrapper of the louvain function from the python package (louvain-igraph, https://github.com/vtraag/louvain-igraph) 
#' The following description is from the original package "This package implements the louvain algorithm in C++ and exposes it to python. It relies on (python-)igraph for it to function. 
#' Besides the relative flexibility of the implementation, it also scales well, and can be run on graphs of millions of nodes (as long as they 
#' can fit in memory). The core function is find_partition which finds the optimal partition using the louvain algorithm [1] for a number of 
#' different methods. The methods currently implemented are (1) modularity [2], (2) Reichardt and Bornholdt's model using the configuration null
#' model and the Erdös-Rényi null model [3], (3) the constant Potts model (CPM) [4], (4) Significance [5], and finally (5) Surprise [6]. In 
#' addition, it supports multiplex partition optimisation allowing community detection on for example negative links [7] or multiple time slices 
#' [8]. It also provides some support for community detection on bipartite graphs. See the documentation for more information." Please see the github 
#' above for the citations. Right now we only support CPMVertexPartition, RBConfigurationVertexPartition, RBERVertexPartition, ModularityVertexPartition
#' SignificanceVertexPartition and SurpriseVertexPartition partition methods. 
#' 
#' @param X the dataset upon which to perform umap dimension reduction
#' @param python_home The python home directory where umap is installed
#' @param partition_method character - either the default "CPMVertexPartition" or  "RBConfigurationVertexPartition" / "RBERVertexPartition". 
#' @param initial_membership (list of int) – Initial membership for the partition. If None then defaults to a singleton partition.
#' @param weights (list of double, or edge attribute) – Weights of edges. Can be either an iterable or an edge attribute.
#' @param res (double) – Resolution parameter.
#' @param node_sizes  (list of int, or vertex attribute) – Sizes of nodes are necessary to know the size of communities in aggregate graphs. Usually this is set to 1 for all nodes, but in specific cases this could be changed.
#' @param random_seed  the seed used by the random number generator in louvain-igraph package  
#' @param verbose bool (optional, default False)
#' @param return_all Whether to return all slots after louvain 
#' @return The cluster id if return_all set to be FALSE, otherwise all slots from the louvain function 
#' @encoding UTF-8
#' @export
#' 
louvain_R <- function(X, python_home = system('which python', intern = TRUE), 
                 partition_method = 'CPMVertexPartition', 
                 initial_membership = NULL, 
                 weights = NULL, 
                 res = 0.6, 
                 node_sizes = NULL, 
                 random_seed = 0L, 
                 verbose = FALSE,
                 return_all = FALSE) {
  
  reticulate::use_python(python_home)
  
  tryCatch({
    reticulate::import("louvain")
  }, warning = function(w) {
  }, error = function(e) {
    print (e)
    stop('please pass the python home directory where louvain is installed with python_home argument!')
  }, finally = {
  })
  
  reticulate::source_python(paste(system.file(package="monocle"), "louvain.py", sep="/"))
  # X <- Matrix::t(X)
  if(length(grep('Matrix', class(X))) == 0){
    X <- as(as.matrix(X), 'TsparseMatrix')
  } else {
    X <- as(X, 'TsparseMatrix')
  }
  
  i <- as.integer(X@i)
  j <- as.integer(X@j)
  val <- X@x
    
  dim <- as.integer(X@Dim)

  if(is.null(partition_method) == F) {
    partition_method <- as.character(partition_method)
  }
  if(!is.null(random_seed)) {
    random_seed <- as.integer(random_seed)
  }
  # if(is.null(initial_membership) == F) { #initial_membership (list of int) 
  #   a <- as.numeric(a)
  # }
  # if(is.null(weights) == F) { # weights (list of double, or edge attribute) 
  #   n_epochs <- as.numeric(b)
  # }
  # if(is.null(res) == F) { # node_sizes (list of int, or vertex attribute)
  #   metric_kwds <- reticulate::dict()
  # } 
  # if(is.null(node_sizes) == F) { # resolution_parameter (double)
  #   metric_kwds <- reticulate::dict(metric_kwds)
  # }
  
  louvain_res <- louvain(i, j, val, dim, 
                   as.character(partition_method), 
                   initial_membership, 
                   weights, 
                   as.numeric(res),
                   node_sizes,
                   random_seed, 
                   as.logical(verbose))
  
  if(return_all) {
    return(louvain_res)
  } else {
    list(membership = louvain_res$membership + 1, modularity = louvain_res$modularity)
  }
}

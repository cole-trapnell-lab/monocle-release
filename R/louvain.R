#' Uniform Manifold Approximation and Projection
#' 
#' @description Finds a low dimensional embedding of the data that approximates an underlying manifold.
#' This functions relies on the python implementation of UMAP (https://github.com/lmcinnes/umap). 
#' The original publication of UMAP can be found here: 
#' https://arxiv.org/abs/1802.03426 (McInnes, L, Healy, J, UMAP: Uniform Manifold Approximation and Projection for Dimension Reduction, ArXiv e-prints 1802.03426, 2018).
#' A useful notebook (written in python) to check for the effects of each parameter in UMAP can be found here: https://nbviewer.jupyter.org/github/CrakeNotSnowman/umapNotebooks/blob/master/UMAP%20Usage.ipynb.
#' 
#' @param X the dataset upon which to perform umap dimension reduction
#' @param python_home The python home directory where umap is installed
#' @param partition_method character - either the default "CPMVertexPartition" or  "RBConfigurationVertexPartition" / "RBERVertexPartition". 
#' @param initial_membership (list of int) – Initial membership for the partition. If None then defaults to a singleton partition.
#' @param weights (list of double, or edge attribute) – Weights of edges. Can be either an iterable or an edge attribute.
#' @param res (double) – Resolution parameter.
#' @param node_sizes  (list of int, or vertex attribute) – Sizes of nodes are necessary to know the size of communities in aggregate graphs. Usually this is set to 1 for all nodes, but in specific cases this could be changed.
#' @param verbose bool (optional, default False)
#' @param return_all Whether to return all slots after louvain 
#' Controls verbosity of logging.
#' @return The cluster id if return_all set to be FALSE, otherwise all slots from the louvain function 
#' @import reticulate
#' @export
#' 
louvain_R <- function(X, python_home = system('which python', intern = TRUE), 
                 partition_method = 'CPMVertexPartition', 
                 initial_membership = NULL, 
                 weights = NULL, 
                 res = 0.6, 
                 node_sizes = NULL, 
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
                   as.logical(verbose))
  
  if(return_all) {
    return(louvain_res)
  } else {
    list(membership = louvain_res$membership + 1, modularity = louvain_res$modularity)
  }
}

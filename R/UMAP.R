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
#' @param log A logic argument to determine whether we need to calculate 
#' log of the input data. Default to be true 
#' @param n_neighbors float (optional, default 15)
#' The size of local neighborhood (in terms of number of neighboring
#' sample points) used for manifold approximation. Larger values
#' result in more global views of the manifold, while smaller
#' values result in more local data being preserved. In general
#' values should be in the range 2 to 100.
#' @param n_component int (optional, default 2)
#' The dimension of the space to embed into. This defaults to 2 to
#' provide easy visualization, but can reasonably be set to any
#' integer value in the range 2 to 100.
#' @param metric string or function (optional, default 'correlation')
#' The metric to use to compute distances in high dimensional space.
#' If a string is passed it must match a valid predefined metric. If
#' a general metric is required a function that takes two 1d arrays and
#' returns a float can be provided. For performance purposes it is
#' required that this be a numba jit'd function. Valid string metrics
#' include:
#'     * euclidean
#'     * manhattan
#'     * chebyshev
#'     * minkowski
#'     * canberra
#'     * braycurtis
#'     * mahalanobis
#'     * wminkowski
#'     * seuclidean
#'     * cosine
#'     * correlation
#'     * haversine
#'     * hamming
#'     * jaccard
#'     * dice
#'     * russelrao
#'     * kulsinski
#'     * rogerstanimoto
#'     * sokalmichener
#'     * sokalsneath
#'     * yule
#'   Metrics that take arguments (such as minkowski, mahalanobis etc.)
#'   can have arguments passed via the metric_kwds dictionary. At this
#'   time care must be taken and dictionary elements must be ordered
#'   appropriately; this will hopefully be fixed in the future.
#' @param n_epochs int The number of training epochs to use in optimization.
#' @param negative_sample_rate int (optional, default 5)
#' The number of negative edge/1-simplex samples to use per positive 
#' edge/1-simplex sample in optimizing the low dimensional embedding. 
#' @param alpha float (optional, default 1.0)
#' The initial learning rate for the embedding optimization.
#' @param init string (optional, default 'spectral')
#' How to initialize the low dimensional embedding. Options are:
#'     * 'spectral': use a spectral embedding of the fuzzy 1-skeleton
#'     * 'random': assign initial embedding positions at random.
#'     * A numpy array of initial embedding positions.
#' @param min_dist float (optional, default 0.1)
#' The effective minimum distance between embedded points. Smaller values
#' will result in a more clustered/clumped embedding where nearby points
#' on the manifold are drawn closer together, while larger values will
#' result on a more even dispersal of points. The value should be set
#' relative to the ``spread`` value, which determines the scale at which
#' embedded points will be spread out.
#' @param spread float (optional, default 1.0)
#' The effective scale of embedded points. In combination with ``min_dist``
#' this determines how clustered/clumped the embedded points are.
#' @param set_op_mix_ratio float (optional, default 1.0)
#' Interpolate between (fuzzy) union and intersection as the set operation
#' used to combine local fuzzy simplicial sets to obtain a global fuzzy
#' simplicial sets. Both fuzzy set operations use the product t-norm.
#' The value of this parameter should be between 0.0 and 1.0; a value of
#' 1.0 will use a pure fuzzy union, while 0.0 will use a pure fuzzy
#' intersection.
#' @param local_connectivity int (optional, default 1)
#' The local connectivity required -- i.e. the number of nearest
#' neighbors that should be assumed to be connected at a local level.
#' The higher this value the more connected the manifold becomes
#' locally. In practice this should be not more than the local intrinsic
#' dimension of the manifold.
#' @param gamma float (optional, default 1.0)
#' Weighting applied to negative samples in low dimensional embedding
#' optimization. Values higher than one will result in greater weight
#' being given to negative samples.
#' @param bandwidth float (optional, default 1.0)
#' The effective bandwidth of the kernel if we view the algorithm as
#' similar to Laplacian eigenmaps. Larger values induce more
#' connectivity and a more global view of the data, smaller values
#' concentrate more locally.
#' @param a float (optional, default None) (not passed in)
#' More specific parameters controlling the embedding. If None these
#' values are set automatically as determined by ``min_dist`` and
#' ``spread``.
#' @param b float (optional, default None) (not passed in)
#' More specific parameters controlling the embedding. If None these
#' values are set automatically as determined by ``min_dist`` and
#' ``spread``.
#' @param random_state int, RandomState instance or None, optional (default: None) (not passed in)
#' If int, random_state is the seed used by the random number generator;
#' If RandomState instance, random_state is the random number generator;
#' If None, the random number generator is the RandomState instance used
#' by `np.random`.
#' @param metric_kwds dict (optional, default {})  (not passed in)
#' Arguments to pass on to the metric, such as the ``p`` value for
#' Minkowski distance. In R, a list should be passed in if you want to use this argument. 
#' The dict function from reticulate package will then convert it into a dictionary for python to use. 
#' @param angular_rp_forest bool (optional, default False)
#' Whether to use an angular random projection forest to initialise
#' the approximate nearest neighbor search. This can be faster, but is
#' mostly on useful for metric that use an angular style distance such
#' as cosine, correlation etc. In the case of those metrics angular forests
#' will be chosen automatically.
#' @param verbose bool (optional, default False)
#' Controls verbosity of logging.
#' @param return_all Whether to return all slots after UMAP 
#' @return Embedding of the training data in low-dimensional space if return_all is set to be FALSE, 
#' otherwise the object returned from umap function, including the following elements: 
#' a, fit_transform, metric, random_state, alpha, gamma, metric_kwds, set_op_mix_ratio, angular_rp_forest, 
#' get_params, min_dist, set_params, b, graph, n_components, spread, bandwidth, init, n_epochs, verbose, 
#' embedding_, initial_alpha, n_neighbors, fit, local_connectivity, negative_sample_rate 
#' @export
UMAP <- function(X, python_home = system('which python', intern = TRUE), 
  log = TRUE, 
  n_neighbors = 15L, 
  n_component = 2L, 
  metric = "correlation", 
  n_epochs = NULL, 
  negative_sample_rate = 5L,
  alpha = 1.0,
  init = 'spectral',
  min_dist = 0.1, 
  spread = 1.0,
  set_op_mix_ratio = 1.0,
  local_connectivity = 1L,
  bandwidth = 1.0, 
  gamma = 1.0,
  a = NULL,
  b = NULL, 
  random_state = 0L,
  metric_kwds = reticulate::dict(), 
  angular_rp_forest = FALSE,
  verbose = FALSE,
  return_all = FALSE) {
  
  reticulate::use_python(python_home)
  
  tryCatch({
    reticulate::import("umap")
  }, warning = function(w) {
  }, error = function(e) {
    print (e)
    stop('please pass the python home directory where umap is installed with python_home argument!')
  }, finally = {
  })
  
  reticulate::source_python(paste(system.file(package="monocle"), "umap.py", sep="/"))
  # X <- Matrix::t(X)
  if(length(grep('Matrix', class(X))) == 0){
    X <- as(as.matrix(X), 'TsparseMatrix')
  } else {
    X <- as(X, 'TsparseMatrix')
  }
  
  i <- as.integer(X@i)
  j <- as.integer(X@j)
  
  if(log) {
  val <- log(X@x + 1)
  } else {
    val <- X@x
  }
  dim <- as.integer(X@Dim)

  if(is.null(n_epochs) == F) {
    n_epochs <- as.integer(n_epochs)
  }
  if(is.null(a) == F) {
    a <- as.numeric(a)
  }
  if(is.null(b) == F) {
    n_epochs <- as.numeric(b)
  }
  if(is.list(metric_kwds) == F) {
    metric_kwds <- reticulate::dict()
  } else {
    metric_kwds <- reticulate::dict(metric_kwds)
  }

  umap_res <- umap(i, j, val, dim, 
                    as.integer(n_neighbors), 
                    as.integer(n_component), 
                    as.character(metric), 
                    n_epochs,
                    as.integer(negative_sample_rate),
                    as.numeric(alpha),
                    as.character(init),
                    as.numeric(min_dist), 
                    as.numeric(spread),
                    as.numeric(set_op_mix_ratio),
                    as.integer(local_connectivity),
                    as.numeric(bandwidth),
                    as.numeric(gamma),
                    a,
                    b,
                    as.integer(random_state),
                    metric_kwds,
                    as.logical(angular_rp_forest),
                    as.logical(verbose))
  
  if(return_all) {
    return(umap_res)
  } else {
    umap_res$embedding_
  }
}

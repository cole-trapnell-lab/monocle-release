# run_dpt <- function(data, branching = T, norm_method = 'log', verbose = F){
#   data <- t(data) 
#   data <- data[!duplicated(data), ]
#   dm <- DiffusionMap(as.matrix(data))
#   return(dm@eigenvectors)
# }
# run_dpt <- function(data, branching = T, norm_method = 'log', root = NULL, verbose = F){
#   if (!requireNamespace("destiny", quietly = TRUE)) {
#     stop("destiny package needed for this function to work. Please install it.",
#          call. = FALSE)
#   }
#   
#   if(verbose)
#     message('root should be the id to the cell not the cell name ....')
#   
#   data <- t(data)
#   data <- data[!duplicated(data), ]
#   dm <- DiffusionMap(as.matrix(data))
#   dpt <- DPT(dm, branching = branching)
# 
#  ts <- dm@transitions
#  M <- destiny::accumulated_transitions(dm)
#
#  branch <- dpt@branch
#  row.names(branch) <- row.names(data[!duplicated(data), ])
#
#  if(is.null(root))
#    root <- random_root(dm)[1]
#  pt <- dpt[root, ]
#  dp_res <- list(dm = dm, pt = pt, ts = ts, M = M, ev = dm@eigenvectors, branch = branch)
# 
#   return(dm@eigenvectors)
# }

#' Marks genes for clustering
#' @description The function marks genes that will be used for clustering in subsequent calls to clusterCells. 
#' The list of selected genes can be altered at any time.
#' 
#' @param cds the CellDataSet upon which to perform this operation
#' @param ordering_genes a vector of feature ids (from the CellDataSet's featureData) used for ordering cells
#' @return an updated CellDataSet object
#' @export
setOrderingFilter <- function(cds, ordering_genes){
  fData(cds)$use_for_ordering <- row.names(fData(cds)) %in% ordering_genes
  cds
}

#' @importFrom igraph V distances
extract_general_graph_ordering <- function(cds, root_cell, verbose=T)
{
  pr_graph <- minSpanningTree(cds)
  
  res <- list(subtree = pr_graph, root = root_cell)

  parents = rep(NA, length(V(pr_graph)))
  states = rep(NA, length(V(pr_graph)))
  
  pr_graph_node_distances = distances(pr_graph, v=root_cell)
  if (length(root_cell) > 1){
    node_names = colnames(pr_graph_node_distances)
    pseudotimes = apply(pr_graph_node_distances, 2, min)
  }else{
    node_names = names(pr_graph_node_distances)
    pseudotimes = pr_graph_node_distances
  }
  
  
  names(pseudotimes) <- node_names
  
  ordering_df <- data.frame(sample_name = V(pr_graph)$name,
                            cell_state = states,
                            pseudo_time = as.vector(pseudotimes),
                            parent = parents)
  row.names(ordering_df) <- ordering_df$sample_name
  return(ordering_df)
}


#' Orders cells according to pseudotime.
#'
#' Learns a "trajectory" describing the biological process the cells are
#' going through, and calculates where each cell falls within that trajectory.
#' Monocle learns trajectories in two steps. The first step is reducing the dimensionality
#' of the data with \code{\link{reduceDimension}()}. The second is this function.
#' function. This function takes as input a CellDataSet and returns it with
#' two new columns: \code{Pseudotime} and \code{State}, which together encode
#' where each cell maps to the trajectory. \code{orderCells()} optionally takes
#' a "root" state, which you can use to specify the start of the trajectory. If
#' you don't provide a root state, one is selected arbitrarily.
#'
#' The \code{reduction_method} argument to \code{\link{reduceDimension}()}
#' determines which algorithm is used by \code{orderCells()} to learn the trajectory.
#' If \code{reduction_method == "ICA"}, this function uses \emph{polygonal reconstruction}
#' to learn the underlying trajectory. If \code{reduction_method == "DDRTree"},
#' the trajectory is specified by the principal graph learned by the
#' \code{\link[DDRTree]{DDRTree}()} function.
#'
#' Whichever algorithm you use, the trajectory will be composed of segments.
#' The cells from a segment will share the same value of \code{State}. One of
#' these segments will be selected as the root of the trajectory arbitrarily.
#' The most distal cell on that segment will be chosen as the "first" cell in the
#' trajectory, and will have a Pseudotime value of zero. \code{orderCells()} will
#' then "walk" along the trajectory, and as it encounters additional cells, it
#' will assign them increasingly large values of Pseudotime.
#'
#' @param cds the CellDataSet upon which to perform this operation
#' @param root_state The state to use as the root of the trajectory.
#' You must already have called orderCells() once to use this argument.
#' @param num_paths the number of end-point cell states to allow in the biological process.
#' @param reverse whether to reverse the beginning and end points of the learned biological process.
#' 
#' @importFrom stats dist
#' @importFrom igraph graph.adjacency V as.undirected
#'
#' @return an updated CellDataSet object, in which phenoData contains values for State and Pseudotime for each cell
#' @export
orderCells <- function(cds,
                       root_pr_nodes=NULL,
                       root_cells=NULL){
  
  if(class(cds)[1] != "CellDataSet") {
    stop("Error cds is not of type 'CellDataSet'")
  }
  
  if (is.null(cds@dim_reduce_type)){
    stop("Error: dimensionality not yet reduced. Please call reduceDimension() before calling this function.")
  }
  # reducedDimA, S, and K are not NULL in the cds
  if (any(c(length(cds@reducedDimS) == 0, length(cds@reducedDimK) == 0))) {
    stop("Error: dimension reduction didn't prodvide correct results. Please check your reduceDimension() step and ensure correct dimension reduction are performed before calling this function.")
  }
  
  if(igraph::vcount(minSpanningTree(cds)) > 50000) {
    stop("orderCells doesn't support more than 50k centroids (cells)")
  }
  if (is.null(root_pr_nodes) == FALSE & is.null(root_cells) == FALSE){
    stop("Error: please specify either root_pr_nodes or root_cells, not both")
  }
  if (is.null(root_pr_nodes) & is.null(root_cells)){
    if (interactive()){
      root_pr_nodes = selectTrajectoryRoots(cds)
    }else{
      stop("Error: You must provide one or more root cells (or principal graph nodes) in non-interactive mode")
    }
  }else if(is.null(root_pr_nodes)){
    valid_root_cells = intersect(root_cells, row.names(pData(cds)))
    if (length(valid_root_cells) == 0){
      stop("Error: no such cell")
    }
    closest_vertex = cds@auxOrderingData[[cds@dim_reduce_type]]$pr_graph_cell_proj_closest_vertex
    root_pr_nodes = closest_vertex[valid_root_cells,]
  }else{
    if (length(intersect(root_pr_nodes, V(minSpanningTree(cds))$name)) == 0){
      stop("Error: no such principal graph node")
    }
  }
  
  if (is.null(root_pr_nodes) || length(root_pr_nodes) == 0){
    stop("Error: no valid root principal graph nodes.")
  }
  
  cds@auxOrderingData[[cds@dim_reduce_type]]$root_pr_nodes <- root_pr_nodes
  
  if (cds@dim_reduce_type == "L1graph"){
    cc_ordering <- extract_general_graph_ordering(cds, root_pr_nodes)
    closest_vertex = cds@auxOrderingData$L1graph$pr_graph_cell_proj_closest_vertex
   
    pData(cds)$Pseudotime = cc_ordering[closest_vertex[row.names(pData(cds)),],]$pseudo_time
    pData(cds)$State = cc_ordering[closest_vertex[row.names(pData(cds)),],]$state
    
    cds@auxOrderingData[[cds@dim_reduce_type]]$root_pr_nodes <- root_pr_nodes
    
    mst_branch_nodes <- NULL
  }else if (cds@dim_reduce_type == "DDRTree"){
    
    cc_ordering <- extract_mst_ordering(cds, root_pr_nodes)
    
    R <- cds@auxOrderingData$DDRTree$R
    edge <- data.frame(start = 1:nrow(R), end = apply(R, 1, which.max), weight = R[cbind(1:nrow(R), apply(R, 1, which.max))])
    
    pData(cds)$Pseudotime <- cc_ordering[edge$end, 'pseudo_time']
    if(is.null(root_state) == TRUE) {
      pData(cds)$State <- cc_ordering[edge$end, 'cell_state']
    }
    
    mst_branch_nodes <- V(minSpanningTree(cds))[which(degree(minSpanningTree(cds)) > 2)]$name
  } else if (cds@dim_reduce_type == "SimplePPT"){
    cc_ordering <- extract_mst_ordering(cds, root_pr_nodes)
    
    pData(cds)$Pseudotime <-  cc_ordering[row.names(pData(cds)),]$pseudo_time
    pData(cds)$State <- cc_ordering[row.names(pData(cds)),]$cell_state
    
    mst_branch_nodes <- V(minSpanningTree(cds))[which(degree(minSpanningTree(cds)) > 2)]$name
  } else if (cds@dim_reduce_type %in% c("UMAP", "UMAPSSE", "SSE")){
    
    ########################################################################################################################################################################
    # downstream pseudotime and branch analysis 
    ########################################################################################################################################################################
    pc_g <- minSpanningTree(cds)
    pData(cds)$Pseudotime <- as.vector(distances(pc_g, v = root_cell, to = as.character(1:igraph::vcount(pc_g))))
    # identify branch
    mst_branch_nodes <- V(minSpanningTree(cds))[which(degree(minSpanningTree(cds)) > 2)]$name
  }
  
  cds@auxOrderingData[[cds@dim_reduce_type]]$branch_points <- mst_branch_nodes
  # FIXME: the scaling code is totally broken after moving to DDRTree. Disabled
  # for now
  #if(scale_pseudotime) {
  #cds <- scale_pseudotime(cds)
  #}
  
  cds
}

# Helper function to normalize the expression data prior to dimensionality
# reduction
normalize_expr_data <- function(cds,
                                norm_method = c("log", "vstExprs", "none"),
                                pseudo_expr = 1,
                                relative_expr = TRUE){
  FM <- exprs(cds)
  use_for_ordering <- NULL
  # If the user has selected a subset of genes for use in ordering the cells
  # via setOrderingFilter(), subset the expression matrix.
  if (is.null(fData(cds)$use_for_ordering) == FALSE &&
      nrow(subset(fData(cds), use_for_ordering == TRUE)) > 0) {
    FM <- FM[fData(cds)$use_for_ordering, ]
  }

  norm_method <- match.arg(norm_method)
  if (cds@expressionFamily@vfamily %in% c("negbinomial", "negbinomial.size")) {

    # If we're going to be using log, and the user hasn't given us a pseudocount
    # set it to 1 by default.
    if (is.null(pseudo_expr)){
      if(norm_method == "log")
        pseudo_expr = 1
      else
        pseudo_expr = 0
    }

    checkSizeFactors(cds)

    if (norm_method == "vstExprs") {
      if (relative_expr == FALSE)
        message("Warning: relative_expr is ignored when using norm_method == 'vstExprs'")

      if (is.null(fData(cds)$use_for_ordering) == FALSE &&
          nrow(subset(fData(cds), use_for_ordering == TRUE)) > 0) {
        VST_FM <- vstExprs(cds[fData(cds)$use_for_ordering,], round_vals = FALSE)
      }else{
        VST_FM <- vstExprs(cds, round_vals = FALSE)
      }

      if (is.null(VST_FM) == FALSE) {
        FM <- VST_FM
      }
      else {
        stop("Error: set the variance-stabilized value matrix with vstExprs(cds) <- computeVarianceStabilizedValues() before calling this function with use_vst=TRUE")
      }
    }else if (norm_method == "log") {
      # If we are using log, normalize by size factor before log-transforming
      
      if (relative_expr)
        FM <- Matrix::t(Matrix::t(FM)/sizeFactors(cds))

      if(is.null(pseudo_expr))
        pseudo_expr <- 1
      if (pseudo_expr != 1 || isSparseMatrix(exprs(cds)) == FALSE){
        FM <- FM + pseudo_expr
        FM <- log2(FM)
      }else{
        FM@x = log2(FM@x + 1)
      }
      

    }else if (norm_method == "none"){
      # If we are using log, normalize by size factor before log-transforming
      FM <- Matrix::t(Matrix::t(FM)/sizeFactors(cds))
      FM <- FM + pseudo_expr
    }
  }else if (cds@expressionFamily@vfamily == "binomialff") {
    if (norm_method == "none"){
      #If this is binomial data, transform expression values into TF-IDF scores.
      ncounts <- FM > 0
      ncounts[ncounts != 0] <- 1
      FM <- Matrix::t(Matrix::t(ncounts) * log(1 + ncol(ncounts)/rowSums(ncounts)))
    }else{
      stop("Error: the only normalization method supported with binomial data is 'none'")
    }
  }else if (cds@expressionFamily@vfamily == "Tobit") {
    FM <- FM + pseudo_expr
    if (norm_method == "none"){

    }else if (norm_method == "log"){
      FM <- log2(FM)
    }else{
      stop("Error: the only normalization methods supported with Tobit-distributed (e.g. FPKM/TPM) data are 'log' (recommended) or 'none'")
    }
  }else if (cds@expressionFamily@vfamily == "gaussianff") {
    if (norm_method == "none"){
      FM <- FM + pseudo_expr
    }else{
      stop("Error: the only normalization method supported with gaussian data is 'none'")
    }
  }
  # if(norm_method != "none")
    #normalize_expr_data
  return (FM)
}

#' project a CellDataSet object into a lower dimensional PCA space
#'
#' @description For most analysis (including trajectory inference, clustering) in Monocle 3, it requires us to to start from a 
#' low dimensional PCA space. projectPCA will be used to first project a CellDataSet object into a lower dimensional PCA space 
#' before we apply clustering with community detection algorithm or other non-linear dimension reduction method, for example 
#' UMAP, tSNE, DDRTree, SSE, L1-graph, SGL-tree, etc.  
#' While tSNE is especially suitable for visualizing clustering result, comparing to UMAP, the global distance in tSNE space is 
#' not meaningful. UMAP can either be used for visualizing clustering result or as a general non-linear dimension reduction method. 
#' It can be used in conjunction with SSE to obtain smooth skeleton representation of the data. DDRTree and L1-graph are two complementary 
#' trajectory inference method where the first one is very great at learning a tree structure but the later is general and can 
#' learn any arbitrary graph structure. Both methods can be applied to the UMAP space or the smoothed SSE space.   
#'
#' @param cds the CellDataSet upon which to perform this operation
#' @param num_dim the dimensionality of the reduced space
#' @param norm_method Determines how to transform expression values prior to reducing dimensionality
#' @param residualModelFormulaStr A model formula specifying the effects to subtract from the data before clustering.
#' @param pseudo_expr amount to increase expression values before dimensionality reduction
#' @param relative_expr When this argument is set to TRUE (default), we intend to convert the expression into a relative expression.
#' @param scaling When this argument is set to TRUE (default), it will scale each gene before running trajectory reconstruction.
#' @param verbose Whether to emit verbose output during dimensionality reduction
#' @param ... additional arguments to pass to the dimensionality reduction function
#' @return an updated CellDataSet object
#' @import methods
#' @importFrom matrixStats rowSds
#' @importFrom limma removeBatchEffect
#' @importFrom fastICA  ica.R.def ica.R.par
#' @import irlba
#' @importFrom stats dist prcomp
#' @export
projectPCA <- function(cds, num_dim=50,
                        norm_method = c("log", "vstExprs", "none"),
                        residualModelFormulaStr=NULL,
                        pseudo_expr=1,
                        relative_expr=TRUE,
                        auto_param_selection = TRUE,
                        verbose=FALSE,
                        scaling = TRUE,
                        ...) {
  extra_arguments <- list(...)
  set.seed(2016) #ensure results from RNG sensitive algorithms are the same on all calls
  
  FM <- normalize_expr_data(cds, norm_method, pseudo_expr, relative_expr)

  # For NB: Var(Y)=mu*(1+mu/k)
  #f_expression_var <- DelayedMatrixStats::rowVars(FM)
  #FM <- FM[f_expression_var > 0,]

  if (is.null(residualModelFormulaStr) == FALSE) {
    if (verbose)
      message("Removing batch effects")
    X.model_mat <- sparse.model.matrix(as.formula(residualModelFormulaStr),
                                       data = pData(cds), drop.unused.levels = TRUE)

    fit <- limma::lmFit(FM, X.model_mat, ...)
    beta <- fit$coefficients[, -1, drop = FALSE]
    beta[is.na(beta)] <- 0
    FM <- as.matrix(FM) - beta %*% t(X.model_mat[, -1])
  }else{
    X.model_mat <- NULL
  }

  if (nrow(FM) == 0) {
    stop("Error: all rows have standard deviation zero")
  }
  
  fm_rowsums = Matrix::rowSums(FM)
  FM <- FM[is.finite(fm_rowsums) & fm_rowsums != 0, ]

  if (verbose)
    message("Remove noise by PCA ...")
  
  irlba_res <- sparse_prcomp_irlba(t(FM), n = min(num_dim, min(dim(FM)) - 1),
                            center = scaling, scale. = scaling)
  irlba_pca_res <- irlba_res$x
  reducedDimA(cds) <- t(irlba_pca_res) # get top 50 PCs, which can be used for louvain clustering later 
  cds@auxOrderingData[["PCA"]]$irlba_pca_res <- irlba_pca_res

  cds
}

#' project a CellDataSet object into a lower dimensional PCA space after normalize the data 
#'
#' @description For most analysis (including trajectory inference, clustering) in Monocle 3, it requires us to to start from a 
#' low dimensional PCA space. preprocessCDS will be used to first project a CellDataSet object into a lower dimensional PCA space 
#' before we apply clustering with community detection algorithm or other non-linear dimension reduction method, for example 
#' UMAP, tSNE, DDRTree, SSE, L1-graph, SGL-tree, etc.  While tSNE is especially suitable for visualizing clustering results, comparing
#' to UMAP, the global distance in tSNE space is not meaningful. UMAP can either be used for visualizing clustering result or as a general 
#' non-linear dimension reduction method. It can be used in conjunction with SSE to obtain smooth skeleton representation of the data. 
#' SimplePPT, DDRTree and L1-graph are two complementary trajectory inference method where the first one is very great at learning a tree structure 
#' but the later is general and can learn any arbitrary graph structure. Both methods can be applied to the UMAP space or the smoothed SSE space.   
#'
#' @details 
#' In Monocle 3, we overhauled the code from Monocle2 so that a standard Monocle 3 workingflow works as following: 
#' 1. run \code{preprocessCDS} to project a CellDataSet object into a lower dimensional PCA space after 
#' normalize the data 
#' 2. run \code{reduceDimension} to further project the PCA space into much lower dimension space with non-linear 
#' dimension reduction techniques, including tSNE, UMAP. 
#' 3. run \code{smoothEmbedding} (optional) to smooth noisy embedding from 2 to facilitate visualization and learning 
#' of the graph structure.
#' 4. run \code{learnGraph} to reconstruct developmental trajectory with reversed graph embedding algorithms. In monocle 3, we enabled the 
#' the capability to learn multiple disjointed trajectory with either tree or loop structure, etc. 
#'
#' Prior to reducing the dimensionality of the data, it usually helps
#' to normalize it so that highly expressed or highly variable genes don't
#' dominate the computation. \code{reduceDimension()} automatically transforms
#' the data in one of several ways depending on the \code{expressionFamily} of
#' the CellDataSet object. If the expressionFamily is \code{negbinomial} or \code{negbinomial.size}, the
#' data are variance-stabilized. If the expressionFamily is \code{Tobit}, the data
#' are adjusted by adding a pseudocount (of 1 by default) and then log-transformed.
#' If you don't want any transformation at all, set norm_method to "none" and
#' pseudo_expr to 0. This maybe useful for single-cell qPCR data, or data you've
#' already transformed yourself in some way.

#' @param cds the CellDataSet upon which to perform this operation
#' @param method the initial dimension method to use, current either PCA or LSI. For LSI (latent semantic indexing), 
#' it converts the (sparse) expression matrix into tf-idf (term-frequency-inverse document frequency) matrix and then performs a 
#' SVD to decompose the gene expression / cells into certain modules / topics. This method can be used to find associated gene modules 
#  and cell clusters at the same time. It removes noise in the data and thus makes the UMAP result even better. 
#' @param use_tf_idf a logic argument to determine whether we should convert the normalized gene expression value into tf-idf value before performing PCA 
#' @param num_dim the dimensionality of the reduced space
#' @param norm_method Determines how to transform expression values prior to reducing dimensionality
#' @param residualModelFormulaStr A model formula specifying the effects to subtract from the data before clustering.
#' @param pseudo_expr amount to increase expression values before dimensionality reduction
#' @param relative_expr When this argument is set to TRUE (default), we intend to convert the expression into a relative expression.
#' @param scaling When this argument is set to TRUE (default), it will scale each gene before running trajectory reconstruction.
#' @param verbose Whether to emit verbose output during dimensionality reduction
#' @param ... additional arguments to pass to the dimensionality reduction function
#' @return an updated CellDataSet object
#' @import methods
#' @importFrom matrixStats rowSds
#' @importFrom limma removeBatchEffect
#' @importFrom fastICA  ica.R.def ica.R.par
#' @import irlba
#' @importFrom stats dist prcomp
#' @export
preprocessCDS <- function(cds, method = c('PCA', 'none'), #, 'LSI' , 'NMF'
                          use_tf_idf = FALSE, 
                          num_dim=50,
                          norm_method = c("log", "vstExprs", "none"),
                          residualModelFormulaStr=NULL,
                          pseudo_expr=1,
                          relative_expr=TRUE,
                          scaling = TRUE,
                          verbose=FALSE,
                          ...) {
  extra_arguments <- list(...)
  set.seed(2016) #ensure results from RNG sensitive algorithms are the same on all calls
  
  FM <- normalize_expr_data(cds, norm_method, pseudo_expr, relative_expr)
  
  # For NB: Var(Y)=mu*(1+mu/k)
  #f_expression_var <- DelayedMatrixStats::rowVars(FM)
  #FM <- FM[f_expression_var > 0,]
  
  if (is.null(residualModelFormulaStr) == FALSE) {
    if (verbose)
      message("Removing batch effects")
    X.model_mat <- sparse.model.matrix(as.formula(residualModelFormulaStr),
                                       data = pData(cds), drop.unused.levels = TRUE)
    
    fit <- limma::lmFit(FM, X.model_mat, ...)
    beta <- fit$coefficients[, -1, drop = FALSE]
    beta[is.na(beta)] <- 0
    FM <- as.matrix(FM) - beta %*% t(X.model_mat[, -1])
  }else{
    X.model_mat <- NULL
  }
  
  if (nrow(FM) == 0) {
    stop("Error: all rows have standard deviation zero")
  }
  
  fm_rowsums = Matrix::rowSums(FM)
  FM <- FM[is.finite(fm_rowsums) & fm_rowsums != 0, ]
  
  cds@auxOrderingData$normalize_expr_data <- FM
  
  if(method == 'PCA') {
    if (verbose)
      message("Remove noise by PCA ...")

    irlba_res <- sparse_prcomp_irlba(t(FM), n = min(num_dim, min(dim(FM)) - 1),
                                     center = scaling, scale. = scaling)
    irlba_pca_res <- irlba_res$x
    reducedDimA(cds) <- t(irlba_pca_res) # get top 50 PCs, which can be used for louvain clustering later 
  } 
  
  cds@normalized_data_projection <- irlba_pca_res
  
  cds
}

#' Compute a projection of a CellDataSet object into a lower dimensional space
#' 
#' @description Monocle aims to learn how cells transition through a biological program of 
#' gene expression changes in an experiment. Each cell can be viewed as a point 
#' in a high-dimensional space, where each dimension describes the expression of 
#' a different gene in the genome. Identifying the program of gene expression 
#' changes is equivalent to learning a \emph{trajectory} that the cells follow
#' through this space. However, the more dimensions there are in the analysis,
#' the harder the trajectory is to learn. Fortunately, many genes typically
#' co-vary with one another, and so the dimensionality of the data can be
#' reduced with a wide variety of different algorithms. Monocle provides two
#' different algorithms for dimensionality reduction via \code{reduceDimension}.
#' Both take a CellDataSet object and a number of dimensions allowed for the
#' reduced space. You can also provide a model formula indicating some variables
#' (e.g. batch ID or other technical factors) to "subtract" from the data so it
#' doesn't contribute to the trajectory. 
#' 
#' @details You can choose a few different reduction algorithms: Independent Component 
#' Analysis (ICA) and Discriminative Dimensionality Reduction with Trees (DDRTree).
#' The choice impacts numerous downstream analysis steps, including \code{\link{orderCells}}.
#' Choosing ICA will execute the ordering procedure described in Trapnell and Cacchiarelli et al.,
#' which was implemented in Monocle version 1. \code{\link[DDRTree]{DDRTree}} is a more recent manifold
#' learning algorithm developed by Qi Mao and colleages. It is substantially more
#' powerful, accurate, and robust for single-cell trajectory analysis than ICA,
#' and is now the default method.
#'
#' Often, experiments include cells from different batches or treatments. You can
#' reduce the effects of these treatments by transforming the data with a linear
#' model prior to dimensionality reduction. To do so, provide a model formula
#' through \code{residualModelFormulaStr}.
#'
#' Prior to reducing the dimensionality of the data, it usually helps
#' to normalize it so that highly expressed or highly variable genes don't
#' dominate the computation. \code{reduceDimension()} automatically transforms
#' the data in one of several ways depending on the \code{expressionFamily} of
#' the CellDataSet object. If the expressionFamily is \code{negbinomial} or \code{negbinomial.size}, the
#' data are variance-stabilized. If the expressionFamily is \code{Tobit}, the data
#' are adjusted by adding a pseudocount (of 1 by default) and then log-transformed.
#' If you don't want any transformation at all, set norm_method to "none" and
#' pseudo_expr to 0. This maybe useful for single-cell qPCR data, or data you've
#' already transformed yourself in some way.
#'
#' @param cds the CellDataSet upon which to perform this operation
#' @param max_components the dimensionality of the reduced space
#' @param reduction_method A character string specifying the algorithm to use for dimensionality reduction.
#' @param auto_param_selection when this argument is set to TRUE (default), it will automatically calculate the proper value for the ncenter (number of centroids) parameters which will be passed into DDRTree call.
#' @param scaling When this argument is set to TRUE (default), it will scale each gene before running trajectory reconstruction.
#' @param verbose Whether to emit verbose output during dimensionality reduction
#' @param ... additional arguments to pass to the dimensionality reduction function
#' @return an updated CellDataSet object
#' @import methods
#' @importFrom matrixStats rowSds
#' @importFrom limma removeBatchEffect
#' @importFrom fastICA  ica.R.def ica.R.par
#' @import irlba
#' @import DDRTree
#' @import Rtsne
#' @importFrom stats dist prcomp
#' @importFrom igraph graph.adjacency
#' @export
reduceDimension <- function(cds,
                            max_components=2,
                            reduction_method=c("DDRTree", "ICA", 'tSNE', "UMAP"), # , "SSE", "SimplePPT", 'L1-graph', 'graphL1'
                            auto_param_selection = TRUE,
                            scaling = TRUE,
                            verbose=FALSE,
                            ...){
  extra_arguments <- list(...)
  set.seed(2016) #ensure results from RNG sensitive algorithms are the same on all calls
  
  if (verbose)
    message("Retrieving normalized and PCA (LSI) reduced data ...")
  
  FM <- cds@auxOrderingData$normalize_expr_data
  irlba_pca_res <- cds@normalized_data_projection
  
  if(is.null(FM)) {
    message('Warning: The cds is not normalized or PCA (LSI) reduced with preprocessCDS function yet, running preprocessCDS with default parameters!')
    cds <- preprocessCDS(cds)
    FM <- cds@auxOrderingData$normalize_expr_data
    irlba_pca_res <- cds@normalized_data_projection
  }
  
  #FM <- FM[apply(FM, 1, function(x) all(is.finite(x))), ] #ensure all the expression values are finite values
  if (is.function(reduction_method)) {
    
    if(scaling){
      FM <- as.matrix(Matrix::t(scale(Matrix::t(FM))))
      FM <- FM[!is.na(row.names(FM)), ]
    } else FM <- as.matrix(FM)
    
    reducedDim <- reduction_method(FM, ...)
    colnames(reducedDim) <- colnames(FM)
    reducedDimW(cds) <- as.matrix(reducedDim)
    reducedDimA(cds) <- as.matrix(reducedDim)
    reducedDimS(cds) <- as.matrix(reducedDim)
    reducedDimK(cds) <- as.matrix(reducedDim)
    dp <- as.matrix(dist(reducedDim))
    cellPairwiseDistances(cds) <- dp
    gp <- graph.adjacency(dp, mode = "undirected", weighted = TRUE)
    dp_mst <- minimum.spanning.tree(gp)
    minSpanningTree(cds) <- dp_mst
    cds@dim_reduce_type <- "function_passed"
  }
  else{
    reduction_method <- match.arg(reduction_method)
    if (reduction_method == "ICA") {
      .Deprecated(msg = 'ICA (used in monocle 1) is not supported anymore, please use RGE (reversed graph embedding) instead!')
      # if (verbose)
      #   message("Reducing to independent components")
      # init_ICA <- ica_helper(Matrix::t(FM), max_components,
      #                        use_irlba = TRUE, ...)
      # x_pca <- Matrix::t(Matrix::t(FM) %*% init_ICA$K)
      # W <- Matrix::t(init_ICA$W)
      # weights <- W
      # A <- Matrix::t(solve(weights) %*% Matrix::t(init_ICA$K))
      # colnames(A) <- colnames(weights)
      # rownames(A) <- rownames(FM)
      # S <- weights %*% x_pca
      # rownames(S) <- colnames(weights)
      # colnames(S) <- colnames(FM)
      # reducedDimW(cds) <- as.matrix(W)
      # reducedDimA(cds) <- as.matrix(A)
      # reducedDimS(cds) <- as.matrix(S)
      # reducedDimK(cds) <- as.matrix(init_ICA$K)
      # adjusted_S <- Matrix::t(reducedDimS(cds))
      # dp <- as.matrix(dist(adjusted_S))
      # cellPairwiseDistances(cds) <- dp
      # gp <- graph.adjacency(dp, mode = "undirected", weighted = TRUE)
      # dp_mst <- minimum.spanning.tree(gp)
      # minSpanningTree(cds) <- dp_mst
      # cds@dim_reduce_type <- "ICA"
    } else if (reduction_method == "tSNE") {
      if("num_dim" %in% names(extra_arguments)){ #when you pass pca_dim to the function, the number of dimension used for tSNE dimension reduction is used
        num_dim <- extra_arguments$num_dim #variance_explained
      }
      else{
        num_dim <- 50
      }
      
      topDim_pca <- irlba_pca_res#[, 1:num_dim]
      
      #then run tSNE
      if (verbose)
        message("Reduce dimension by tSNE ...")
      
      tsne_res <- Rtsne::Rtsne(as.matrix(topDim_pca), dims = max_components, pca = F, check_duplicates=FALSE, ...)
      
      tsne_data <- tsne_res$Y[, 1:max_components]
      row.names(tsne_data) <- colnames(tsne_data)
      
      reducedDimA(cds) <- t(tsne_data) #this may move to the auxClusteringData environment
      
      #set the important information from densityClust to certain part of the cds object:
      cds@auxClusteringData[["tSNE"]]$pca_components_used <- num_dim

      cds@dim_reduce_type <- "tSNE"
      
      pData(cds)$tsne_1 = reducedDimA(cds)[1,]
      pData(cds)$tsne_2 = reducedDimA(cds)[2,]
    }
    else if (reduction_method %in% c("DDRTree")) {
      
      message('DDRTree will be eventually deprecated in reduceDimension call and be used in RGE function instead. We are calling RGE for you now.')
      cds <- learnGraph(cds, method = 'DDRTree', ...)
      
    }else if (reduction_method %in% c("UMAP") ) {  
      if (verbose)
        message("Running Uniform Manifold Approximation and Projection")
      
      umap_args <- c(list(X = irlba_pca_res, log = F, n_component = as.integer(max_components), verbose = verbose, return_all = T),
                     extra_arguments[names(extra_arguments) %in% 
                                       c("python_home", "n_neighbors", "metric", "n_epochs", "negative_sample_rate", "alpha", "init", "min_dist", "spread", 
                                         'set_op_mix_ratio', 'local_connectivity', 'bandwidth', 'gamma', 'a', 'b', 'random_state', 'metric_kwds', 'angular_rp_forest', 'verbose')])
      tmp <- do.call(UMAP, umap_args)
      tmp$embedding_ <- (tmp$embedding_ - min(tmp$embedding_)) / max(tmp$embedding_) # normalize UMAP space
      umap_res <- tmp$embedding_; 
      
      adj_mat <- Matrix::sparseMatrix(i = tmp$graph$indices, p = tmp$graph$indptr, 
                                      x = -as.numeric(tmp$graph$data), dims = c(ncol(cds), ncol(cds)), index1 = F, 
                                      dimnames = list(colnames(cds), colnames(cds)))
      
      S <- t(umap_res)
      
      Y <- S
      W <- t(irlba_pca_res)
      
      if(verbose)
        message("Running louvain clustering algorithm ...")
      row.names(umap_res) <- colnames(FM)
      louvain_clustering_args <- c(list(data = umap_res, pd = pData(cds)[colnames(FM), ], verbose = verbose),
                                   extra_arguments[names(extra_arguments) %in% c("k", "weight", "louvain_iter")])
      louvain_res <- do.call(louvain_clustering, louvain_clustering_args)
      
      if("louvain_qval" %in% names(extra_arguments)){ 
        louvain_qval <- extra_arguments$louvain_qval 
      }
      else{
        louvain_qval <- 0.05
      }
      
      cluster_graph_res <- compute_louvain_connected_components(louvain_res$g, louvain_res$optim_res, louvain_qval, verbose)
      louvain_component = components(cluster_graph_res$cluster_g)$membership[louvain_res$optim_res$membership]
      cds@auxOrderingData[["L1graph"]]$louvain_component = louvain_component
      names(louvain_component) = colnames(FM)
      louvain_component = as.factor(louvain_component)
      pData(cds)$louvain_component <- louvain_component
      
      minSpanningTree(cds) <- louvain_res$g
      
      A <- S
      colnames(A) <- colnames(FM)
      reducedDimA(cds) <- A
      
      colnames(S) <- colnames(FM)
      colnames(Y) <- colnames(FM)
      reducedDimW(cds) <- W 
      reducedDimS(cds) <- as.matrix(Y)
      reducedDimK(cds) <- S
      
      cds@auxOrderingData$UMAP <- list(umap_res = umap_res, louvain_res = louvain_res, adj_mat = adj_mat, cluster_graph_res = cluster_graph_res)
      cds@dim_reduce_type <- reduction_method
    } else {
      stop("Error: unrecognized dimensionality reduction method")
    }
  }
  cds
}

#' This function tries to learn a smooth embedding from the noisy reduced dimension using different techniques 
#' @description The function relies on smooth skeleton learning or force direct layout function to learn a smoothier
#' representation of the data. It can be used to facilitate visualization of the data or downstream graph learning 
#' 
#' @param cds the CellDataSet upon which to perform this operation
#' @param max_components the dimensionality of the smoothed reduced space
#' @param smooth_method A character string specifying the algorithm to use for smoothing the embedding
#' @param verbose Whether to emit verbose output during dimensionality reduction
#' @param ... additional arguments to pass to the smoothEmbedding function
#' @return an updated CellDataSet object
smoothEmbedding <- function(cds,
                            max_components = 2, 
                            smooth_method = c('SSE'), #, 'FDL'
                            verbose = FALSE, 
                            ...){
  extra_arguments <- list(...)
  set.seed(2016) #ensure results from RNG sensitive algorithms are the same on all calls
  
  if (verbose)
    message("Retrieving normalized and PCA (LSI) reduced data ...")
  
  FM <- cds@auxOrderingData$normalize_expr_data
  irlba_pca_res <- cds@normalized_data_projection
  
  landmark_id <- NULL
  
  if(smooth_method == 'SSE') {
    if(c("UMAP") %in% names(cds@auxOrderingData)) {
      # if UMAP or L1graph has already done, use those information 
      irlba_pca_res <- cds@normalized_data_projection
      S <- t(cds@auxOrderingData$UMAP$umap_res)
      umap_res <- cds@auxOrderingData$UMAP$umap_res
      data <- umap_res
      adj_mat <- cds@auxOrderingData$UMAP$adj_mat
    } else {
      irlba_pca_res <- cds@normalized_data_projection
      data = irlba_pca_res
      S <- t(irlba_pca_res)
    }
    
    # if number of cells is larger than 20 k, peform landmark selection and do SSE, L1 on the landmarks. We will project others cells on the learn SSE space 
    louvain_res_ori <- NULL
    if(ncol(cds) > 5000) {
      data_ori <- data 
      adj_mat_ori <- adj_mat
      FM_ori <- FM
      
      if("landmark_num" %in% names(extra_arguments)) {
        landmark_num <- extra_arguments$landmark_num
      } else {
        landmark_num <- 2000
      }
      
      centers <- data_ori[seq(1, nrow(data_ori), length.out=landmark_num), ]
      centers <- centers + matrix(rnorm(length(centers), sd = 1e-10), nrow = nrow(centers)) # add random noise 
      kmean_res <- kmeans(data_ori, landmark_num, centers=centers, iter.max = 100)
      landmark_id <- sort(unique(apply(as.matrix(proxy::dist(data_ori, kmean_res$centers)), 2, which.min))) # avoid duplicated points 

      data <- data_ori[landmark_id, ]
      FM <- FM_ori[, landmark_id]
      
      # run UMAP to get the adjacency graph for downstream SSE learning 
      umap_args <- c(list(X = irlba_pca_res[landmark_id, ], log = F, n_component = as.integer(max_components), verbose = verbose, return_all = T),
                     extra_arguments[names(extra_arguments) %in% 
                                       c("python_home", "n_neighbors", "metric", "n_epochs", "negative_sample_rate", "alpha", "init", "min_dist", "spread", 
                                         'set_op_mix_ratio', 'local_connectivity', 'bandwidth', 'gamma', 'a', 'b', 'random_state', 'metric_kwds', 'angular_rp_forest', 'verbose')])
      tmp <- do.call(UMAP, umap_args)
      adj_mat <- Matrix::sparseMatrix(i = tmp$graph$indices, p = tmp$graph$indptr, 
                                      x = -as.numeric(tmp$graph$data), dims = c(length(landmark_id), length(landmark_id)), index1 = F, 
                                      dimnames = list(colnames(cds)[landmark_id], colnames(cds)[landmark_id]))
      
      if(verbose)
        message("Running louvain clustering algorithm (UMAP space) ...")
      row.names(umap_res) <- colnames(FM_ori)
      louvain_clustering_args <- c(list(data = umap_res, pd = pData(cds)[colnames(FM_ori), ], verbose = verbose),
                                   extra_arguments[names(extra_arguments) %in% c("k", "weight", "louvain_iter")])
      louvain_res <- do.call(louvain_clustering, louvain_clustering_args)
      
      louvain_res_ori <- louvain_res 
    }
    
    if (verbose)
      message("Running Smooth Skeleton Embedding ...")
    
    # we are gonna ignore the method arugment here 
    sse_args <- c(list(data=data, dist_mat = adj_mat, embeding_dim=max_components, d = max_components, verbose = verbose),
                  extra_arguments[names(extra_arguments) %in% c("method", "para.gamma", "knn", "C", "maxiter", "beta")])
    SSE_res <- do.call(SSE, sse_args)
    
    # "Y" "K" "W" "U" "V"
    if(verbose)
      message("Running louvain clustering algorithm ...")
    
    SSE_res$Y <- SSE_res$Y * 100 
    K <- SSE_res$Y
    row.names(K) <- paste('cell_', 1:nrow(K), sep = '')
    S <- t(K)
    
    Y <- S
    W <- as.matrix(t(SSE_res$W))
    
    # SSE_res$Y <- (SSE_res$Y - min(SSE_res$Y)) / max(SSE_res$Y) # normalize the SSE space to avoid space shrinking in L1graph step 
    row.names(SSE_res$Y) <- colnames(FM)
    louvain_clustering_args <- c(list(data = SSE_res$Y, pd = pData(cds)[colnames(FM), ], verbose = verbose),
                                 extra_arguments[names(extra_arguments) %in% c("k", "weight", "louvain_iter")])
    louvain_res <- do.call(louvain_clustering, louvain_clustering_args)
    
    # now let us project other cells back to the landmark space 
    projection_res <- NULL
    if(ncol(cds) > 5000) {
      projection_res <- matrix(0, nrow = max_components, ncol = ncol(cds))
      
      # get a graph with distance between cells as the weight 
      relations <- louvain_res_ori$relations
      relations$weight <- reshape2::melt(t(louvain_res_ori$distMatrix))[, 3]
      g <- igraph::graph.data.frame(relations, directed = FALSE)
      
      # iterate each non-landmark cells and project it to the SSE space
      if("cores" %in% names(extra_arguments)) {
        cores <- extra_arguments$landmark_num
      } else {
        cores <- detectCores() 
      }
      
      if(verbose) 
        message('project other non-landmark cells to the landmark SSE space ...')
      
      # using matrix multiplication to accelerate the process (optimize it to handle millions of points)
      block_size <- 50000
      num_blocks = ceiling(nrow(data_ori) / block_size)
      weight_mat <- NULL
      for (i in 1:num_blocks){
        if (i < num_blocks){
          block = data_ori[((((i-1) * block_size)+1):(i*block_size)), ]
        }else{
          block = data_ori[((((i-1) * block_size)+1):(nrow(data_ori))),]
        }
        distances_Z_to_Y <- proxy::dist(block, data_ori[landmark_id, ])
        
        tmp <- as(t(apply(distances_Z_to_Y , 1, function(x) {
          tmp <- sort(x)[6] 
          x[x > tmp] <- 0; p <- rep(0, length(x))
          bandwidth <- mean(range(x[x > 0])) # half of the range of the nearest neighbors as bindwidth 
          p[x > 0] <- exp(-x[x > 0]/bandwidth) # Gaussian kernel 
          p / sum(p) 
        })), 'sparseMatrix')
        
        weight_mat <- rBind(weight_mat, tmp)
      }
      
      projection_res <- t(weight_mat %*% as(t(Y), 'sparseMatrix'))
      
      projection_res[, landmark_id] <- Y
      Y <- as.matrix(projection_res)
      
      FM <- FM_ori
    }
    
    colnames(Y) <- colnames(FM)
    reducedDimA(cds) <- Y
    # colnames(S) <- colnames(FM)
    colnames(Y) <- colnames(FM)
    reducedDimW(cds) <- W # update this !!! 
    reducedDimS(cds) <- as.matrix(Y)
    reducedDimK(cds) <- t(K)
    
    dimnames(SSE_res$W) <- list(paste('cell_', 1:nrow(W), sep = ''), paste('cell_', 1:nrow(W), sep = ''))
    gp <- graph_from_adjacency_matrix(SSE_res$W, weighted=TRUE, add.rownames="code")
    minSpanningTree(cds) <- gp
    
    # pData(cds)$louvain_cluster <- as.character(igraph::membership(louvain_res$optim_res)) 
    cds@auxOrderingData[[smooth_method]] <- list(SSE_res = SSE_res, projection_res = projection_res, 
                                                 louvain_module = as.factor(igraph::membership(louvain_res$optim_res)), 
                                                 louvain_res_ori = louvain_res_ori, louvain_res = louvain_res, adj_mat = adj_mat, landmark_id = landmark_id)
    
    cds@dim_reduce_type <- smooth_method
  } else if(smooth_method == 'FDL') {
    FDL_args <- c(list(cds=cds, verbose = verbose, cell_num_threshold = 0),
                  extra_arguments[names(extra_arguments) %in% c("separate_group", "method", "merge_coords_method", "start.temp", "k", "cell_num_threshold")])
    cds <- do.call(FDL, FDL_args)
  }
  
  cds
}

#' Learn principal graph from the reduced space using reversed graph embedding 
#' 
#' @description Monocle aims to learn how cells transition through a biological program of 
#' gene expression changes in an experiment. Each cell can be viewed as a point 
#' in a high-dimensional space, where each dimension describes the expression of 
#' a different gene in the genome. Identifying the program of gene expression 
#' changes is equivalent to learning a \emph{trajectory} that the cells follow
#' through this space. However, the more dimensions there are in the analysis,
#' the harder the trajectory is to learn. Fortunately, many genes typically
#' co-vary with one another, and so the dimensionality of the data can be
#' reduced with a wide variety of different algorithms. Monocle provides two
#' different algorithms for dimensionality reduction via \code{reduceDimension}.
#' Both take a CellDataSet object and a number of dimensions allowed for the
#' reduced space. You can also provide a model formula indicating some variables
#' (e.g. batch ID or other technical factors) to "subtract" from the data so it
#' doesn't contribute to the trajectory. 
#' 
#' @details You can choose two different reduction algorithms: Independent Component 
#' Analysis (ICA) and Discriminative Dimensionality Reduction with Trees (DDRTree).
#' The choice impacts numerous downstream analysis steps, including \code{\link{orderCells}}.
#' Choosing ICA will execute the ordering procedure described in Trapnell and Cacchiarelli et al.,
#' which was implemented in Monocle version 1. \code{\link[DDRTree]{DDRTree}} is a more recent manifold
#' learning algorithm developed by Qi Mao and colleages. It is substantially more
#' powerful, accurate, and robust for single-cell trajectory analysis than ICA,
#' and is now the default method.
#'
#' Often, experiments include cells from different batches or treatments. You can
#' reduce the effects of these treatments by transforming the data with a linear
#' model prior to dimensionality reduction. To do so, provide a model formula
#' through \code{residualModelFormulaStr}.
#'
#'
#' @param cds the CellDataSet upon which to perform this operation
#' @param max_components the dimensionality of the reduced space
#' @param RGE_method Determines how to transform expression values prior to reducing dimensionality
#' @param auto_param_selection when this argument is set to TRUE (default), it will automatically calculate the proper value for the ncenter (number of centroids) parameters which will be passed into DDRTree call.
#' @param partition_component When this argument is set to TRUE (default to be FALSE), we will learn a tree structure for each separate over-connected louvain component. 
#' @param scaling When this argument is set to TRUE (default), it will scale each gene before running trajectory reconstruction.
#' @param verbose Whether to emit verbose output during dimensionality reduction
#' @param ... additional arguments to pass to the dimensionality reduction function
#' @return an updated CellDataSet object
#' @import methods
#' @importFrom matrixStats rowSds
#' @importFrom limma removeBatchEffect
#' @importFrom fastICA  ica.R.def ica.R.par
#' @import irlba
#' @import DDRTree
#' @import Rtsne
#' @importFrom stats dist prcomp
#' @importFrom igraph graph.adjacency
#' @export
learnGraph <- function(cds,
                       max_components=2,
                       RGE_method = c('L1graph', 'SimplePPT', 'DDRTree'), 
                       auto_param_selection = TRUE, 
                       partition_component = TRUE, 
                       scale = FALSE, 
                       verbose = FALSE, 
                       ...){
  
  extra_arguments <- list(...)
  FM <- cds@auxOrderingData$normalize_expr_data
  irlba_pca_res <- cds@normalized_data_projection
  
  if(RGE_method == 'L1graph') { 
    if(cds@dim_reduce_type == "SSE") {
      Y <- cds@auxOrderingData[['SSE']]$SSE_res$Y
      Y <- Y 
      # SSE_res$Y <- (SSE_res$Y - min(SSE_res$Y)) / max(SSE_res$Y) # normalize the SSE space to avoid space shrinking in L1graph step 
      
      louvain_res <- cds@auxOrderingData[["SSE"]]$louvain_res
      louvain_module_length = length(levels(cds@auxOrderingData[['SSE']]$louvain_module))
      
      landmark_id <- cds@auxOrderingData[['SSE']]$landmark_id
      if(is.null(landmark_id)) {
        reduced_dim_res = t(Y) 
        row.names(Y) <- colnames(FM)
        landmark_id <- 1:ncol(cds)
      } else {
        reduced_dim_res = t(Y)           
        row.names(Y) <- colnames(FM[, landmark_id])
      }
    } else if(cds@dim_reduce_type == "UMAP") {
      Y <- cds@auxOrderingData[['UMAP']]$umap_res
      row.names(Y) <- colnames(FM)
      reduced_dim_res = t(Y) 
      
      louvain_res <- cds@auxOrderingData[["UMAP"]]$louvain_res
      louvain_module_length = length(unique(sort(louvain_res$optim_res$membership)))
      
      if(ncol(cds) < 5000) {
        landmark_id <- 1:ncol(cds)
      } else {
        if("landmark_num" %in% names(extra_arguments)) {
          landmark_num <- extra_arguments$landmark_num
        } else {
          landmark_num <- 2000
        }
        
        centers <- Y[seq(1, nrow(Y), length.out=landmark_num), ]
        centers <- centers + matrix(rnorm(length(centers), sd = 1e-10), nrow = nrow(centers)) # add random noise 
        kmean_res <- kmeans(Y, landmark_num, centers=centers, iter.max = 100)
        landmark_id <- sort(unique(apply(as.matrix(proxy::dist(Y, kmean_res$centers)), 2, which.min))) # avoid duplicated points 
        reduced_dim_res <- reduced_dim_res[, landmark_id]
      }
      
    } else {
      stop('L1graph can be only applied to either the MAP or SSE reduced space, please first apply those dimension reduction techniques!')
    }
    
    if("ncenter" %in% names(extra_arguments)){ #avoid overwrite the ncenter parameter
      ncenter <- extra_arguments$ncenter
    }else{
      if("L1.pr_graph_vertex_per_louvain_module" %in% names(extra_arguments)){ #avoid overwrite the ncenter parameter
        L1.pr_graph_vertex_per_louvain_module <- extra_arguments$L1.pr_graph_vertex_per_louvain_module
      }else{
        L1.pr_graph_vertex_per_louvain_module = 3
      }
      ncenter = L1.pr_graph_vertex_per_louvain_module * louvain_module_length
      ncenter = min(ncol(FM) / 2, ncenter)
    }
    
    if (ncenter > ncol(reduced_dim_res))
      stop("Error: ncenters must be less than or equal to ncol(X)")
    
    centers <- t(reduced_dim_res)[seq(1, ncol(reduced_dim_res), length.out=ncenter),]
    centers <- centers + matrix(rnorm(length(centers), sd = 1e-10), nrow = nrow(centers)) # add random noise 
    
    kmean_res <- kmeans(t(reduced_dim_res), ncenter, centers=centers, iter.max = 100)
    if (kmean_res$ifault != 0){
      message(paste("Warning: kmeans returned ifault =", kmean_res$ifault))
    }
    nearest_center = findNearestVertex(t(kmean_res$centers), reduced_dim_res, process_targets_in_blocks=TRUE)
    medioids = reduced_dim_res[,unique(nearest_center)]
    reduced_dim_res <- t(medioids)

    if(verbose)
      message('running L1-graph ...')
    
    X <- t(reduced_dim_res)

    if('C0' %in% names(extra_arguments)){
      C0 <- extra_arguments$C0
    }
    else
      C0 <- X
    Nz <- ncol(C0)
    
    if('nn' %in% names(extra_arguments))
      G_T = get_mst_with_shortcuts(C0, K = extra_arguments$nn)
    else
      G_T = get_mst_with_shortcuts(C0, K = 5)
    
    
    G = G_T$G #+ G_knn$G
    
    G[G > 0] = 1
    W = G_T$W #+ G_knn$W
    
    if("louvain_qval" %in% names(extra_arguments)){ 
      louvain_qval <- extra_arguments$louvain_qval 
    }
    else{
      louvain_qval <- 0.05
    }
    
    cluster_graph_res <- compute_louvain_connected_components(louvain_res$g, louvain_res$optim_res, louvain_qval, verbose)
    louvain_component = components(cluster_graph_res$cluster_g)$membership[louvain_res$optim_res$membership]
    cds@auxOrderingData[["L1graph"]]$louvain_component = louvain_component
    names(louvain_component) <- colnames(cds[, landmark_id])
    louvain_component <- louvain_component[rownames(reduced_dim_res)]
    louvain_component <- as.factor(louvain_component)
    if (length(levels(louvain_component)) > 1){
      louvain_component_mask = as.matrix(tcrossprod(sparse.model.matrix( ~ louvain_component + 0)))
      
      G = G * louvain_component_mask
      W = W * louvain_component_mask
      rownames(G) = rownames(W)
      colnames(G) = colnames(W)
    }

    l1graph_args <- c(list(X = t(reduced_dim_res), C0 = C0, G = G, gstruct = 'l1-graph', verbose = verbose),
                      extra_arguments[names(extra_arguments) %in% c('maxiter', 'eps', 'L1.lambda', 'L1.gamma', 'L1.sigma', 'nn')])
    
    
    l1_graph_res <- do.call(principal_graph, l1graph_args)
    
    colnames(l1_graph_res$C) <-  rownames(reduced_dim_res)
    DCs <- t(reduced_dim_res) #FM
    
    colnames(l1_graph_res$W) <- rownames(reduced_dim_res)
    rownames(l1_graph_res$W) <- rownames(reduced_dim_res)
    
    
    # row.names(l1_graph_res$X) <- colnames(cds)
    reducedDimW(cds) <- l1_graph_res$W
    # reducedDimS(cds) <- DCs
    reducedDimK(cds) <- l1_graph_res$C
    cds@auxOrderingData[["L1graph"]]$objective_vals <- tail(l1_graph_res$objs, 1)
    cds@auxOrderingData[["L1graph"]]$W <- l1_graph_res$W
    cds@auxOrderingData[["L1graph"]]$P <- l1_graph_res$P
    
    adjusted_K <- Matrix::t(reducedDimK(cds))
    dp <- as.matrix(dist(adjusted_K))
    cellPairwiseDistances(cds) <- dp
    
    W <- l1_graph_res$W
    dimnames(l1_graph_res$W) <- list(paste('cell_', 1:nrow(W), sep = ''), paste('cell_', 1:nrow(W), sep = ''))
    W[W < 1e-5] <- 0
    gp <- graph.adjacency(W, mode = "undirected", weighted = TRUE)
    # dp_mst <- minimum.spanning.tree(gp)
    minSpanningTree(cds) <- gp
    cds@dim_reduce_type <- "L1graph"
    cds <- findNearestPointOnMST(cds)
  } else if(RGE_method == 'SimplePPT') {
    if(ncol(cds@reducedDimS) > 1) {
      irlba_pca_res <- t(cds@reducedDimS)
    }
    
    if(auto_param_selection & ncol(cds) >= 100){
      if("ncenter" %in% names(extra_arguments)) #avoid overwrite the ncenter parameter
        ncenter <- extra_arguments$ncenter
      else
        ncenter <- cal_ncenter(nrow(irlba_pca_res))
      #add other parameters...
      if(scale) 
        X <- as.matrix(scale(t(irlba_pca_res)))
      else 
        X <- t(irlba_pca_res)
      
      ddr_args <- c(list(X=X, dimensions=ncol(X), ncenter=ncenter, no_reduction = T, verbose = verbose),
                    extra_arguments[names(extra_arguments) %in% c("initial_method", "maxIter", "sigma", "lambda", "param.gamma", "tol")])
      #browser()
      ddrtree_res <- do.call(DDRTree, ddr_args)
    } else{
      if(scale) 
        X <- as.matrix(scale(t(irlba_pca_res)))
      else 
        X <- t(irlba_pca_res)
      
      ddrtree_res <- DDRTree(X, dimensions=ncol(X), no_reduction = T, verbose = verbose, ...)
    }
    
    # ddrtree_res <- DDRTree(as.matrix(scale(t(irlba_pca_res))), max_components, no_reduction = T, verbose = verbose, ...)
    if(ncol(ddrtree_res$Y) == ncol(cds))
      colnames(ddrtree_res$Y) <- colnames(FM) #paste("Y_", 1:ncol(ddrtree_res$Y), sep = "")
    else
      colnames(ddrtree_res$Y) <- paste("Y_", 1:ncol(ddrtree_res$Y), sep = "")
    
    colnames(ddrtree_res$Z) <- colnames(FM)
    reducedDimW(cds) <- ddrtree_res$W
    reducedDimS(cds) <- ddrtree_res$Z
    reducedDimK(cds) <- ddrtree_res$Y
    cds@auxOrderingData[["DDRTree"]] <- ddrtree_res[c('stree', 'Q', 'R', 'objective_vals', 'history')]
    
    adjusted_K <- Matrix::t(reducedDimK(cds))
    dp <- as.matrix(dist(adjusted_K))
    cellPairwiseDistances(cds) <- dp
    gp <- graph.adjacency(dp, mode = "undirected", weighted = TRUE)
    dp_mst <- minimum.spanning.tree(gp)
    minSpanningTree(cds) <- dp_mst
    
    cds@dim_reduce_type <- "DDRTree"
    
    if(ncol(cds) < 100) { 
      cds <- findNearestPointOnMST(cds)
    } else {
      tmp <- matrix(apply(cds@auxOrderingData$DDRTree$R, 1, which.max))
      row.names(tmp) <- colnames(cds)
      cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex <- tmp
    }
    # if(cds@dim_reduce_type == "SSE") {
    #   Y <- cds@auxOrderingData[['SSE']]$SSE_res$Y
    #   Y <- Y 
    #   # SSE_res$Y <- (SSE_res$Y - min(SSE_res$Y)) / max(SSE_res$Y) # normalize the SSE space to avoid space shrinking in L1graph step 
    
    #   louvain_res <- cds@auxOrderingData[["SSE"]]$louvain_res
    #   louvain_module_length = length(levels(cds@auxOrderingData[['SSE']]$louvain_module))
    
    #   landmark_id <- cds@auxOrderingData[['SSE']]$landmark_id
    #   if(is.null(landmark_id)) {
    #     reduced_dim_res = t(Y) 
    #     row.names(Y) <- colnames(FM)
    #     landmark_id <- 1:ncol(cds)
    #   } else {
    #     reduced_dim_res = t(Y)           
    #     row.names(Y) <- colnames(FM[, landmark_id])
    #   }
    # } else if(cds@dim_reduce_type == "UMAP") {
    #   Y <- cds@auxOrderingData[['UMAP']]$umap_res
    #   row.names(Y) <- colnames(FM)
    #   reduced_dim_res = t(Y) 
    
    #   louvain_res <- cds@auxOrderingData[["UMAP"]]$louvain_res
    #   louvain_module_length = length(unique(sort(louvain_res$optim_res$membership)))
    
    #   if(ncol(cds) < 5000) {
    #     landmark_id <- 1:ncol(cds)
    #   } else {
    #     if("landmark_num" %in% names(extra_arguments)) {
    #       landmark_num <- extra_arguments$landmark_num
    #     } else {
    #       landmark_num <- 2000
    #     }
    
    #     centers <- Y[seq(1, nrow(Y), length.out=landmark_num), ]
    #     centers <- centers + matrix(rnorm(length(centers), sd = 1e-10), nrow = nrow(centers)) # add random noise 
    #     kmean_res <- kmeans(Y, landmark_num, centers=centers, iter.max = 100)
    #     landmark_id <- sort(unique(apply(as.matrix(proxy::dist(Y, kmean_res$centers)), 2, which.min))) # avoid duplicated points 
    #     reduced_dim_res <- reduced_dim_res[, landmark_id]
    #   }
    
    # } else {
    #   stop('L1graph can be only applied to either the MAP or SSE reduced space, please first apply those dimension reduction techniques!')
    # }
    
    # if("ncenter" %in% names(extra_arguments)){ #avoid overwrite the ncenter parameter
    #   ncenter <- extra_arguments$ncenter
    # }else{
    #   if("L1.pr_graph_vertex_per_louvain_module" %in% names(extra_arguments)){ #avoid overwrite the ncenter parameter
    #     L1.pr_graph_vertex_per_louvain_module <- extra_arguments$L1.pr_graph_vertex_per_louvain_module
    #   }else{
    #     L1.pr_graph_vertex_per_louvain_module = 3
    #   }
    #   ncenter = L1.pr_graph_vertex_per_louvain_module * louvain_module_length
    #   ncenter = min(ncol(FM) / 2, ncenter)
    # }
    
    # if (ncenter > ncol(reduced_dim_res))
    #   stop("Error: ncenters must be less than or equal to ncol(X)")
    
    # centers <- t(reduced_dim_res)[seq(1, ncol(reduced_dim_res), length.out=ncenter),]
    # centers <- centers + matrix(rnorm(length(centers), sd = 1e-10), nrow = nrow(centers)) # add random noise 
    
    # kmean_res <- kmeans(t(reduced_dim_res), ncenter, centers=centers, iter.max = 100)
    # if (kmean_res$ifault != 0){
    #   message(paste("Warning: kmeans returned ifault =", kmean_res$ifault))
    # }
    # nearest_center = findNearestVertex(t(kmean_res$centers), reduced_dim_res, process_targets_in_blocks=TRUE)
    # medioids = reduced_dim_res[,unique(nearest_center)]
    # reduced_dim_res <- t(medioids)
    # #reduced_dim_res = t(reduced_dim_res)[centers,]
    # #reduced_dim_res <- t(reduced_dim_res)
    # #rownames(reduced_dim_res) = paste("Y_", 1:nrow(reduced_dim_res), sep = "")
    
    # if(verbose)
    #   message('running L1-graph ...')
    
    # #X <- t(reduced_dim_res)
    # X <- t(reduced_dim_res)
    # # D <- nrow(X); N <- ncol(X)
    # # Z <- X
    
    # if('C0' %in% names(extra_arguments)){
    #   C0 <- extra_arguments$C0
    # }
    # else
    #   C0 <- X
    # Nz <- ncol(C0)
    
    
    # #G_T = get_mst(C0)
    
    # # # print(extra_arguments)
    # # if('nn' %in% names(extra_arguments))
    # #   G_knn <- get_knn(C0, K = extra_arguments$nn)
    # # else
    # #   G_knn <- get_knn(C0, K = 5)
    
    # # print(extra_arguments)
    # if('nn' %in% names(extra_arguments))
    #   G_T = get_mst_with_shortcuts(C0, K = extra_arguments$nn)
    # else
    #   G_T = get_mst_with_shortcuts(C0, K = 5)
    
    
    # G = G_T$G #+ G_knn$G
    
    # G[G > 0] = 1
    # W = G_T$W #+ G_knn$W
    
    # if("louvain_qval" %in% names(extra_arguments)){ 
    #   louvain_qval <- extra_arguments$louvain_qval 
    # }
    # else{
    #   louvain_qval <- 0.05
    # }
    
    # cluster_graph_res <- compute_louvain_connected_components(louvain_res$g, louvain_res$optim_res, louvain_qval, verbose)
    # louvain_component = components(cluster_graph_res$cluster_g)$membership[louvain_res$optim_res$membership]
    # cds@auxOrderingData[["L1graph"]]$louvain_component = louvain_component
    # # pData(cds)$louvain_component <- as.factor(louvain_component)
    # names(louvain_component) <- colnames(cds[, landmark_id])
    # louvain_component <- louvain_component[rownames(reduced_dim_res)]
    # louvain_component <- as.factor(louvain_component)
    # if (length(levels(louvain_component)) > 1){
    #   louvain_component_mask = as.matrix(tcrossprod(sparse.model.matrix( ~ louvain_component + 0)))
    
    #   G = G * louvain_component_mask
    #   W = W * louvain_component_mask
    #   rownames(G) = rownames(W)
    #   colnames(G) = colnames(W)
    # }
    # #louvain_component = data.frame(component=louvain_component)
    
    # l1graph_args <- c(list(X = t(reduced_dim_res), C0 = C0, G = G, gstruct = 'span-tree', verbose = verbose),
    #                   extra_arguments[names(extra_arguments) %in% c('maxiter', 'eps', 'L1.lambda', 'L1.gamma', 'L1.sigma', 'nn')])
    
    
    # l1_graph_res <- do.call(principal_graph, l1graph_args)
    
    # colnames(l1_graph_res$C) <-  rownames(reduced_dim_res)
    # DCs <- t(reduced_dim_res) #FM
    
    # colnames(l1_graph_res$W) <- rownames(reduced_dim_res)
    # rownames(l1_graph_res$W) <- rownames(reduced_dim_res)
    
    
    # # row.names(l1_graph_res$X) <- colnames(cds)
    # reducedDimW(cds) <- l1_graph_res$W
    # reducedDimS(cds) <- DCs
    # reducedDimK(cds) <- l1_graph_res$C
    # cds@auxOrderingData[["SimplePPT"]]$objective_vals <- tail(l1_graph_res$objs, 1)
    # cds@auxOrderingData[["SimplePPT"]]$W <- l1_graph_res$W
    # cds@auxOrderingData[["SimplePPT"]]$P <- l1_graph_res$P
    
    # adjusted_K <- Matrix::t(reducedDimK(cds))
    # dp <- as.matrix(dist(adjusted_K))
    # cellPairwiseDistances(cds) <- dp
    
    # W <- l1_graph_res$W
    # dimnames(l1_graph_res$W) <- list(paste('cell_', 1:nrow(W), sep = ''), paste('cell_', 1:nrow(W), sep = ''))
    # W[W < 1e-5] <- 0
    # gp <- graph.adjacency(W, mode = "undirected", weighted = TRUE)
    # # dp_mst <- minimum.spanning.tree(gp)
    # minSpanningTree(cds) <- gp
    # cds@dim_reduce_type <- "SimplePPT"
    # cds <- findNearestPointOnMST(cds)
  } else if(RGE_method == 'DDRTree') {
    if(ncol(cds@reducedDimS) > 1) {
      irlba_pca_res <- t(cds@reducedDimS)
    }
    
    row.names(irlba_pca_res) <- colnames(FM)
    
    if (verbose)
      message("Learning principal graph with DDRTree")
    
    # TODO: DDRTree should really work with sparse matrices.
    louvain_component <- pData(cds)$louvain_component
    
    if(length(louvain_component) == ncol(cds) & partition_component) {
      X <- t(irlba_pca_res)
      
      reducedDimK_coord <- NULL  
      dp_mst <- NULL 
      pr_graph_cell_proj_closest_vertex <- NULL 
      cell_name_vec <- NULL
      
      for(cur_comp in unique(louvain_component)) {
        X_subset <- X[, louvain_component == cur_comp]
        
        #add other parameters...
        if(scale) 
          X_subset <- t(as.matrix(scale(t(X_subset))))
        
        ncenter <- cal_ncenter(ncol(X_subset))
        
        ddr_args <- c(list(X=X_subset, dimensions=max_components, ncenter=ncenter, no_reduction = T, verbose = verbose),
                      extra_arguments[names(extra_arguments) %in% c("initial_method", "maxIter", "sigma", "lambda", "param.gamma", "tol")])
        #browser()
        ddrtree_res <- do.call(DDRTree, ddr_args)
        
        if(is.null(reducedDimK_coord)) {
          curr_cell_names <- paste("Y_", 1:ncol(ddrtree_res$Y), sep = "")
          pr_graph_cell_proj_closest_vertex <- matrix(apply(ddrtree_res$R, 1, which.max))
          cell_name_vec <- colnames(X_subset)
        } else {
          curr_cell_names <- paste("Y_", ncol(reducedDimK_coord) + 1:ncol(ddrtree_res$Y), sep = "")
          pr_graph_cell_proj_closest_vertex <- rbind(pr_graph_cell_proj_closest_vertex, matrix(apply(ddrtree_res$R, 1, which.max) + ncol(reducedDimK_coord)))
          cell_name_vec <- c(cell_name_vec, colnames(X_subset))
        }
        
        tmp <- ddrtree_res$Y
        
        reducedDimK_coord <- cbind(reducedDimK_coord, tmp)
        
        
        dp <- ddrtree_res$stree[1:ncol(ddrtree_res$Y), 1:ncol(ddrtree_res$Y)]
        dimnames(dp) <- list(curr_cell_names, curr_cell_names)
        
        dp_mst <- graph.union(dp_mst, graph.adjacency(dp, mode = "undirected", weighted = TRUE))
        
        tmp <- matrix(apply(ddrtree_res$R, 1, which.max))
        
      }
      
      row.names(pr_graph_cell_proj_closest_vertex) <- cell_name_vec
      
      ddrtree_res_W <- ddrtree_res$W
      ddrtree_res_Z <- cds@reducedDimS
      ddrtree_res_Y <- reducedDimK_coord
      
      colnames(ddrtree_res_Y) <- paste0("Y_", 1:ncol(ddrtree_res_Y), sep = "")
      
      cds@auxOrderingData[["DDRTree"]] <- ddrtree_res[c('stree', 'Q', 'R', 'objective_vals', 'history')]
      cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex <- pr_graph_cell_proj_closest_vertex
      
    } else {
      ncenter <- NULL
      if(auto_param_selection & ncol(cds) >= 100){
        if("ncenter" %in% names(extra_arguments)) #avoid overwrite the ncenter parameter
          ncenter <- extra_arguments$ncenter
        else
          ncenter <- cal_ncenter(nrow(irlba_pca_res))
        
      } 
      
      if(scale) 
        X <- as.matrix(scale(t(irlba_pca_res)))
      else 
        X <- t(irlba_pca_res)
      
      #add other parameters...
      ddr_args <- c(list(X=X, dimensions=max_components, ncenter=ncenter, verbose = verbose),
                    extra_arguments[names(extra_arguments) %in% c("initial_method", "maxIter", "sigma", "lambda", "param.gamma", "tol", "no_reduction")])
      #browser()
      ddrtree_res <- do.call(DDRTree, ddr_args)
      
      if(ncol(ddrtree_res$Y) == ncol(cds))
        colnames(ddrtree_res$Y) <- colnames(FM) #paste("Y_", 1:ncol(ddrtree_res$Y), sep = "")
      else
        colnames(ddrtree_res$Y) <- paste("Y_", 1:ncol(ddrtree_res$Y), sep = "")
      
      colnames(ddrtree_res$Z) <- colnames(FM)
      
      ddrtree_res_W <- ddrtree_res$W
      ddrtree_res_Z <- ddrtree_res$Z
      ddrtree_res_Y <- ddrtree_res$Y
      
      adjusted_K <- Matrix::t(reducedDimK(cds))
      dp <- as.matrix(dist(adjusted_K))
      cellPairwiseDistances(cds) <- dp
      gp <- graph.adjacency(dp, mode = "undirected", weighted = TRUE)
      dp_mst <- minimum.spanning.tree(gp)
      
      cds@auxOrderingData[["DDRTree"]] <- ddrtree_res[c('stree', 'Q', 'R', 'objective_vals', 'history')]
      
      if(ncol(cds) < 100) { 
        cds <- findNearestPointOnMST(cds)
      } else {
        tmp <- matrix(apply(cds@auxOrderingData$DDRTree$R, 1, which.max))
        row.names(tmp) <- colnames(cds)
        cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex <- tmp
      }
    }
    
    reducedDimW(cds) <- ddrtree_res_W
    reducedDimS(cds) <- ddrtree_res_Z
    reducedDimK(cds) <- ddrtree_res_Y
    
    minSpanningTree(cds) <- dp_mst
    
    cds@dim_reduce_type <- "DDRTree"
    
  }
  cds 
}

#' Finds the nearest principal graph node
#' @param data_matrix the input matrix
#' @param target_points the target points
#' @param block_size the number of input matrix rows to process per blocl
#' @param process_targets_in_blocks whether to process the targets points in blocks instead
findNearestVertex = function(data_matrix, target_points, block_size=50000, process_targets_in_blocks=FALSE){
  closest_vertex = c()
  if (process_targets_in_blocks == FALSE){
    num_blocks = ceiling(ncol(data_matrix) / block_size)
    for (i in 1:num_blocks){
      if (i < num_blocks){
        block = data_matrix[,((((i-1) * block_size)+1):(i*block_size))]
      }else{
        block = data_matrix[,((((i-1) * block_size)+1):(ncol(data_matrix)))]
      }
      distances_Z_to_Y <- proxy::dist(t(block), t(target_points))
      closest_vertex_for_block <- apply(distances_Z_to_Y, 1, function(z) { which.min(z) } )
      closest_vertex = append(closest_vertex, closest_vertex_for_block)
    }
  }else{
    num_blocks = ceiling(ncol(target_points) / block_size)
    dist_to_closest_vertex = rep(Inf, length(ncol(data_matrix)))
    closest_vertex = rep(NA, length(ncol(data_matrix)))
    for (i in 1:num_blocks){
      if (i < num_blocks){
        block = target_points[,((((i-1) * block_size)+1):(i*block_size))]
      }else{
        block = target_points[,((((i-1) * block_size)+1):(ncol(target_points)))]
      }
      distances_Z_to_Y <- proxy::dist(t(data_matrix), t(block))
      closest_vertex_for_block <- apply(distances_Z_to_Y, 1, function(z) { which.min(z) } )
      new_block_distances = distances_Z_to_Y[cbind(1:nrow(distances_Z_to_Y), closest_vertex_for_block)]
      updated_nearest_idx = which(new_block_distances < dist_to_closest_vertex)
      closest_vertex[updated_nearest_idx] = closest_vertex_for_block[updated_nearest_idx] + (i-1) * block_size
      dist_to_closest_vertex[updated_nearest_idx] = new_block_distances[updated_nearest_idx]
      #closest_vertex = append(closest_vertex, closest_vertex_for_block)
    }
  }
  stopifnot(length(closest_vertex) == ncol(data_matrix))
  #closest_vertex <- which(distance_to_closest == min(distance_to_closest))
  return (closest_vertex)
}

# Project each point to the nearest on the MST:
findNearestPointOnMST <- function(cds){
  dp_mst <- minSpanningTree(cds)
  Z <- reducedDimS(cds)
  Y <- reducedDimK(cds)

  tip_leaves <- names(which(degree(dp_mst) == 1))

  #distances_Z_to_Y <- proxy::dist(t(Z), t(Y))
  #closest_vertex <- apply(distances_Z_to_Y, 1, function(z) { which ( z == min(z) )[1] } )
  closest_vertex = findNearestVertex(Z, Y)
  
  #closest_vertex <- as.vector(closest_vertex)
  closest_vertex_names <- colnames(Y)[closest_vertex]
  closest_vertex_df <- as.matrix(closest_vertex) #index on Z
  row.names(closest_vertex_df) <- names(closest_vertex) #original cell names for projection

  cds@auxOrderingData[[cds@dim_reduce_type]]$pr_graph_cell_proj_closest_vertex <- closest_vertex_df #as.matrix(closest_vertex)
  cds
}

#' @importFrom igraph graph.adjacency V
#' @importFrom stats dist
project2MST <- function(cds, Projection_Method){
  dp_mst <- minSpanningTree(cds)
  Z <- reducedDimS(cds)
  Y <- reducedDimK(cds)

  cds <- findNearestPointOnMST(cds)
  closest_vertex <- cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex

  #closest_vertex <- as.vector(closest_vertex)
  closest_vertex_names <- colnames(Y)[closest_vertex]
  closest_vertex_df <- as.matrix(closest_vertex)
  row.names(closest_vertex_df) <- row.names(closest_vertex)
  #closest_vertex_names <- as.vector(closest_vertex)

  tip_leaves <- names(which(degree(dp_mst) == 1))

  if(!is.function(Projection_Method)) {
    P <- Y[, closest_vertex]
  }
  else{
    P <- matrix(rep(0, length(Z)), nrow = nrow(Z)) #Y
    for(i in 1:length(closest_vertex)) {
      neighbors <- names(V(dp_mst) [ suppressWarnings(nei(closest_vertex_names[i], mode="all")) ])
      projection <- NULL
      distance <- NULL
      Z_i <- Z[, i]

      for(neighbor in neighbors) {
        if(closest_vertex_names[i] %in% tip_leaves) {
          tmp <- projPointOnLine(Z_i, Y[, c(closest_vertex_names[i], neighbor)]) #projPointOnLine: always perform orthogonal projection to the line
        }
        else {
          tmp <- Projection_Method(Z_i, Y[, c(closest_vertex_names[i], neighbor)])
        }
        projection <- rbind(projection, tmp)
        distance <- c(distance, dist(rbind(Z_i, tmp)))
      }
      if(class(projection) != 'matrix')
        projection <- as.matrix(projection)
      P[, i] <- projection[which(distance == min(distance))[1], ] #use only the first index to avoid assignment error
    }
  }
    # tip_leaves <- names(which(degree(minSpanningTree(cds)) == 1))

  colnames(P) <- colnames(Z)

  #reducedDimK(cds) <- P
  dp <- as.matrix(dist(t(P)))
  #dp <- as.matrix(dist(t(reducedDimS(cds))))

  min_dist = min(dp[dp!=0])
  #dp[dp == 0] <- min_dist
  dp <- dp + min_dist #to avoid exact Pseudotime for a lot cells
  diag(dp) <- 0

  cellPairwiseDistances(cds) <- dp
  gp <- graph.adjacency(dp, mode = "undirected", weighted = TRUE)
  dp_mst <- minimum.spanning.tree(gp)

  cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_tree <- dp_mst
  cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_dist <- P #dp, P projection point not output
  cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex <- closest_vertex_df #as.matrix(closest_vertex)

  cds
}

#project points to a line
projPointOnLine <- function(point, line) {
  ap <- point - line[, 1]
  ab <- line[, 2] - line[, 1]

  res <- line[, 1] + c((ap %*% ab) / (ab %*% ab)) * ab
  return(res)
}

# projPointOnLine <- function(point, line){
#   vx = line[1, 2]
#   vy = line[2, 2]

#   # difference of point with line origin
#   dx = point[1] - line[1,1]
#   dy = point[2] - line[2,1]

#   # Position of projection on line, using dot product
#   tp = (dx * vx + dy * vy ) / (vx * vx + vy * vy)

#   # convert position on line to cartesian coordinates
#   point = c(line[1,1] + tp * vx, line[2,1] + tp * vy)

#   return(point)
# }

# Project point to line segment
project_point_to_line_segment <- function(p, df){
  # returns q the closest point to p on the line segment from A to B
  A <- df[, 1]
  B <- df[, 2]
  # vector from A to B
  AB <- (B-A)
  # squared distance from A to B
  AB_squared = sum(AB^2)
  if(AB_squared == 0) {
    # A and B are the same point
    q <- A
  }
  else {
    # vector from A to p
    Ap <- (p-A)
    # from http://stackoverflow.com/questions/849211/
    # Consider the line extending the segment, parameterized as A + t (B - A)
    # We find projection of point p onto the line.
    # It falls where t = [(p-A) . (B-A)] / |B-A|^2
    # t <- max(0, min(1, sum(Ap * AB) / AB_squared))
    t <- sum(Ap * AB) / AB_squared

    if (t < 0.0) {
      # "Before" A on the line, just return A
      q <- A
    }
    else if (t > 1.0) {
      # "After" B on the line, just return B
      q <- B
    }
    else {
      # projection lines "inbetween" A and B on the line
      q <- A + t * AB#
    }
  }
  return(q)
}

# #' traverse from one cell to another cell
# #'
# #' @param g the tree graph learned from monocle 2 during trajectory reconstruction
# #' @param starting_cell the initial vertex for traversing on the graph
# #' @param end_cells the terminal vertex for traversing on the graph
# #' @return a list of shortest path from the initial cell and terminal cell, geodestic distance between initial cell and terminal cells and branch point passes through the shortest path
#' @importFrom igraph shortest.paths shortest_paths degree
traverseTree <- function(g, starting_cell, end_cells){
  distance <- shortest.paths(g, v=starting_cell, to=end_cells)
  branchPoints <- which(degree(g) == 3)
  path <- shortest_paths(g, from = starting_cell, end_cells)

  return(list(shortest_path = path$vpath, distance = distance, branch_points = intersect(branchPoints, unlist(path$vpath))))
}

# #' Make a cds by traversing from one cell to another cell
# #'
# #' @param cds a cell dataset after trajectory reconstruction
# #' @param starting_cell the initial vertex for traversing on the graph
# #' @param end_cells the terminal vertex for traversing on the graph
# #' @return a new cds containing only the cells traversed from the intial cell to the end cell
traverseTreeCDS <- function(cds, starting_cell, end_cells){
  subset_cell <- c()
  dp_mst <- cds@minSpanningTree

  for(end_cell in end_cells) {
    traverse_res <- traverseTree(dp_mst, starting_cell, end_cell)
    path_cells <- names(traverse_res$shortest_path[[1]])

    subset_cell <- c(subset_cell, path_cells)
  }

  subset_cell <- unique(subset_cell)
  cds_subset <- SubSet_cds(cds, subset_cell)

  root_state <- pData(cds_subset[, starting_cell])[, 'State']
  cds_subset <- orderCells(cds_subset, root_state = as.numeric(root_state))

  return(cds_subset)
}

# #' Subset a cds which only includes cells provided with the argument cells
# #'
# #' @param cds a cell dataset after trajectory reconstruction
# #' @param cells a vector contains all the cells you want to subset
# #' @return a new cds containing only the cells from the cells argument
#' @importFrom igraph graph.adjacency
SubSet_cds <- function(cds, cells){
  cells <- unique(cells)
  if(ncol(reducedDimK(cds)) != ncol(cds))
    stop("SubSet_cds doesn't support cds with ncenter run for now. You can try to subset the data and do the construction of trajectory on the subset cds")

  exprs_mat <- as(as.matrix(exprs(cds[, cells])), "sparseMatrix")
  cds_subset <- newCellDataSet(exprs_mat,
                                 phenoData = new("AnnotatedDataFrame", data = pData(cds)[colnames(exprs_mat), ]),
                                 featureData = new("AnnotatedDataFrame", data = fData(cds)),
                                 expressionFamily=negbinomial.size(),
                                 lowerDetectionLimit=1)
  sizeFactors(cds_subset) <- sizeFactors(cds[, cells])
  cds_subset@dispFitInfo <- cds@dispFitInfo

  cds_subset@reducedDimW <- cds@reducedDimW
  cds_subset@reducedDimS <- cds@reducedDimS[, cells]
  cds_subset@reducedDimK <- cds@reducedDimK[, cells]

  cds_subset@cellPairwiseDistances <- cds@cellPairwiseDistances[cells, cells]

  adjusted_K <- Matrix::t(reducedDimK(cds_subset))
  dp <- as.matrix(dist(adjusted_K))
  cellPairwiseDistances(cds_subset) <- dp
  gp <- graph.adjacency(dp, mode = "undirected", weighted = TRUE)
  dp_mst <- minimum.spanning.tree(gp)
  minSpanningTree(cds_subset) <- dp_mst
  cds_subset@dim_reduce_type <- "DDRTree"
  cds_subset <- findNearestPointOnMST(cds_subset)

  cds_subset <- orderCells(cds_subset)
}

# #' Reverse embedding latent graph coordinates back to the high dimension
# #'
# #' @param cds a cell dataset after trajectory reconstruction
# #' @return a new cds containing only the genes used in reducing dimension. Expression values are reverse embedded and rescaled.

reverseEmbeddingCDS <- function(cds) {
  if(nrow(cds@reducedDimW) < 1)
    stop('You need to first apply reduceDimension function on your cds before the reverse embedding')
  
  FM <- normalize_expr_data(cds, norm_method = 'log')
  
  #FM <- FM[unlist(sparseApply(FM, 1, sd, convert_to_dense=TRUE)) > 0, ]
  xm <- Matrix::rowMeans(FM)
  xsd <- sqrt(Matrix::rowMeans((FM - xm)^2))
  FM <- FM[xsd > 0,]
  
  reverse_embedding_data <- reducedDimW(cds) %*% reducedDimS(cds)
  row.names(reverse_embedding_data) <- row.names(FM)
  
  cds_subset <- cds[row.names(FM), ]
  
  #make every value larger than 1: 
  reverse_embedding_data <- t(apply(reverse_embedding_data, 1, function(x) x + abs(min(x))))
  
  #rescale to the original scale: 
  raw_data <- as.matrix(exprs(cds)[row.names(FM), ]) 
  reverse_embedding_data <- reverse_embedding_data * (apply(raw_data, 1, function(x) quantile(x, 0.99)) ) / apply(reverse_embedding_data, 1, max)

  Biobase::exprs(cds_subset) <- reverse_embedding_data
  
  return(cds_subset)
}

# Function to decide a good number of centers for running DDRTree on big datasets
cal_ncenter <- function(ncells, ncells_limit = 100){
  round(2 * ncells_limit * log(ncells)/ (log(ncells) + log(ncells_limit)))
}

#' Select the roots of the principal graph 
#' 
#' 
selectTrajectoryRoots <- function(cds, x=1, y=2, num_roots = NULL, pch = 19, ...)
{
  #TODO: need to validate cds as ready for this plot (need mst, pseudotime, etc)
  lib_info_with_pseudo <- pData(cds)
  
  if (is.null(cds@dim_reduce_type)){
    stop("Error: dimensionality not yet reduced. Please call reduceDimension() before calling this function.")
  }
  
  if (cds@dim_reduce_type == "ICA"){
    reduced_dim_coords <- reducedDimS(cds)
  }else if (cds@dim_reduce_type %in% c("simplePPT", "DDRTree", "SSE", "UMAPSSE", "UMAP", 'L1graph') ){
    reduced_dim_coords <- reducedDimK(cds)
  } else {
    stop("Error: unrecognized dimensionality reduction method.")
  }
  
  ica_space_df <- Matrix::t(reduced_dim_coords) %>%
    as.data.frame() %>%
    select_(prin_graph_dim_1 = x, prin_graph_dim_2 = y) %>%
    mutate(sample_name = rownames(.), sample_state = rownames(.))
  
  dp_mst <- minSpanningTree(cds)
  
  if (is.null(dp_mst)){
    stop("You must first call orderCells() before using this function")
  }
  
  edge_df <- dp_mst %>%
    igraph::as_data_frame() %>%
    select_(source = "from", target = "to") %>%
    left_join(ica_space_df %>% select_(source="sample_name", source_prin_graph_dim_1="prin_graph_dim_1", source_prin_graph_dim_2="prin_graph_dim_2"), by = "source") %>%
    left_join(ica_space_df %>% select_(target="sample_name", target_prin_graph_dim_1="prin_graph_dim_1", target_prin_graph_dim_2="prin_graph_dim_2"), by = "target")
  
  if (is.null(num_roots)){
    num_roots = nrow(ica_space_df)
  }
  #xy <- xy.coords(x, y); x <- xy$x; y <- xy$y
  sel <- rep(FALSE, nrow(ica_space_df))
  
  plot(ica_space_df$prin_graph_dim_1[!sel], ica_space_df$prin_graph_dim_2[!sel]);
  segments(edge_df$source_prin_graph_dim_1, edge_df$source_prin_graph_dim_2, edge_df$target_prin_graph_dim_1, edge_df$target_prin_graph_dim_2)
  
  while(sum(sel) < num_roots) {
    ans <- identify(ica_space_df$prin_graph_dim_1[!sel], ica_space_df$prin_graph_dim_2[!sel], labels = which(!sel), n = 1, ...)
    if(!length(ans)) break
    ans <- which(!sel)[ans]
    points(ica_space_df$prin_graph_dim_1[ans], ica_space_df$prin_graph_dim_2[ans], pch = pch)
    sel[ans] <- TRUE
  }
  ## return indices of selected points
  as.character(ica_space_df$sample_name[which(sel)])
}


#' Run improved force directed layout for cells in low dimensional space.
#'
#' @param cds CellDataSet for the experiment
#' @param separate_group Whether or not to separate participation groups and apply FDL in each partition or directly select equal number of representatives in each louvain groups    
#' @param method The force directed layout method to use  
#' @param merge_coords_method The method used to patch different disconnected component into a single graph  
#' @param start.temp argument passed into layout_with_fr function 
#' @param k Number of nearest neighbors used in calculating the force direct layout 
#' @param cell_num_threshold The minimum number of cells to be ignore for each partition component  
#' @param verbose Wheter to print all running details 
#' @param ... additional arguments passed to functions (louvain_clustering) called by this function. 
#' @return a ggplot2 plot object or a list of all computated information from this function 
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom viridis scale_color_viridis
#' @examples
#' \dontrun{
#' library(HSMMSingleCell)
#' HSMM <- load_HSMM()
#' HSMM <- reduceDimension(HSMM, reduction_method = 'UMAP')
#' HSMM <- smoothEmbedding(HSMM, smooth_method = 'FDL', verbose = T)
#' }
FDL <- function(cds, 
                separate_group = TRUE, 
                method = c('drl', 'fr', 'kk'), 
                merge_coords_method = c('procrutes', 'dla'), 
                start.temp = NULL, 
                color_by,
                markers = NULL,
                cell_size = 1,
                cell_link_size = 1,
                k = 25, 
                cell_num_threshold = 100, 
                verbose = FALSE,
                ...) {
  extra_arguments <- list(...)
  
  # for each louvain cluster, identify equal number of representative cells in each cluster, up to 2000 cells in total 
  if(separate_group == FALSE) {
    data_ori <- t(cds@reducedDimS)
    
    if(ncol(cds) > 2000) {
      cell_cluster <- cds@auxOrderingData$UMAP$louvain_res$optim_res$membership
      cluster_ids <- unique(cell_cluster) 
      landmark_ratio <- 2000 / length(cell_cluster) # determine number of representatives in each louvain cluster 
      
      landmark_id <- c()
      for(current_cluster in cluster_ids) {
        current_cell_ids <- which(cell_cluster == current_cluster)
        data <- cds@auxOrderingData$UMAP$umap_res[current_cell_ids, ] 
        
        cell_num_in_cluster <- round(landmark_ratio * length(current_cell_ids))
        centers <- data[seq(1, nrow(data), length.out=cell_num_in_cluster), ]
        centers <- centers + matrix(rnorm(length(centers), sd = 1e-10), nrow = nrow(centers)) # add random noise 
        kmean_res <- kmeans(data, cell_num_in_cluster, centers=centers, iter.max = 100)
        landmark_id_tmp <- unique(findNearestVertex(t(kmean_res$centers), t(data), process_targets_in_blocks=TRUE))
        
        landmark_id <- c(landmark_id, landmark_id_tmp)
      }
      
      data <- data_ori[landmark_id, ]
    } else {
      data <- t(cds@reducedDimS)
      landmark_id <- 1:ncol(cds)
    }
    
    res <- project_to_representatives(data, 
                                      data_ori, 
                                      landmark_id, 
                                      cds@auxOrderingData$UMAP$louvain_res, 
                                      pd = pData(cds)[landmark_id, ], 
                                      method, #(fr, kk) SSE or other methods? 
                                      start.temp, 
                                      verbose) 
    
    coord <- res$sub_coord_mat
    g <- res$sub_g
    row.names(coord) <- V(g)$name  
    
  } else {
    if(is.null(pData(cds)$louvain_component))
      stop('please first run UMAP or SSE and calculate assign louvain_component for each cell before running this FDL function')
    
    group_stat <- table(pData(cds)$louvain_component)
    cell_num_threshold <- min(max(group_stat) / 2, cell_num_threshold)
    
    valid_groups <- which(group_stat > cell_num_threshold)
    sub_cds_list <- vector('list', length = length(valid_groups))
    sub_g_list <- vector('list', length = length(valid_groups))
    sub_coord_mat_list <- vector('list', length = length(valid_groups))
    
    if(verbose)
      message(paste("total number of valid components is", length(valid_groups)))
    
    for (i in 1:length(valid_groups)){
      if(verbose)
        message(paste("Processing subtrajectory", valid_groups[i]))
      
      t_cds <- cds[, pData(cds)$louvain_component == valid_groups[i]]
      
      if(verbose)
        message(paste("t_cds has cell number:", ncol(t_cds)))
      
      irlba_pca_res <- cds@auxOrderingData$PCA$irlba_pca_res[pData(cds)$louvain_component == valid_groups[i], ]
      umap_res <- cds@auxOrderingData$UMAP$umap_res[pData(cds)$louvain_component == valid_groups[i], ]
      adj_mat <- cds@auxOrderingData$UMAP$adj_mat[pData(cds)$louvain_component == valid_groups[i], pData(cds)$louvain_component == valid_groups[i]]
      
      t_cds@auxOrderingData$PCA$irlba_pca_res <- irlba_pca_res
      t_cds@auxOrderingData$UMAP$umap_res <- umap_res
      t_cds@auxOrderingData$UMAP$adj_mat <- adj_mat
      t_cds@reducedDimS <- t_cds@reducedDimS[, colnames(t_cds)]
      
      sub_cds_list[[i]] <- t_cds
      
      # let us downsample the data 
      if(ncol(t_cds) > 2000) {
        if("landmark_num" %in% names(extra_arguments)) {
          landmark_num <- extra_arguments$landmark_num
        } else {
          landmark_num <- 2000
        }
        
        data_ori <- t(t_cds@reducedDimS)
        
        centers <- data_ori[seq(1, nrow(data_ori), length.out=landmark_num), ]
        kmean_res <- kmeans(data_ori, landmark_num, centers=centers, iter.max = 100)
        landmark_id <- unique(findNearestVertex(t(kmean_res$centers), t(data_ori), process_targets_in_blocks=TRUE))
        
        data <- data_ori[landmark_id, ]
        
      } else {
        data_ori <- t(t_cds@reducedDimS)
        data <- data_ori
        landmark_id <- 1:ncol(t_cds)
      } 
      
      res <- project_to_representatives(data, 
                                        data_ori, 
                                        landmark_id, 
                                        cds@auxOrderingData$UMAP$louvain_res, 
                                        pd = pData(t_cds)[landmark_id, ], 
                                        method, #(fr, kk) SSE or other methods? 
                                        start.temp, 
                                        verbose) 
      
      sub_coord_mat_list[[i]] <- res$sub_coord_mat
      sub_g_list[[i]] <- res$sub_g 
    }
    
    if(merge_coords_method == 'procrutes' & length(valid_groups) > 1) {
      set.seed(2018)
      
      # get high and low dimension landmark data point 
      landmark_num <- 5
      landmark_res <- lapply(1:length(sub_cds_list), function(id) {
        x <- sub_cds_list[[id]]
        landmark_ids <- which(landmark_selection(x, landmark_num = landmark_num)$flag == 1)
        cell_names <- colnames(x)[landmark_ids]
        coords <- sub_coord_mat_list[[id]][landmark_ids, ]
        
        return(list(cell_names = cell_names, coords = coords))
      })
      
      cell_names <- lapply(landmark_res, function(x) x$cell_names)
      high_landmakrs_coord <- t(cds@reducedDimS[, unlist(cell_names)])
      
      # run SSE (MVU version)
      SSE_res <- SSE(high_landmakrs_coord, embeding_dim = 2, C = Inf, knn = 3)
      # SSE_res <- SSE(high_landmakrs_coord, embeding_dim = 2)
      low_landmakrs_coord <- lapply(landmark_res, function(x) x$coords)
      
      mvu_res <- as.data.frame(SSE_res$Y) 
      mvu_res$component <- rep(as.character(1:length(sub_cds_list)), each = landmark_num)
      sd_df <- mvu_res %>% dplyr::group_by(component) %>% dplyr::summarize(sd_d1 = sd(V1), sd_d2 = sd(V2))
      
      # avoid any singular results from SSE 
      if(any(sd_df[, 2:3] < 1e-5)) {
        mvu_res[, 1:2] <- high_landmakrs_coord[, 1:2]
      }
      mean_df <- mvu_res %>% dplyr::group_by(component) %>% dplyr::summarize(mean_d1 = mean(V1), mean_d2 = mean(V2))
      
      # apply rotation, scaling and translation function to all other reduced data point 
      transform_data_by_procrutes <- function(coord, procrutes_res) {
        # ((coord %*% procrutes_res$rotation) * procrutes_res$scale) + matrix(rep(procrutes_res$translation, each = nrow(coord)), nrow = nrow(coord), byrow = F)
        (coord %*% procrutes_res$rotation) + matrix(rep(procrutes_res$translation, each = nrow(coord)), nrow = nrow(coord), byrow = F)
      }
      
      # increase the distance between centroids accordinly to avoid any overlapping 
      centroid_dist <- as.vector(dist(mean_df[, 2:3]))
      
      radius_vec <- unlist(lapply(low_landmakrs_coord, function(x) {
        sqrt((diff(range(x[, 1])) / 2)^2 + (diff(range(x[, 2])) / 2)^2)
      })) 
      
      coord_scale <- min(1000, max(combn(radius_vec, 2, FUN = sum) / centroid_dist))
      
      # procrutes analysis on the landmark points 
      coord_after_procrutes <- lapply(1:length(sub_cds_list), function(x, cd_scale = coord_scale) {
        procrutes_res <- procrustes(subset(mvu_res, component == x)[, 1:2] * cd_scale * 1.1, low_landmakrs_coord[[x]])
        transform_data_by_procrutes(sub_coord_mat_list[[x]], procrutes_res)
      })
      
      coord <- do.call(rbind.data.frame, coord_after_procrutes)
      
    } else {
      coord <- igraph::merge_coords(sub_g_list, sub_coord_mat_list)
    }
    g <- disjoint_union(sub_g_list)
    row.names(coord) <- V(g)$name
  }
  
  cds@reducedDimS <- t(coord)
  cds@reducedDimK <- t(coord)
  
  minSpanningTree(cds) <- g
  
  return(cds)
}

#' apply force-directed layout on the kNN graph built from the downsampled representive cells and then project other 
#' non-representative cells to the low dimensional space with five nearest (on original UMAP space) representive cells' coordinates 
#'
#' @param data (dowmsampled) data to perform Louvain clustering and kNN graph construction 
#' @param data_ori All data points without downsampling. This data space will be used to find the nearest five points for projecting non-landmark points  
#' @param landmark_id The index of cells   
#' @param louvain_res The result of louvain clustering clustering from the original space 
#' @param pd the data frame of phenotype information 
#' @param method The force directed layout function to use  
#' @param start.temp argument passed into layout_with_fr function 
#' @param verbose Wheter to print all running details 
#' @param ... additional arguments passed to functions (louvain_clustering) called by this function. 
project_to_representatives <- function(data, 
                                       data_ori, 
                                       landmark_id, 
                                       louvain_res, 
                                       pd, 
                                       method = 'drl', #(fr, kk) SSE or other methods? 
                                       start.temp = NULL, 
                                       verbose, 
                                       ...) {
  extra_arguments <- list(...)
  
  # build kNN graph 
  if(verbose) 
    message("Running louvain clustering algorithm ...")
  
  louvain_clustering_args <- c(list(data = data, pd = pd, verbose = verbose),
                               extra_arguments[names(extra_arguments) %in% c("k", "weight", "louvain_iter")])
  data_louvain_res <- do.call(louvain_clustering, louvain_clustering_args)
  
  # layout kNN with force direct layout 
  if (method=="fr") coord <- layout_with_fr(data_louvain_res$g, dim=2, coords=data[, 1:2], start.temp=start.temp)
  if (method=="drl") coord <- layout_with_drl(data_louvain_res$g, dim=2, options=list(edge.cut=0))
  if (method=="kk") coord <- layout_with_kk(data_louvain_res$g, dim=2, coords=data[, 1:2])
  
  # this function can be integrated with SSE too? 
  
  sub_g <- igraph::induced_subgraph(louvain_res$g, row.names(data_ori))
  
  # project other cells to the current location 
  if(nrow(data_ori) > 2000) {
    block_size <- 50000
    num_blocks = ceiling(nrow(data_ori) / block_size)
    weight_mat <- NULL
    
    for (j in 1:num_blocks){
      if (j < num_blocks){
        block <- data_ori[((((j-1) * block_size)+1):(j*block_size)), ]
      }else{
        block <- data_ori[((((j-1) * block_size)+1):(nrow(data_ori))), ]
      }
      distances_Z_to_Y <- proxy::dist(block, data_ori[landmark_id, ])
      
      tmp <- as(t(apply(distances_Z_to_Y , 1, function(x) {
        tmp <- sort(x)[6] 
        x[x > tmp] <- 0; p <- rep(0, length(x))
        bandwidth <- mean(range(x[x > 0])) # half of the range of the nearest neighbors as bindwidth 
        p[x > 0] <- exp(-x[x > 0]/bandwidth) # Gaussian kernel 
        p / sum(p) 
      })), 'sparseMatrix')
      
      weight_mat <- rBind(weight_mat, tmp)
    }
    
    projection_res <- as.matrix(weight_mat %*% as(as.matrix(coord), 'sparseMatrix'))
    projection_res[landmark_id, ] <- as.matrix(coord)
    
    coord <- as.matrix(projection_res)
  }
  
  return(list(sub_coord_mat = coord, sub_g = sub_g))
}

# the following functioin is taken from vegan package for performing procrustes analysis 
# Function procrustes rotates a configuration to maximum similarity with another configuration. Function protest tests the non-randomness (significance) between two configurations.
#' @param X Target matrix
#' @param Y Matrix to be rotated.
#' @param scale Allow scaling of axes of Y.
#' @param symmetric Use symmetric Procrustes statistic (the rotation will still be non-symmetric).
procrustes <- function (X, Y, scale = TRUE, symmetric = FALSE) 
{
  # X <- scores(X, display = scores, ...)
  # Y <- scores(Y, display = scores, ...)
  if (nrow(X) != nrow(Y)) 
    stop(gettextf("matrices have different number of rows: %d and %d", 
                  nrow(X), nrow(Y)))
  if (ncol(X) < ncol(Y)) {
    warning("X has fewer axes than Y: X adjusted to comform Y\n")
    addcols <- ncol(Y) - ncol(X)
    for (i in 1:addcols) X <- cbind(X, 0)
  }
  ctrace <- function(MAT) sum(MAT^2)
  c <- 1
  if (symmetric) {
    X <- scale(X, scale = FALSE)
    Y <- scale(Y, scale = FALSE)
    X <- X/sqrt(ctrace(X))
    Y <- Y/sqrt(ctrace(Y))
  }
  xmean <- apply(X, 2, mean)
  ymean <- apply(Y, 2, mean)
  if (!symmetric) {
    X <- scale(X, scale = FALSE)
    Y <- scale(Y, scale = FALSE)
  }
  XY <- crossprod(X, Y)
  sol <- svd(XY)
  A <- sol$v %*% t(sol$u)
  if (scale) {
    c <- sum(sol$d)/ctrace(Y)
  }
  Yrot <- c * Y %*% A
  b <- xmean - c * ymean %*% A
  R2 <- ctrace(X) + c * c * ctrace(Y) - 2 * c * sum(sol$d)
  reslt <- list(Yrot = Yrot, X = X, ss = R2, rotation = A, 
                translation = b, scale = c, xmean = xmean, symmetric = symmetric, 
                call = match.call())
  reslt$svd <- sol
  class(reslt) <- "procrustes"
  reslt
}

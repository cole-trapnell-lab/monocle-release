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

  if(any(is.na(E(pr_graph)$weight))) {
    E(pr_graph)$weight <- 1
  }  
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
  
  cc_ordering <- extract_general_graph_ordering(cds, root_pr_nodes)
  closest_vertex = cds@auxOrderingData[[cds@dim_reduce_type]]$pr_graph_cell_proj_closest_vertex
  pData(cds)$Pseudotime = cc_ordering[closest_vertex[row.names(pData(cds)),],]$pseudo_time
  cds@auxOrderingData[[cds@dim_reduce_type]]$root_pr_nodes <- root_pr_nodes

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

#' project a CellDataSet object into a lower dimensional PCA (or ISI) space after normalize the data 
#'
#' @description For most analysis (including trajectory inference, clustering) in Monocle 3, it requires us to to start from a 
#' low dimensional PCA space. preprocessCDS will be used to first project a CellDataSet object into a lower dimensional PCA space 
#' before we apply clustering with community detection algorithm or other non-linear dimension reduction method, for example 
#' UMAP, tSNE, DDRTree, L1-graph, etc.  While tSNE is especially suitable for visualizing clustering results, comparing
#' to UMAP, the global distance in tSNE space is not meaningful. UMAP can either be used for visualizing clustering result or as a general 
#' non-linear dimension reduction method. 
#' SimplePPT, DDRTree and L1-graph are two complementary trajectory inference method where the first one is very great at learning a tree structure 
#' but the later is general and can learn any arbitrary graph structure. Both methods can be applied to the UMAP space.   
#'
#' @details 
#' In Monocle 3, we overhauled the code from Monocle2 so that a standard Monocle 3 workingflow works as following: 
#' 1. run \code{preprocessCDS} to project a CellDataSet object into a lower dimensional PCA space after 
#' normalize the data 
#' 2. run \code{reduceDimension} to further project the PCA space into much lower dimension space with non-linear 
#' dimension reduction techniques, including tSNE, UMAP. 
#' 3. run \code{smoothEmbedding} (optional) to smooth noisy embedding from 2 to facilitate visualization and learning 
#' of the graph structure.
#' 4. run \code{partitionCells} to partition cells into different graphs based on a similar approach proposed by Alex Wolf and colleagues. 
#' We then reconstruct the trajectory in each partition with the \code{learnGraph} function. 
#' 5. run \code{learnGraph} to reconstruct developmental trajectory with reversed graph embedding algorithms. In monocle 3, we enabled the 
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
#' it converts the (sparse) expression matrix into tf-idf (term-frequency-inverse document frequency 
#' which increases proportionally to the gene expression value appears in the cell and is offset by the frequency of 
#' the gene in the entire dataset, which helps to adjust for the fact that some gene appear more frequently across cell in general.) matrix and then performs a 
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
preprocessCDS <- function(cds, method = c('PCA', 'LSI', 'none'), #, 'LSI' , 'NMF'
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
    
    if(use_tf_idf == TRUE) {
      FM <- as(FM, "dgCMatrix")
      cds_dfm <- new("dfmSparse", FM)
      cds_dfm <- dfm_tfidf(cds_dfm)
      FM <- sparseMatrix(i = cds_dfm@i, p = cds_dfm@p, x = cds_dfm@x, dimnames = cds_dfm@Dimnames, dims = cds_dfm@Dim, index1 = F)
    }
    
    irlba_res <- sparse_prcomp_irlba(t(FM), n = min(num_dim, min(dim(FM)) - 1),
                                     center = scaling, scale. = scaling)
    irlba_pca_res <- irlba_res$x
    # reducedDimA(cds) <- t(irlba_pca_res) # get top 50 PCs, which can be used for louvain clustering later 
  } else if(method == 'LSI') {
    FM <- as(FM, "dgCMatrix")
    cds_dfm <- new("dfmSparse", FM)
    cds_dfm <- dfm_tfidf(cds_dfm)
    cds_dfm_lsa <- textmodel_lsa(cds_dfm, nd = num_dim, margin = c("both"))
    irlba_pca_res <- cds_dfm_lsa$features
    
  } else if(method == 'none') {
    irlba_pca_res <- t(FM)
  } else {
    stop('unknown preprocessing method, stop!')
  }
  row.names(irlba_pca_res) <- colnames(cds)
  cds@normalized_data_projection <- irlba_pca_res
  
  cds
}

#' Compute a projection of a CellDataSet object into a lower dimensional space with non-linear dimension reduction methods 
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
#' @references DDRTree: Qi Mao, Li Wang, Steve Goodison, and Yijun Sun. Dimensionality reduction via graph structure learning. In Proceedings of the 21th ACM SIGKDD International Conference on Knowledge Discovery and Data Mining, pages 765–774. ACM, 2015.
#' @references UMAP: McInnes, L, Healy, J, UMAP: Uniform Manifold Approximation and Projection for Dimension Reduction, ArXiv e-prints 1802.03426, 2018
#' @references tSNE: Laurens van der Maaten and Geoffrey Hinton. Visualizing data using t-SNE. J. Mach. Learn. Res., 9(Nov):2579– 2605, 2008.
#' @export
reduceDimension <- function(cds,
                            max_components=2,
                            reduction_method=c("DDRTree", "ICA", 'tSNE', "UMAP", 'none'),
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
    else if (reduction_method == c("DDRTree")) {
      
      message('DDRTree will be eventually deprecated in reduceDimension call and be used in RGE function instead. We are calling RGE for you now.')
      cds@reducedDimS <- t(cds@normalized_data_projection)
      cds <- partitionCells(cds)
      cds <- learnGraph(cds, RGE_method = 'DDRTree', do_partition = F, ...)
      
    }else if (reduction_method == c("UMAP") ) {  
      if (verbose)
        message("Running Uniform Manifold Approximation and Projection")
      
      umap_args <- c(list(X = irlba_pca_res, log = F, n_component = as.integer(max_components), verbose = verbose, return_all = T),
                     extra_arguments[names(extra_arguments) %in% 
                                       c("python_home", "n_neighbors", "metric", "n_epochs", "negative_sample_rate", "learning_rate", "init", "min_dist", "spread", 
                                         'set_op_mix_ratio', 'local_connectivity', 'repulsion_strength', 'a', 'b', 'random_state', 'metric_kwds', 'angular_rp_forest', 'verbose')])
      tmp <- do.call(UMAP, umap_args)
      tmp$embedding_ <- (tmp$embedding_ - min(tmp$embedding_)) / max(tmp$embedding_) # normalize UMAP space
      umap_res <- tmp$embedding_; 
      
      adj_mat <- Matrix::sparseMatrix(i = tmp$graph_$indices, p = tmp$graph_$indptr, 
                                      x = -as.numeric(tmp$graph_$data), dims = c(ncol(cds), ncol(cds)), index1 = F, 
                                      dimnames = list(colnames(cds), colnames(cds)))
      
      S <- t(umap_res)
      
      Y <- S
      W <- t(irlba_pca_res)
      
      minSpanningTree(cds) <- graph_from_adjacency_matrix(adj_mat, weighted=TRUE)
      
      A <- S
      colnames(A) <- colnames(FM)
      reducedDimA(cds) <- A
      
      colnames(S) <- colnames(FM)
      colnames(Y) <- colnames(FM)
      reducedDimW(cds) <- W 
      reducedDimS(cds) <- as.matrix(Y)
      reducedDimK(cds) <- S
      
      #cds@auxOrderingData$UMAP <- list(umap_res = umap_res, adj_mat = adj_mat)
      cds@dim_reduce_type <- reduction_method
    } else if(reduction_method == 'none') {
      irlba_pca_res <- t(irlba_pca_res)
      colnames(irlba_pca_res) <- colnames(FM)
      reducedDimS(cds) <- irlba_pca_res
      reducedDimK(cds) <- irlba_pca_res
    }else {
      stop("Error: unrecognized dimensionality reduction method")
    }
  }
  cds
}

#' This function tries to learn a smooth embedding from the noisy reduced dimension using different techniques 
#' @description The function relies on smooth skeleton learning or force direct layout function to learn a smoothier
#' representation of the data. It can be used to facilitate visualization of the data or downstream graph learning 
#' 
#' @param cds CellDataSet for the experiment
#' @param max_components the dimensionality of the reduced space
#' @param do_partition Whether or not to separate participation groups and apply FDL in each partition or directly select equal number of representatives in each louvain groups.     
#' @param use_pca Whether or not to cluster cells based on top PCA component. Default to be FALSE. 
#' @param method The embedding smoothing technique to use, including force directed layout (which includes drl, fr, kk three methods), our new method PSL and SSE  
#' @param merge_coords_method The method used to patch different smoothed embedding from disconnected component into a coordinate system
#' @param start.temp argument passed into layout_with_fr function 
#' @param k Number of nearest neighbors used in calculating the force direct layout 
#' @param landmark_num Number of landmark used for downsampling the data. Currently landmarks are defined as the medoids from kmean clustering.  
#' @param cell_num_threshold The minimum number of cells in each partition component required for running embedding smooth  
#' @param verbose Wheter to print all running details 
#' @param ... additional arguments passed to functions (louvain_clustering) called by this function. 
#' @return a cds with smoothed embedding for each disconnected component learned (stored in reduceDimS, reduceDimK) and the corresponding graph stored in minSpanningTree
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom viridis scale_color_viridis
#' @references PSL: Li Wang, Qi Mao (2018). Probabilistic Dimensionality Reduction via Structure Learning. IEEE Transactions on Pattern Analysis and Machine Intelligence
#' @references SSE: Li Wang, Qi Mao, Ivor W. Tsang (2017). Latent Smooth Skeleton Embedding. Proceedings of the 31th AAAI Conference on Artificial Intelligence. 2017.
#' @seealso \code{\link[monocle]{patchEmbedding}}
#' @export 
smoothEmbedding <- function(cds,
                           max_components = 2, 
                           do_partition = FALSE, 
                           use_pca = FALSE, 
                           method = c('PSL', 'drl', 'fr', 'kk', 'SSE'), 
                           merge_coords_method = c('dla', 'procrutes'), 
                           start.temp = NULL, 
                           k = 20, 
                           landmark_num = 2000, 
                           cell_num_threshold = 0, 
                           verbose = FALSE,
                            ...){
  cds <- patchEmbedding(cds = cds, 
                        max_components = max_components, 
                        do_partition = do_partition, 
                        use_pca = use_pca, 
                        method = method, 
                        merge_coords_method = merge_coords_method, 
                        start.temp = start.temp,
                        k = k, 
                        landmark_num = landmark_num,
                        cell_num_threshold = cell_num_threshold, 
                        verbose = verbose, 
                        ...)

  cds
}

#' This function tries to partition cells into different graphs based on a similar approach proposed by Alex Wolf and colleagues 
#' @description Recently Alex Wolf and colleague first proposed the idea to represent the data with an “abstract partition graph”
#' of clusters identified by Louvain clustering by simply connecting significantly overlapping Louvain clusters (Wolf et al. 2017). 
#' Similar methods for “abstract partition graph” are also recently developed and applied in analyzing the zebrafish / frog cell 
#' atlas datasets (Wagner et al. 2018; Briggs et al. 2018). This coarse-graining representation of the data address a few limitations 
#' of tree-based trajectory inference algorithms, for example, the default principal tree learning algorithm (DDRTree) in Monocle 2. 
#' Although the particion graph doesn’t learn an explicit simplified principal tree as DDRTree in Monocle 2, it can naturally separate 
#' outlier cell groups and potentially also parallel trajectories while DDRTree often requires pre-processing before hand to robustly 
#' reconstruct trajectory and cannot handle non-tree like structure. Instead of directly learn a coarse-graining graph of clusters, 
#' we instead take advantage of the participation graph and use it as merely a heuristic initial condition for L1-graph algorithm 
#' (as explained in the next section) to learn the principal points and principal graph at the same time directly from the reduced 
#' UMAP data space. In contrast to the cluster participation method, the principal graph learnt provides an abstraction of the data 
#' manifold while also preserves the local information from the original data space as it is directly embedded in the original data space. 
#' In Monocle 3, we uses the clustering_louvain function from the igraph package to perform community detection and implemented an efficient
#' version of “abstract partition graph” from Alex Wolf. Basically, we first create a design matrix $$X$$ representing the allocation of 
#' each cell to a particular louvain cluster. The column of $$X$$ represents a louvain cluster while the row of $$X$$ a particular cell. 
#' $$X_{ij} = 1$$ if cell $$i$$ belongs to cluster $$j$$, otherwise 0. We can further obtain the adjacency matrix $$A$$ of the kNN graph 
#' used to perform the louvain clustering where $$A_{ij} = 1$$ if cell $$i$$ connects to $$j$$ in the kNN graph. Then the connection 
#' matrix $$M$$ between each cluster is calculated as, $$M = X‘ x A x X$$. Once $$M$$ is constructed, we can then follow 
#' Supplemental Note 3.1 from (Wolf et al. 2017) to calculate the significance of the connection between each louvain clustering and 
#' consider any clusters with p-value larger than 0.05 by default as not disconnected. 

#' 
#' @param cds the CellDataSet upon which to perform this operation
#' @param k number of nearest neighbors used for Louvain clustering (pass to louvain_clustering function)
#' @param weight whether or not to calculate the weight for each edge in the kNN graph (pass to louvain_clustering function)
#' @param louvain_iter the number of iteraction for louvain clustering (pass to louvain_clustering function)
#' @param resolution resolution of clustering result, specifiying the granularity of clusters. 
#' Default to not use resolution and the standard igraph louvain clustering algorithm will be used. 
#' @param louvain_qval The q-val threshold used to determine the partition of cells (pass to compute_louvain_connected_components)
#' @param return_all Whether to return all saved objects from compute_louvain_connected_components function. 
#' @param verbose Whether to emit verbose output during louvain clustering
#' @param ... additional arguments to pass to the smoothEmbedding function
#' @return an updated CellDataSet object
#' @export
partitionCells <- function(cds,
                           k = 20, 
                           weight = F, 
                           louvain_iter = 1, 
                           resolution = NULL,
                           louvain_qval = 0.05, 
                           return_all = FALSE, 
                           verbose = FALSE, ...){
  extra_arguments <- list(...)
  FM <- cds@auxOrderingData$normalize_expr_data
  irlba_pca_res <- cds@normalized_data_projection
  
  Y <- reducedDimS(cds)
  reduced_dim_res = Y 
  
  if(verbose)
    message("Running louvain clustering algorithm ...")
  #row.names(umap_res) <- colnames(FM)
  if(nrow(Y) == 0) {
    reduced_dim_res <- t(irlba_pca_res)
  }
  louvain_clustering_args <- c(list(data = t(reduced_dim_res), pd = pData(cds)[colnames(FM), ], k = k, 
                                    resolution = resolution, weight = weight, louvain_iter = louvain_iter, verbose = verbose)) # , extra_arguments[names(extra_arguments) %in% c("k", "weight", "louvain_iter")]
  louvain_res <- do.call(louvain_clustering, louvain_clustering_args)
  
  if(length(unique(louvain_res$optim_res$membership)) == 1) {
    pData(cds)$louvain_component <- 1
    cds@auxClusteringData$partitionCells <- louvain_res
    
    return(cds)
  }
  
  cluster_graph_res <- compute_louvain_connected_components(louvain_res$g, louvain_res$optim_res, louvain_qval, verbose)
  louvain_component = components(cluster_graph_res$cluster_g)$membership[louvain_res$optim_res$membership]
  names(louvain_component) = colnames(FM)
  louvain_component = as.factor(louvain_component)
  pData(cds)$louvain_component <- louvain_component
  
  cds@auxClusteringData$partitionCells <- louvain_res
  
  if(return_all) {
    return(list(cds = cds, cluster_graph_res = cluster_graph_res))
  } else {
    return(cds)
  }
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
#' @param partition_group When this argument is set to TRUE (default to be FALSE), we will learn a tree structure for each separate over-connected louvain component. 
#' @param do_partition When this argument is set to TRUE (default to be FALSE), we will learn a tree structure for each separate over-connected louvain component. 
#' @param scale When this argument is set to TRUE (default), it will scale each gene before running trajectory reconstruction.
#' @param euclidean_distance_ratio The maximal ratio between the euclidean distance of two tip nodes in the spanning tree inferred from SimplePPT algorithm and 
#' that of the maximum distance between any connecting points on the spanning tree allowed to be connected during the loop closure procedure .   
#' @param geodestic_distance_ratio  The minimal ratio between the geodestic distance of two tip nodes in the spanning tree inferred from SimplePPT algorithm and 
#' that of the length of the diameter path on the spanning tree allowed to be connected during the loop closure procedure. (Both euclidean_distance_ratio and geodestic_distance_ratio 
#' need to be satisfied to introduce the edge for loop closure.)    
#' @param close_loop Whether or not to perform an additional run of loop closing after running DDRTree or SimplePPT to identify potential loop structure in the data space
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
#' @references DDRTree: Qi Mao, Li Wang, Steve Goodison, and Yijun Sun. Dimensionality reduction via graph structure learning. In Proceedings of the 21th ACM SIGKDD International Conference on Knowledge Discovery and Data Mining, pages 765–774. ACM, 2015.
#' @references L1graph (generalized SimplePPT): Qi Mao, Li Wang, Ivor Tsang, and Yijun Sun. Principal graph and structure learning based on reversed graph embedding . IEEE Trans. Pattern Anal. Mach. Intell., 5 December 2016.
#' @references Original SimplePPT: Qi Mao, Le Yang, Li Wang, Steve Goodison, Yijun Sun. SimplePPT: A Simple Principal Tree Algorithm https://epubs.siam.org/doi/10.1137/1.9781611974010.89
#' @export
learnGraph <- function(cds,
                       max_components=2,
                       RGE_method = c('SimplePPT', 'L1graph', 'DDRTree'), 
                       auto_param_selection = TRUE, 
                       partition_group = 'louvain_component', 
                       do_partition = TRUE, 
                       scale = FALSE, 
                       close_loop = FALSE, 
                       euclidean_distance_ratio = 1, 
                       geodestic_distance_ratio = 1/3, 
                       verbose = FALSE, 
                       ...){
  RGE_method <- RGE_method[1]
  extra_arguments <- list(...)
  FM <- cds@auxOrderingData$normalize_expr_data
  irlba_pca_res <- cds@normalized_data_projection
  
  Y <- reducedDimS(cds)
  reduced_dim_res = Y 
  
  # 
  if(do_partition && !(partition_group %in% colnames(pData(cds))))
    stop('Please make sure the partition_group you want to partition the dataset based on is included in the pData of the cds!')
  
  if(length(unique(pData(cds)[, partition_group])) <= 1) {
    do_partition <- FALSE 
  }
  louvain_res <- cds@auxClusteringData$partitionCells
  
  if(is.null(louvain_res))
    stop('Please run partitionCells function before run learnGraph!')
  
  louvain_module_length = length(unique(sort(louvain_res$optim_res$membership)))
  louvain_component <- pData(cds)$louvain_component
  names(louvain_component) <- colnames(cds)
  
  if(RGE_method == 'L1graph') { 
    # FIXME: This case is broken, because I didn't have time to update the landmark
    # stuff during the refactor.
    if(cds@dim_reduce_type == "UMAP") {
      
    } else {
      stop('L1graph can be only applied to the UMAP space, please first call reduceDimension() using UMAP!')
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
    
    if (ncenter > ncol(Y))
      stop("Error: ncenters must be less than or equal to ncol(X)")
    
    centers <- reduced_dim_res[,seq(1, ncol(reduced_dim_res), length.out=ncenter)]
    #centers <- centers + matrix(rnorm(length(centers), sd = 1e-10), nrow = nrow(centers)) # add random noise 
    
    kmean_res <- kmeans(t(reduced_dim_res), ncenter, centers=t(centers), iter.max = 100)
    if (kmean_res$ifault != 0){
      message(paste("Warning: kmeans returned ifault =", kmean_res$ifault))
    }
    nearest_center = findNearestVertex(t(kmean_res$centers), reduced_dim_res, process_targets_in_blocks=TRUE)
    medioids = reduced_dim_res[,unique(nearest_center)]
    reduced_dim_res <- medioids
    
    if(verbose)
      message('running L1-graph ...')
    
    #X <- t(reduced_dim_res)
    
    if('C0' %in% names(extra_arguments)){
      C0 <- extra_arguments$C0
    }
    else
      C0 <- reduced_dim_res
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
    
    # cluster_graph_res <- compute_louvain_connected_components(louvain_res$g, louvain_res$optim_res, louvain_qval, verbose)
    # louvain_component = components(cluster_graph_res$cluster_g)$membership[louvain_res$optim_res$membership]
    # cds@auxOrderingData[["L1graph"]]$louvain_component = louvain_component
    # names(louvain_component) <- colnames(cds)
    louvain_component_for_medioids <- louvain_component[colnames(reduced_dim_res)]
    #louvain_component_for_medioids <- as.factor(louvain_component_for_medioids)
    if (do_partition && length(levels(louvain_component_for_medioids)) > 1){
      louvain_component_mask = as.matrix(tcrossprod(sparse.model.matrix( ~ louvain_component_for_medioids + 0)))
      
      G = G * louvain_component_mask
      W = W * louvain_component_mask
      rownames(G) = rownames(W)
      colnames(G) = colnames(W)
    }
    
    l1graph_args <- c(list(X = reduced_dim_res, C0 = C0, G = G, gstruct = 'l1-graph', verbose = verbose),
                      extra_arguments[names(extra_arguments) %in% c('maxiter', 'eps', 'L1.lambda', 'L1.gamma', 'L1.sigma', 'nn')])
    
    
    l1_graph_res <- do.call(principal_graph, l1graph_args)
    
    colnames(l1_graph_res$C) <-  colnames(reduced_dim_res)
    #DCs <- reduced_dim_res #FM
    
    colnames(l1_graph_res$W) <- colnames(reduced_dim_res)
    rownames(l1_graph_res$W) <- colnames(reduced_dim_res)
    
    
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
    #cds@dim_reduce_type <- "L1graph"
    cds@dim_reduce_type <- RGE_method
    cds <- findNearestPointOnMST(cds)
  } else if(RGE_method == 'SimplePPT') {
    if(ncol(cds@reducedDimS) > 1) {
      irlba_pca_res <- t(cds@reducedDimS)
    }
    
    #louvain_component <- pData(cds)[, partition_group]
    if(do_partition && length(louvain_component) == ncol(cds)) {
      multi_tree_DDRTree_res <- multi_component_RGE(cds, scale = scale, 
        RGE_method = RGE_method, 
        partition_group = partition_group, 
        irlba_pca_res = irlba_pca_res, 
        max_components = max_components, 
        extra_arguments = extra_arguments, 
        close_loop = close_loop, 
        euclidean_distance_ratio = euclidean_distance_ratio, 
        geodestic_distance_ratio = geodestic_distance_ratio, 
        verbose = verbose)
      
      ddrtree_res_W <- multi_tree_DDRTree_res$ddrtree_res_W
      ddrtree_res_Z <- multi_tree_DDRTree_res$ddrtree_res_Z
      ddrtree_res_Y <- multi_tree_DDRTree_res$ddrtree_res_Y
      cds <- multi_tree_DDRTree_res$cds
      dp_mst <- multi_tree_DDRTree_res$dp_mst
    } else {  ## need to change the following to SimplePPT soon! 
      ncenter <- NULL
      if(auto_param_selection & ncol(cds) >= 100) {
        if("ncenter" %in% names(extra_arguments)) #avoid overwrite the ncenter parameter
          ncenter <- extra_arguments$ncenter
        else
          ncenter <- cal_ncenter(nrow(irlba_pca_res))
        
      } 
      
      if(scale) {
        X <- as.matrix(scale(t(irlba_pca_res)))
      }
      else {
        X <- t(irlba_pca_res)
      }
      
      ddr_args <- c(list(X=X, dimensions=ncol(X), ncenter=ncenter, no_reduction = T, verbose = verbose),
                    extra_arguments[names(extra_arguments) %in% c("initial_method", "maxIter", "sigma", "lambda", "param.gamma", "tol")])
      
      ddrtree_res <- do.call(DDRTree, ddr_args)
      
      # ddrtree_res <- DDRTree(as.matrix(scale(t(irlba_pca_res))), max_components, no_reduction = T, verbose = verbose, ...)
      if(ncol(ddrtree_res$Y) == ncol(cds))
        colnames(ddrtree_res$Y) <- colnames(FM) #paste("Y_", 1:ncol(ddrtree_res$Y), sep = "")
      else
        colnames(ddrtree_res$Y) <- paste("Y_", 1:ncol(ddrtree_res$Y), sep = "")
      
      colnames(ddrtree_res$Z) <- colnames(FM)
      
      ddrtree_res_W <- ddrtree_res$W
      ddrtree_res_Z <- ddrtree_res$Z
      ddrtree_res_Y <- ddrtree_res$Y
      
      adjusted_K <- t(ddrtree_res_Y)
      dp <- as.matrix(dist(adjusted_K))
      cellPairwiseDistances(cds) <- dp
      gp <- graph.adjacency(dp, mode = "undirected", weighted = TRUE)
      dp_mst <- minimum.spanning.tree(gp)
      
      ddrtree_res$stree <- ddrtree_res$stree[1:ncol(ddrtree_res$Y), 1:ncol(ddrtree_res$Y)]
      dimnames(ddrtree_res$stree) <- list(paste("Y_", 1:ncol(ddrtree_res$Y), sep = ""), paste("Y_", 1:ncol(ddrtree_res$Y), sep = ""))
      row.names(ddrtree_res$R) <- colnames(cds); colnames(ddrtree_res$R) <- paste("Y_", 1:ncol(ddrtree_res$Y), sep = "")
      colnames(ddrtree_res$Q) <- colnames(cds)
      cds@auxOrderingData[["SimplePPT"]] <- ddrtree_res[c('stree', 'Q', 'R', 'objective_vals', 'history')]
      
      if(ncol(cds) < 100) { 
        cds <- findNearestPointOnMST(cds)
      } else {
        tmp <- matrix(apply(ddrtree_res$R, 1, which.max))
        row.names(tmp) <- colnames(cds)
        cds@auxOrderingData[["SimplePPT"]]$pr_graph_cell_proj_closest_vertex <- tmp
      }
    }
    
    reducedDimW(cds) <- ddrtree_res_W
    reducedDimS(cds) <- ddrtree_res_Z
    reducedDimK(cds) <- ddrtree_res_Y
    
    minSpanningTree(cds) <- dp_mst
    
    #cds@dim_reduce_type <- "SimplePPT"
    
  } else if(RGE_method == 'L1_SimplePPT') {

  } else if(RGE_method == 'DDRTree') {
    if(ncol(cds@reducedDimS) > 1) {
      irlba_pca_res <- t(cds@reducedDimS)
    }
    
    row.names(irlba_pca_res) <- colnames(FM)
    
    if (verbose)
      message("Learning principal graph with DDRTree")
    
    # TODO: DDRTree should really work with sparse matrices.
    #louvain_component <- pData(cds)$louvain_component
    
    if(do_partition && length(louvain_component) == ncol(cds)) {
      X <- t(irlba_pca_res)
      
      reducedDimK_coord <- NULL  
      dp_mst <- NULL 
      pr_graph_cell_proj_closest_vertex <- NULL 
      cell_name_vec <- NULL
      
      multi_tree_DDRTree_res <- multi_component_RGE(cds, scale = scale, 
        RGE_method = RGE_method, 
        partition_group = partition_group, 
        irlba_pca_res = irlba_pca_res, 
        max_components = max_components, 
        extra_arguments = extra_arguments, 
        close_loop = close_loop, 
        euclidean_distance_ratio = euclidean_distance_ratio, 
        geodestic_distance_ratio = geodestic_distance_ratio, 
        verbose = verbose)

      ddrtree_res_W <- multi_tree_DDRTree_res$ddrtree_res_W
      ddrtree_res_Z <- multi_tree_DDRTree_res$ddrtree_res_Z
      ddrtree_res_Y <- multi_tree_DDRTree_res$ddrtree_res_Y
      cds <- multi_tree_DDRTree_res$cds
      dp_mst <- multi_tree_DDRTree_res$dp_mst
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
      
      adjusted_K <- t(ddrtree_res_Y)
      dp <- as.matrix(dist(adjusted_K))
      cellPairwiseDistances(cds) <- dp
      gp <- graph.adjacency(dp, mode = "undirected", weighted = TRUE)
      dp_mst <- minimum.spanning.tree(gp)
      
      ddrtree_res$stree <- ddrtree_res$stree[1:ncol(ddrtree_res$Y), 1:ncol(ddrtree_res$Y)]
      dimnames(ddrtree_res$stree) <- list(paste("Y_", 1:ncol(ddrtree_res$Y), sep = ""), paste("Y_", 1:ncol(ddrtree_res$Y), sep = ""))
      row.names(ddrtree_res$R) <- colnames(cds); colnames(ddrtree_res$R) <- paste("Y_", 1:ncol(ddrtree_res$Y), sep = "")
      colnames(ddrtree_res$Q) <- colnames(cds)

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
    
    #cds@dim_reduce_type <- "DDRTree"
    
  }
  
  cds@dim_reduce_type <- RGE_method
  
  cds 
}

#' Finds the nearest principal graph node
#' @param data_matrix the input matrix
#' @param target_points the target points
#' @param block_size the number of input matrix rows to process per bloclk
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
  cds_subset <- subset_cds(cds, subset_cell)

  root_state <- pData(cds_subset[, starting_cell])[, 'State']
  cds_subset <- orderCells(cds_subset, root_state = as.numeric(root_state))

  return(cds_subset)
}

#' Subset a cds which only includes cells provided with the argument cells
#'
#' @param cds a cell dataset after trajectory reconstruction
#' @param cells a vector contains all the cells you want to subset
#' @return a new cds containing only the cells from the cells argument
#' @importFrom igraph graph.adjacency
#' @export
#' @examples 
#' \dontrun{
#' lung <- load_lung()
#' tmp <- subset_cds(lung, cells = row.names(subset(pData(lung), State == 1)))
#' plot_cell_trajectory(tmp)
#' }
subset_cds <- function(cds, cells){
  cells <- unique(intersect(cells, colnames(cds)))
  if(length(cells) == 0) {
    stop("Cannot find any cell from the cds matches with the cell name from the cells argument! Please make sure the cell name you input is correct.")
  }
  
  exprs_mat <- exprs(cds[, cells])
  cds_subset <- newCellDataSet(exprs_mat,
                                 phenoData = new("AnnotatedDataFrame", data = pData(cds)[cells, ]),
                                 featureData = new("AnnotatedDataFrame", data = fData(cds)),
                                 lowerDetectionLimit=cds@lowerDetectionLimit, 
                                 expressionFamily=cds@expressionFamily)
  
  cds_subset@dispFitInfo <- cds@dispFitInfo
  
  if(ncol(cds@reducedDimS) == ncol(cds)) {
    cds_subset@reducedDimS <- cds@reducedDimS[, cells]
  } else {
    cds_subset@reducedDimS <- cds@reducedDimS
  }
  if(ncol(cds@reducedDimW) == ncol(cds)) {
    cds_subset@reducedDimW <- cds@reducedDimW[, cells]
  } else {
    cds_subset@reducedDimW <- cds@reducedDimW
  }
  if(ncol(cds@reducedDimA) == ncol(cds)) {
    cds_subset@reducedDimA <- cds@reducedDimA[, cells]
  } else {
    cds_subset@reducedDimA <- cds@reducedDimA
  }

  if(nrow(cds@normalized_data_projection) == ncol(cds)) {
    cds_subset@normalized_data_projection <- cds@normalized_data_projection[cells, ]
  } else {
    cds_subset@normalized_data_projection <- cds@normalized_data_projection
  }
  
  # we may also subset results from any RGE methods 
  if('DDRTree' %in% names(cds@auxOrderingData)) {
    cell_ids <- cds@auxOrderingData$DDRTree$pr_graph_cell_proj_closest_vertex[cells, 1]
    cds_subset@auxOrderingData$DDRTree$stree <- cds@auxOrderingData$DDRTree$stree[cell_ids, cell_ids]
    cds_subset@auxOrderingData$DDRTree$R <- cds@auxOrderingData$DDRTree$R[cells, cell_ids]
    cds_subset@auxOrderingData$DDRTree$pr_graph_cell_proj_closest_vertex <- cds@auxOrderingData$DDRTree$pr_graph_cell_proj_closest_vertex[cells, , drop = F]
  }
  if('SimplePPT' %in% names(cds@auxOrderingData)) {
    cell_ids <- cds@auxOrderingData$SimplePPT$pr_graph_cell_proj_closest_vertex[cells, 1]
    cds_subset@auxOrderingData$SimplePPT$stree <- cds@auxOrderingData$SimplePPT$stree[cell_ids, cell_ids]
    cds_subset@auxOrderingData$SimplePPT$R <- cds@auxOrderingData$SimplePPT$R[cells, cell_ids]
    cds_subset@auxOrderingData$SimplePPT$pr_graph_cell_proj_closest_vertex <- cds@auxOrderingData$SimplePPT$pr_graph_cell_proj_closest_vertex[cells, , drop = F]
  }
  if('L1Graph' %in% names(cds@auxOrderingData)) {
    cell_ids <- cds@auxOrderingData$L1graph$pr_graph_cell_proj_closest_vertex[cells, 1]
    cds_subset@auxOrderingData$L1graph$stree <- cds@auxOrderingData$L1graph$stree[cell_ids, cell_ids]
    cds_subset@auxOrderingData$L1graph$R <- cds@auxOrderingData$L1graph$R[cells, cell_ids]
    cds_subset@auxOrderingData$L1graph$pr_graph_cell_proj_closest_vertex <- cds@auxOrderingData$L1graph$pr_graph_cell_proj_closest_vertex[cells, , drop = F]
  }
  
  # find the corresponding principal graph nodes for those selected cells, followed by subseting the trajectories 
  principal_graph_points <- cds@auxOrderingData[[cds@dim_reduce_type]]$pr_graph_cell_proj_closest_vertex[cells, 1]
  cds_subset@minSpanningTree <- induced_subgraph(cds@minSpanningTree, paste0('Y_', principal_graph_points))
  
  cds_subset@reducedDimK <- cds@reducedDimK[, principal_graph_points]
  cds_subset@dim_reduce_type <- cds@dim_reduce_type
  cds_subset 
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
  if(ncells <= ncells_limit) {
    return(NULL)
  }
  
  round(2 * ncells_limit * log(ncells)/ (log(ncells) + log(ncells_limit)))
}

#' Select the roots of the principal graph 
#' @param cds CellDataSet where roots will be selected from
#' @param x The first dimension to plot 
#' @param y The number of dimension to plot 
#' @param num_roots Number of roots for the trajectory 
#' @param pch Size of the principal graph node
#' 
selectTrajectoryRoots <- function(cds, x=1, y=2, num_roots = NULL, pch = 19, ...)
{
  #TODO: need to validate cds as ready for this plot (need mst, pseudotime, etc)
  lib_info_with_pseudo <- pData(cds)
  
  if (is.null(cds@dim_reduce_type)){
    stop("Error: dimensionality not yet reduced. Please call reduceDimension() before calling this function.")
  }
  
  if (cds@dim_reduce_type %in% c("SimplePPT", "DDRTree", "UMAP") ){
    reduced_dim_coords <- reducedDimK(cds)
  } else {
    stop("Error: unrecognized dimensionality reduction method.")
  }
  
  ica_space_df <- Matrix::t(reduced_dim_coords) %>%
    as.data.frame() 
  use_3d = ncol(ica_space_df) >= 3
  if (use_3d){
    colnames(ica_space_df) = c("prin_graph_dim_1", "prin_graph_dim_2", "prin_graph_dim_3")
  }
  else{
    colnames(ica_space_df) = c("prin_graph_dim_1", "prin_graph_dim_2")
  }
  ica_space_df = ica_space_df %>% mutate(sample_name = rownames(.), sample_state = rownames(.))
  
  dp_mst <- minSpanningTree(cds)
  
  if (is.null(dp_mst)){
    stop("You must first call orderCells() before using this function")
  }
  
  if (use_3d){
    edge_df <- dp_mst %>%
      igraph::as_data_frame() %>%
      select_(source = "from", target = "to") %>%
      left_join(ica_space_df %>% select_(source="sample_name", source_prin_graph_dim_1="prin_graph_dim_1", source_prin_graph_dim_2="prin_graph_dim_2", source_prin_graph_dim_3="prin_graph_dim_3"), by = "source") %>%
      left_join(ica_space_df %>% select_(target="sample_name", target_prin_graph_dim_1="prin_graph_dim_1", target_prin_graph_dim_2="prin_graph_dim_2", target_prin_graph_dim_3="prin_graph_dim_3"), by = "target")
  }else{
    edge_df <- dp_mst %>%
      igraph::as_data_frame() %>%
      select_(source = "from", target = "to") %>%
      left_join(ica_space_df %>% select_(source="sample_name", source_prin_graph_dim_1="prin_graph_dim_1", source_prin_graph_dim_2="prin_graph_dim_2"), by = "source") %>%
      left_join(ica_space_df %>% select_(target="sample_name", target_prin_graph_dim_1="prin_graph_dim_1", target_prin_graph_dim_2="prin_graph_dim_2"), by = "target")
  }

  if (is.null(num_roots)){
    num_roots = nrow(ica_space_df)
  }
  #xy <- xy.coords(x, y); x <- xy$x; y <- xy$y
  sel <- rep(FALSE, nrow(ica_space_df))
  
  if (use_3d){
    open3d(windowRect=c(0,0,1024,1024))
    segments3d(matrix(as.matrix(t(edge_df[,c(3,4,5,6,7,8)])), ncol=3, byrow=T), lwd=2, 
               col="black",
               line_antialias=TRUE)
    points3d(Matrix::t(reduced_dim_coords[1:3,]), col="black")
    while(sum(sel) < num_roots) {
      ans <- identify3d(Matrix::t(reduced_dim_coords[1:3,!sel]), labels = which(!sel), n = 1, buttons = c("left", "right"), ...)  
      if(!length(ans)) break
      ans <- which(!sel)[ans]
      #points3d(Matrix::t(reduced_dim_coords[1:3,ans]), col="red")
      sel[ans] <- TRUE
    }
  }else{
    plot(ica_space_df$prin_graph_dim_1[!sel], ica_space_df$prin_graph_dim_2[!sel]);
    segments(edge_df$source_prin_graph_dim_1, edge_df$source_prin_graph_dim_2, edge_df$target_prin_graph_dim_1, edge_df$target_prin_graph_dim_2)
  
    while(sum(sel) < num_roots) {
      ans <- identify(ica_space_df$prin_graph_dim_1[!sel], ica_space_df$prin_graph_dim_2[!sel], labels = which(!sel), n = 1, ...)
      if(!length(ans)) break
      ans <- which(!sel)[ans]
      points(ica_space_df$prin_graph_dim_1[ans], ica_space_df$prin_graph_dim_2[ans], pch = pch)
      sel[ans] <- TRUE
    }
  }
  ## return indices of selected points
  as.character(ica_space_df$sample_name[which(sel)])
}

#' 
#' #' the following functioin is used to learn trajectory on each disjointed components 
#' #' @param cds CellDataSet  The CellDataSet upon which to perform this operation
#' #' @param scale A logical argument to determine whether or not we should scale the data before constructing trajectory (default to be FALSE)
#' #' @param RGE_method The method for reversed graph embedding  
#' #' @param partition_group The column name in the pData used to partition cells 
#' #' @param irlba_pca_res The matrix for PCA top components (retrieved with irlba by default)
#' #' @param max_components Number of maximum component 
#' #' @param extra_arguments Extra arguments passed into learnGraph (which calls this function) 
#' #' @param close_loop A logical argument to determine whether or not we should close loop for the trajectory we learned (default to be FALSE)
#' #' @param verbose Whether to emit verbose output when running this function 
#' multi_tree_DDRTree <- function(cds, scale = FALSE, RGE_method, partition_group = 'louvain_component', irlba_pca_res, max_components, extra_arguments, close_loop = FALSE, verbose = FALSE) {
#'   louvain_component <- pData(cds)[, partition_group]
#'   
#'   X <- t(irlba_pca_res)
#'   
#'   reducedDimK_coord <- NULL  
#'   dp_mst <- NULL 
#'   pr_graph_cell_proj_closest_vertex <- NULL 
#'   cell_name_vec <- NULL
#'   
#'   for(cur_comp in unique(louvain_component)) {
#'     X_subset <- X[, louvain_component == cur_comp]
#'     
#'     #add other parameters...
#'     if(scale) 
#'       X_subset <- t(as.matrix(scale(t(X_subset))))
#'     
#'     ncenter <- cal_ncenter(ncol(X_subset))
#'     
#'     ddr_args <- c(list(X=X_subset, dimensions=max_components, ncenter=ncenter, no_reduction = T, verbose = verbose),
#'                   extra_arguments[names(extra_arguments) %in% c("initial_method", "maxIter", "sigma", "lambda", "param.gamma", "tol")])
#'     #browser()
#'     ddrtree_res <- do.call(DDRTree, ddr_args)
#'     
#'     if(is.null(reducedDimK_coord)) {
#'       curr_cell_names <- paste("Y_", 1:ncol(ddrtree_res$Y), sep = "")
#'       pr_graph_cell_proj_closest_vertex <- matrix(apply(ddrtree_res$R, 1, which.max))
#'       cell_name_vec <- colnames(X_subset)
#'     } else {
#'       curr_cell_names <- paste("Y_", ncol(reducedDimK_coord) + 1:ncol(ddrtree_res$Y), sep = "")
#'       pr_graph_cell_proj_closest_vertex <- rbind(pr_graph_cell_proj_closest_vertex, matrix(apply(ddrtree_res$R, 1, which.max) + ncol(reducedDimK_coord)))
#'       cell_name_vec <- c(cell_name_vec, colnames(X_subset))
#'     }
#'     
#'     curr_reducedDimK_coord <- ddrtree_res$Y
#'     
#'     dp <- as.matrix(dist(t(curr_reducedDimK_coord))) #ddrtree_res$stree[1:ncol(ddrtree_res$Y), 1:ncol(ddrtree_res$Y)]
#'     dimnames(dp) <- list(curr_cell_names, curr_cell_names)
#'     
#'     cur_dp_mst <- mst(graph.adjacency(dp, mode = "undirected", weighted = TRUE))
#'     
#'     tmp <- matrix(apply(ddrtree_res$R, 1, which.max))
#'     
#'     if(length(close_loop) == length(unique(louvain_component)))
#'       curr_close_loop <- close_loop[which(unique(louvain_component) %in% cur_comp)]
#'     else 
#'       curr_close_loop <- close_loop[1]
#'     
#'     if(curr_close_loop == TRUE) {
#'       colnames(curr_reducedDimK_coord) <- curr_cell_names
#'       connectTips_res <- connectTips(pData(cds)[louvain_component == cur_comp, ], ddrtree_res$R, cur_dp_mst, 
#'                                      curr_reducedDimK_coord, cds@reducedDimS[, louvain_component == cur_comp])
#'       
#'       curr_reducedDimK_coord <- connectTips_res$reducedDimK_df
#'       cur_dp_mst <- connectTips_res$mst_g
#'       
#'       current_W <- as.matrix(get.adjacency(cur_dp_mst))
#'       message('current_W dim is ', nrow(current_W), ncol(current_W))
#'       l1graph_args <- c(list(X = curr_reducedDimK_coord, C0 = curr_reducedDimK_coord, G = current_W, gstruct = 'l1-graph', verbose = verbose),
#'                         extra_arguments[names(extra_arguments) %in% c('maxiter', 'eps', 'L1.lambda', 'L1.gamma', 'L1.sigma', 'nn')])
#'       
#'       l1_graph_res <- do.call(principal_graph, l1graph_args)
#'       
#'       W <- l1_graph_res$W
#'       message('W dim is ', nrow(W), ncol(W))
#'       message('nrow(ddrtree_res$R) is ', nrow(ddrtree_res$R))
#'       start_id <- min(nrow(ddrtree_res$R), nrow(current_W))
#'       
#'       dimnames(W) <- list(V(cur_dp_mst)$name, V(cur_dp_mst)$name)
#'       current_W[, start_id:nrow(current_W)] <- W[, start_id:nrow(current_W)]
#'       current_W[start_id:nrow(current_W), ] <- W[start_id:nrow(current_W), ]
#'       
#'       # W[W < 1e-5] <- 0
#'       cur_dp_mst <- graph.adjacency(current_W, mode = "undirected", weighted = TRUE)
#'     }
#'     
#'     dp_mst <- graph.union(dp_mst, cur_dp_mst)
#'     reducedDimK_coord <- cbind(reducedDimK_coord, curr_reducedDimK_coord)
#'     
#'   }
#'   
#'   row.names(pr_graph_cell_proj_closest_vertex) <- cell_name_vec
#'   
#'   ddrtree_res_W <- ddrtree_res$W
#'   ddrtree_res_Z <- cds@reducedDimS
#'   ddrtree_res_Y <- reducedDimK_coord
#'   
#'   cds@auxOrderingData[[RGE_method]] <- ddrtree_res[c('stree', 'Q', 'R', 'objective_vals', 'history')]
#'   cds@auxOrderingData[[RGE_method]]$pr_graph_cell_proj_closest_vertex <- pr_graph_cell_proj_closest_vertex
#'   
#'   colnames(ddrtree_res_Y) <- paste0("Y_", 1:ncol(ddrtree_res_Y), sep = "")
#'   
#'   return(list(cds = cds, 
#'               ddrtree_res_W = ddrtree_res_W, 
#'               ddrtree_res_Z = ddrtree_res_Z, 
#'               ddrtree_res_Y = ddrtree_res_Y, 
#'               dp_mst = dp_mst))
#' }

multi_component_RGE <- function(cds, 
                                scale = FALSE, 
                                RGE_method, 
                                partition_group = 'louvain_component', 
                                irlba_pca_res, 
                                max_components, 
                                extra_arguments, 
                                close_loop = FALSE, 
                                euclidean_distance_ratio = 1, 
                                geodestic_distance_ratio = 1/3, 
                                verbose = FALSE) {
  louvain_component <- pData(cds)[, partition_group]
  
  X <- t(irlba_pca_res)
  
  reducedDimK_coord <- NULL  
  dp_mst <- NULL 
  pr_graph_cell_proj_closest_vertex <- NULL 
  cell_name_vec <- NULL
  
  merge_rge_res <- NULL
  max_ncenter <- 0
  for(cur_comp in unique(louvain_component)) {
    X_subset <- X[, louvain_component == cur_comp]
    if(verbose) message('Current louvain_component is ', cur_comp)
    if(ncol(X_subset) < 10) {
      message('Louvain component with less than 10 cells will be ignored!')
      next; 
    }
    #add other parameters...
    if(scale) {
      X_subset <- t(as.matrix(scale(t(X_subset))))
    }

    if(!("ncenter" %in% names(extra_arguments))) {
      ncenter <- cal_ncenter(ncol(X_subset))
      if(is.null(ncenter)) {
        ncenter <- ncol(X_subset) - 1
      }
    } else {
      ncenter <- min(ncol(X_subset) - 1, extra_arguments$ncenter)
    }
    
    kmean_res <- NULL 
    if(RGE_method == 'DDRTree') {
      ddr_args <- c(list(X=X_subset, dimensions=max_components, ncenter=ncenter, no_reduction = T, verbose = verbose),
                    extra_arguments[names(extra_arguments) %in% c("initial_method", "maxIter", "sigma", "lambda", "param.gamma", "tol")])
      #browser()
      rge_res <- do.call(DDRTree, ddr_args)
      medioids <- rge_res$Y
      stree <- rge_res$stree[1:ncol(medioids), 1:ncol(medioids)]
      
      if(!close_loop) {
        if(is.null(merge_rge_res)) {
          colnames(rge_res$Y) <- paste0('Y_', 1:ncol(rge_res$Y))
          merge_rge_res <- rge_res
          colnames(merge_rge_res$X) <- colnames(X_subset)
          colnames(merge_rge_res$Z) <- colnames(X_subset)
          colnames(merge_rge_res$Q) <- colnames(X_subset)
          row.names(merge_rge_res$R) <- colnames(X_subset); colnames(merge_rge_res$R) <- paste0('Y_', 1:ncol(merge_rge_res$Y))
          merge_rge_res$R <- list(merge_rge_res$R)
          merge_rge_res$stree <- list(stree)
          merge_rge_res$objective_vals <- list(merge_rge_res$objective_vals)
        } else {
          colnames(rge_res$X) <- colnames(X_subset)
          # assign R column names first
          row.names(rge_res$R) <- colnames(X_subset); colnames(rge_res$R) <- paste0('Y_', (ncol(merge_rge_res$Y) + 1):(ncol(merge_rge_res$Y) + ncol(rge_res$Y)), sep = "")
          colnames(rge_res$Y) <- paste("Y_", (ncol(merge_rge_res$Y) + 1):(ncol(merge_rge_res$Y) + ncol(rge_res$Y)), sep = "")
          merge_rge_res$Y <- cbind(merge_rge_res$Y, rge_res$Y)
          colnames(rge_res$Z) <- colnames(X_subset)
          merge_rge_res$R <- c(merge_rge_res$R, list(rge_res$R))
          colnames(rge_res$Q) <- colnames(X_subset)
          rge_res$Q <- cbind(merge_rge_res$Q, rge_res$Q) 
          merge_rge_res$stree <- c(merge_rge_res$stree, list(stree))
          merge_rge_res$objective_vals <- c(merge_rge_res$objective_vals, list(rge_res$objective_vals))   
        }
      }
    } else if(RGE_method == 'SimplePPT') {
      centers <- t(X_subset)[seq(1, ncol(X_subset), length.out=ncenter), ]
      centers <- centers + matrix(rnorm(length(centers), sd = 1e-10), nrow = nrow(centers)) # add random noise 
      
      kmean_res <- tryCatch({
        kmeans(t(X_subset), ncenter, centers=centers, iter.max = 100)
        }, error = function() {
          kmeans(t(X_subset), ncenter, iter.max = 100)  
        })
      
      if (kmean_res$ifault != 0){
        message(paste("Warning: kmeans returned ifault =", kmean_res$ifault))
      }
      nearest_center <- findNearestVertex(t(kmean_res$centers), X_subset, process_targets_in_blocks=TRUE)
      medioids <- X_subset[, unique(nearest_center)]
      reduced_dim_res <- t(medioids)
      k <- 25
      mat <- t(X_subset)
      if (is.null(k)) {
        k <- round(sqrt(nrow(mat))/2)
        k <- max(10, k)
      }
      if (verbose)
        message("Finding kNN using RANN with ", k, " neighbors")
      dx <- RANN::nn2(mat, k = min(k, nrow(mat) - 1))
      nn.index <- dx$nn.idx[, -1]
      nn.dist <- dx$nn.dists[, -1]

      if (verbose) 
        message("Calculating the local density for each sample based on kNNs ...")

      rho <- exp(-rowMeans(nn.dist))
      mat_df <- as.data.frame(mat)
      tmp <- mat_df %>% dplyr::add_rownames() %>% dplyr::mutate(cluster = kmean_res$cluster, density = rho) %>% dplyr::group_by(cluster) %>% dplyr::top_n(n = 1, wt = density) %>% arrange(-desc(cluster))
      medioids <- X_subset[, tmp$rowname] # select representative cells by highest density
      
      reduced_dim_res <- t(medioids)
      l1graph_args <- c(list(X = X_subset, C0 = medioids, G = NULL, gstruct = 'span-tree', verbose = verbose),
                        extra_arguments[names(extra_arguments) %in% c('maxiter', 'eps', 'L1.lambda', 'L1.gamma', 'L1.sigma', 'nn')])
      
      rge_res <- do.call(principal_graph, l1graph_args)
      
      names(rge_res)[c(2, 4, 5)] <- c('Y', 'R','objective_vals')
      stree <- rge_res$W
      
      if(!close_loop) {
        if(is.null(merge_rge_res)) {
          colnames(rge_res$Y) <- paste0('Y_', 1:ncol(rge_res$Y))
          merge_rge_res <- rge_res
          colnames(merge_rge_res$X) <- colnames(X_subset)
          # colnames(merge_rge_res$Z) <- colnames(X_subset)
          # colnames(merge_rge_res$Q) <- colnames(X_subset) # no Q for principal_graph
          row.names(merge_rge_res$R) <- colnames(X_subset); colnames(merge_rge_res$R) <- paste0('Y_', 1:ncol(merge_rge_res$Y))
          merge_rge_res$R <- list(merge_rge_res$R)
          merge_rge_res$stree <- list(stree)
          merge_rge_res$objective_vals <- list(merge_rge_res$objective_vals)
        } else {
          colnames(rge_res$X) <- colnames(X_subset)
          row.names(rge_res$R) <- colnames(X_subset); colnames(rge_res$R) <- paste0('Y_', (ncol(merge_rge_res$Y) + 1):(ncol(merge_rge_res$Y) + ncol(rge_res$Y)), sep = "")
          colnames(rge_res$Y) <- paste("Y_", (ncol(merge_rge_res$Y) + 1):(ncol(merge_rge_res$Y) + ncol(rge_res$Y)), sep = "")
          merge_rge_res$Y <- cbind(merge_rge_res$Y, rge_res$Y)
          # colnames(rge_res$Z) <- colnames(X_subset)
          merge_rge_res$R <- c(merge_rge_res$R, list(rge_res$R))
          # colnames(rge_res$Q) <- colnames(X_subset)
          # rge_res$Q <- cbind(merge_rge_res$Q, rge_res$Q) 
          merge_rge_res$stree <- c(merge_rge_res$stree, list(stree))
          merge_rge_res$objective_vals <- c(merge_rge_res$objective_vals, list(rge_res$objective_vals))  
        }
      }
    }
    
    if(close_loop) {
      connectTips_res <- connectTips(pData(cds)[louvain_component == cur_comp, ], 
                                           R = rge_res$R, 
                                           stree = stree, 
                                           reducedDimK_old = rge_res$Y, 
                                           reducedDimS_old = cds@reducedDimS[, louvain_component == cur_comp],
                                           kmean_res = kmean_res, 
                                           euclidean_distance_ratio = 1, 
                                           geodestic_distance_ratio = 1/3, 
                                           verbose = verbose)
      stree <- connectTips_res$stree    
      # use louvain clustering method to get a better initial graph? 
      G <- connectTips_res$G

      l1graph_args <- c(list(X = X_subset, G = G, C0 = medioids, stree = as.matrix(stree), gstruct = 'l1-graph', verbose = verbose),
                        extra_arguments[names(extra_arguments) %in% c('eps', 'L1.lambda', 'L1.gamma', 'L1.sigma', 'nn', "maxiter")])
      
      rge_res <- do.call(principal_graph, l1graph_args)
      names(rge_res)[c(2, 4, 5)] <- c('Y', 'R','objective_vals')
      stree <- as(rge_res$W, 'sparseMatrix')

      if(is.null(merge_rge_res)) {
        colnames(rge_res$Y) <- paste0('Y_', 1:ncol(rge_res$Y))
        merge_rge_res <- rge_res
        colnames(merge_rge_res$X) <- colnames(X_subset)
        # colnames(merge_rge_res$Z) <- colnames(X_subset)
        # colnames(merge_rge_res$Q) <- colnames(X_subset) # no Q for principal_graph
        row.names(merge_rge_res$R) <- colnames(X_subset); colnames(merge_rge_res$R) <- paste0('Y_', 1:ncol(merge_rge_res$Y))
        merge_rge_res$R <- list(merge_rge_res$R)
        merge_rge_res$stree <- list(stree)
        merge_rge_res$objective_vals <- list(merge_rge_res$objective_vals)
      } else {
        colnames(rge_res$X) <- colnames(X_subset)
        row.names(rge_res$R) <- colnames(X_subset); colnames(rge_res$R) <- paste0('Y_', (ncol(merge_rge_res$Y) + 1):(ncol(merge_rge_res$Y) + ncol(rge_res$Y)), sep = "")
        colnames(rge_res$Y) <- paste("Y_", (ncol(merge_rge_res$Y) + 1):(ncol(merge_rge_res$Y) + ncol(rge_res$Y)), sep = "")
        merge_rge_res$Y <- cbind(merge_rge_res$Y, rge_res$Y)
        # colnames(rge_res$Z) <- colnames(X_subset)
        merge_rge_res$R <- c(merge_rge_res$R, list(rge_res$R))
        # colnames(rge_res$Q) <- colnames(X_subset)
        # rge_res$Q <- cbind(merge_rge_res$Q, rge_res$Q) 
        merge_rge_res$stree <- c(merge_rge_res$stree, list(stree))
        merge_rge_res$objective_vals <- c(merge_rge_res$objective_vals, list(rge_res$objective_vals))  
      }
    }
    
    if(is.null(reducedDimK_coord)) {
      curr_cell_names <- paste("Y_", 1:ncol(rge_res$Y), sep = "")
      pr_graph_cell_proj_closest_vertex <- matrix(apply(rge_res$R, 1, which.max))
      cell_name_vec <- colnames(X_subset)
    } else {
      curr_cell_names <- paste("Y_", (ncol(reducedDimK_coord) + 1):(ncol(reducedDimK_coord) + ncol(rge_res$Y)), sep = "")
      pr_graph_cell_proj_closest_vertex <- rbind(pr_graph_cell_proj_closest_vertex, matrix(apply(rge_res$R, 1, which.max) + ncol(reducedDimK_coord)))
      cell_name_vec <- c(cell_name_vec, colnames(X_subset))
    }
    
    curr_reducedDimK_coord <- rge_res$Y
    
    dimnames(stree) <- list(curr_cell_names, curr_cell_names)
    cur_dp_mst <- graph.adjacency(stree, mode = "undirected", weighted = TRUE)
    
    dp_mst <- graph.union(dp_mst, cur_dp_mst)
    reducedDimK_coord <- cbind(reducedDimK_coord, curr_reducedDimK_coord)
  }
  
  row.names(pr_graph_cell_proj_closest_vertex) <- cell_name_vec
  
  ddrtree_res_W <- as.matrix(rge_res$W)
  ddrtree_res_Z <- cds@reducedDimS
  ddrtree_res_Y <- reducedDimK_coord # ensure the order of column names matches that of the original name ids 
  # correctly set up R, stree -- the mapping from each cell to the principal graph points 
  R <- sparseMatrix(i = 1, j = 1, x = 0, dims = c(ncol(cds), ncol(merge_rge_res$Y))) # use sparse matrix for large datasets 
  stree <- sparseMatrix(i = 1, j = 1, x = 0, dims = c(ncol(merge_rge_res$Y), ncol(merge_rge_res$Y)))
  curr_row_id <- 1
  curr_col_id <- 1
  R_row_names <- NULL 
  for(i in 1:length(merge_rge_res$R)) {
    current_R <- merge_rge_res$R[[i]]
    # R[curr_row_id:(curr_row_id + nrow(current_R) - 1), curr_col_id:(curr_col_id + ncol(current_R) - 1)] <- current_R # this is why learnGraph is very slow ...
    
    stree[curr_col_id:(curr_col_id + ncol(current_R) - 1), curr_col_id:(curr_col_id + ncol(current_R) - 1)] <- merge_rge_res$stree[[i]]

    curr_row_id <- curr_row_id + nrow(current_R)
    curr_col_id <- curr_col_id + ncol(current_R)
    R_row_names <- c(R_row_names, row.names(current_R))
  }

  row.names(R) <- R_row_names
  R <- R[colnames(cds), ] # reorder the colnames 
  pr_graph_cell_proj_closest_vertex <- slam::rowapply_simple_triplet_matrix(slam::as.simple_triplet_matrix(R), function(x) {
    which.max(x)
  })

  cds@auxOrderingData[[RGE_method]] <- list(stree = stree, Q = merge_rge_res$Q, R = R, objective_vals = merge_rge_res$objective_vals, history = merge_rge_res$history) # rge_res[c('stree', 'Q', 'R', 'objective_vals', 'history')] # 
  cds@auxOrderingData[[RGE_method]]$pr_graph_cell_proj_closest_vertex <- as.data.frame(pr_graph_cell_proj_closest_vertex)[colnames(cds), , drop = F] # Ensure the row order matches up that of the column order of the cds 
  
  colnames(ddrtree_res_Y) <- paste0("Y_", 1:ncol(ddrtree_res_Y), sep = "")
  
  return(list(cds = cds, 
              ddrtree_res_W = ddrtree_res_W, 
              ddrtree_res_Z = ddrtree_res_Z, 
              ddrtree_res_Y = ddrtree_res_Y, 
              dp_mst = dp_mst))
}

#' functions to connect the tip points after learning the DDRTree or simplePPT tree
connectTips <- function(pd,
                        R, # kmean cluster
                        stree, 
                        reducedDimK_old, 
                        reducedDimS_old, 
                        k = 25, 
                        weight = F,
                        qval_thresh = 0.05, 
                        kmean_res, 
                        euclidean_distance_ratio = 1, 
                        geodestic_distance_ratio = 1/3, 
                        verbose = FALSE,
                        ...) {
  mst_g_old <- igraph::graph_from_adjacency_matrix(stree, mode = 'undirected')
  if(is.null(kmean_res)) {
    tmp <- matrix(apply(R, 1, which.max))
    
    row.names(tmp) <- colnames(reducedDimS_old)
    
    tip_pc_points <- which(igraph::degree(mst_g_old) == 1)

    data <- t(reducedDimS_old[, ])
    
    louvain_res <- louvain_clustering(data, pd[, ], k = k, weight = weight, verbose = verbose)
    
    # louvain_res$optim_res$memberships[1, ] <-  tmp[raw_data_tip_pc_points, 1]
    louvain_res$optim_res$membership <- tmp[, 1]
  } else { # use kmean clustering result 
    tip_pc_points <- which(igraph::degree(mst_g_old) == 1)
    tip_pc_points_kmean_clusters <- sort(kmean_res$cluster[names(tip_pc_points)])
    # raw_data_tip_pc_points <- which(kmean_res$cluster %in% tip_pc_points_kmean_clusters) # raw_data_tip_pc_points <- which(kmean_res$cluster %in% tip_pc_points)
    
    data <- t(reducedDimS_old[, ]) # raw_data_tip_pc_points
    
    louvain_res <- louvain_clustering(data, pd[row.names(data), ], k = k, weight = weight, verbose = verbose)
    
    # louvain_res$optim_res$memberships[4, ] <-  kmean_res$cluster #[raw_data_tip_pc_points]
    louvain_res$optim_res$membership <- kmean_res$cluster #[raw_data_tip_pc_points]   
  }
  
  # identify edges between only tip cells  
  cluster_graph_res <- compute_louvain_connected_components(louvain_res$g, louvain_res$optim_res, qval_thresh=qval_thresh, verbose = verbose)
  dimnames(cluster_graph_res$cluster_mat) <- dimnames(cluster_graph_res$num_links)
  valid_connection <- which(cluster_graph_res$cluster_mat < qval_thresh, arr.ind = T) 
  valid_connection <- valid_connection[apply(valid_connection, 1, function(x) all(x %in% tip_pc_points)), ] # only the tip cells 

  # prepare the PAGA graph 
  G <- cluster_graph_res$cluster_mat
  G[cluster_graph_res$cluster_mat < qval_thresh] <- -1
  G[cluster_graph_res$cluster_mat > 0] <- 0
  G <- - G
  
  if(nrow(valid_connection) == 0) {
    return(list(stree = igraph::get.adjacency(mst_g_old), Y = reducedDimK_old, G = G))
  }
  
  # calculate length of the MST diameter path   
  mst_g <- mst_g_old
  diameter_dis <- igraph::diameter(mst_g_old)
  reducedDimK_df <- reducedDimK_old
  
  # find the maximum distance between nodes from the MST   
  res <- dist(t(reducedDimK_old))
  g <- igraph::graph_from_adjacency_matrix(as.matrix(res), weighted = T, mode = 'undirected')
  mst <- igraph::minimum.spanning.tree(g)
  max_node_dist <- max(igraph::E(mst)$weight)

  # append new edges to close loops in the spanning tree returned from SimplePPT   
  for(i in 1:nrow(valid_connection)) {
    # cluster id for the tip point; if kmean_res return valid_connection[i, ] is itself; otherwise the id identified in the tmp file 
    edge_vec <- sort(unique(louvain_res$optim_res$membership))[valid_connection[i, ]]
    edge_vec_in_tip_pc_point <- igraph::V(mst_g_old)$name[edge_vec]
    
    if(length(edge_vec_in_tip_pc_point) == 1) next; 
    
    if(all(edge_vec %in% tip_pc_points) & (igraph::distances(mst_g_old, edge_vec_in_tip_pc_point[1], edge_vec_in_tip_pc_point[2]) > geodestic_distance_ratio * diameter_dis) & 
                                           (euclidean_distance_ratio * max_node_dist > dist(t(reducedDimK_old[, edge_vec]))) ) {
      if(verbose) message('edge_vec is ', edge_vec[1], '\t', edge_vec[2])
      if(verbose) message('edge_vec_in_tip_pc_point is ', edge_vec_in_tip_pc_point[1], '\t', edge_vec_in_tip_pc_point[2])
      
      mst_g <- igraph::add_edges(mst_g, edge_vec_in_tip_pc_point) 
    }
  }

  list(stree = igraph::get.adjacency(mst_g), Y = reducedDimK_df, G = G)
}

#' Run embedding smoothing techniques (drl, PSL and SSE) for each isolated group and then put different embedding in the same coordinate system 
#'
#' @param cds CellDataSet for the experiment
#' @param max_components the dimensionality of the reduced space
#' @param do_partition Whether or not to separate participation groups and apply FDL in each partition or directly select equal number of representatives in each louvain groups    
#' @param use_pca Whether or not to cluster cells based on top PCA component. Default to be FALSE. 
#' @param method The embedding smoothing technique to use, including force directed layout (which includes drl, fr, kk three methods), our new method PSL and SSE  
#' @param merge_coords_method The method used to patch different smoothed embedding from disconnected component into a coordinate system
#' @param start.temp argument passed into layout_with_fr function 
#' @param k Number of nearest neighbors used in calculating the force direct layout 
#' @param landmark_num Number of landmark used for downsampling the data. Currently landmarks are defined as the medoids from kmean clustering.  
#' @param cell_num_threshold The minimum number of cells in each partition component required for running embedding smooth  
#' @param verbose Wheter to print all running details 
#' @param ... additional arguments passed to functions (louvain_clustering) called by this function. 
#' @return a cds with smoothed embedding for each disconnected component learned (stored in reduceDimS, reduceDimK) and the corresponding graph stored in minSpanningTree
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom viridis scale_color_viridis
#' @references PSL: Li Wang, Qi Mao (2018). Probabilistic Dimensionality Reduction via Structure Learning. IEEE Transactions on Pattern Analysis and Machine Intelligence
#' @references SSE: Li Wang, Qi Mao, Ivor W. Tsang (2017). Latent Smooth Skeleton Embedding. Proceedings of the 31th AAAI Conference on Artificial Intelligence. 2017.
#' @export
patchEmbedding <- function(cds, 
                           max_components = 2, 
                           do_partition = FALSE, 
                           use_pca = FALSE, 
                           method = c('PSL', 'drl', 'fr', 'kk', 'SSE'), 
                           merge_coords_method = c('dla', 'procrutes'), 
                           start.temp = NULL, 
                           k = 20, 
                           landmark_num = 2000, 
                           cell_num_threshold = 0, 
                           verbose = FALSE,
                           ...) {
  extra_arguments <- list(...)
  
  if(method == 'SSE' & do_partition == FALSE) {
    message('Note that if your data includes separate groups, you should set do_partition to be TRUE when using SSE method!')
  }

  if(any(nrow(cds@normalized_data_projection) == 0 | nrow(cds@reducedDimS) == 0 | is.null(pData(cds)$louvain_component))) {
    stop('Please first run preprocessCDS, reduceDimension, partitionCells (in order) before running this function!')
  }
  # for each louvain cluster, identify equal number of representative cells in each cluster, up to 2000 cells in total 
  if(do_partition == FALSE) {
    if(use_pca) {
      data_ori <- cds@normalized_data_projection
    } else {
      data_ori <- t(cds@reducedDimS)
    }

    if(ncol(cds) > landmark_num) {
      cell_cluster <- cds$louvain_component
      cluster_ids <- unique(cell_cluster) 
      landmark_ratio <- landmark_num / ncol(cds) # determine number of representatives in each louvain cluster 
      
      landmark_id <- c()
      for(current_cluster in cluster_ids) {
        current_cell_ids <- which(cell_cluster == current_cluster)
        data <- data_ori[current_cell_ids, ] 
        
        cell_num_in_cluster <- round(landmark_ratio * length(current_cell_ids))
        if(cell_num_in_cluster == 0) cell_num_in_cluster <- 1  # avoid the case where cell_num_in_cluster = 0
        centers <- data[seq(1, nrow(data), length.out=cell_num_in_cluster), , drop = F] # avoid the case where cell_num_in_cluster = 0
        centers <- centers + matrix(rnorm(length(centers), sd = 1e-10), nrow = nrow(centers)) # add random noise 
        kmean_res <- kmeans(data, cell_num_in_cluster, centers=centers, iter.max = 100)
        landmark_id_tmp <- unique(findNearestVertex(t(kmean_res$centers), t(data), process_targets_in_blocks=TRUE))
        
        landmark_id <- c(landmark_id, current_cell_ids[landmark_id_tmp]) # landmark index from the original data 
      }
      
      data <- data_ori[landmark_id, ]
    } else {
      data <- data_ori
      landmark_id <- 1:ncol(cds)
    }
    
    res <- project_to_representatives(data, 
                                      data_ori, 
                                      landmark_id, 
                                      cds@auxClusteringData$partitionCells, 
                                      pd = pData(cds)[landmark_id, ], 
                                      method, 
                                      start.temp, 
                                      k, 
                                      do_partition = do_partition, 
                                      max_components = max_components, 
                                      verbose,
                                      ...) 
    
    coord <- res$sub_coord_mat
    row.names(coord) <- colnames(cds) 
    g <- res$sub_g   
  } else {    
    group_stat <- table(pData(cds)$louvain_component)
    cell_num_threshold <- min(max(group_stat) / 2, cell_num_threshold)
    
    valid_groups <- which(group_stat > cell_num_threshold)
    sub_cds_list <- vector('list', length = length(valid_groups))
    sub_g_list <- vector('list', length = length(valid_groups))
    sub_coord_mat_list <- vector('list', length = length(valid_groups))
    
    if(verbose) {
      message(paste("total number of valid components is", length(valid_groups)))
    }
    
    for (i in 1:length(valid_groups)){
      if(verbose) {
        message(paste("Processing subtrajectory", valid_groups[i]))
      }
      
      t_cds <- cds[, pData(cds)$louvain_component == valid_groups[i]]
      
      if(verbose) {
        message(paste("t_cds has cell number:", ncol(t_cds)))
      }
      
      irlba_pca_res <- cds@normalized_data_projection[pData(cds)$louvain_component == valid_groups[i], ]
      # umap_res <- cds@reducedDimS[, pData(cds)$louvain_component == valid_groups[i]]
      # adj_mat <- cds@minSpanningTree[pData(cds)$louvain_component == valid_groups[i], pData(cds)$louvain_component == valid_groups[i]]
      
      t_cds@normalized_data_projection <- irlba_pca_res
      # t_cds@minSpanningTree <- adj_mat
      t_cds@reducedDimS <- t_cds@reducedDimS[, colnames(t_cds)]
      
      sub_cds_list[[i]] <- t_cds
      
      # let us downsample the data 
      if(use_pca) {
        data_ori <- cds@normalized_data_projection[pData(cds)$louvain_component == valid_groups[i], ]
      } else {
        data_ori <- t(cds@reducedDimS[, pData(cds)$louvain_component == valid_groups[i]])
      }
      if(ncol(t_cds) > landmark_num) {      
        centers <- data_ori[seq(1, nrow(data_ori), length.out=landmark_num), ]
        kmean_res <- kmeans(data_ori, landmark_num, centers=centers, iter.max = 100)
        landmark_id <- unique(findNearestVertex(t(kmean_res$centers), t(data_ori), process_targets_in_blocks=TRUE))
        
        data <- data_ori[landmark_id, ]
        
      } else {
        data <- data_ori
        landmark_id <- 1:ncol(t_cds)
      } 
      
      res <- project_to_representatives(data, 
                                        data_ori, 
                                        landmark_id, 
                                        cds@auxClusteringData$partitionCells,
                                        pd = pData(t_cds)[landmark_id, ], 
                                        method, 
                                        start.temp, 
                                        k, 
                                        do_partition = do_partition, 
                                        max_components = max_components, 
                                        verbose,
                                        ...) 
      
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
      row.names(coord) <- unlist(lapply(sub_cds_list, colnames))
    } else {
      coord <- igraph::merge_coords(sub_g_list, sub_coord_mat_list)
      row.names(coord) <- unlist(lapply(sub_cds_list, colnames))
    }
    g <- disjoint_union(sub_g_list)
    # browser()
    # row.names(coord) <- V(g)$name
  }
  
  cds@reducedDimS <- t(coord[colnames(cds), ])
  cds@reducedDimK <- t(coord[colnames(cds), ])
  
  minSpanningTree(cds) <- g
  
  return(cds)
}

#' Apply embedding smoothing techniques, inlcuding PSL, SSE, force-directed layout, etc. on the kNN graph built from all cells or the downsampled representive 
#' cells (if number of samples is large than 2000) and then project other non-representative cells to the low dimensional space with five nearest (on original 
#' UMAP or PCA space) representive cells' coordinates.
#'
#' @param data (dowmsampled) data to perform Louvain clustering and kNN graph construction 
#' @param data_ori All data points without downsampling. This data space will be used to find the nearest five points for projecting non-landmark points  
#' @param landmark_id The index of cells   
#' @param louvain_res The result of louvain clustering from the original space 
#' @param pd the data frame of phenotype information 
#' @param method The force directed layout function to use. Only relevant to when drl, fr or kk are used as fdl methods 
#' @param start.temp argument passed into layout_with_fr function 
#' @param do_partition Whether or not to separate participation groups and apply FDL in each partition or directly select equal number of representatives in each louvain groups.     
#' @param max_components the dimensionality of the reduced space
#' @param verbose Wheter to print all running details 
#' @param ... additional arguments passed to functions (louvain_clustering) called by this function. 
project_to_representatives <- function(data, 
                                       data_ori, 
                                       landmark_id, 
                                       louvain_res, 
                                       pd, 
                                       method = c('PSL', 'drl', 'fr', 'kk', 'SSE'),  
                                       start.temp = NULL, 
                                       k = 20,
                                       do_partition = F, 
                                       max_components = 2, 
                                       verbose, 
                                       ...) {
  extra_arguments <- list(...)
  
  if(verbose) {
    message("Running louvain clustering algorithm ...")
  }

  if(do_partition == FALSE) {
    d <- max_components
  } else {
    message('if do_partition is TRUE, we can only set the dimensionality of reduced space to be 2 because the drl coordinates-merging algorithm also support for two dimensions')
    d <- 2
  }
  
  # build an asymmetric kNN graph -- replace with the louvain_clustering one 
  adj_mat <- build_asym_kNN_graph(data, k, ...)
  # this one gives an symmetric matrix
  # louvain_clustering_args <- c(list(data = data, pd = pd, k = k, verbose = verbose),
  #                              extra_arguments[names(extra_arguments) %in% 
  #                              c("weight", "louvain_iter", "resolution", "random_seed")])
  # data_louvain_res <- do.call(louvain_clustering, louvain_clustering_args)
  # 
  # adj_mat <- igraph::get.adjacency(data_louvain_res$g)
  # adj_mat@x[adj_mat@x > 0] <- 1 # Ensure a connectivity graph (only 1 represents connectivity)
  if(method %in% c('drl', 'fr', 'kk')) {  
    # layout kNN with force direct layout: https://github.com/TypeFox/R-Examples/blob/d0917dbaf698cb8bc0789db0c3ab07453016eab9/igraph/R/layout_drl.R
    sub_g <- igraph::graph_from_adjacency_matrix(adj_mat, mode = 'direct', weighted = T) # 
    if (method=="fr") coord <- igraph::layout_with_fr(sub_g, dim=d, coords=data[, 1:d], weights = NULL, start.temp=start.temp)
    if (method=="drl") coord <- igraph::layout_with_drl(sub_g, dim=d, weights = NULL, options=drl_defaults$final) # list(edge.cut=0, simmer.attraction=0)
    if (method=="kk") coord <- igraph::layout_with_kk(sub_g, dim=d, coords=data[, 1:d]) # weights = E(sub_g)$weight / max(E(sub_g)$weight)
    
    colnames(coord) <- paste0(method, 1:d)
    row.names(coord) <- row.names(data)
  } else if(method == 'SSE') {
    # d <- 2 # d can only be 2 for dla coordinates merging method 
    sse_args <- c(list(data=data, dist_mat = adj_mat, verbose = verbose, embeding_dim = d),
                  extra_arguments[names(extra_arguments) %in% c("method", "para.gamma", "knn", "C", "maxiter", "beta")])
    SSE_res <- do.call(SSE, sse_args)    

    colnames(SSE_res$Y) <- paste0('SSE_', 1:d)
    rownames(SSE_res$Y) <- rownames(data)
    dimnames(SSE_res$W) <- list(rownames(data), rownames(data))
    
    coord <- SSE_res$Y 
    sub_g <- igraph::graph.adjacency(SSE_res$W, mode = "undirected", weighted = TRUE)
  } else if(method == 'PSL') {
    # d <- 2 # d can only be 2 for dla coordinates merging method 
    PSL_args <- c(list(Y = data, d = d, K = k, sG = adj_mat), # as.matrix(adj_mat)), #, # add verbose ncol(data)
                  extra_arguments[names(extra_arguments) %in% c('C', 'param.gamma', 'maxIter')])
    
    PSL_res <- do.call(psl, PSL_args)
    
    colnames(PSL_res$Z) <- paste0('PSL_', 1:d)
    rownames(PSL_res$Z) <- rownames(data)
    dimnames(PSL_res$S) <- list(rownames(data), rownames(data))
    
    coord <- as.matrix(PSL_res$Z) 
    sub_g <- igraph::graph.adjacency(PSL_res$S, mode = "undirected", weighted = TRUE)
  } else {
    stop("Unknown method. Please ensure method to be any one from 'drl', 'fr', 'kk', 'SSE' or 'PSL'")
  }
    
  # project other cells to the current location 
  if(nrow(data_ori) > nrow(data)) {
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
        tmp <- sort(x)[6] # choose 5 nearest neighbors for projection 
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
    row.names(coord) <- row.names(data_ori)
  }
  
  return(list(sub_coord_mat = coord, sub_g = sub_g))
}


#' The following functioin is taken from vegan package for performing procrustes analysis 
#'
#' @description Function procrustes rotates a configuration to maximum similarity with another configuration. Function protest tests the non-randomness (significance) between two configurations.
#' 
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
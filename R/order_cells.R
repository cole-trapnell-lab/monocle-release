#' Scale pseudotime to be in the range from 0 to 100 (it works both for situations involving only one state and complex states)
#' @param cds the CellDataSet upon which to perform this operation
#' @param ordering_genes a vector of feature ids (from the CellDataSet's featureData) used for ordering cells
#' @return an updated CellDataSet object which an
#' @export
scale_pseudotime <- function(cds) {
    pd <- pData(cds)
    pd$Cell_name <- row.names(pd)
    range_df <- plyr::ddply(pd, .(State), function(x) {
        min_max <- range(x$Pseudotime)
        min_cell <- subset(x, Pseudotime %in% min_max[1])
        max_cell <- subset(x, Pseudotime %in% min_max[2])
        min_cell$fate_type <- 'Start'
        max_cell$fate_type <- 'End'
        rbind(min_cell, max_cell)
    }) #pseudotime range for each state
    
    #1. construct a tree of the selected cells
    adj_list <- data.frame(Source = subset(range_df, length(Parent) > 0 & !is.na(Parent))[, 'Parent'], Target = subset(range_df, length(Parent) > 0 & !is.na(Parent))[, 'Cell_name'])
    #convert to cell fate adjancency list:
    adj_list$Source <- pd[as.character(adj_list$Source), 'State']
    adj_list$Target <- pd[as.character(adj_list$Target), 'State']
    
    uniq_cell_list <- unique(c(as.character(adj_list$Source), as.character(adj_list$Target)))
    adj_mat <- matrix(rep(0, length(uniq_cell_list)^2), nrow = length(uniq_cell_list), ncol = length(uniq_cell_list), dimnames = list(uniq_cell_list, uniq_cell_list))
    adj_mat[as.matrix(adj_list)] <- 1
    net <- graph.adjacency(as.matrix(adj_mat), mode = 'directed', weighted = NULL, diag = FALSE)
    
    # plot(net, layout = layout.fruchterman.reingold,
    #        vertex.size = 25,
    #        vertex.color="red",
    #        vertex.frame.color= "white",
    #        vertex.label.color = "white",
    #        vertex.label.cex = .5,
    #        vertex.label.family = "sans",
    #        edge.width=2,
    #        edge.color="black")
    
    #dfs_net <- graph.dfs(net, root = 1, order.out = T) #DFS search for the cell fate tree
    #net_order_out <- as.vector(dfs_net$order.out)
    net_leaves <- which(degree(net, v = V(net), mode = "out")==0, useNames = T)
    
    pd$scale_pseudotime <- NA
    
    #2. scale the psedutime:
    for(i in net_leaves) {
        path_vertex <- as.vector(get.all.shortest.paths(net, from = 1, to = i, mode="out")$res[[1]])
        pd_subset <- subset(pd, State %in% path_vertex & is.na(scale_pseudotime))
        
        #scale the pseudotime between the parent cell of the first cell to the last cell on the remaing branch
        # print(path_vertex)
        min_cell_name <- row.names(subset(pd_subset, Pseudotime == min(Pseudotime)))
        print(min_cell_name)
        if(!is.na(pd[min_cell_name, 'Parent'])) {
            parent_min_cell <- as.character(pd[min_cell_name, 'Parent'])
            subset_min_pseudo <- pd[parent_min_cell, 'Pseudotime']
            scale_pseudotime_ini <- pd[parent_min_cell, 'scale_pseudotime']
            scaling_factor <- (100 - pd[parent_min_cell, 'scale_pseudotime']) / c(max(pd_subset$Pseudotime) - subset_min_pseudo)
        }
        else {
            subset_min_pseudo <- min(pd[, 'Pseudotime'])
            scale_pseudotime_ini <- 0
            scaling_factor <- 100 / c(max(pd_subset$Pseudotime) - min(pd_subset$Pseudotime))
        }
        
        scale_pseudotime <- (pd_subset$Pseudotime - subset_min_pseudo) * scaling_factor + scale_pseudotime_ini
        message(i, '\t', range(scale_pseudotime)[1],'\t', range(scale_pseudotime)[2])
        pd[row.names(pd_subset), 'scale_pseudotime'] <- scale_pseudotime
    }	
    
    pData(cds) <- pd
    
    return(cds)
}

# Methods for PQ-tree based ordering

get_next_node_id <- function()
{
  next_node <<- next_node + 1
  return (next_node) 
}

#' Recursively builds and returns a PQ tree for the MST
#' @import igraph
pq_helper<-function(mst, use_weights=TRUE, root_node=NULL)
{
  new_subtree <- graph.empty()
  
  root_node_id <- paste("Q_", get_next_node_id(), sep="")
  
  new_subtree <- new_subtree + vertex(root_node_id, type="Q", color="black")
  
  if (is.null(root_node) == FALSE){
    sp <- get.all.shortest.paths(mst, from=V(mst)[root_node])
    #print(sp)
    sp_lengths <- sapply(sp$res, length)
    target_node_idx <- which(sp_lengths == max(sp_lengths))[1]
    #print(unlist(sp$res[target_node_idx]))
    diam <- V(mst)[unlist(sp$res[target_node_idx])]
    #print(diam)
  }else{
    if (use_weights){
      diam <- V(mst)[get.diameter(mst)]
    }else{
      diam <- V(mst)[get.diameter(mst, weights=NA)]
    }
  }
  
  
  #print (diam)
  
  V(new_subtree)[root_node_id]$diam_path_len = length(diam)
  
  diam_decisiveness <- igraph::degree(mst, v=diam) > 2
  ind_nodes <- diam_decisiveness[diam_decisiveness == TRUE]
  
  first_diam_path_node_idx <- head(as.vector(diam), n=1)
  last_diam_path_node_idx <- tail(as.vector(diam), n=1)
  if (sum(ind_nodes) == 0 || 
        (igraph::degree(mst, first_diam_path_node_idx) == 1 && 
           igraph::degree(mst, last_diam_path_node_idx) == 1))
  {
    ind_backbone <- diam
  }
  else 
  {
    last_bb_point <- names(tail(ind_nodes, n=1))[[1]]
    first_bb_point <- names(head(ind_nodes, n=1))[[1]]  
    #diam_path_vertex_names <- as.vector()
    #print (last_bb_point)
    #print (first_bb_point)
    diam_path_names <- V(mst)[as.vector(diam)]$name
    last_bb_point_idx <- which(diam_path_names == last_bb_point)[1]
    first_bb_point_idx <- which(diam_path_names == first_bb_point)[1]
    ind_backbone_idxs <- as.vector(diam)[first_bb_point_idx:last_bb_point_idx]
    #print (ind_backbone_idxs)
    ind_backbone <- V(mst)[ind_backbone_idxs]
    
    #ind_backbone <- diam[first_bb_point:last_bb_point]
  }
  
  
  
  mst_no_backbone <- mst - ind_backbone
  #print (V(mst_no_backbone)$name)
  
  for (backbone_n in ind_backbone)
  {
    #print (n)
    #backbone_n <- ind_backbone[[i]]
    
    if (igraph::degree(mst, v=backbone_n) > 2)
    {
      new_p_id <- paste("P_", get_next_node_id(), sep="")
      #print(new_p_id)
      new_subtree <- new_subtree + vertex(new_p_id, type="P", color="grey")
      new_subtree <- new_subtree + vertex(V(mst)[backbone_n]$name, type="leaf", color="white")
      new_subtree <- new_subtree + edge(new_p_id, V(mst)[backbone_n]$name)
      new_subtree <- new_subtree + edge(root_node_id, new_p_id)
      
      nb <- graph.neighborhood(mst, 1, nodes=backbone_n)[[1]]
      
      #print (E(nb))
      #print (V(nb))
      for (n_i in V(nb))
      {
        n <- V(nb)[n_i]$name			
        if (n %in% V(mst_no_backbone)$name)
        {	
          #print (n)
          
          sc <- subcomponent(mst_no_backbone, n)
          
          sg <- induced.subgraph(mst_no_backbone, sc, impl="copy_and_delete")
          
          
          if (ecount(sg) > 0)
          {
            #print (E(sg))	
            sub_pq <- pq_helper(sg, use_weights)
            
            
            # Works, but slow:
            for (v in V(sub_pq$subtree))
            {
              new_subtree <- new_subtree + vertex(V(sub_pq$subtree)[v]$name, type=V(sub_pq$subtree)[v]$type, color=V(sub_pq$subtree)[v]$color, diam_path_len=V(sub_pq$subtree)[v]$diam_path_len)
            }
            
            edge_list <- get.edgelist(sub_pq$subtree)
            for (i in 1:nrow(edge_list))
            {
              new_subtree <- new_subtree + edge(V(sub_pq$subtree)[edge_list[i, 1]]$name, V(sub_pq$subtree)[edge_list[i, 2]]$name)
            }   					
            #plot (new_subtree)
            
            new_subtree <- new_subtree + edge(new_p_id, V(sub_pq$subtree)[sub_pq$root]$name)  
          }
          else
          {
            new_subtree <- new_subtree + vertex(n, type="leaf", color="white")
            new_subtree <- new_subtree + edge(new_p_id, n)
          }
        }
        
      }
      #print ("##########################")
    }
    else
    {
      new_subtree <- new_subtree + vertex(V(mst)[backbone_n]$name, type="leaf", color="white")
      new_subtree <- new_subtree + edge(root_node_id, V(mst)[backbone_n]$name)
    }
  }
  # else
  # {
  #     for (backbone_n in diam)
  #     {
  #           new_subtree <- new_subtree + vertex(backbone_n, type="leaf")
  #           new_subtree <- new_subtree + edge(root_node_id, backbone_n)
  #     }
  # }
  
  return (list(root=root_node_id, subtree=new_subtree))
}

make_canonical <-function(pq_tree)
{
  nei <- NULL
  
  canonical_pq <- pq_tree
  
  V(canonical_pq)[type == "P" & igraph::degree(canonical_pq, mode="out") == 2]$color="black"
  V(canonical_pq)[type == "P" &  igraph::degree(canonical_pq, mode="out")== 2]$type="Q"
  
  single_child_p <- V(canonical_pq)[type == "P" & igraph::degree(canonical_pq, mode="out")== 1]
  V(canonical_pq)[type == "P" & igraph::degree(canonical_pq, mode="out")==1]$color="blue"
  
  for (p_node in single_child_p)
  {
    child_of_p_node <- V(canonical_pq) [ .nei(p_node, mode="out") ]
    parent_of_p_node <- V(canonical_pq) [ .nei(p_node, mode="in") ]
    
    for (child_of_p in child_of_p_node)
    {
      canonical_pq[parent_of_p_node, child_of_p] <- TRUE
      # print (p_node)
      # print (child_of_p)
      # print (parent_of_p_node)
      # print ("*********")
    }
  }
  
  canonical_pq <- delete.vertices(canonical_pq, V(canonical_pq)[type == "P" & igraph::degree(canonical_pq, mode="out")==1])
  #print (V(canonical_pq)[type == "Q" & igraph::degree(canonical_pq, mode="in")==0])
  return (canonical_pq)
}

count_leaf_descendents <- function(pq_tree, curr_node, children_counts)
{
  nei <- NULL
  
  if (V(pq_tree)[curr_node]$type == "leaf")
  {
    #V(pq_tree)[curr_node]$num_desc = 0
    children_counts[curr_node] = 0
    return(children_counts)
  } else {
    children_count = 0
    for (child in V(pq_tree) [ .nei(curr_node, mode="out") ])
    {
      children_counts <- count_leaf_descendents(pq_tree, child, children_counts)
      if (V(pq_tree)[child]$type == "leaf")
      {
        children_count <- children_count + 1
      }
      else
      {
        children_count <- children_count + children_counts[child]
      }
    }
    #print (curr_node)
    children_counts[curr_node] = children_count
    return(children_counts)
  }
}

measure_diameter_path <- function(pq_tree, curr_node, path_lengths)
{
  nei <- NULL
  
  if (V(pq_tree)[curr_node]$type != "Q")
  {
    #V(pq_tree)[curr_node]$num_desc = 0
    path_lengths[curr_node] = 0
    return(path_lengths)
  } else {
    
    children_count = 0
    for (child in V(pq_tree) [ .nei(curr_node, mode="out") ])
    {
      children_counts <- count_leaf_descendents(pq_tree, child, children_counts)
      if (V(pq_tree)[child]$type == "leaf")
      {
        children_count <- children_count + 1
      }
      else
      {
        children_count <- children_count + children_counts[child]
      }
    }
    
    
    path_lengths[curr_node] = children_count
    return(children_counts)
  }
}

# Assign leaf nodes reachable in pq_tree from curr_node to assigned_state
assign_cell_lineage <- function(pq_tree, curr_node, assigned_state, node_states)
{
  nei <- NULL
  
  if (V(pq_tree)[curr_node]$type == "leaf")
  {
    #V(pq_tree)[curr_node]$num_desc = 0
    #print (curr_node)
    node_states[V(pq_tree)[curr_node]$name] = assigned_state
    return(node_states)
  } else {
    for (child in V(pq_tree) [ .nei(curr_node, mode="out") ])
    {
      node_states <- assign_cell_lineage(pq_tree, child, assigned_state, node_states)
    }
    return(node_states)
  }
}

reverse_ordering <- function(pseudo_time_ordering)
{
  pt <- pseudo_time_ordering$pseudo_time
  names(pt) <- pseudo_time_ordering$sample_name
  rev_pt <- -((pt - max(pt)))
  rev_df <- pseudo_time_ordering
  rev_df$pseudo_time <- rev_pt
  return(rev_df)
}


weight_of_ordering <- function(ordering, dist_matrix)
{
  time_delta <- c(0)
  curr_weight <- 0
  ep <- 0.01
  for (i in 2:length(ordering))
  {
    d <- dist_matrix[ordering[[i]], ordering[[i-1]]]
    curr_weight <- curr_weight + d + ep
    time_delta <- c(time_delta, curr_weight)
  }
  
  return(time_delta)
}


#' Sets the features (e.g. genes) to be used for ordering cells in pseudotime.
#' @param cds the CellDataSet upon which to perform this operation
#' @param ordering_genes a vector of feature ids (from the CellDataSet's featureData) used for ordering cells
#' @return an updated CellDataSet object
#' @export
setOrderingFilter <- function(cds, ordering_genes){
  fData(cds)$use_for_ordering <- row.names(fData(cds)) %in% ordering_genes
  cds
}

#' Run the fastICA algorithm on a numeric matrix.
#' @importFrom fastICA  ica.R.def ica.R.par
#' @importFrom irlba irlba
ica_helper <- function(X, n.comp, alg.typ = c("parallel", "deflation"), fun = c("logcosh", "exp"), alpha = 1, 
                       row.norm = TRUE, maxit = 200, tol = 1e-4, verbose = FALSE, w.init = NULL, use_irlba=TRUE){
  dd <- dim(X) 
  
  d <- dd[dd != 1L]
  if (length(d) != 2L) 
    stop("data must be matrix-conformal")
  X <- if (length(d) != length(dd)) 
    matrix(X, d[1L], d[2L])
  else as.matrix(X)
  if (alpha < 1 || alpha > 2) 
    stop("alpha must be in range [1,2]")
  alg.typ <- match.arg(alg.typ)
  fun <- match.arg(fun)
  n <- nrow(X)
  p <- ncol(X)
  if (n.comp > min(n, p)) {
    message("'n.comp' is too large: reset to ", min(n, p))
    n.comp <- min(n, p)
  }
  if (is.null(w.init)) 
    w.init <- matrix(rnorm(n.comp^2), n.comp, n.comp)
  else {
    if (!is.matrix(w.init) || length(w.init) != (n.comp^2)) 
      stop("w.init is not a matrix or is the wrong size")
  }
  
  if (verbose) 
    message("Centering")
  X <- scale(X, scale = FALSE)
  X <- if (row.norm) 
    t(scale(X, scale = row.norm))
  else t(X)
  if (verbose) 
    message("Whitening")
  V <- X %*% t(X)/n
  
  if (verbose) 
    message("Finding SVD")
  if (use_irlba)
  {
    s <- irlba::irlba(V, n.comp, n.comp)  
    svs <- s$d  
  }
  else
  {
    s <- La.svd(V)
    svs <- s$d  
  }
  
  D <- diag(c(1/sqrt(s$d)))
  K <- D %*% t(s$u)
  K <- matrix(K[1:n.comp, ], n.comp, p)
  X1 <- K %*% X
  
  if (verbose) 
    message("Running ICA")
  if (alg.typ == "deflation") {
    a <- fastICA::ica.R.def(X1, n.comp, tol = tol, fun = fun, 
                   alpha = alpha, maxit = maxit, verbose = verbose, 
                   w.init = w.init)
  }
  else if (alg.typ == "parallel") {
    a <- fastICA::ica.R.par(X1, n.comp, tol = tol, fun = fun, 
                   alpha = alpha, maxit = maxit, verbose = verbose, 
                   w.init = w.init)
  }
  w <- a %*% K
  S <- w %*% X
  A <- t(w) %*% solve(w %*% t(w))
  return(list(X = t(X), K = t(K), W = t(a), A = t(A), S = t(S), svs=svs))
}

#
extract_ddrtree_ordering <- function(cds, root_cell, verbose=T)
{
  nei <- NULL
  type <- NULL
  pseudo_time <- NULL
   
  dp_mst <- minSpanningTree(cds) 
  
  # terminal_cell_ids <- V(dp_mst)[which(degree(dp_mst, mode = 'total') == 1)]
  # if(verbose) {
  #   print('the cells on the end of MST: ')
  #   print((degree(dp_mst, mode = 'total') == 1)[terminal_cell_ids])
  # }
  
  # if(is.null(root_cell))
  #   root_cell = terminal_cell_ids[1]
  
  
  Pseudotime <- shortest.paths(dp_mst, v=root_cell, to=V(dp_mst))
 
  curr_state <- 1
  #' a function to assign pseudotime for the MST
  assign_cell_state_helper <- function(ordering_tree_res, curr_cell, visited_node = curr_cell) {
    nei <- NULL

    cell_tree <- ordering_tree_res$subtree
    V(cell_tree)[curr_cell]$cell_state = curr_state

    children <- V(cell_tree) [ .nei(curr_cell, mode="all") ]$name
    children <- setdiff(children, V(cell_tree)[visited_node]$name)

    ordering_tree_res$subtree <- cell_tree
#     message('curr_cell: ', curr_cell)
#     message('children: ', children)

    if (length(children) == 0){
      return (ordering_tree_res)
    }else if (length(children) == 1){
        #visited_node <- union(children, visited_node)
        V(cell_tree)[children]$parent = rep(V(cell_tree)[visited_node]$name, length(children))
        ordering_tree_res <- assign_cell_state_helper(ordering_tree_res, V(cell_tree)[children]$name, curr_cell)
    }else{
        for (child in children) {
            #visited_node <- union(child, visited_node)
            V(cell_tree)[children]$parent = rep(V(cell_tree)[visited_node]$name, length(children))
            curr_state <<- curr_state + 1
            ordering_tree_res <- assign_cell_state_helper(ordering_tree_res, V(cell_tree)[child]$name, curr_cell)
        }
    }
    return (ordering_tree_res)
  }

  
  res <- list(subtree = dp_mst, root = root_cell)
  V(res$subtree)$parent <- rep(NA, nrow(pData(cds)))
  res <- assign_cell_state_helper(res, res$root)
  states <- V(res$subtree)[colnames(cds)]$cell_state
  
  cell_names <-  colnames(Pseudotime)
  cell_states <- states
  cell_pseudotime <- Pseudotime
  
  #add parents
  cell_parents <- V(res$subtree)$parent

  # ordering_df <- data.frame(sample_name = cell_names,
  #                           cell_state = factor(cell_states),
  #                           pseudo_time = cell_pseudotime,
  #                           parent = cell_parents)
  
  ordering_df <- data.frame(sample_name = cell_names,
                            cell_state = factor(cell_states),
                            pseudo_time = as.vector(cell_pseudotime),
                            parent = cell_parents)

  # ordering_df <- plyr::arrange(ordering_df, pseudo_time)
  return(ordering_df)
}


#' Orders cells according to progress through a learned biological process.
#' @param cds the CellDataSet upon which to perform this operation
#' @param num_paths the number of end-point cell states to allow in the biological process.
#' @param reverse whether to reverse the beginning and end points of the learned biological process.
#' @param root_cell the name of a cell to use as the root of the ordering tree.
#' @param scale_pseudotime a logic flag to determine whether or not to scale the pseudotime. If this is set to be true, then the pData dataframe of the returned CDS included an additional column called scale_pseudotime which store the scaled pseudotime values.
#' @return an updated CellDataSet object, in which phenoData contains values for State and Pseudotime for each cell
#' @export
orderCells <- function(cds, num_paths=1, root_state=NULL, scale_pseudotime = F){

  if (is.null(root_state)){
    diameter <- get.diameter(minSpanningTree(cds))
    root_cell = diameter[1]
  }else{
    # FIXME: Need to gaurd against the case when the supplied root state isn't actually a terminal state in the tree.
    root_cell_candidates <- subset(pData(cds), State = root_state)
    root_cell <-row.names(root_cell_candidates)[which(root_cell_candidates$Pseudotime == max(root_cell_candidates$Pseudotime))]
  }
  
  cc_ordering <- extract_ddrtree_ordering(cds, root_cell)
  row.names(cc_ordering) <- cc_ordering$sample_name
  
  pData(cds)$Pseudotime <-  cc_ordering[row.names(pData(cds)),]$pseudo_time
  pData(cds)$State <-  cc_ordering[row.names(pData(cds)),]$cell_state
  pData(cds)$Parent <-  cc_ordering[row.names(pData(cds)),]$parent
  
  if(scale_pseudotime) {
      cds <- scale_pseudotime(cds)
  }
  
  cds
}


#' Computes a projection of a CellDataSet object into a lower dimensional space
#' @param cds the CellDataSet upon which to perform this operation
#' @param max_components the dimensionality of the reduced space
#' @param use_irlba Whether to use the IRLBA package for ICA reduction.
#' @param pseudo_expr amount to increase expression values before dimensionality reduction
#' @param batch a vector of labels specifying batch for each cell, the effects of which will be removed prior to dimensionality reduction.
#' @param covariates a numeric vector or matrix specifying continuous effects to be removed prior to dimensionality reduction
#' @param ... additional arguments to pass to the dimensionality reduction function
#' @return an updated CellDataSet object
#' @details Currently, Monocle supports dimensionality reduction with Independent Component Analysis (ICA).
#' @importFrom matrixStats rowSds
#' @importFrom limma removeBatchEffect
#' @export
reduceDimension <- function(cds, 
                            max_components=2, 
                            use_irlba=TRUE, 
                            pseudo_expr=1, 
                            batch=NULL, 
                            batch2=NULL, 
                            covariates=NULL, 
                            use_vst=FALSE,
                            verbose=FALSE,
                            ...){
  
  FM <- exprs(cds)
  
  if (is.null(use_vst) && cds@expressionFamily@vfamily == "negbinomial"){
    use_vst = TRUE
    pseudo_expr = 0
  }
  
  # If we aren't using VST, then normalize the expression values by size factor
  if (use_vst == FALSE && cds@expressionFamily@vfamily == "negbinomial")
  {
    checkSizeFactors(cds)
    size_factors <- sizeFactors(cds)
    #print (size_factors)
    FM <- t(t(FM) / size_factors)
    #FM <- log2(FM)
  }
  
  if (is.null(fData(cds)$use_for_ordering) == FALSE && nrow(subset(fData(cds), use_for_ordering == TRUE)) > 0)
    FM <- FM[fData(cds)$use_for_ordering,]
  
  if (cds@expressionFamily@vfamily == "binomialff"){
    ncounts <- FM
    ncounts[ncounts != 0] <- 1
    FM <- t(t(ncounts) * log(1 + ncol(ncounts) / rowSums(ncounts)))
  }
  
  if (cds@expressionFamily@vfamily != "binomialff"){
    FM <- FM + pseudo_expr
  }
  
  FM <- FM[matrixStats::rowSds(FM) > 0,]
  

  if (cds@expressionFamily@vfamily != "binomialff"){
    if (use_vst){
      VST_FM <- vstExprs(cds, round_vals=FALSE)
      if (is.null(VST_FM) == FALSE){
        FM <- VST_FM
        FM <- FM[fData(cds)$use_for_ordering,]
      }else{
        stop("Error: set the variance-stabilized value matrix with vstExprs(cds) <- computeVarianceStabilizedValues() before calling this function with use_vst=TRUE")
      }
      #
    }else{
      FM <- log2(FM)
    }
  }

  # TODO: get rid of this if possible.
  if (is.null(batch) == FALSE || is.null(batch2) == FALSE|| is.null(covariates) == FALSE)
  {
    if (verbose)
      message("Removing batch effects")
    #FM <- log2(FM)
    FM <- limma::removeBatchEffect(FM, batch=batch, batch2=batch2, covariates=covariates)
    if (cds@expressionFamily@vfamily != "binomialff"){
      if (use_vst == FALSE) {
        FM <- 2^FM
      }
    }
  }
  
  #FM <- log2(FM)
  if (verbose)
    message("Reducing to independent components")
  
  ddrtree_res <- DDRTree_cpp(FM, max_components, verbose=F)
  
#   x_pca <- t(t(FM) %*% init_ICA$K)
#   W <- t(init_ICA$W)
#   
#   weights <- W
# 
#   A <- t(solve(weights) %*% t(init_ICA$K))
#   
#   colnames(A) <- colnames(weights)
#   rownames(A) <- rownames(FM)
#   
#   S <- weights %*% x_pca
#   
   #rownames(ddrtree_res$Y) <- colnames(weights)
  colnames(ddrtree_res$Y) <- colnames(FM) 
  colnames(ddrtree_res$Z) <- colnames(FM) 
  
  #reducedDimA(cds) <- A
  reducedDimS(cds) <- ddrtree_res$Z
  reducedDimK(cds) <- ddrtree_res$Y
  #reducedDimK(cds) <- init_ICA$K
  
  # Xiaojie, why did you change this? We should preserve this way.
  #stree <- ddrtree_res$stree
  #stree[lower.tri(stree)] = Matrix::t(stree)[lower.tri(stree)]

  #tmp <- as.matrix(ddrtree_res$stree)
  #stree <- tmp +  t(tmp) 
  
  # stree <- as.matrix(stree)
  # stree[stree == T] <- 1
  # stree[stree == F] <- 0

  #stree[upper.tri(stree)] <- 0
  #dimnames(stree) <- list(colnames(FM), colnames(FM))
  #gp <- graph.adjacency(stree, mode = "undirected", diag = F, weighted = TRUE)

  adjusted_K <- t(reducedDimK(cds))
  dp <- as.matrix(dist(adjusted_K))
  cellPairwiseDistances(cds) <- as.matrix(dist(adjusted_K))
  gp <- graph.adjacency(dp, mode = "undirected", weighted = TRUE)
  dp_mst <- minimum.spanning.tree(gp)
  minSpanningTree(cds) <- dp_mst
  
  # gp <- graph.adjacency(ddrtree_res$stree, mode="undirected", weighted=TRUE)
  # minSpanningTree(cds) <- gp
  
  cds
}


#' a function to assign pseudotime and states based on the projected coordinates from the DDRTree algorithm, the stree matrix maybe used later
#' also: fix the bug when the scale_pseudotime can be used
#' @param cds a matrix with N x N dimension
#' @param root_cell a dataframe used to generate new data for interpolation of time points
#' @param scale_pseudotime a matrix with N x N dimension
#' @param verbose a matrix with N x N dimension
#' @return a cds object with the states and the branch assigned correctly
#' @export
#' 
assignPseudotimeBranchPT <- function(cds, root_cell = NULL, scale_pseudotime = F, verbose = F) {

}



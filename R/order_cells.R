run_pca <- function(data, ...) {
  res <- prcomp(data, center = F, scale = F)
  res$x
}
# function to run umap with cds as input 
umapM <- function(cds, n_neighbors = 10, n_component = 2, min_dist = 0.1, metric = "euclidean") {
  x <- t(cds@assayData$exprs)
  if(length(grep('Matrix', class(x))) == 0){
    x <- as(as.matrix(x), 'TsparseMatrix')
  } else {
    x <- as(x, 'TsparseMatrix')
  }
  
  i <- x@i
  j <- x@j
  val <- log(x@x + 1)
  dim <- x@Dim
  
  rPython::python.exec( c( "def umap(i, j, val, dim, n, n_c,mdist, metric):",
                           "\timport umap" ,
                           "\timport numpy",
                           "\tfrom scipy.sparse import csc_matrix",
                           "\tdata = csc_matrix((val, (i, j)), shape = dim)",
                           "\tembedding = umap.UMAP(n_neighbors = n, n_components = n_c, min_dist = mdist, metric = metric).fit_transform(data)",
                           "\tres = embedding.tolist()",
                           "\treturn res"))
  
  res <- rPython::python.call( "umap", i, j, val, dim, n_neighbors, n_component, min_dist, metric, simplify = F)
  do.call("rbind",res)
}

umap <- function(x,n_neighbors=10,n_component=2,min_dist=0.1,metric="euclidean"){
  x <- as.matrix(x)
  colnames(x) <- NULL
  rPython::python.exec( c( "def umap(data,n,n_c,mdist,metric):",
                  "\timport umap" ,
                  "\timport numpy",
                  "\tembedding = umap.UMAP(n_neighbors=n,n_components=n_c,min_dist=mdist,metric=metric).fit_transform(data)",
                  "\tres = embedding.tolist()",
                  "\treturn res"))
  
  res <- rPython::python.call( "umap", x,n_neighbors,n_component,min_dist,metric)
  do.call("rbind",res)
}

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

# #' Scale pseudotime to be in the range from 0 to 100
# #'
# #' This function transforms the pseudotime scale so that it ranges from 0 to 100. If there are multiple branches, each leaf is set to be 100, with branches stretched accordingly.
# #'
# #' @param cds the CellDataSet upon which to perform this operation
# #' @param verbose Whether to emit verbose output
# #' @return an updated CellDataSet object which an
#' @importFrom igraph graph.adjacency V graph.dfs get.all.shortest.paths
scale_pseudotime <- function(cds, verbose = F) {
  Parent <- NA
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
    # print(min_cell_name)

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

    pseudotime_scaled <- (pd_subset$Pseudotime - subset_min_pseudo) * scaling_factor + scale_pseudotime_ini

    if(verbose)
      message(i, '\t', range(pseudotime_scaled)[1],'\t', range(pseudotime_scaled)[2])

    pd[row.names(pd_subset), 'ori_pseudotime'] <- pd[row.names(pd_subset), 'Pseudotime']
    pd[row.names(pd_subset), 'Pseudotime'] <- pseudotime_scaled
  }
  scale_pseudotime <- (pd_subset$Pseudotime - subset_min_pseudo) * scaling_factor + scale_pseudotime_ini
  message(i, '\t', range(scale_pseudotime)[1],'\t', range(scale_pseudotime)[2])
  pd[row.names(pd_subset), 'scale_pseudotime'] <- scale_pseudotime

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
#' @param mst The minimum spanning tree, as an igraph object.
#' @param use_weights Whether to use edge weights when finding the diameter path of the tree.
#' @param root_node The name of the root node to use for starting the path finding.
#' @importFrom igraph V vertex degree get.diameter edge graph.neighborhood 
#' @importFrom igraph graph.empty get.edgelist get.all.shortest.paths
#' @importFrom igraph subcomponent induced.subgraph ecount
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

#' @importFrom igraph V degree V<- delete.vertices
make_canonical <-function(pq_tree)
{
  type <- NA
  canonical_pq <- pq_tree

  V(canonical_pq)[type == "P" & igraph::degree(canonical_pq, mode="out") == 2]$color="black"
  V(canonical_pq)[type == "P" &  igraph::degree(canonical_pq, mode="out")== 2]$type="Q"

  single_child_p <- V(canonical_pq)[type == "P" & igraph::degree(canonical_pq, mode="out")== 1]
  V(canonical_pq)[type == "P" & igraph::degree(canonical_pq, mode="out")==1]$color="blue"

  for (p_node in single_child_p)
  {
    child_of_p_node <- V(canonical_pq) [ suppressWarnings(nei(p_node, mode="out")) ]
    parent_of_p_node <- V(canonical_pq) [ suppressWarnings(nei(p_node, mode="in")) ]

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

#' @importFrom igraph V
count_leaf_descendents <- function(pq_tree, curr_node, children_counts)
{


  if (V(pq_tree)[curr_node]$type == "leaf")
  {
    #V(pq_tree)[curr_node]$num_desc = 0
    children_counts[curr_node] = 0
    return(children_counts)
  } else {
    children_count = 0
    for (child in V(pq_tree) [ suppressWarnings(nei(curr_node, mode="out")) ])
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


#' Return an ordering for a P node in the PQ tree
#' @param q_level_list A list of Q nodes in the PQ tree
#' @param dist_matrix A symmetric matrix of pairwise distances between cells
#' @importFrom combinat permn
order_p_node <- function(q_level_list, dist_matrix)
{
  q_order_res <- combinat::permn(q_level_list, fun=order_q_node, dist_matrix)
  #print (q_order_res)
  all_perms <- lapply(q_order_res, function(x) { x$ql } )
  #print ("perm ql:")
  #print(all_perms)
  all_perms_weights <- unlist(lapply(q_order_res, function(x) { x$wt }))
  #print ("perm weights:")
  #print (all_perms_weights)

  opt_perm_idx <- head((which(all_perms_weights == min(all_perms_weights))), 1)
  opt_perm <- all_perms[[opt_perm_idx]]

  #print ("opt_path:")
  #print (opt_perm)
  #print ("opt_all_weight:")
  #print (min(all_perms_weights))
  #print ("weights:")
  #print (all_perms_weights)
  # print ("q_level_list:")
  # print (q_level_list)
  stopifnot (length(opt_perm) == length(q_level_list))

  return(opt_perm)
}

#' @importFrom igraph V vertex edge graph.empty get.shortest.paths E
order_q_node <- function(q_level_list, dist_matrix)
{
  new_subtree <- graph.empty()

  if (length(q_level_list) == 1)
  {
    return (list(ql=q_level_list, wt=0))
  }
  for (i in 1:length(q_level_list))
  {
    new_subtree <- new_subtree + vertex(paste(i,"F"), type="forward")
    new_subtree <- new_subtree + vertex(paste(i,"R"), type="reverse")
  }

  for (i in (1:(length(q_level_list)-1)))
  {
    cost <- dist_matrix[q_level_list[[i]][length(q_level_list[[i]])], q_level_list[[i+1]][1]]
    new_subtree <- new_subtree + edge(paste(i,"F"), paste(i+1,"F"), weight=cost)

    cost <- dist_matrix[q_level_list[[i]][length(q_level_list[[i]])], q_level_list[[i+1]][length(q_level_list[[i+1]])]]
    new_subtree <- new_subtree + edge(paste(i,"F"), paste(i+1,"R"), weight=cost)

    cost <- dist_matrix[q_level_list[[i]][1], q_level_list[[i+1]][1]]
    new_subtree <- new_subtree + edge(paste(i,"R"), paste(i+1,"F"), weight=cost)

    cost <- dist_matrix[q_level_list[[i]][1], q_level_list[[i+1]][length(q_level_list[[i+1]])]]
    new_subtree <- new_subtree + edge(paste(i,"R"), paste(i+1,"R"), weight=cost)
  }

  first_fwd = V(new_subtree)[paste(1,"F")]
  first_rev = V(new_subtree)[paste(1,"R")]
  last_fwd = V(new_subtree)[paste(length(q_level_list),"F")]
  last_rev = V(new_subtree)[paste(length(q_level_list),"R")]

  FF_path <- unlist(get.shortest.paths(new_subtree, from=as.vector(first_fwd), to=as.vector(last_fwd), mode="out", output="vpath")$vpath)
  FR_path <- unlist(get.shortest.paths(new_subtree, from=as.vector(first_fwd), to=as.vector(last_rev), mode="out", output="vpath")$vpath)
  RF_path <- unlist(get.shortest.paths(new_subtree, from=as.vector(first_rev), to=as.vector(last_fwd), mode="out", output="vpath")$vpath)
  RR_path <- unlist(get.shortest.paths(new_subtree, from=as.vector(first_rev), to=as.vector(last_rev), mode="out", output="vpath")$vpath)

  # print (FF_path)
  # print (FR_path)
  # print (RF_path)
  # print (RR_path)

  FF_weight <- sum(E(new_subtree, path=FF_path)$weight)
  FR_weight <- sum(E(new_subtree, path=FR_path)$weight)
  RF_weight <- sum(E(new_subtree, path=RF_path)$weight)
  RR_weight <- sum(E(new_subtree, path=RR_path)$weight)

  # print (FF_weight)
  # print (FR_weight)
  # print (RF_weight)
  # print (RR_weight)

  paths <- list(FF_path, FR_path, RF_path, RR_path)
  path_weights <- c(FF_weight, FR_weight, RF_weight, RR_weight)
  opt_path_idx <- head((which(path_weights == min(path_weights))), 1)
  opt_path <- paths[[opt_path_idx]]

  # print ("opt_path:")
  # print (opt_path)
  # print ("q_level_list:")
  # print (q_level_list)
  stopifnot (length(opt_path) == length(q_level_list))

  directions <- V(new_subtree)[opt_path]$type
  #print (directions)
  q_levels <- list()
  for (i in 1:length(directions))
  {
    if (directions[[i]] == "forward"){
      q_levels[[length(q_levels)+1]] <- q_level_list[[i]]
    }else{
      q_levels[[length(q_levels)+1]] <- rev(q_level_list[[i]])
    }
  }

  return(list(ql=q_levels, wt=min(path_weights)))
}

#' @importFrom igraph V
measure_diameter_path <- function(pq_tree, curr_node, path_lengths)
{

  if (V(pq_tree)[curr_node]$type != "Q")
  {
    #V(pq_tree)[curr_node]$num_desc = 0
    path_lengths[curr_node] = 0
    return(path_lengths)
  } else {

    children_count = 0
    for (child in V(pq_tree) [ suppressWarnings(nei(curr_node, mode="out")) ])
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
#' @importFrom igraph V
assign_cell_lineage <- function(pq_tree, curr_node, assigned_state, node_states)
{


  if (V(pq_tree)[curr_node]$type == "leaf")
  {
    #V(pq_tree)[curr_node]$num_desc = 0
    #print (curr_node)
    node_states[V(pq_tree)[curr_node]$name] = assigned_state
    return(node_states)
  } else {
    for (child in V(pq_tree) [ suppressWarnings(nei(curr_node, mode="out")) ])
    {
      node_states <- assign_cell_lineage(pq_tree, child, assigned_state, node_states)
    }
    return(node_states)
  }
}

#' @importFrom igraph V
extract_good_ordering <- function(pq_tree, curr_node, dist_matrix)
{


  if (V(pq_tree)[curr_node]$type == "leaf")
  {
    #print ("ordering leaf node")
    return (V(pq_tree)[curr_node]$name)
  }else if (V(pq_tree)[curr_node]$type == "P"){
    #print ("ordering P node")
    p_level <- list()
    for (child in V(pq_tree) [ suppressWarnings(nei(curr_node, mode="out")) ])
    {
      p_level[[length(p_level)+1]] <- extract_good_ordering(pq_tree, child, dist_matrix)
    }
    p_level <- order_p_node(p_level, dist_matrix)
    p_level <- unlist(p_level)
    #print (p_level)
    return (p_level)
  }else if(V(pq_tree)[curr_node]$type == "Q"){
    #print ("ordering Q node")
    q_level <- list()
    for (child in V(pq_tree) [ suppressWarnings(nei(curr_node, mode="out")) ])
    {
      q_level[[length(q_level)+1]] <- extract_good_ordering(pq_tree, child, dist_matrix)
    }
    q_level <- order_q_node(q_level, dist_matrix)
    q_level <- q_level$ql
    q_level <- unlist(q_level)
    #print (q_level)
    return (q_level)
  }
}



#' Extract a linear ordering of cells from a PQ tree
#'
#' @param orig_pq_tree The PQ object to use for ordering
#' @param curr_node The node in the PQ tree to use as the start of ordering
#' @param dist_matrix A symmetric matrix containing pairwise distances between cells
#' @param num_branches The number of outcomes allowed in the trajectory.
#' @param reverse_main_path Whether to reverse the direction of the trajectory
#' 
#' @importFrom igraph V vertex edge graph.empty get.edgelist
extract_good_branched_ordering <- function(orig_pq_tree, curr_node, dist_matrix, num_branches, reverse_main_path=FALSE)
{
  requireNamespace("plyr")
  nei <- NULL
  type <- NA
  pseudo_time <- NA

  pq_tree <- orig_pq_tree

  # children_counts <- rep(0, length(as.vector(V(pq_tree))))
  #     names(children_counts) <- V(pq_tree)$name
  # children_counts <- count_leaf_descendents(pq_tree, curr_node, children_counts)
  #
  # branch_node_counts <- children_counts[V(res$subtree)[type == "P"]]
  # branch_node_counts <- sort(branch_node_counts, decreasing=TRUE)
  # print (branch_node_counts)


  branch_node_counts <- V(pq_tree)[type == "Q"]$diam_path_len
  names(branch_node_counts) <- V(pq_tree)[type == "Q"]$name
  if(length(names(branch_node_counts)) < num_branches)
    stop('Number of branches attempted is larger than the branches constructed from pq_tree algorithm')

  branch_node_counts <- sort(branch_node_counts, decreasing=TRUE)
  #print (branch_node_counts)


  cell_states <- rep(NA, length(as.vector(V(pq_tree)[type=="leaf"])))
  names(cell_states) <- V(pq_tree)[type=="leaf"]$name

  cell_states <- assign_cell_lineage(pq_tree, curr_node, 1, cell_states)

  branch_point_roots <- list()

  # Start building the ordering tree. Each pseudo-time segment will be a node.
  branch_tree <- graph.empty()
  #root_branch_id <- "Q_1"
  #branch_tree <- branch_tree + vertex(root_branch_id)

  for (i in 1:num_branches)
  {
    #cell_states <- assign_cell_lineage(pq_tree, names(branch_node_counts)[i], i+1, cell_states)
    #print (head(cell_states))
    #print(names(branch_node_counts)[i])

    branch_point_roots[[length(branch_point_roots) + 1]] <- names(branch_node_counts)[i]
    branch_id <- names(branch_node_counts)[i]
    #print (branch_id)
    branch_tree <- branch_tree + vertex(branch_id)
    parents <- V(pq_tree)[suppressWarnings(nei(names(branch_node_counts)[i], mode="in"))]
    if (length(parents) > 0 && parents$type == "P")
    {
      p_node_parent <- V(pq_tree)[suppressWarnings(nei(names(branch_node_counts)[i], mode="in"))]
      parent_branch_id <- V(pq_tree)[suppressWarnings(nei(p_node_parent, mode="in"))]$name
      #print (parent_branch_id)
      #print (branch_id)
      branch_tree <- branch_tree + edge(parent_branch_id, branch_id)
    }
    pq_tree[V(pq_tree) [ suppressWarnings(nei(names(branch_node_counts)[i], mode="in")) ], names(branch_node_counts)[i] ] <- FALSE
  }

  #branch_point_roots[[length(branch_point_roots) + 1]] <- curr_node
  #branch_point_roots <- rev(branch_point_roots)
  branch_pseudotimes <- list()

  for (i in 1:length(branch_point_roots))
  {
    branch_ordering <- extract_good_ordering(pq_tree, branch_point_roots[[i]], dist_matrix)
    branch_ordering_time <- weight_of_ordering(branch_ordering, dist_matrix)
    names(branch_ordering_time) <- branch_ordering
    branch_pseudotimes[[length(branch_pseudotimes) + 1]] = branch_ordering_time
    names(branch_pseudotimes)[length(branch_pseudotimes)] = branch_point_roots[[i]]
  }

  cell_ordering_tree <- graph.empty()
  curr_branch <- "Q_1"

  extract_branched_ordering_helper <- function(branch_tree, curr_branch, cell_ordering_tree, branch_pseudotimes, dist_matrix, reverse_ordering=FALSE)
  {
    nei <- NULL

    curr_branch_pseudotimes <- branch_pseudotimes[[curr_branch]]
    #print (curr_branch_pseudotimes)
    curr_branch_root_cell <- NA
    for (i in 1:length(curr_branch_pseudotimes))
    {
      cell_ordering_tree <- cell_ordering_tree + vertex(names(curr_branch_pseudotimes)[i])
      if (i > 1)
      {
        if (reverse_ordering == FALSE){
          cell_ordering_tree <- cell_ordering_tree + edge(names(curr_branch_pseudotimes)[i-1], names(curr_branch_pseudotimes)[i])
        }else{
          cell_ordering_tree <- cell_ordering_tree + edge(names(curr_branch_pseudotimes)[i], names(curr_branch_pseudotimes)[i-1])
        }
      }
    }

    if (reverse_ordering == FALSE)
    {
      curr_branch_root_cell <- names(curr_branch_pseudotimes)[1]
    }else{
      curr_branch_root_cell <- names(curr_branch_pseudotimes)[length(curr_branch_pseudotimes)]
    }

    for (child in V(branch_tree) [ suppressWarnings(nei(curr_branch, mode="out")) ])
    {
      child_cell_ordering_subtree <- graph.empty()

      child_head <- names(branch_pseudotimes[[child]])[1]
      child_tail <- names(branch_pseudotimes[[child]])[length(branch_pseudotimes[[child]])]

      # find the closest cell in the parent branch for each of the head and the tail

      curr_branch_cell_names <- names(branch_pseudotimes[[curr_branch]])
      head_dist_to_curr <- dist_matrix[child_head, curr_branch_cell_names]
      closest_to_head <- names(head_dist_to_curr)[which(head_dist_to_curr == min(head_dist_to_curr))]

      head_dist_to_anchored_branch = NA
      branch_index_for_head <- NA

      head_dist_to_anchored_branch <- dist_matrix[closest_to_head, child_head]

      tail_dist_to_curr <- dist_matrix[child_tail, curr_branch_cell_names]
      closest_to_tail <- names(tail_dist_to_curr)[which(tail_dist_to_curr == min(tail_dist_to_curr))]

      tail_dist_to_anchored_branch = NA
      branch_index_for_tail <- NA

      tail_dist_to_anchored_branch <- dist_matrix[closest_to_tail, child_tail]

      if (tail_dist_to_anchored_branch < head_dist_to_anchored_branch)
      {
        reverse_child <- TRUE
      }else{
        reverse_child <- FALSE
      }

      res <- extract_branched_ordering_helper(branch_tree, child, child_cell_ordering_subtree, branch_pseudotimes, dist_matrix, reverse_child)
      child_cell_ordering_subtree <- res$subtree
      child_subtree_root <- res$root

      # Works, but slow:
      for (v in V(child_cell_ordering_subtree))
      {
        cell_ordering_tree <- cell_ordering_tree + vertex(V(child_cell_ordering_subtree)[v]$name)
      }

      edge_list <- get.edgelist(child_cell_ordering_subtree)
      for (i in 1:nrow(edge_list))
      {
        cell_ordering_tree <- cell_ordering_tree + edge(V(cell_ordering_tree)[edge_list[i, 1]]$name, V(cell_ordering_tree)[edge_list[i, 2]]$name)
      }

      if (tail_dist_to_anchored_branch < head_dist_to_anchored_branch)
      {
        cell_ordering_tree <- cell_ordering_tree + edge(closest_to_tail, child_subtree_root)
      }else{
        cell_ordering_tree <- cell_ordering_tree + edge(closest_to_head, child_subtree_root)
      }

    }

    return (list(subtree=cell_ordering_tree, root=curr_branch_root_cell, last_cell_state=1, last_cell_pseudotime=0.0))
  }

  res <- extract_branched_ordering_helper(branch_tree, curr_branch, cell_ordering_tree, branch_pseudotimes, dist_matrix, reverse_main_path)
  cell_ordering_tree <- res$subtree

  curr_state <- 1

  assign_cell_state_helper <- function(ordering_tree_res, curr_cell)
  {
    nei <- NULL

    cell_tree <- ordering_tree_res$subtree
    V(cell_tree)[curr_cell]$cell_state = curr_state

    children <- V(cell_tree) [ suppressWarnings(nei(curr_cell, mode="out")) ]
    ordering_tree_res$subtree <- cell_tree

    if (length(children) == 1){
      ordering_tree_res <- assign_cell_state_helper(ordering_tree_res, V(cell_tree)[children]$name)
    }else{
      for (child in children) {
        curr_state <<- curr_state + 1
        ordering_tree_res <- assign_cell_state_helper(ordering_tree_res, V(cell_tree)[child]$name)
      }
    }
    return (ordering_tree_res)
  }

  res <- assign_cell_state_helper(res, res$root)

  assign_pseudotime_helper <- function(ordering_tree_res, dist_matrix, last_pseudotime, curr_cell)
  {
    nei <- NULL

    cell_tree <- ordering_tree_res$subtree
    curr_cell_pseudotime <- last_pseudotime
    V(cell_tree)[curr_cell]$pseudotime = curr_cell_pseudotime
    V(cell_tree)[curr_cell]$parent =  V(cell_tree)[ suppressWarnings(nei(curr_cell, mode="in")) ]$name
    #print (curr_cell_pseudotime)

    ordering_tree_res$subtree <- cell_tree
    children <- V(cell_tree) [ suppressWarnings(nei(curr_cell, mode="out")) ]

    for (child in children) {
      next_node <- V(cell_tree)[child]$name
      delta_pseudotime <- dist_matrix[curr_cell, next_node]
      ordering_tree_res <- assign_pseudotime_helper(ordering_tree_res, dist_matrix, last_pseudotime + delta_pseudotime, next_node)
    }

    return (ordering_tree_res)
  }

  res <- assign_pseudotime_helper(res, dist_matrix, 0.0, res$root)

  cell_names <- V(res$subtree)$name
  cell_states <- V(res$subtree)$cell_state
  cell_pseudotime <- V(res$subtree)$pseudotime
  cell_parents <- V(res$subtree)$parent
  # print (cell_names)
  # print (cell_states)
  # print (cell_pseudotime)
  ordering_df <- data.frame(sample_name = cell_names,
                            cell_state = factor(cell_states),
                            pseudo_time = cell_pseudotime,
                            parent = cell_parents)

  ordering_df <- plyr::arrange(ordering_df, pseudo_time)
  return(list("ordering_df"=ordering_df, "cell_ordering_tree"=cell_ordering_tree))
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

# Run the fastICA algorithm on a numeric matrix.
#' @importFrom stats rnorm qnorm
ica_helper <- function(X, n.comp, alg.typ = c("parallel", "deflation"), fun = c("logcosh", "exp"), alpha = 1,
                       row.norm = TRUE, maxit = 200, tol = 1e-4, verbose = FALSE, w.init = NULL, use_irlba=TRUE){
  dd <- dim(X)
  #FIXME: This will internally convert to a dense matrix
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

  initial_v <- as.matrix(qnorm(1:(ncol(V) + 1)/(ncol(V) + 1))[1:ncol(V)])
  s <- irlba::irlba(V, n.comp, n.comp, v = initial_v)
  svs <- s$d

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

#' @importFrom igraph V minimum.spanning.tree graph.adjacency degree get.diameter graph.dfs
extract_ddrtree_ordering <- function(cds, root_cell, verbose=T)
{

  dp <- cellPairwiseDistances(cds)
  dp_mst <- minSpanningTree(cds)

  curr_state <- 1

  res <- list(subtree = dp_mst, root = root_cell)

  states = rep(1, ncol(dp))
  names(states) <- V(dp_mst)$name

  pseudotimes = rep(0, ncol(dp))
  names(pseudotimes) <- V(dp_mst)$name

  parents = rep(NA, ncol(dp))
  names(parents) <- V(dp_mst)$name

  mst_traversal <- graph.dfs(dp_mst,
                             root=root_cell,
                             neimode = "all",
                             unreachable=FALSE,
                             father=TRUE)
  mst_traversal$father <- as.numeric(mst_traversal$father)
  curr_state <- 1

  for (i in 1:length(mst_traversal$order)){
    curr_node = mst_traversal$order[i]
    curr_node_name = V(dp_mst)[curr_node]$name

    if (is.na(mst_traversal$father[curr_node]) == FALSE){
      parent_node = mst_traversal$father[curr_node]
      parent_node_name = V(dp_mst)[parent_node]$name
      parent_node_pseudotime = pseudotimes[parent_node_name]
      parent_node_state = states[parent_node_name]
      curr_node_pseudotime = parent_node_pseudotime + dp[curr_node_name, parent_node_name]
      if (degree(dp_mst, v=parent_node_name) > 2){
        curr_state <- curr_state + 1
      }
    }else{
      parent_node = NA
      parent_node_name = NA
      curr_node_pseudotime = 0
    }

    curr_node_state = curr_state
    pseudotimes[curr_node_name] <- curr_node_pseudotime
    states[curr_node_name] <- curr_node_state
    parents[curr_node_name] <- parent_node_name
  }

  ordering_df <- data.frame(sample_name = names(states),
                            cell_state = factor(states),
                            pseudo_time = as.vector(pseudotimes),
                            parent = parents)
  row.names(ordering_df) <- ordering_df$sample_name
  # ordering_df <- plyr::arrange(ordering_df, pseudo_time)
  return(ordering_df)
}

#' @importFrom stats dist
#' @importFrom igraph graph.adjacency minimum.spanning.tree V
select_root_cell <- function(cds, root_state=NULL, reverse=FALSE){
  if (is.null(root_state) == FALSE & vcount(minSpanningTree(cds)) == ncol(cds)) {
    if (is.null(pData(cds)$State)){
      stop("Error: State has not yet been set. Please call orderCells() without specifying root_state, then try this call again.")
    }
    # FIXME: Need to gaurd against the case when the supplied root state isn't actually a terminal state in the tree.
    root_cell_candidates <- subset(pData(cds), State == root_state)
    if (nrow(root_cell_candidates) == 0){
      stop(paste("Error: no cells for State =", root_state))
    }

    # build a local MST to find a good root cell for this state
    dp <- as.matrix(dist(t(reducedDimS(cds)[,row.names(root_cell_candidates)])))
    gp <- graph.adjacency(dp, mode = "undirected", weighted = TRUE)
    dp_mst <- minimum.spanning.tree(gp)

    # Make sure to use the real MST here
    tip_leaves <- names(which(degree(minSpanningTree(cds)) == 1))
    #root_cell_candidates <- root_cell_candidates[row.names(root_cell_candidates) %in% tip_leaves,]
    #sg <- make_ego_graph(dp_mst, nodes=row.names(root_cell_candidates))[[1]]

    # if(length(intersect(tip_leaves, row.names(root_cell_candidates))) == 0){
    #   stop(paste("Error: no valid root cells for State =", root_state))
    # }
      
    diameter <- get.diameter(dp_mst)

    if (length(diameter) == 0){
      stop(paste("Error: no valid root cells for State =", root_state))
    }

    #root_cell = names(diameter)[tip_leaves %in% names(diameter)]
    root_cell_candidates <- root_cell_candidates[names(diameter),]
    if (is.null(cds@auxOrderingData[[cds@dim_reduce_type]]$root_cell) == FALSE &&
        pData(cds)[cds@auxOrderingData[[cds@dim_reduce_type]]$root_cell,]$State == root_state){
        root_cell <- row.names(root_cell_candidates)[which(root_cell_candidates$Pseudotime == min(root_cell_candidates$Pseudotime))]
    }else{
      root_cell <- row.names(root_cell_candidates)[which(root_cell_candidates$Pseudotime == max(root_cell_candidates$Pseudotime))]
    }
    if (length(root_cell) > 1)
      root_cell <- root_cell[1]

    # If we used DDRTree, we need to go from this actual cell to the nearst
    # point on the principal graph
    if (cds@dim_reduce_type == "DDRTree"){
      #root_cell_idx <- which(V(minSpanningTree(cds))$name == root_cell, arr.ind=T)
      graph_point_for_root_cell <- cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex[root_cell,]
      root_cell = V(minSpanningTree(cds))[graph_point_for_root_cell]$name
    }

  }else{ 
    if(is.null(root_state)) {
      if (is.null(minSpanningTree(cds))){
        stop("Error: no spanning tree found for CellDataSet object. Please call reduceDimension before calling orderCells()")
      }
      diameter <- get.diameter(minSpanningTree(cds))
      if (is.null(reverse) == FALSE && reverse == TRUE){
        root_cell = names(diameter[length(diameter)])
      } else {
        root_cell = names(diameter[1])
      }
    } else {
      root_cell_candidates <- which(pData(cds)$State == root_state)
      leaf_cells <- which(degree(minSpanningTree(cds)) == 1)
      R <- cds@auxOrderingData$DDRTree$R
      edge <- data.frame(start = 1:nrow(R), end = apply(R, 1, which.max), weight = R[cbind(1:nrow(R), apply(R, 1, which.max))])
      
      root_cell <- V(minSpanningTree(cds))$name[intersect(edge[root_cell_candidates, 'end'], leaf_cells)]

    }
  }
  return(root_cell)
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
                       root_state=NULL,
                       num_paths = NULL,
                       reverse=NULL){
  
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

  root_cell <- select_root_cell(cds, root_state, reverse)

  if (cds@dim_reduce_type == "ICA"){
    cds@auxOrderingData <- new.env( hash=TRUE )
    if (is.null(num_paths)){
      num_paths = 1
    }
    adjusted_S <- t(cds@reducedDimS)

    dp <- as.matrix(dist(adjusted_S))

    cellPairwiseDistances(cds) <- as.matrix(dist(adjusted_S))
    # Build an MST of the cells in ICA space.
    gp <- graph.adjacency(dp, mode="undirected", weighted=TRUE)
    dp_mst <- minimum.spanning.tree(gp)
    minSpanningTree(cds) <- dp_mst
    # Build the PQ tree
    next_node <<- 0
    res <- pq_helper(dp_mst, use_weights=FALSE, root_node=root_cell)

    cds@auxOrderingData[[cds@dim_reduce_type]]$root_cell <- root_cell

    order_list <- extract_good_branched_ordering(res$subtree, res$root, cellPairwiseDistances(cds), num_paths, FALSE)
    cc_ordering <- order_list$ordering_df
    row.names(cc_ordering) <- cc_ordering$sample_name

    minSpanningTree(cds)  <- as.undirected(order_list$cell_ordering_tree)

    pData(cds)$Pseudotime <-  cc_ordering[row.names(pData(cds)),]$pseudo_time
    pData(cds)$State <-  cc_ordering[row.names(pData(cds)),]$cell_state
    #pData(cds)$Parent <-  cc_ordering[row.names(pData(cds)),]$parent

    mst_branch_nodes <- V(minSpanningTree(cds))[which(degree(minSpanningTree(cds)) > 2)]$name

    minSpanningTree(cds) <- dp_mst
    cds@auxOrderingData[[cds@dim_reduce_type]]$cell_ordering_tree <- as.undirected(order_list$cell_ordering_tree)

  } else if (cds@dim_reduce_type == "DDRTree"){
    if (is.null(num_paths) == FALSE){
      message("Warning: num_paths only valid for method 'ICA' in reduceDimension()")
    }
    cc_ordering <- extract_ddrtree_ordering(cds, root_cell)

    if(ncol(cds) > 100) {
      R <- cds@auxOrderingData$DDRTree$R
      edge <- data.frame(start = 1:nrow(R), end = apply(R, 1, which.max), weight = R[cbind(1:nrow(R), apply(R, 1, which.max))])

      pData(cds)$Pseudotime <- cc_ordering[edge$end, 'pseudo_time']
      if(is.null(root_state) == TRUE) {
        pData(cds)$State <- cc_ordering[edge$end, 'cell_state']
      }
      
      cds@auxOrderingData[[cds@dim_reduce_type]]$root_cell <- root_cell
    } else {
      cds@auxOrderingData <- new.env( hash=TRUE )
      pData(cds)$Pseudotime <-  cc_ordering[row.names(pData(cds)),]$pseudo_time

      K_old <- reducedDimK(cds)
      old_dp <- cellPairwiseDistances(cds)
      old_mst <- minSpanningTree(cds)
      old_A <- reducedDimA(cds)
      old_W <- reducedDimW(cds)

      cds <- project2MST(cds, project_point_to_line_segment) #project_point_to_line_segment can be changed into other states
      minSpanningTree(cds) <- cds@auxOrderingData[[cds@dim_reduce_type]]$pr_graph_cell_proj_tree

      root_cell_idx <- which(V(old_mst)$name == root_cell, arr.ind=T)
      cells_mapped_to_graph_root <- which(cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex == root_cell_idx)
      if(length(cells_mapped_to_graph_root) == 0) { #avoid the issue of multiple cells projected to the same point on the principal graph
        cells_mapped_to_graph_root <- root_cell_idx
      }

      cells_mapped_to_graph_root <- V(minSpanningTree(cds))[cells_mapped_to_graph_root]$name

      tip_leaves <- names(which(degree(minSpanningTree(cds)) == 1))
      root_cell <- cells_mapped_to_graph_root[cells_mapped_to_graph_root %in% tip_leaves][1]
      if(is.na(root_cell)) {
        root_cell <- select_root_cell(cds, root_state, reverse)
      }

      cds@auxOrderingData[[cds@dim_reduce_type]]$root_cell <- root_cell

      cc_ordering_new_pseudotime <- extract_ddrtree_ordering(cds, root_cell) #re-calculate the pseudotime again

      pData(cds)$Pseudotime <- cc_ordering_new_pseudotime[row.names(pData(cds)),]$pseudo_time
      if (is.null(root_state) == TRUE) {
        closest_vertex <- cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex
        pData(cds)$State <- cc_ordering[closest_vertex[, 1],]$cell_state #assign the state to the states from the closet vertex
      }

      reducedDimK(cds) <-  K_old
      cellPairwiseDistances(cds) <- old_dp
      minSpanningTree(cds) <- old_mst
      reducedDimA(cds) <- old_A
      reducedDimW(cds) <- old_W
    }

    mst_branch_nodes <- V(minSpanningTree(cds))[which(degree(minSpanningTree(cds)) > 2)]$name
  } else if (cds@dim_reduce_type == "SimplePPT"){
    if (is.null(num_paths) == FALSE){
      message("Warning: num_paths only valid for method 'ICA' in reduceDimension()")
    }
    cc_ordering <- extract_ddrtree_ordering(cds, root_cell)

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
  FM <- FM[is.finite(fm_rowsums) | fm_rowsums != 0, ]

  if (verbose)
    message("Remove noise by PCA ...")
  
  irlba_res <- sparse_prcomp_irlba(t(FM), n = min(num_dim, min(dim(FM)) - 1),
                            center = scaling, scale. = scaling)
  irlba_pca_res <- irlba_res$x
  reducedDimA(cds) <- t(irlba_pca_res) # get top 50 PCs, which can be used for louvain clustering later 
  cds@auxOrderingData[["PCA"]]$irlba_pca_res <- irlba_pca_res

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
#' @param norm_method Determines how to transform expression values prior to reducing dimensionality
#' @param residualModelFormulaStr A model formula specifying the effects to subtract from the data before clustering.
#' @param pseudo_expr amount to increase expression values before dimensionality reduction
#' @param relative_expr When this argument is set to TRUE (default), we intend to convert the expression into a relative expression.
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
                            reduction_method=c("DDRTree", "ICA", 'tSNE', "UMAP", "UMAPDDRTree", "UMAPSSE", "SSE", "SimplePPT", 'L1-graph', 'L1graph'),
                            norm_method = c("log", "vstExprs", "none"),
                            residualModelFormulaStr=NULL,
                            pseudo_expr=1,
                            relative_expr=TRUE,
                            auto_param_selection = TRUE,
                            verbose=FALSE,
                            scaling = TRUE,
                            ...){
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
  FM = FM[is.finite(fm_rowsums) & fm_rowsums != 0,]
  
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
    if (reduction_method == "tSNE") {
      #first perform PCA
      if (verbose)
          message("Remove noise by PCA ...")

      # # Calculate the variance across genes without converting to a dense
      # # matrix:
      # FM_t <- Matrix::t(FM)
      # cell_means <- Matrix::rowMeans(FM_t)
      # cell_vars <- Matrix::rowMeans((FM_t - cell_means)^2)
      # Filter out genes that are constant across all cells:
      #genes_to_keep <- expression_vars > 0
      #FM <- FM[genes_to_keep,]
      #expression_means <- expression_means[genes_to_keep]
      #expression_vars <- expression_vars[genes_to_keep]
      # Here✬s how to take the top PCA loading genes, but using
      # sparseMatrix operations the whole time, using irlba.


      if("num_dim" %in% names(extra_arguments)){ #when you pass pca_dim to the function, the number of dimension used for tSNE dimension reduction is used
        num_dim <- extra_arguments$num_dim #variance_explained
      }
      else{
        num_dim <- 50
      }
      
      irlba_res <- sparse_prcomp_irlba(t(FM), n = min(num_dim, min(dim(FM)) - 1),
                                center = scaling, scale. = scaling)
      irlba_pca_res <- irlba_res$x

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
      # cds@auxClusteringData[["tSNE"]]$reduced_dimension <- t(topDim_pca)
      #cds@auxClusteringData[["tSNE"]]$variance_explained <- prop_varex

      cds@dim_reduce_type <- "tSNE"

      pData(cds)$tsne_1 = reducedDimA(cds)[1,]
      pData(cds)$tsne_2 = reducedDimA(cds)[2,]
    }
    else if (reduction_method %in% c("DDRTree", "UMAPDDRTree")) {
      # FM <- as.matrix(Matrix::t(scale(Matrix::t(FM))))
      # FM <- FM[!is.na(row.names(FM)), ]
      
      # when num_dim is passed or the number of cells is more than 5 k cells and the feature number is large than 50 (implying the feature is not PCA space), do an intial PCA 
      if("num_dim" %in% names(extra_arguments) | (ncol(FM) > 5000 & nrow(FM) > 50)) { 
        if("num_dim" %in% names(extra_arguments)){ #when you pass pca_dim to the function, the number of dimension used for tSNE dimension reduction is used
          num_dim <- extra_arguments$num_dim #variance_explained
        }
        else{
          num_dim <- 50
        }

        if (verbose)
          message("Remove noise by PCA ...")
        
        irlba_res <- sparse_prcomp_irlba(t(FM), n = min(num_dim, min(dim(FM)) - 1),
                                  center = scaling, scale. = scaling)
        irlba_pca_res <- t(irlba_res$x)
        colnames(irlba_pca_res) <- colnames(FM)
        FM <- irlba_pca_res #[, 1:num_dim]
      }else{
        if(scaling){
          FM <- as.matrix(Matrix::t(scale(Matrix::t(FM))))
          FM <- FM[!is.na(row.names(FM)), ]
        } else FM <- as.matrix(FM)
      }

      if(reduction_method == "UMAPDDRTree") {
        umap_args <- c(list(X = t(FM), log = F, n_component = as.integer(max_components), verbose = verbose, return_all = T),
                               extra_arguments[names(extra_arguments) %in% 
                               c("python_home", "n_neighbors", "metric", "n_epochs", "negative_sample_rate", "alpha", "init", "min_dist", "spread", 
                                'set_op_mix_ratio', 'local_connectivity', 'bandwidth', 'gamma', 'a', 'b', 'random_state', 'metric_kwds', 'angular_rp_forest', 'verbose')])
        tmp <- do.call(UMAP, umap_args)
        umap_res <- t(tmp$embedding_); 

        adj_mat <- Matrix::sparseMatrix(i = tmp$graph$indices, p = tmp$graph$indptr, 
                        x = -as.numeric(tmp$graph$data), dims = c(ncol(cds), ncol(cds)), index1 = F, 
                        dimnames = list(colnames(cds), colnames(cds)))

        colnames(umap_res) <- colnames(FM)
        FM <- umap_res
        reducedDimA(cds) <- FM

        cds@auxOrderingData[["DDRTree"]]$adj_mat <- adj_mat
      }


      if (verbose)
        message("Learning principal graph with DDRTree")
      
      if(reduction_method == "UMAPDDRTree"){
        no_reduction=TRUE
      }else{
        no_reduction=FALSE
      }
      
      # TODO: DDRTree should really work with sparse matrices.
      if(auto_param_selection & ncol(cds) >= 100){
        if("ncenter" %in% names(extra_arguments)) #avoid overwrite the ncenter parameter
          ncenter <- extra_arguments$ncenter
        else
          ncenter <- cal_ncenter(ncol(FM))
        #add other parameters...
        
        ddr_args <- c(list(X=FM, dimensions=max_components, ncenter=ncenter, no_reduction=no_reduction, verbose = verbose),
                      extra_arguments[names(extra_arguments) %in% c("initial_method", "maxIter", "sigma", "lambda", "param.gamma", "tol")])
        #browser()
        ddrtree_res <- do.call(DDRTree, ddr_args)
      } else{
        ddrtree_res <- DDRTree(FM, max_components, no_reduction=no_reduction, verbose = verbose, ...)
      }
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
    }else if(reduction_method == "L1-graph") {
        if("num_dim" %in% names(extra_arguments)){ #when you pass pca_dim to the function, the number of dimension used for tSNE dimension reduction is used
          num_dim <- extra_arguments$num_dim #variance_explained
        }
        else{
          num_dim <- 50
        }
        
        if (verbose)
          message("Remove noise by PCA ...")
        
        if(verbose)
          message('running PCA (no further scaling or center) ...')
        irlba_res <- sparse_prcomp_irlba(t(FM), n = min(num_dim, min(dim(FM)) - 1),
                                  center = scaling, scale. = scaling)
        irlba_pca_res <- as.matrix(t(irlba_res$x))
        colnames(irlba_pca_res) <- colnames(FM)
        FM <- irlba_pca_res
        
        umap_args <- c(list(X = t(FM), log = F, n_component = as.integer(max_components), verbose = verbose, return_all = T),
                       extra_arguments[names(extra_arguments) %in% 
                                         c("python_home", "n_neighbors", "metric", "n_epochs", "negative_sample_rate", "alpha", "init", "min_dist", "spread", 
                                'set_op_mix_ratio', 'local_connectivity', 'bandwidth', 'gamma', 'a', 'b', 'random_state', 'metric_kwds', 'angular_rp_forest', 'verbose')])
        tmp <- do.call(UMAP, umap_args)
        umap_res <- t(tmp$embedding_); 
        
        adj_mat <- Matrix::sparseMatrix(i = tmp$graph$indices, p = tmp$graph$indptr, 
                                        x = -as.numeric(tmp$graph$data), dims = c(ncol(cds), ncol(cds)), index1 = F, 
                                        dimnames = list(colnames(cds), colnames(cds)))
        
        colnames(umap_res) <- colnames(FM)
        FM <- umap_res
        reducedDimA(cds) <- FM
        
        louvain_clustering_args <- c(list(data = t(umap_res), pd = pData(cds)[colnames(FM), ], verbose = verbose),
                                     extra_arguments[names(extra_arguments) %in% c("k", "weight", "louvain_iter")])
        louvain_res <- do.call(louvain_clustering, louvain_clustering_args)
        cds@auxOrderingData[["L1graph"]]$louvain_module <- as.factor(igraph::membership(louvain_res$optim_res))
        #cds@auxOrderingData[["L1graph"]]$louvain_res <- louvain_res
        #cds@auxOrderingData[["L1graph"]]$adj_mat <- adj_mat
        
        reduced_dim_res = FM 
        
        if("ncenter" %in% names(extra_arguments)){ #avoid overwrite the ncenter parameter
          ncenter <- extra_arguments$ncenter
        }else{
          if("L1.pr_graph_vertex_per_louvain_module" %in% names(extra_arguments)){ #avoid overwrite the ncenter parameter
            L1.pr_graph_vertex_per_louvain_module <- extra_arguments$L1.pr_graph_vertex_per_louvain_module
          }else{
            L1.pr_graph_vertex_per_louvain_module = 3
          }
          #ncenter <- cal_ncenter(ncol(FM))
          ncenter = L1.pr_graph_vertex_per_louvain_module * length(levels(cds@auxOrderingData[["L1graph"]]$louvain_module))
          ncenter = min(ncol(FM) / 2, ncenter)
        }
         
        
        if (ncenter > ncol(reduced_dim_res))
          stop("Error: ncenters must be less than or equal to ncol(X)")
        
        centers = t(reduced_dim_res)[seq(1, ncol(reduced_dim_res), length.out=ncenter),]
        
        kmean_res <- kmeans(t(reduced_dim_res), ncenter, centers=centers, iter.max = 100)
        if (kmean_res$ifault != 0){
          message(paste("Warning: kmeans returned ifault =", kmean_res$ifault))
        }
        nearest_center = findNearestVertex(t(kmean_res$centers), FM, process_targets_in_blocks=TRUE)
        medioids = reduced_dim_res[,unique(nearest_center)]
        reduced_dim_res <- t(medioids)
        #reduced_dim_res = t(reduced_dim_res)[centers,]
        #reduced_dim_res <- t(reduced_dim_res)
        #rownames(reduced_dim_res) = paste("Y_", 1:nrow(reduced_dim_res), sep = "")
      
      if(verbose)
        message('running L1-graph ...')
      
      #X <- t(reduced_dim_res)
      X <- t(reduced_dim_res)
      # D <- nrow(X); N <- ncol(X)
      # Z <- X
      
      if('C0' %in% names(extra_arguments)){
        C0 <- extra_arguments$C0
      }
      else
        C0 <- X
      Nz <- ncol(C0)
      
      
      #G_T = get_mst(C0)
      
      # # print(extra_arguments)
      # if('nn' %in% names(extra_arguments))
      #   G_knn <- get_knn(C0, K = extra_arguments$nn)
      # else
      #   G_knn <- get_knn(C0, K = 5)
      
      # print(extra_arguments)
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
      names(louvain_component) = colnames(FM)
      louvain_component = louvain_component[rownames(reduced_dim_res)]
      louvain_component = as.factor(louvain_component)
      if (length(levels(louvain_component)) > 1){
        louvain_component_mask = as.matrix(tcrossprod(sparse.model.matrix( ~ louvain_component + 0)))
        
        G = G * louvain_component_mask
        W = W * louvain_component_mask
        rownames(G) = rownames(W)
        colnames(G) = colnames(W)
      }
      #louvain_component = data.frame(component=louvain_component)
     
      l1graph_args <- c(list(X = t(reduced_dim_res), C0 = C0, G = G, gstruct = 'l1-graph', verbose = verbose),
                        extra_arguments[names(extra_arguments) %in% c('maxiter', 'eps', 'L1.lambda', 'L1.gamma', 'L1.sigma', 'nn')])
      
      
      l1_graph_res <- do.call(principal_graph, l1graph_args)
      
      colnames(l1_graph_res$C) <-  rownames(reduced_dim_res)
      DCs <- FM
      
      colnames(l1_graph_res$W) <- rownames(reduced_dim_res)
      rownames(l1_graph_res$W) <- rownames(reduced_dim_res)
      
      
      # row.names(l1_graph_res$X) <- colnames(cds)
      reducedDimW(cds) <- l1_graph_res$W
      reducedDimS(cds) <- DCs
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
    }else if (reduction_method %in% c("UMAP") ) {
      # FM <- as.matrix(Matrix::t(scale(Matrix::t(FM))))
      # FM <- FM[!is.na(row.names(FM)), ]
      if(c("PCA") %in% names(cds@auxOrderingData)) {
        # if UMAP or L1graph has already done, use those information 
        irlba_pca_res <- cds@auxOrderingData$PCA$irlba_pca_res
      } else {
        stop('Please first project the high dimension data to PCA space with projectPCA function before applying UMAP dimension reduction!')
      }
     
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
      
      minSpanningTree(cds) <- louvain_res$g

      A <- S
      colnames(A) <- colnames(FM)
      reducedDimA(cds) <- A

      colnames(S) <- colnames(FM)
      colnames(Y) <- colnames(FM)
      reducedDimW(cds) <- W 
      reducedDimS(cds) <- as.matrix(Y)
      reducedDimK(cds) <- S

      cds@auxOrderingData$UMAP <- list(umap_res = umap_res, louvain_res = louvain_res, adj_mat = adj_mat)
      cds@dim_reduce_type <- reduction_method
    } else if(reduction_method %in% c("SSE")) {
      landmark_id <- NULL

      if(c("UMAP") %in% names(cds@auxOrderingData)) {
        # if UMAP or L1graph has already done, use those information 
        irlba_pca_res <- cds@auxOrderingData$PCA$irlba_pca_res
        S <- t(cds@auxOrderingData$UMAP$umap_res)
        umap_res <- cds@auxOrderingData$UMAP$umap_res
        data <- umap_res
        adj_mat <- cds@auxOrderingData$UMAP$adj_mat
      } else {
         irlba_pca_res <- cds@auxOrderingData$PCA$irlba_pca_res
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
        kmean_res <- kmeans(data_ori, landmark_num, centers=centers, iter.max = 100)
        landmark_id <- sort(unique(apply(as.matrix(proxy::dist(data_ori, kmean_res$centers)), 2, which.min))) # avoid duplicated points 

        # data_ori <- as(data_ori, 'sparseMatrix')  
        # landmark_res <- monocle:::select_landmarks(data_ori@x, data_ori@i, data_ori@p, data_ori@Dim[2], data_ori@Dim[1], landmark_num)

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

        # adj_mat <- NULL # adj_mat_ori[landmark_id, landmark_id]  # # adj_mat_ori[landmark_id, landmark_id] 
        # adj_mat <- adj_mat_ori[landmark_id, landmark_id]  # # adj_mat_ori[landmark_id, landmark_id] 

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

      # reduced_dim_res = t(SSE_res$Y) 
      
      # if("ncenter" %in% names(extra_arguments)){ #avoid overwrite the ncenter parameter
      #   ncenter <- extra_arguments$ncenter
      # }else{
      #   #ncenter <- cal_ncenter(ncol(FM))
      #   ncenter = 3 * length(levels(cds@auxOrderingData[[reduction_method]]$louvain_module))
      #   ncenter = min(ncol(FM) / 2, ncenter)
      # }
               
      # if (ncenter > ncol(reduced_dim_res))
      #   stop("Error: ncenters must be less than or equal to ncol(X)")
      
      # centers = t(reduced_dim_res)[seq(1, ncol(reduced_dim_res), length.out=ncenter),]

      # kmean_res <- kmeans(t(reduced_dim_res), ncenter, centers=centers, iter.max = 100)
      # medioids = reduced_dim_res[, unique(apply(as.matrix(proxy::dist(t(reduced_dim_res), kmean_res$centers)), 2, which.min))] # avoid duplicated points 
      # reduced_dim_res <- t(medioids)

      # if(verbose)
      #   message("Running generalized SimplePPT algorithm ...")

      # if('C0' %in% names(extra_arguments)){
      #   C0 <- extra_arguments$C0
      # }
      # else
      #   C0 <- medioids
      # Nz <- ncol(C0)
      
      # # print(extra_arguments)
      # if('nn' %in% names(extra_arguments))
      #   G <- get_knn(C0, K = extra_arguments$nn)
      # else
      #   G <- get_knn(C0, K = 5)

      # if("louvain_qval" %in% names(extra_arguments)){ 
      #   louvain_qval <- extra_arguments$louvain_qval 
      # }
      # else{
      #   louvain_qval <- 0.05
      # }

      # # reduced_dim_res <- t(SSE_res$Y) 
      # cluster_graph_res <- compute_louvain_connected_components(louvain_res$g, louvain_res$optim_res, louvain_qval, verbose)
      # louvain_component = components(cluster_graph_res$cluster_g)$membership[louvain_res$optim_res$membership]
      # cds@auxOrderingData[[reduction_method]]$louvain_component = louvain_component
      # names(louvain_component) = colnames(FM)
      # louvain_component = louvain_component[rownames(reduced_dim_res)]
      # louvain_component = as.factor(louvain_component)
      # if (length(levels(louvain_component)) > 1){
      #   louvain_component_mask = as.matrix(tcrossprod(sparse.model.matrix( ~ louvain_component + 0)))
        
      #   G$G = G$G * louvain_component_mask
      #   G$W = G$W * louvain_component_mask
      #   rownames(G$G) = rownames(G$W)
      #   colnames(G$G) = colnames(G$W)
      # }

      # pg_args <- c(list(X = t(reduced_dim_res), C0 = C0, G = G$G, gstruct = c("l1-graph"),  verbose = verbose), #span-tree
      #               extra_arguments[names(extra_arguments) %in% c("maxiter", "eps", "L1.lambda", "L1.gamma", "L1.sigma", "nn")]) # "l1-graph", 
      # l1_graph_res <- do.call(principal_graph, pg_args)
   
      # colnames(l1_graph_res$C) <-  rownames(reduced_dim_res)
      # Y <- t(SSE_res$Y)

      # # now let us project other cells back to the landmark space 
      projection_res1 <- NULL
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

        # tmp <- mclapply(setdiff(1:ncol(cds), landmark_id), function(x) {
        #     dist_mat <- igraph::distances(g, v = colnames(cds)[x], to = colnames(cds)[landmark_id])

        #     # rank top 5 landmark cells
        #     top_5 <- sort(dist_mat, index.return = T, decreasing = FALSE) # find the closest cells 
        #     valid_id <- which(is.finite(top_5$x[1:5]))
        #     bandwidth <- mean(range(top_5$x[valid_id])) # half of the range of the 5 nearest landmark as bindwidth 
        #     p <- exp(-top_5$x[valid_id]/bandwidth) # Gaussian kernel 
        #     weight <- p / sum(p)  
        #     Y[, top_5$ix[valid_id]] %*% matrix(weight, ncol = 1)            
              
        #   }, mc.cores = cores)

        # tmp <- do.call(cbind, tmp)
        # projection_res[, setdiff(1:ncol(cds), landmark_id)] <- tmp
        # rm(g, dist_mat, links, relations)

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
      # minSpanningTree(cds) <- gp
      # colnames(l1_graph_res$W) <- rownames(reduced_dim_res)
      # rownames(l1_graph_res$W) <- rownames(reduced_dim_res)
      # reducedDimW(cds) <- l1_graph_res$W

      # reducedDimS(cds) <- Y
      # reducedDimK(cds) <- l1_graph_res$C
      # K <- l1_graph_res$C

      # cds@auxOrderingData[[reduction_method]]$objective_vals <- tail(l1_graph_res$objs, 1)
      # cds@auxOrderingData[[reduction_method]]$W <- l1_graph_res$W
      # cds@auxOrderingData[[reduction_method]]$P <- l1_graph_res$P
      
      # adjusted_K <- Matrix::t(reducedDimK(cds))
      # dp <- as.matrix(dist(adjusted_K))
      # cellPairwiseDistances(cds) <- dp
      
      # W <- l1_graph_res$W
      # dimnames(l1_graph_res$W) <- list(paste('cell_', 1:nrow(W), sep = ''), paste('cell_', 1:nrow(W), sep = ''))
      # W[W < 1e-5] <- 0
      # gp <- graph.adjacency(W, mode = "undirected", weighted = TRUE)
      # # dp_mst <- minimum.spanning.tree(gp)
      # minSpanningTree(cds) <- gp
      # cds@dim_reduce_type <- reduction_method
      # cds <- findNearestPointOnMST(cds)

      # cluster_graph_res <- cluster_graph(minSpanningTree(cds), louvain_res$g, louvain_res$optim_res, t(as.matrix(Y)))
   
      # pData(cds)$louvain_cluster <- as.character(igraph::membership(louvain_res$optim_res)) 
      cds@auxOrderingData[[reduction_method]] <- list(SSE_res = SSE_res, projection_res1 = projection_res1, 
        louvain_module = as.factor(igraph::membership(louvain_res$optim_res)), 
        louvain_res_ori = louvain_res_ori, louvain_res = louvain_res, adj_mat = adj_mat, landmark_id = landmark_id)

      cds@dim_reduce_type <- reduction_method
    } else if(reduction_method == "L1graph") {
      if(c("SSE") %in% names(cds@auxOrderingData)) {
        Y <- cds@auxOrderingData[['SSE']]$SSE_res$Y
        Y <- Y 
        # SSE_res$Y <- (SSE_res$Y - min(SSE_res$Y)) / max(SSE_res$Y) # normalize the SSE space to avoid space shrinking in L1graph step 
        
        louvain_res <- cds@auxOrderingData[["SSE"]]$louvain_res
        louvain_module_length = length(levels(cds@auxOrderingData[['SSE']]$louvain_module))

        landmark_id <- cds@auxOrderingData[['SSE']]$landmark_id
        if(is.null(landmark_id)) {
          reduced_dim_res = t(Y) 
          row.names(Y) <- colnames(FM)
        } else {
          reduced_dim_res = t(Y)           
          row.names(Y) <- colnames(FM[, landmark_id])
      }
      } else if(c("UMAP") %in% names(cds@auxOrderingData)) {
        Y <- cds@auxOrderingData[['UMAP']]$umap_res
        row.names(Y) <- colnames(FM)
        reduced_dim_res = t(Y) 

        louvain_res <- cds@auxOrderingData[["UMAP"]]$louvain_res
        louvain_module_length = length(levels(cds@auxOrderingData[['UMAP']]$louvain_module))

        landmark_id <- cds@auxOrderingData[['UMAP']]$landmark_id
      } else {
        stop('L1graph can be only applied to either the MAP or SSE Ureduced space, please first apply those dimension reduction techniques!')
      }
       
      if("ncenter" %in% names(extra_arguments)){ #avoid overwrite the ncenter parameter
        ncenter <- extra_arguments$ncenter
      }else{
        if("L1.pr_graph_vertex_per_louvain_module" %in% names(extra_arguments)){ #avoid overwrite the ncenter parameter
          L1.pr_graph_vertex_per_louvain_module <- extra_arguments$L1.pr_graph_vertex_per_louvain_module
        }else{
          L1.pr_graph_vertex_per_louvain_module <- 3
        }
        #ncenter <- cal_ncenter(ncol(FM))
        ncenter <- L1.pr_graph_vertex_per_louvain_module * louvain_module_length
        ncenter <- min(ncol(FM) / 2, ncenter)
      }

      if (ncenter > ncol(reduced_dim_res))
        stop("Error: ncenters must be less than or equal to ncol(X)")
      
      centers = t(reduced_dim_res)[seq(1, ncol(reduced_dim_res), length.out=ncenter),]

      kmean_res <- kmeans(t(reduced_dim_res), ncenter, centers=centers, iter.max = 100)
      medioids = reduced_dim_res[, unique(apply(as.matrix(proxy::dist(t(reduced_dim_res), kmean_res$centers)), 2, which.min))] # avoid duplicated points 
      reduced_dim_res <- t(medioids)

      if(verbose)
        message("Running generalized SimplePPT algorithm ...")

      if('C0' %in% names(extra_arguments)){
        C0 <- extra_arguments$C0
      }
      else
        C0 <- medioids
      Nz <- ncol(C0)
      
      # print(extra_arguments)
      if('nn' %in% names(extra_arguments))
        G <- get_knn(C0, K = extra_arguments$nn)
      else
        G <- get_knn(C0, K = 5)

      if("louvain_qval" %in% names(extra_arguments)){ 
        louvain_qval <- extra_arguments$louvain_qval 
      }
      else{
        louvain_qval <- 0.05
      }

      # reduced_dim_res <- t(SSE_res$Y) 
      cluster_graph_res <- compute_louvain_connected_components(louvain_res$g, louvain_res$optim_res, louvain_qval, verbose)
      louvain_component = components(cluster_graph_res$cluster_g)$membership[louvain_res$optim_res$membership]
      cds@auxOrderingData[[reduction_method]]$louvain_component = louvain_component
      names(louvain_component) = colnames(FM[, landmark_id])
      louvain_component = louvain_component[rownames(reduced_dim_res)]
      louvain_component = as.factor(louvain_component)
      if (length(levels(louvain_component)) > 1){
        louvain_component_mask = as.matrix(tcrossprod(sparse.model.matrix( ~ louvain_component + 0)))
        
        G$G = G$G * louvain_component_mask
        G$W = G$W * louvain_component_mask
        rownames(G$G) = rownames(G$W)
        colnames(G$G) = colnames(G$W)
      }

      pg_args <- c(list(X = t(reduced_dim_res), C0 = C0, G = G$G, gstruct = c("l1-graph"),  verbose = verbose), #span-tree
                    extra_arguments[names(extra_arguments) %in% c("maxiter", "eps", "L1.lambda", "L1.gamma", "L1.sigma", "nn")]) # "l1-graph", 
      l1_graph_res <- do.call(principal_graph, pg_args)
   
      colnames(l1_graph_res$C) <-  rownames(reduced_dim_res)
      Y <- t(Y)

      # # now let us project other cells back to the landmark space 
      # if(ncol(cds) > 5000) {
      #   projection_res <- matrix(0, nrow = max_components, ncol = ncol(cds))

      #   # get a graph with distance between cells as the weight 
      #   relations <- louvain_res_ori$relations
      #   relations$weight <- reshape2::melt(t(louvain_res_ori$distMatrix))[, 3]
      #   g <- igraph::graph.data.frame(relations, directed = FALSE)

      #   # iterate each non-landmark cells and project it to the SSE space
      #   if("cores" %in% names(extra_arguments)) {
      #     cores <- extra_arguments$landmark_num
      #   } else {
      #     cores <- detectCores() 
      #   }

      #   if(verbose) 
      #     message('project other non-landmark cells to the landmark SSE space ...')
        
      #   # using matrix multiplication to accelerate the process (optimize it to handle millions of points)
      #   block_size <- 50000
      #   num_blocks = ceiling(nrow(data_ori) / block_size)
      #   weight_mat <- NULL
      #   for (i in 1:num_blocks){
      #     if (i < num_blocks){
      #       block = data_ori[((((i-1) * block_size)+1):(i*block_size)), ]
      #     }else{
      #       block = data_ori[((((i-1) * block_size)+1):(nrow(data_ori))),]
      #     }
      #     distances_Z_to_Y <- proxy::dist(block, data_ori[landmark_id, ])

      #     tmp <- as(t(apply(distances_Z_to_Y , 1, function(x) {
      #         tmp <- sort(x)[6] 
      #         x[x > tmp] <- 0; p <- rep(0, length(x))
      #         bandwidth <- mean(range(x[x > 0])) # half of the range of the nearest neighbors as bindwidth 
      #         p[x > 0] <- exp(-x[x > 0]/bandwidth) # Gaussian kernel 
      #         p / sum(p) 
      #     })), 'sparseMatrix')

      #     weight_mat <- rBind(weight_mat, tmp)
      #   }

      #   projection_res1 <- t(weight_mat %*% as(t(Y), 'sparseMatrix'))

      #   tmp <- mclapply(setdiff(1:ncol(cds), landmark_id), function(x) {
      #       dist_mat <- igraph::distances(g, v = colnames(cds)[x], to = colnames(cds)[landmark_id])

      #       # rank top 5 landmark cells
      #       top_5 <- sort(dist_mat, index.return = T, decreasing = FALSE) # find the closest cells 
      #       valid_id <- which(is.finite(top_5$x[1:5]))
      #       bandwidth <- mean(range(top_5$x[valid_id])) # half of the range of the 5 nearest landmark as bindwidth 
      #       p <- exp(-top_5$x[valid_id]/bandwidth) # Gaussian kernel 
      #       weight <- p / sum(p)  
      #       Y[, top_5$ix[valid_id]] %*% matrix(weight, ncol = 1)            
              
      #     }, mc.cores = cores)

      #   tmp <- do.call(cbind, tmp)
      #   projection_res[, setdiff(1:ncol(cds), landmark_id)] <- tmp
      #   rm(g, dist_mat, links, relations)

      #   projection_res[, landmark_id] <- Y
      #   Y <- as.matrix(projection_res)

      #   FM <- FM_ori
      # }

      # colnames(Y) <- colnames(FM)[, ]
      # reducedDimA(cds) <- Y

      colnames(l1_graph_res$W) <- rownames(reduced_dim_res)
      rownames(l1_graph_res$W) <- rownames(reduced_dim_res)
      reducedDimW(cds) <- l1_graph_res$W

      # reducedDimS(cds) <- Y
      colnames(l1_graph_res$C) <- paste('cell_', 1:ncol(l1_graph_res$C), sep = '')
      reducedDimK(cds) <- l1_graph_res$C
      K <- l1_graph_res$C

      cds@auxOrderingData[[reduction_method]]$objective_vals <- tail(l1_graph_res$objs, 1)
      cds@auxOrderingData[[reduction_method]]$W <- l1_graph_res$W
      cds@auxOrderingData[[reduction_method]]$P <- l1_graph_res$P
      
      adjusted_K <- Matrix::t(reducedDimK(cds))
      dp <- as.matrix(dist(adjusted_K))
      cellPairwiseDistances(cds) <- dp
      
      dimnames(l1_graph_res$W) <- list(paste('cell_', 1:nrow(l1_graph_res$W), sep = ''), paste('cell_', 1:nrow(l1_graph_res$W), sep = ''))
      W <- l1_graph_res$W
      W[W < 1e-5] <- 0
      gp <- graph.adjacency(W, mode = "undirected", weighted = TRUE)
      # dp_mst <- minimum.spanning.tree(gp)
      minSpanningTree(cds) <- gp
      cds@dim_reduce_type <- reduction_method
      cds <- findNearestPointOnMST(cds)

      # cluster_graph_res <- cluster_graph(minSpanningTree(cds), louvain_res$g, louvain_res$optim_res, t(as.matrix(Y)))
   
      # pData(cds)$louvain_cluster <- as.character(igraph::membership(louvain_res$optim_res)) 
      cds@auxOrderingData[[reduction_method]] <- list(PG_res = l1_graph_res, landmark_id = landmark_id)
    } else {
      stop("Error: unrecognized dimensionality reduction method")
    }
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

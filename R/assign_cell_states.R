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

#' Assign a state for each cell.
#'
#' @description This function first decomposes the learned principal graph for the data into separate graph components. 
#' It then applies depth first traversal (DFS) for each graph component sequentially and assigns each principal point 
#' a State number. The State number increases every time we find a node has degree larger than 2 or when we 
#' move to the next graph component. Thus the state assignment are more meaningful if we run learnGraph with 
#' DDRTree, SimplePPT (and without loop clousre by setting close_loop = T) but not with L1graph. The cells inherit
#' the corresponding principal points' state (pseudotime, see below) value. The first cell 
#' with degree equal to 2 (the first cell if there is no cells has degree equal to 2) identified is used as the 
#' root cell. Note that this function also assigns pseudotime by default. But the orderCells function should be 
#' used in order to obtain correct ordering of the developmental trajectory under study. 
#'
#' @param cds the CellDataSet upon which to perform this operation
#' 
#' @importFrom stats dist
#' @importFrom igraph decompose.graph V graph.dfs
#'
#' @return an updated CellDataSet object, in which phenoData contains values for State and Pseudotime for each cell
#' @export
#' @examples
#' \dontrun{
#' lung <- load_lung()
#' lung <- assign_cell_states(lung)
#' }
assign_cell_states <- function(cds){
  
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
  
  # iterate over each graph component and assign branches, pseduotime for each component 
  if(is.null(cds@minSpanningTree)) {
    stop('Error: please run partitionCells and learnGraph before assign cell states')
  }

  dp_mst <- cds@minSpanningTree
  dp <- as.matrix(dist(t(cds@reducedDimK))) 
  g_list <- decompose.graph(dp_mst) 
  curr_state <- 1

  states = rep(1, ncol(dp))
  names(states) <- V(dp_mst)$name

  pseudotimes = rep(0, ncol(dp))
  names(pseudotimes) <- V(dp_mst)$name

  parents = rep(NA, ncol(dp))
  names(parents) <- V(dp_mst)$name

  for(cur_g in g_list) {
    root_cell <- which(neighborhood.size(cur_g) == 2)[1] 
    if(is.na(root_cell)) {
      root_cell <- V(cur_g)$name[1]
    } 

    mst_traversal <- graph.dfs(cur_g,
                               root = root_cell,
                               neimode = "all",
                               unreachable=FALSE,
                               father=TRUE)
    mst_traversal$father <- as.numeric(mst_traversal$father)

    for (i in 1:length(mst_traversal$order)){
      curr_node <- mst_traversal$order[i]
      curr_node_name <- V(cur_g)[curr_node]$name

      if (is.na(mst_traversal$father[curr_node]) == FALSE){
        parent_node <- mst_traversal$father[curr_node]
        parent_node_name <- V(cur_g)[parent_node]$name
        parent_node_pseudotime <- pseudotimes[parent_node_name]
        parent_node_state <- states[parent_node_name]
        curr_node_pseudotime <- parent_node_pseudotime + dp[curr_node_name, parent_node_name]
        if (degree(cur_g, v=parent_node_name) > 2){
          curr_state <- curr_state + 1
        }
      }else{
        parent_node = NA
        parent_node_name = NA
        curr_node_pseudotime = 0
      }

      curr_node_state <- curr_state
      pseudotimes[curr_node_name] <- curr_node_pseudotime
      states[curr_node_name] <- curr_node_state
      parents[curr_node_name] <- parent_node_name
    }
    
    curr_state <- curr_state + 1 # update after each graph component 
    
  }

  ordering_df <- data.frame(sample_name = names(states),
                            cell_state = as.character(states),
                            pseudo_time = as.vector(pseudotimes),
                            parent = parents, stringsAsFactors = F)
  row.names(ordering_df) <- ordering_df$sample_name
  
  pData(cds)$State <- NULL # reset state 
  pr_graph_cell_proj_closest_vertexordering_df <- cds@auxOrderingData[[cds@dim_reduce_type]]$pr_graph_cell_proj_closest_vertex
  pData(cds)[row.names(pr_graph_cell_proj_closest_vertexordering_df), 'State'] <- ordering_df[paste0('Y_', pr_graph_cell_proj_closest_vertexordering_df[, 1]), 'cell_state']
  pData(cds)[row.names(pr_graph_cell_proj_closest_vertexordering_df), 'Pseudotime'] <- ordering_df[paste0('Y_', pr_graph_cell_proj_closest_vertexordering_df[, 1]), 'pseudo_time']
 
  # Ensure states follows a consectutive sequence 
  tmp <- 1:length(sort(unique(cds$State)))
  names(tmp) <- sort(unique(cds$State))
  cds$State <- tmp[cds$State]
  pData(cds)$State <- as.factor(pData(cds)$State)
  
  cds@auxOrderingData[[cds@dim_reduce_type]]$branch_points <- V(dp_mst)[which(degree(dp_mst) > 2)]$name

  cds
}




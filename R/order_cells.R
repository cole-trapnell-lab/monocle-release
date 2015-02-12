
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
    child_of_p_node <- V(canonical_pq) [ nei(p_node, mode="out") ]
    parent_of_p_node <- V(canonical_pq) [ nei(p_node, mode="in") ]
    
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

extract_ordering <- function(pq_tree, curr_node)
{
  nei <- NULL
  
  if (V(pq_tree)[curr_node]$type == "leaf")
  {
    return (V(pq_tree)[curr_node]$name)
  }
  else if (V(pq_tree)[curr_node]$type == "P")
  {
    p_level <- list()
    for (child in V(pq_tree) [ nei(curr_node, mode="out") ])
    {
      p_level[[length(p_level)+1]] <- extract_ordering(pq_tree, child)
    }
    p_level <- p_level[sample(length(p_level))]
    p_level <- unlist(p_level)
    #print (p_level)
    return (p_level)
  }
  else if(V(pq_tree)[curr_node]$type == "Q")
  {
    q_level <- list()
    for (child in V(pq_tree) [ nei(curr_node, mode="out") ])
    {
      q_level[[length(q_level)+1]] <- extract_ordering(pq_tree, child)
    }
    if (runif(1) >= 0.5)
    {
      q_level <- rev(q_level)
    }
    q_level <- unlist(q_level)
    #print (q_level)
    return (q_level)
  }
}

extract_fixed_ordering <- function(pq_tree, curr_node)
{
  nei <- NULL
  
  if (V(pq_tree)[curr_node]$type == "leaf")
  {
    return (V(pq_tree)[curr_node]$name)
  }
  else if (V(pq_tree)[curr_node]$type == "P")
  {
    p_level <- list()
    for (child in V(pq_tree) [ nei(curr_node, mode="out") ])
    {
      p_level[[length(p_level)+1]] <- extract_ordering(pq_tree, child)
    }
    #p_level <- p_level[sample(length(p_level))]
    p_level <- unlist(p_level)
    #print (p_level)
    return (p_level)
  }
  else if(V(pq_tree)[curr_node]$type == "Q")
  {
    q_level <- list()
    for (child in V(pq_tree) [ nei(curr_node, mode="out") ])
    {
      q_level[[length(q_level)+1]] <- extract_ordering(pq_tree, child)
    }
    # if (runif(1) >= 0.5)
    # {
    #     q_level <- rev(q_level)
    # }
    q_level <- unlist(q_level)
    #print (q_level)
    return (q_level)
  }
}

weight_of_p_node_order<-function(p_level_list, dist_matrix)
{
  cost <- 0
  if (length(p_level_list) <= 1)
  {
    return (0)
  }
  #print (p_level_list)
  #print (length(p_level_list))
  for (i in (1:(length(p_level_list)-1)))
  {
    #print (i)
    #print ("***")
    #print (p_level_list[[i]])
    #print (paste("...", p_level_list[[i]][length(p_level_list[[i]])]))
    #print (p_level_list[[i+1]])
    #print (paste("...", p_level_list[[i+1]][1]))
    #print (dist_matrix[p_level_list[[i]][length(p_level_list[[i]])], p_level_list[[i+1]][1]])
    cost <- cost + dist_matrix[p_level_list[[i]][length(p_level_list[[i]])], p_level_list[[i+1]][1]]
  }
  return(cost)
}

#' Return an ordering for a P node in the PQ tree
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

#order_q_node(q_level, dp)

extract_good_ordering <- function(pq_tree, curr_node, dist_matrix)
{
  nei <- NULL
  
  if (V(pq_tree)[curr_node]$type == "leaf")
  {
    #print ("ordering leaf node")
    return (V(pq_tree)[curr_node]$name)
  }else if (V(pq_tree)[curr_node]$type == "P"){
    #print ("ordering P node")
    p_level <- list()
    for (child in V(pq_tree) [ nei(curr_node, mode="out") ])
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
    for (child in V(pq_tree) [ nei(curr_node, mode="out") ])
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
    for (child in V(pq_tree) [ nei(curr_node, mode="out") ])
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
    for (child in V(pq_tree) [ nei(curr_node, mode="out") ])
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
    for (child in V(pq_tree) [ nei(curr_node, mode="out") ])
    {
      node_states <- assign_cell_lineage(pq_tree, child, assigned_state, node_states)
    }
    return(node_states)
  }
}

#' Extract a linear ordering of cells from a PQ tree
#' @importFrom plyr arrange
extract_good_branched_ordering <- function(orig_pq_tree, curr_node, dist_matrix, num_branches, reverse_main_path=FALSE)
{
  nei <- NULL
  type <- NULL
  pseudo_time <- NULL
  
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
    parents <- V(pq_tree)[nei(names(branch_node_counts)[i], mode="in")]
    if (length(parents) > 0 && parents$type == "P")
    {
      p_node_parent <- V(pq_tree)[nei(names(branch_node_counts)[i], mode="in")]
      parent_branch_id <- V(pq_tree)[nei(p_node_parent, mode="in")]$name
      #print (parent_branch_id)
      #print (branch_id)
      branch_tree <- branch_tree + edge(parent_branch_id, branch_id)
    }
    pq_tree[V(pq_tree) [ nei(names(branch_node_counts)[i], mode="in") ], names(branch_node_counts)[i] ] <- FALSE
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
    
    for (child in V(branch_tree) [ nei(curr_branch, mode="out") ])
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
    
    children <- V(cell_tree) [ nei(curr_cell, mode="out") ]
    ordering_tree_res$subtree <- cell_tree
    
    if (length(children) == 1){
      ordering_tree_res <- assign_cell_state_helper(ordering_tree_res, V(cell_tree)[children]$name)
    }else{
      for (child in children)	{
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
    V(cell_tree)[curr_cell]$parent =  V(cell_tree)[ nei(curr_cell, mode="in") ]$name
    #print (curr_cell_pseudotime)
    
    ordering_tree_res$subtree <- cell_tree
    children <- V(cell_tree) [ nei(curr_cell, mode="out") ]
    
    for (child in children)	{
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
  return(ordering_df)
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
    s <- irlba::irlba(V, min(n,p), min(n,p))  
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


#' Orders cells according to progress through a learned biological process.
#' @param cds the CellDataSet upon which to perform this operation
#' @param num_paths the number of end-point cell states to allow in the biological process.
#' @param reverse whether to reverse the beginning and end points of the learned biological process.
#' @param root_cell the name of a cell to use as the root of the ordering tree.
#' @return an updated CellDataSet object, in which phenoData contains values for State and Pseudotime for each cell
#' @export
orderCells <- function(cds, num_paths=1, reverse=FALSE, root_cell=NULL){
  
  adjusted_S <- t(cds@reducedDimS)
  
  dp <- as.matrix(dist(adjusted_S))
  
  #print (sum(rowSums(dp)))
  #dp <- as.matrix(dist(dp))
  #dp <- as.matrix(dist(adjusted_S))
  cellPairwiseDistances(cds) <- as.matrix(dist(adjusted_S))
  # Build an MST of the cells in ICA space.
  gp <- graph.adjacency(dp, mode="undirected", weighted=TRUE)
  dp_mst <- minimum.spanning.tree(gp)
  minSpanningTree(cds) <- dp_mst
  # Build the PQ tree
  next_node <<- 0
  res <- pq_helper(dp_mst, use_weights=FALSE, root_node=root_cell)
  #stopifnot(length(V(res$subtree)[type == "leaf"]) == nrow(pData(cds)))
  
  cc_ordering <- extract_good_branched_ordering(res$subtree, res$root, cellPairwiseDistances(cds), num_paths, reverse)
  row.names(cc_ordering) <- cc_ordering$sample_name
  
  pData(cds)$Pseudotime <-  cc_ordering[row.names(pData(cds)),]$pseudo_time
  pData(cds)$State <-  cc_ordering[row.names(pData(cds)),]$cell_state
  pData(cds)$Parent <-  cc_ordering[row.names(pData(cds)),]$parent
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
                            use_vst=F, ...){
  
  FM <- exprs(cds)
  
  # If we aren't using VST, then normalize the expression values by size factor
  if (use_vst == FALSE && cds@expressionFamily@vfamily == "negbinomial")
  {
    checkSizeFactors(cds)
    size_factors <- sizeFactors(cds)
    #print (size_factors)
    FM <- t(t(FM) / size_factors)
    #FM <- log2(FM)
  }
  
  if (is.null(fData(cds)$use_for_ordering) == FALSE)
    FM <- FM[fData(cds)$use_for_ordering,]
  
  FM <- FM + pseudo_expr
  FM <- FM[matrixStats::rowSds(FM) > 0,]
  
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
  
  # TODO: get rid of this if possible.
  if (is.null(batch) == FALSE || is.null(batch2) == FALSE|| is.null(covariates) == FALSE)
  {
    message("Removing batch effects")
    #FM <- log2(FM)
    FM <- limma::removeBatchEffect(FM, batch=batch, batch2=batch2, covariates=covariates)
    if (use_vst == FALSE) {
      FM <- 2^FM
    }
  }
  
  #FM <- log2(FM)
  
  message("Reducing to independent components")
  
  #FM <- t(scale(t(FM)))
  #FM <- FM[rowSds(FM) > 0,]
  init_ICA <- ica_helper(t(FM), max_components, use_irlba=use_irlba, ...)
  
  x_pca <- t(t(FM) %*% init_ICA$K)
  W <- t(init_ICA$W)
  
  weights <- W
  
  # print(dim (init_ICA$K))
  # print(dim (solve(weights)))
  
  A <- t(solve(weights) %*% t(init_ICA$K))
  
  colnames(A) <- colnames(weights)
  rownames(A) <- rownames(FM)
  
  S <- weights %*% x_pca
  
  rownames(S) <- colnames(weights)
  colnames(S) <- colnames(FM) 
  
  reducedDimW(cds) <- W
  reducedDimA(cds) <- A
  reducedDimS(cds) <- S
  reducedDimK(cds) <- init_ICA$K
  
  cds
}




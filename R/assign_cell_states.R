
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
#' @param minimal_branch_len the minimal length of the principal tree segment to be treated as a true branch. Default to be 10. 
#' @importFrom stats dist
#' @importFrom igraph graph.empty neighborhood.size decompose.graph V graph.dfs E add_edges add_vertices degree delete_vertices
#'
#' @return an updated CellDataSet object, in which phenoData contains values for State and Pseudotime for each cell
#' @export
#' @examples
#' \dontrun{
#' lung <- load_lung()
#' lung <- assign_cell_states(lung)
#' }
setStates <- function(cds, minimal_branch_len = 10){
  
  if(class(cds)[1] != "CellDataSet") {
    stop("Error cds is not of type 'CellDataSet'")
  }
  
  if (is.null(cds@dim_reduce_type)){
    stop("Error: dimensionality not yet reduced. Please call reduceDimension() and learnGraph() (for learning principal graph) before calling this function.")
  }
  if (is.null(cds@rge_method)){
    stop("Error: principal graph has not learned yet. Please call learnGraph() before calling this function.")
  }
  # reducedDimA, S, and K are not NULL in the cds
  if (length(cds@reducedDimS) == 0) {
    stop("Error: dimension reduction didn't prodvide correct results. Please check your reduceDimension() step and ensure correct dimension reduction are performed before calling this function.")
  }
  if (length(cds@reducedDimK) == 0) {
    stop("Error: principal graph learning didn't prodvide correct results. Please check your learnGraph() step and ensure correct principal graph learning are performed before calling this function.")
  }
  if(igraph::vcount(minSpanningTree(cds)) > 10000) {
    stop("orderCells doesn't support more than 10k centroids (cells)")
  }

  # iterate over each graph component and assign branches, pseduotime for each component 
  if(is.null(cds@minSpanningTree)) {
    stop('Error: please run partitionCells and learnGraph before assign cell states')
  }

  principal_points_coord <- cds@reducedDimK
  dp_mst <- cds@minSpanningTree
  dp <- as.matrix(dist(t(principal_points_coord))) 
  g_list <- decompose.graph(dp_mst) 
  curr_state <- 1

  states = rep(1, ncol(dp))
  names(states) <- V(dp_mst)$name

  pseudotimes = rep(0, ncol(dp))
  names(pseudotimes) <- V(dp_mst)$name

  parents = rep(NA, ncol(dp))
  names(parents) <- V(dp_mst)$name

  # create a state graph: 
  # vertex: - name: cell state; properties: (1) cells in the state 
  # edge: - name: branch principal point:  properties: (1) true branch point? 
  state_graph <- graph.empty(n=0, directed=TRUE)

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

          state_graph <- add_vertices(state_graph, nv = 1, attr = list(name = as.character(curr_state), cell_names = list(curr_node_name) ))
          state_graph <- add_edges(state_graph, c(curr_state - 1, curr_state), attr = list(name = parent_node_name)) # , branch == TRUE
        } else {
          # add names for new cell - Note that we need to use list so that the length of the value for 'cell_names' is 1 
          state_graph <- set.vertex.attribute(state_graph, 'cell_names', curr_state, list(c(V(state_graph)[curr_state]$cell_names[[1]], curr_node_name) )) 
        }
      }else{
        state_graph <- add_vertices(state_graph, nv = 1, attr = list(name = as.character(curr_state), cell_names = list(curr_node_name) ))
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

  # prune the graph with small insignificant branches   
  if(minimal_branch_len != 1) {
    branch_points <- NULL 

    g_list <- decompose.graph(state_graph) 
    curr_state <- 1

    for(cur_g in g_list) {
      state_degree <- degree(cur_g)
      
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

          # if the current graph segment is less than minimal_branch_len, use the parent cell's state 
          if(length(V(cur_g)[curr_node]$cell_names[[1]]) < minimal_branch_len) {
            states[V(cur_g)[curr_node]$cell_names[[1]]] <- V(cur_g)[parent_node_name]$cell_state
            V(cur_g)[curr_node]$cell_state <- V(cur_g)[parent_node_name]$cell_state

            # if(state_degree[curr_node] == 1) {
            #   dp_mst <- delete_vertices(dp_mst, V(cur_g)[curr_node]$cell_names[[1]]) 
            #   principal_points_coord <- principal_points_coord[, setdiff(colnames(principal_points_coord), V(cur_g)[curr_node]$cell_names[[1]])]
            # }
          } else {
            curr_state <- curr_state + 1
            states[V(cur_g)[curr_node]$cell_names[[1]]] <- curr_state
            V(cur_g)[curr_node]$cell_state <- curr_state

            # only if the current graph segment is no less than minimal_branch_len, we will append the current node to branch_points set
            branch_points <- c(branch_points, E(cur_g)[parent_node_name %--% curr_node_name]$name)
          }

        }else{
          parent_node = NA

          states[V(cur_g)[curr_node]$cell_names[[1]]] <- curr_state
          V(cur_g)[curr_node]$cell_state <- curr_state
        }

        curr_node_state <- curr_state
      }
      
      curr_state <- curr_state + 1 # update after each graph component 
    }
  } else {
    branch_points <- V(dp_mst)[which(degree(dp_mst) > 2)]$name
  }

  ordering_df <- data.frame(sample_name = names(states),
                            cell_state = as.character(states),
                            pseudo_time = as.vector(pseudotimes),
                            parent = parents, stringsAsFactors = F)
  row.names(ordering_df) <- ordering_df$sample_name
  
  pData(cds)$State <- NULL # reset state 
  pr_graph_cell_proj_closest_vertexordering_df <- cds@auxOrderingData[[cds@rge_method]]$pr_graph_cell_proj_closest_vertex
  pData(cds)[row.names(pr_graph_cell_proj_closest_vertexordering_df), 'State'] <- ordering_df[paste0('Y_', pr_graph_cell_proj_closest_vertexordering_df[, 1]), 'cell_state']
  # pData(cds)[row.names(pr_graph_cell_proj_closest_vertexordering_df), 'Pseudotime'] <- ordering_df[paste0('Y_', pr_graph_cell_proj_closest_vertexordering_df[, 1]), 'pseudo_time']
 
  # Ensure states follows a consectutive sequence 
  tmp <- 1:length(sort(unique(cds$State)))
  names(tmp) <- sort(unique(cds$State))
  cds$State <- tmp[cds$State]
  pData(cds)$State <- as.factor(pData(cds)$State)
  
  cds@auxOrderingData[[cds@rge_method]]$branch_points <- unique(branch_points) 

  # if(state_degree[curr_node] == 1) {
  #   dp_mst <- delete_vertices(dp_mst, V(cur_g)[curr_node]$cell_names[[1]]) 
  #   principal_points_coord <- principal_points_coord[, setdiff(colnames(principal_points_coord), V(cur_g)[curr_node]$cell_names[[1]])]
  # }
  # cds@minSpanningTree <- dp_mst
  # cds@reducedDimK <- principal_points_coord
  cds
}

pruneTree <- function(cds, minimal_branch_len = 10){
  
  if(class(cds)[1] != "CellDataSet") {
    stop("Error cds is not of type 'CellDataSet'")
  }
  
  if (is.null(cds@dim_reduce_type)){
    stop("Error: dimensionality not yet reduced. Please call reduceDimension() and learnGraph() (for learning principal graph) before calling this function.")
  }
  if (is.null(cds@rge_method)){
    stop("Error: principal graph has not learned yet. Please call learnGraph() before calling this function.")
  }
  # reducedDimA, S, and K are not NULL in the cds
  if (length(cds@reducedDimS) == 0) {
    stop("Error: dimension reduction didn't prodvide correct results. Please check your reduceDimension() step and ensure correct dimension reduction are performed before calling this function.")
  }
  if (length(cds@reducedDimK) == 0) {
    stop("Error: principal graph learning didn't prodvide correct results. Please check your learnGraph() step and ensure correct principal graph learning are performed before calling this function.")
  }
  if(igraph::vcount(minSpanningTree(cds)) > 10000) {
    stop("orderCells doesn't support more than 10k centroids (cells)")
  }

  # iterate over each graph component and assign branches, pseduotime for each component 
  if(is.null(cds@minSpanningTree)) {
    stop('Error: please run partitionCells and learnGraph before assign cell states')
  }

  principal_points_coord <- cds@reducedDimK
  dp_mst <- cds@minSpanningTree
  dp <- as.matrix(dist(t(principal_points_coord))) 
  g_list <- decompose.graph(dp_mst) 
  curr_state <- 1

  states = rep(1, ncol(dp))
  names(states) <- V(dp_mst)$name

  pseudotimes = rep(0, ncol(dp))
  names(pseudotimes) <- V(dp_mst)$name

  parents = rep(NA, ncol(dp))
  names(parents) <- V(dp_mst)$name

  # create a state graph: 
  # vertex: - name: cell state; properties: (1) cells in the state 
  # edge: - name: branch principal point:  properties: (1) true branch point? 
  state_graph <- graph.empty(n=0, directed=TRUE)

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
          parent_neighbors <- neighbors(cur_g, v=parent_node_name, mode = 'all')
          # let us first assume the data is a tree for now (latter for loops and combination of loop and tree structure)
          # browser()
          # first identify if there is are loop the loop
          # the edges and vertex need to be removed if they doesn't satisfy the thresold
          tmp <- mst_traversal$order$name[(match(parent_node_name, mst_traversal$order$name) + 1):vcount(cur_g)]
          parent_node_name_other_neighbor <- setdiff(intersect(parent_neighbors$name, tmp), curr_node$name) # the other node of the parent node (not the current point)
          
          parent_neighbors_index <- sort(match(parent_neighbors$name, mst_traversal$order$name)) # follow the order of gene expression 
          parent_neighbors <- mst_traversal$order$name[parent_neighbors_index]
          vertex_name_on_branch_a <- mst_traversal$order$name[parent_neighbors_index[2]:(parent_neighbors_index[3] - 1)]
          
          loop_node <- intersect(neighbors(cur_g, v=parent_node_name_other_neighbor, mode = 'all')$name, vertex_name_on_branch_a) # the other points' other neighbor (not the parent_node_name)
          if(length(loop_node)) { # if it is loop 
            browser()
            tmp <- mst_traversal$order$name[(match(parent_node_name, mst_traversal$order$name)):(match(loop_node, mst_traversal$order$name))]
            degree_ <- degree(cur_g, tmp)

            loop <- induced.subgraph(cur_g, tmp[1:which.min(degree_ > 2)])
            
            if(diameter(loop) > minimal_branch_len) {
              curr_state <- curr_state + 1
              
              state_graph <- add_vertices(state_graph, nv = 1, attr = list(name = as.character(curr_state), cell_names = list(curr_node_name) ))
              state_graph <- add_edges(state_graph, c(curr_state - 1, curr_state), attr = list(name = parent_node_name)) # , branch == TRUE
              
              # branch_points <- c(branch_points, E(cur_g)[parent_node_name %--% curr_node_name]$name)
            }
          }
          else { # not loop 
            
            sub_cur_g_a <- induced_subgraph(cur_g, vertex_name_on_branch_a)
            diameter_len_a <- diameter(sub_cur_g_a) 
            
            # find the previous branch point 
            tmp <- mst_traversal$order$name[1:parent_neighbors_index[1]]
            previous_branch_point <- which(degree(cur_g)[tmp] > 2)
            
            if(length(previous_branch_point) == 0) { # if no previous branch point, this must be the first branch point, so take all other cells as the second branch 
              vertex_name_on_one_branch_b <- mst_traversal$order$name[parent_neighbors_index[3]:(vcount(cur_g))]
            } else {
              previous_branch_point_name <- names(which.max(previous_branch_point)) # [-length(previous_branch_point)] get the name for the previous branch point 
              tmp <- mst_traversal$order$name[parent_neighbors_index[3]:(vcount(cur_g))] # all nodes 
              vertex_name_on_one_branch_b <-  tryCatch({
                vertex_name_on_one_branch_b <- tmp[1:(match(intersect(tmp, neighbors(cur_g, v=previous_branch_point_name, mode = 'all')$name), tmp) - 1)] # restrict to only the current branch
              }, error = function(err) {
                browser()
              })
            }
            
            sub_cur_g_b <- induced_subgraph(cur_g, vertex_name_on_one_branch_b)
            diameter_len_b <- diameter(sub_cur_g_b) 
            
            if(diameter_len_a > minimal_branch_len & diameter_len_b > minimal_branch_len) {
              curr_state <- curr_state + 1
              # browser()
              state_graph <- add_vertices(state_graph, nv = 1, attr = list(name = as.character(curr_state), cell_names = list(curr_node_name) ))
              state_graph <- add_edges(state_graph, c(curr_state - 1, curr_state), attr = list(name = parent_node_name)) # , branch == TRUE
              # branch_points <- c(branch_points, E(cur_g)[parent_node_name %--% curr_node_name]$name)
            }
          }

        } else {
          # add names for new cell - Note that we need to use list so that the length of the value for 'cell_names' is 1 
          state_graph <- set.vertex.attribute(state_graph, 'cell_names', curr_state, list(c(V(state_graph)[curr_state]$cell_names[[1]], curr_node_name) )) 
        }
      }else{
        state_graph <- add_vertices(state_graph, nv = 1, attr = list(name = as.character(curr_state), cell_names = list(curr_node_name) ))
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
  branch_points <- NULL
# 
#   # prune the graph with small insignificant branches   
#   if(minimal_branch_len != 1) {
#     branch_points <- NULL 
# 
#     g_list <- decompose.graph(state_graph) 
#     curr_state <- 1
# 
#     for(cur_g in g_list) {
#       state_degree <- degree(cur_g)
#       
#       root_cell <- which(neighborhood.size(cur_g) == 2)[1] 
#       if(is.na(root_cell)) {
#         root_cell <- V(cur_g)$name[1]
#       } 
#       mst_traversal <- graph.dfs(cur_g,
#                                  root = root_cell,
#                                  neimode = "all",
#                                  unreachable=FALSE,
#                                  father=TRUE)
#       mst_traversal$father <- as.numeric(mst_traversal$father)
# 
#       for (i in 1:length(mst_traversal$order)){
#         curr_node <- mst_traversal$order[i]
#         curr_node_name <- V(cur_g)[curr_node]$name
# 
#         if (is.na(mst_traversal$father[curr_node]) == FALSE){
#           parent_node <- mst_traversal$father[curr_node]
#           parent_node_name <- V(cur_g)[parent_node]$name
# 
#           # if the current graph segment is less than minimal_branch_len, use the parent cell's state 
#           if(length(V(cur_g)[curr_node]$cell_names[[1]]) < minimal_branch_len) {
#             states[V(cur_g)[curr_node]$cell_names[[1]]] <- V(cur_g)[parent_node_name]$cell_state
#             V(cur_g)[curr_node]$cell_state <- V(cur_g)[parent_node_name]$cell_state
# 
#             # if(state_degree[curr_node] == 1) {
#             #   dp_mst <- delete_vertices(dp_mst, V(cur_g)[curr_node]$cell_names[[1]]) 
#             #   principal_points_coord <- principal_points_coord[, setdiff(colnames(principal_points_coord), V(cur_g)[curr_node]$cell_names[[1]])]
#             # }
#           } else {
#             curr_state <- curr_state + 1
#             states[V(cur_g)[curr_node]$cell_names[[1]]] <- curr_state
#             V(cur_g)[curr_node]$cell_state <- curr_state
# 
#             # only if the current graph segment is no less than minimal_branch_len, we will append the current node to branch_points set
#             branch_points <- c(branch_points, E(cur_g)[parent_node_name %--% curr_node_name]$name)
#           }
# 
#         }else{
#           parent_node = NA
# 
#           states[V(cur_g)[curr_node]$cell_names[[1]]] <- curr_state
#           V(cur_g)[curr_node]$cell_state <- curr_state
#         }
# 
#         curr_node_state <- curr_state
#       }
#       
#       curr_state <- curr_state + 1 # update after each graph component 
#     }
#   } else {
#     branch_points <- V(dp_mst)[which(degree(dp_mst) > 2)]$name
#   }

  ordering_df <- data.frame(sample_name = names(states),
                            cell_state = as.character(states),
                            pseudo_time = as.vector(pseudotimes),
                            parent = parents, stringsAsFactors = F)
  row.names(ordering_df) <- ordering_df$sample_name
  
  pData(cds)$State <- NULL # reset state 
  pr_graph_cell_proj_closest_vertexordering_df <- cds@auxOrderingData[[cds@rge_method]]$pr_graph_cell_proj_closest_vertex
  pData(cds)[row.names(pr_graph_cell_proj_closest_vertexordering_df), 'State'] <- ordering_df[paste0('Y_', pr_graph_cell_proj_closest_vertexordering_df[, 1]), 'cell_state']
  # pData(cds)[row.names(pr_graph_cell_proj_closest_vertexordering_df), 'Pseudotime'] <- ordering_df[paste0('Y_', pr_graph_cell_proj_closest_vertexordering_df[, 1]), 'pseudo_time']
 
  # Ensure states follows a consectutive sequence 
  tmp <- 1:length(sort(unique(cds$State)))
  names(tmp) <- sort(unique(cds$State))
  cds$State <- tmp[cds$State]
  pData(cds)$State <- as.factor(pData(cds)$State)
  
  cds@auxOrderingData[[cds@rge_method]]$branch_points <- branch_points 

  # min
  # if(state_degree[curr_node] == 1) {
  #   dp_mst <- delete_vertices(dp_mst, V(cur_g)[curr_node]$cell_names[[1]]) 
  #   principal_points_coord <- principal_points_coord[, setdiff(colnames(principal_points_coord), V(cur_g)[curr_node]$cell_names[[1]])]
  # }
  cds@minSpanningTree <- dp_mst
  cds@reducedDimK <- principal_points_coord
  cds
}

#' Function to prune the graph 
pruneTree_in_learnGraph <- function(stree_ori, stree_loop_clousre, minimal_branch_len = 10){
  if (ncol(stree_ori) < minimal_branch_len)
    return(stree_loop_clousre);
  dimnames(stree_loop_clousre) <- dimnames(stree_ori)
  stree_ori[stree_ori != 0] <- 1
  stree_ori <- graph_from_adjacency_matrix(stree_ori, mode = 'undirected', weight = NULL)
  stree_loop_clousre[stree_loop_clousre != 0] <- 1
  stree_loop_clousre <- graph_from_adjacency_matrix(stree_loop_clousre, mode = 'undirected', weight = NULL)
  
  # get closed loops: 
  added_edges <- get.edgelist(stree_loop_clousre - stree_ori)
  valid_edges <- matrix(ncol = 2, nrow = 0)
  edges_to_remove_df <- matrix(ncol = 2, nrow = 0)
  vertex_top_keep <- NULL
  
  if(nrow(added_edges) > 0) {
    edge_dists <- apply(added_edges, 1, function(x) distances(stree_ori, x[1], x[2]))
    valid_edges <- added_edges[which(edge_dists >= minimal_branch_len), , drop = F]
    edges_to_remove_df <- added_edges[which(edge_dists < minimal_branch_len), , drop = F]
  }
  if(nrow(valid_edges) > 0) {
    vertex_top_keep <- as.character(unlist(apply(valid_edges, 1, function(x) shortest_paths(stree_ori, x[1], x[2])$vpath[[1]]$name )))
  }
  
  root_cell <- which(neighborhood.size(stree_ori) == 2)[1] 
  if(is.na(root_cell)) {
    root_cell <- V(stree_ori)$name[1]
  } 

  mst_traversal <- graph.dfs(stree_ori,
                             root = root_cell,
                             neimode = "all",
                             unreachable=FALSE,
                             father=TRUE)
  mst_traversal$father <- as.numeric(mst_traversal$father)
  vertex_to_be_deleted <- c()

  # further remove other small branches? 
  if (length(mst_traversal$order) > 0){
    for (i in 1:length(mst_traversal$order)){
      curr_node <- tryCatch({mst_traversal$order[i]}, error = function(e) { NA })
      if (is.na(curr_node))
        next;
      curr_node_name <- V(stree_ori)[curr_node]$name
      
      if (is.na(mst_traversal$father[curr_node]) == FALSE){
        parent_node <- mst_traversal$father[curr_node]
        parent_node_name <- V(stree_ori)[parent_node]$name
        
        if (degree(stree_ori, v=parent_node_name) > 2){
          parent_neighbors <- neighbors(stree_ori, v=parent_node_name, mode = 'all')
          
          parent_neighbors_index <- sort(match(parent_neighbors$name, mst_traversal$order$name)) # follow the order of gene expression 
          parent_neighbors <- mst_traversal$order$name[parent_neighbors_index]
          
          tmp <- delete.edges(stree_ori, paste0(parent_node_name, "|", parent_neighbors))
          tmp_decomposed <- decompose.graph(tmp)
          
          comp_a <- tmp_decomposed[unlist(lapply(tmp_decomposed, function(x) {
            parent_neighbors[2] %in% V(x)$name
          }))][[1]]
          
          comp_b <- tmp_decomposed[unlist(lapply(tmp_decomposed, function(x) {
            parent_neighbors[3] %in% V(x)$name
          }))][[1]]
          
          diameter_len_a <- diameter(comp_a) + 1
          diameter_len_b <- diameter(comp_b) + 1
          
          if(diameter_len_a < minimal_branch_len) {# if loop closure is not applied to cells on this branch 
            vertex_to_be_deleted <- c(vertex_to_be_deleted, V(comp_a)$name)
          } else {
            # browser()
            # vertex_top_keep <- c(vertex_top_keep, get_diameter(comp_a)$name)
          }
          if(diameter_len_b < minimal_branch_len) {# if loop closure is not applied to cells on this branch 
            # browser()
            
            vertex_to_be_deleted <- c(vertex_to_be_deleted, V(comp_b)$name)
          } else {
            # vertex_top_keep <- c(vertex_top_keep, get_diameter(comp_b)$name)
          }
          
        }
      }
    }
  }
  
  # browser()
  valid_vertex_to_be_deleted <- setdiff(vertex_to_be_deleted, vertex_top_keep)
  stree_loop_clousre <- delete_vertices(stree_loop_clousre, valid_vertex_to_be_deleted)

  tmp <- edges_to_remove_df[edges_to_remove_df[, 1] %in% V(stree_loop_clousre)$name & edges_to_remove_df[, 2] %in% V(stree_loop_clousre)$name, , drop = FALSE]
  if(nrow(tmp) > 0) {
    edges_to_remove <- paste0(tmp[, 1], '|', tmp[, 2])
    stree_loop_clousre <- delete.edges(stree_loop_clousre, edges_to_remove)
  }

  return(get.adjacency(stree_loop_clousre))
}


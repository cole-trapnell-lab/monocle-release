#' Build a CellDataSet with appropriate duplication along two lineages
#' @param cds CellDataSet for the experiment
#' @param progenitor_method The method to treat the progenitor cells, there are three options: sequential_split, random_split, duplicate. 
#' They represent that the progenitors can be either splitted one by one (sequential_split) or randomly (random_split) or duplicate the 
#' progenitors (duplicate). For random_split and sequential_split, the cells will be evenly split (when number of progenitors are odd, 
#' there will be one more cell in one lineage)
#' @param branch_point The name of a branch point that specifies the lineages. Provide either this argument or lineage_states.
#' @param lineage_states The states for two branching lineages
#' @param lineage_labels The names for each branching lineage
#' @param stretch A logic flag to determine whether or not the pseudotime trajectory for each lineage should be stretched to the same range or not 
#' @return a CellDataSet with the duplicated cells and stretched lineages
#' @export
#'
buildLineageBranchCellDataSet <- function(cds, 
                                          progenitor_method = c('sequential_split','random_split', 'duplicate'), 
                                          lineage_states = NULL, 
                                          branch_point = 1,
                                          lineage_labels = NULL, 
                                          stretch = FALSE)
{
  # TODO: check that lineages are on the same paths
  
  if(is.null(pData(cds)$State) | is.null(pData(cds)$Pseudotime)) 
    stop('Please first order the cells in pseudotime using orderCells()')
  if(!is.null(branch_point) & is.null(lineage_states)) 
    stop('Please either specify the branch_point or lineage_states to select subset of cells')

  if (!is.null(lineage_labels) & !is.null(lineage_states)) {
    if(length(lineage_labels) != length(lineage_states))
      stop("length of lineage_labels doesn't match with that of lineage_states")
    lineage_map <- setNames(lineage_labels, as.character(lineage_states))
  }
  
  if(cds@dim_reduce_type == "DDRTree") {
    pr_graph_cell_proj_mst <- minSpanningTree(cds)
  }
  else {
    pr_graph_cell_proj_mst <- cds@auxOrderingData[[cds@dim_reduce_type]]$cell_ordering_tree
  }

  root_cell <- cds@auxOrderingData[[cds@dim_reduce_type]]$root_cell
  root_state <- pData(cds)[root_cell,]$State
  #root_state <- V(pr_graph_cell_proj_mst)[root_cell,]$State
  
  pr_graph_root <- subset(pData(cds), State == root_state)
  
  if (cds@dim_reduce_type == "DDRTree"){
    closest_vertex <- cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex
    root_cell_point_in_Y <- closest_vertex[row.names(pr_graph_root),]
  }else{
    root_cell_point_in_Y <- row.names(pr_graph_root)
  }
  
  root_cell <- names(which(degree(pr_graph_cell_proj_mst, v = root_cell_point_in_Y, mode = "all")==1, useNames = T))[1]
  
  paths_to_root <- list()
  if (is.null(branch_point)){
    
    # If the user didn't specify a branch point,
    # let's walk back from the lineage states
    for (leaf_state in lineage_states){
      
      curr_cell <- subset(pData(cds), State == leaf_state)
      #Get all the nearest cells in Y for curr_cells:
      
      if (cds@dim_reduce_type == "DDRTree"){
        closest_vertex <- cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex
        curr_cell_point_in_Y <- closest_vertex[row.names(curr_cell),] 
      }else{
        curr_cell_point_in_Y <- row.names(curr_cell)
      }
      
      # Narrow down to a single tip cell in Y:
      curr_cell <- names(which(degree(pr_graph_cell_proj_mst, v = curr_cell_point_in_Y, mode = "all")==1, useNames = T))[1]

      path_to_ancestor <- shortest_paths(pr_graph_cell_proj_mst,curr_cell, root_cell)
      path_to_ancestor <- names(unlist(path_to_ancestor$vpath))
      
      if (cds@dim_reduce_type == "DDRTree"){
        closest_vertex <- cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex
        ancestor_cells_for_lineage <- row.names(closest_vertex)[which(V(pr_graph_cell_proj_mst)[closest_vertex]$name %in% path_to_ancestor)]
      }else if (cds@dim_reduce_type == "ICA"){
        ancestor_cells_for_lineage <- path_to_ancestor
      }
      paths_to_root[[as.character(leaf_state)]] <- ancestor_cells_for_lineage
    }
  }else{
    pr_graph_cell_proj_mst <- minSpanningTree(cds)
    mst_branch_nodes <- cds@auxOrderingData[[cds@dim_reduce_type]]$branch_points
    branch_cell <- mst_branch_nodes[branch_point]
    mst_no_branch_point <- pr_graph_cell_proj_mst - V(pr_graph_cell_proj_mst)[branch_cell]
    
    path_to_ancestor <- shortest_paths(pr_graph_cell_proj_mst, branch_cell, root_cell)
    path_to_ancestor <- names(unlist(path_to_ancestor$vpath))
    
    for (backbone_nei in V(pr_graph_cell_proj_mst)[suppressWarnings(nei(branch_cell))]$name){
      descendents <- bfs(mst_no_branch_point, V(mst_no_branch_point)[backbone_nei], unreachable=FALSE)
      descendents <- descendents$order[!is.na(descendents$order)]
      descendents <- V(mst_no_branch_point)[descendents]$name
      if (root_cell %in% descendents == FALSE){
        path_to_root <- unique(c(path_to_ancestor, branch_cell, descendents))
        
        if (cds@dim_reduce_type == "DDRTree"){
          closest_vertex <- cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex
          path_to_root <- row.names(closest_vertex)[which(V(pr_graph_cell_proj_mst)[closest_vertex]$name %in% path_to_root)]
        }else{
          path_to_root <- path_to_root
        }
        
        paths_to_root[[length(paths_to_root) + 1]] <-  path_to_root
      }
    }
    
    leaf_cells <- setdiff(names(which(degree(pr_graph_cell_proj_mst, v = union(paths_to_root[[1]], paths_to_root[[2]]), mode = "all")==1, useNames = T)), root_cell)
    lineage_states <- pData(cds)[leaf_cells, ]$State
  }
  all_cells_in_subset <- c()
  
  for (path_to_ancestor in paths_to_root){
    if (length(path_to_ancestor) == 0){
      stop("Error: common ancestors between selected State values on path to root State")
    }
    all_cells_in_subset <- c(all_cells_in_subset, path_to_ancestor)
  }
  all_cells_in_subset <- unique(all_cells_in_subset)
  
  # FIXME: This is a slow, terrible way of doing things.
  common_ancestor_cells <- intersect(paths_to_root[[1]], paths_to_root[[2]])
  if (length(paths_to_root) > 2){
    for (i in seq(3,length(paths_to_root))){
      common_ancestor_cells <- intersect(intersect(paths_to_root[i], paths_to_root[i-1]), common_ancestor_cells)
    }
  }
  
  cds <- cds[, row.names(pData(cds[,all_cells_in_subset]))] #or just union(ancestor_cells, lineage_cells)
  
  #State <- pData(cds)$State 
  Pseudotime <- pData(cds)$Pseudotime 
  
  pData <- pData(cds)

  if(stretch) {
    max_pseudotime <- -1
    for (path_to_ancestor in paths_to_root){
      max_pseudotime_on_path <- max(pData[path_to_ancestor,]$Pseudotime)  
      if (max_pseudotime < max_pseudotime_on_path){
        max_pseudotime <- max_pseudotime_on_path
      }
    }
    
    branch_pseudotime <- max(pData[common_ancestor_cells,]$Pseudotime)
    #ancestor_scaling_factor <- branch_pseudotime / max_pseudotime
    
    for (path_to_ancestor in paths_to_root){
      max_pseudotime_on_path <- max(pData[path_to_ancestor,]$Pseudotime) 
      path_scaling_factor <-(max_pseudotime - branch_pseudotime) / (max_pseudotime_on_path - branch_pseudotime)
      if (is.finite(path_scaling_factor)){
        lineage_cells <- setdiff(path_to_ancestor, common_ancestor_cells)
        pData[lineage_cells,]$Pseudotime <- ((pData[lineage_cells,]$Pseudotime - branch_pseudotime) * path_scaling_factor + branch_pseudotime)
      }
    }
    #pData[common_ancestor_cells,]$Pseudotime <- pData[common_ancestor_cells,]$Pseudotime / max_pseudotime
    
    pData$Pseudotime <- 100 * pData$Pseudotime / max_pseudotime
  }
  pData$original_cell_id <- row.names(pData)
  
  pData$original_cell_id <- row.names(pData)
  pData[common_ancestor_cells, "State"] <- lineage_states[1] #set progenitors to the lineage 1

  if (progenitor_method == 'duplicate') {
    # for (i in 1:(length(lineage_states) - 1)) { #duplicate progenitors for multiple branches
    #   if (nrow(exprs_data) == 1)
    #       exprs_data <- cbind(exprs_data, t(as.matrix(exprs_data[,
    #           progenitor_ind])))
    #   else exprs_data <- cbind(exprs_data, exprs_data[, progenitor_ind])
    #   weight_vec <- c(weight_vec, rep(weight_constant, length(progenitor_ind)))

    #   colnames(exprs_data)[(ncol(exprs_data) - length(progenitor_ind) + 1):ncol(exprs_data)] <- 
    #     paste('duplicate', lineage_states[i], 1:length(progenitor_ind), sep = '_')
    #   pData <- rbind(pData, pData[progenitor_ind, ])
      
    #   pData$State[(length(pData$State) - length(progenitor_ind) + 1):length(pData$State)] <- lineage_states[i + 1]
    #   row.names(pData)[(length(pData$State) - length(progenitor_ind) + 1):length(pData$State)] <- 
    #     paste('duplicate', lineage_states[i], 1:length(progenitor_ind), sep = '_')
    # }
    ancestor_exprs <- exprs(cds)[,common_ancestor_cells]
    expr_blocks <- list()
  
    # Duplicate the expression data
    for (i in 1:length(paths_to_root)) { #duplicate progenitors for multiple branches
      if (nrow(ancestor_exprs) == 1)
        exprs_data <- t(as.matrix(ancestor_exprs))
      else exprs_data <- ancestor_exprs
      
      colnames(exprs_data) <- paste('duplicate', i, 1:length(common_ancestor_cells), sep = '_')
      expr_lineage_data <- exprs(cds)[,setdiff(paths_to_root[[i]], common_ancestor_cells)]
      exprs_data <- cbind(exprs_data, expr_lineage_data)
      expr_blocks[[i]] <- exprs_data
    }
    
    # Make a bunch of copies of the pData entries from the common ancestors
    ancestor_pData_block <- pData[common_ancestor_cells,]
    
    pData_blocks <- list()
    
    weight_vec <- c()
    for (i in 1:length(paths_to_root)) {
      #weight_vec <- c(weight_vec, rep(weight_constant, length(common_ancestor_cells)))
      weight_vec_block <- rep(weight_constant, length(common_ancestor_cells))
      
      #pData <- rbind(pData, pData[common_ancestor_cells, ])
      new_pData_block <- ancestor_pData_block
      new_pData_block$Lineage <- lineage_states[i]
      new_pData_block$State <- lineage_states[i]
      
      row.names(new_pData_block) <- paste('duplicate', i, 1:length(common_ancestor_cells), sep = '_')
      
      pData_lineage_cells <- pData[setdiff(paths_to_root[[i]], common_ancestor_cells),]
      pData_lineage_cells$Lineage <- lineage_states[i]
      pData_lineage_cells$State <- lineage_states[i]
      
      weight_vec_block <- c(weight_vec_block, rep(1, nrow(pData_lineage_cells)))
      
      weight_vec <- c(weight_vec, weight_vec_block)
      
      new_pData_block <- rbind(new_pData_block, pData_lineage_cells)
      pData_blocks[[i]] <- new_pData_block
    }
    pData <- do.call(rbind, pData_blocks)
    exprs_data <- do.call(cbind, expr_blocks)
  }
  else if(progenitor_method == 'random_split') {
    if(length(lineage_states) != 2)
      stop('more than 2 lineage states are used!')

    lineageA <- sample(progenitor_ind, round(length(progenitor_ind) / 2))
    pData[lineageA, 'State'] <- lineage_states[1]
    lineageB <- setdiff(progenitor_ind, lineageA)
    pData[lineageB, 'State'] <- lineage_states[2]    
  }
  else if(progenitor_method == 'sequential_split') {
    if(length(lineage_states) != 2)
      stop('more than 2 lineage states are used!')

    # progenitor_pseudotime_order <- order(pData[progenitor_ind, 'Pseudotime'])

    # lineageA <- progenitor_pseudotime_order[seq(1, length(progenitor_ind), by = 2)]
    # pData[progenitor_ind[lineageA], 'State'] <- lineage_states[1]
    # lineageB <- progenitor_pseudotime_order[seq(2, length(progenitor_ind), by = 2)]
    # pData[progenitor_ind[lineageB], 'State'] <- lineage_states[2]    

    progenitor_pseudotime_order <- order(pData[common_ancestor_cells, 'Pseudotime'])

    lineageA <- progenitor_pseudotime_order[seq(1, length(common_ancestor_cells), by = 2)]
    pData[common_ancestor_cells[lineageA], 'State'] <- lineage_states[1]
    lineageB <- progenitor_pseudotime_order[seq(2, length(common_ancestor_cells), by = 2)]
    pData[common_ancestor_cells[lineageB], 'State'] <- lineage_states[2]   

    zero_pseudotime_root_cell <- common_ancestor_cells[progenitor_pseudotime_order[1]]
    exprs_data <- cbind(exprs(cds), 'duplicate_root' = exprs(cds)[, zero_pseudotime_root_cell])
    pData <- rbind(pData, pData[zero_pseudotime_root_cell, ])
    row.names(pData)[nrow(pData)] <- 'duplicate_root'
    pData[nrow(pData), 'State'] <- lineage_states[2]
    weight_vec <- rep(1, nrow(pData))
  }

  if (!is.null(lineage_labels))
    pData$Lineage <- as.factor(lineage_map[as.character(pData$State)])
  else
    pData$Lineage <- as.factor(pData$State)
  
  pData$State <- factor(pData(cds)[as.character(pData$original_cell_id),]$State, 
                        levels =levels(cds$State))
  pData$weight <- weight_vec
  Size_Factor <- pData$Size_Factor
  
  fData <- fData(cds)
  
  colnames(exprs_data) <- row.names(pData) #check this 
  cds_subset <- newCellDataSet(as.matrix(exprs_data),
                               phenoData = new("AnnotatedDataFrame", data = pData),
                               featureData = new("AnnotatedDataFrame", data = fData),
                               expressionFamily=cds@expressionFamily,
                               lowerDetectionLimit=cds@lowerDetectionLimit)
  pData(cds_subset)$State <- as.factor(pData(cds_subset)$State)
  pData(cds_subset)$Size_Factor <- Size_Factor
  
  cds_subset@dispFitInfo <- cds@dispFitInfo
  
  return (cds_subset)
}

#' Peform the branching test
#'
#' This function is used to perform the branching test to ask the question about whether or not the genes under tests are significant 
#' lineage dependent genes. If stretch equal to TRUE, each lineage is firstly stretched into maturation level 0-100 and the progenitor 
#' cells are duplicated and assigned to each lineage. This test can be used to detect lineage dependent genes. 
#'
#' @param cds a CellDataSet object upon which to perform this operation
#' @param fullModelFormulaStr a formula string specifying the full model in differential expression tests (i.e. likelihood ratio tests) for each gene/feature.
#' @param reducedModelFormulaStr a formula string specifying the reduced model in differential expression tests (i.e. likelihood ratio tests) for each gene/feature.
#' @param lineage_states  states corresponding to two branches 
#' @param branch_point the id for the branch point choosed to select cell branches and progenitors 
#' @param relative_expr a logic flag to determine whether or not the relative gene expression should be used
#' @param stretch  a logic flag to determine whether or not each lineage should be stretched
#' @param cores the number of cores to be used while testing each gene for differential expression
#' @param lineage_labels the name for each lineage, for example, AT1 or AT2  
#' @param exprs_thrsld_percentage For the gene under test, the percentage of cells expressed across all cells. Default is 0.05
#' @param verbose Whether to show VGAM errors and warnings. Only valid for cores = 1. 
#' @param ... Additional arguments passed to differentialGeneTest
#' @return a data frame containing the p values and q-values from the likelihood ratio tests on the parallel arrays of models.
#' @export
#'
branchTest <- function(cds, fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3)*Lineage",
                       reducedModelFormulaStr = "~sm.ns(Pseudotime, df = 3)", 
                       lineage_states = NULL, 
                       branch_point=1,
                       relative_expr = TRUE,
                       stretch = TRUE,
                       cores = 1, 
                       lineage_labels = NULL, 
                       exprs_thrsld_percentage = NULL,
                       verbose = F,
                      #  backup_method = c('nb1', 'nb2'), 
                      #  use_epislon = F,
                      # stepsize = NULL,

                        ...) {
  
  if("Lineage" %in% all.vars(terms(as.formula(fullModelFormulaStr)))) {
    cds_subset <- buildLineageBranchCellDataSet(cds = cds, 
                                                lineage_states = lineage_states,
                                                branch_point=branch_point,
                                                lineage_labels = lineage_labels,
                                                stretch = stretch)
  }
  else
    cds_subset <- cds
  
  branchTest_res <- differentialGeneTest(cds_subset, 
                                         fullModelFormulaStr = fullModelFormulaStr, 
                                         reducedModelFormulaStr = reducedModelFormulaStr, 
                                         cores = cores, 
                                         relative_expr = relative_expr, 
                                         exprs_thrsld_percentage = exprs_thrsld_percentage,
                                         verbose=verbose
                       #                   ,

                       # backup_method = backup_method, 
                       # use_epislon = use_epislon,
                       # stepsize = stepsize

                                         )
  
  return(branchTest_res)
}


#to do: 
# 1. think about how calculating ABCs for multiple lineages (use a common reference as implemented?)
# 2. how to store ABCs? 
# 3. how to use fit_model_helper directly for calculating the model fitting?
#' Calculate the area between TWO fitted lineage trajectories  
#' 
#' This function is used to calculate the ABC score based on the the nature spline curves fitted for each lineage. ABC score is used to 
#' quantify the total magnitude of divergence between two lineages. By default, the ABC score is the area between two fitted spline curves. 
#' The ABC score can be used to rank gene divergence. When coupled with p-val calculated from the branchTest, it can be used to identify
#' potential major regulators for lineage bifurcation. 
#'
#' @param cds a CellDataSet object upon which to perform this operation
#' @param trend_formula a formula string specifying the full model in differential expression tests (i.e. likelihood ratio tests) for each gene/feature.
#' @param trajectory_states States corresponding to two branches 
#' @param relative_expr a logic flag to determine whether or not the relative gene expression should be used
#' @param stretch a logic flag to determine whether or not each lineage should be stretched
#' @param cores the number of cores to be used while testing each gene for differential expression
#' @param verbose a logic flag to determine whether or not we should output detailed running information 
#' @param min_expr the lower limit for the expressed gene
#' @param integer_expression the logic flag to determine whether or not the integer numbers are used for calculating the ABCs. Default is False. 
#' @param num number of points on the fitted lineage trajectories used for calculating the ABCs. Default is 5000. 
#' @param lineage_labels the name for each lineage, for example, AT1 or AT2  
#' @param ... Additional arguments passed to buildLineageBranchCellDataSet
#' @return a data frame containing the ABCs (Area under curves) score as the first column and other meta information from fData
#' @export
calABCs <- function(cds,
                    trend_formula = "~sm.ns(Pseudotime, df = 3)*Lineage",
                    trajectory_states = c(2, 3),
                    relative_expr = TRUE, 
                    stretch = TRUE, 
                    cores = 1, 
                    verbose = F,
                    min_expr = 0.5, 
                    integer_expression = FALSE, 
                    num = 5000, 
                    lineage_labels = NULL,
                    ...){
  
  ABC_method = "integral"
  if (length(trajectory_states) != 2)
    stop("Sorry, this function only supports the calculation of ABCs between TWO lineage trajectories")
  
  
    cds_subset <- buildLineageBranchCellDataSet(cds = cds, #lineage_states = trajectory_states,
                                                lineage_labels = lineage_labels, stretch = stretch)
    overlap_rng <- c(0, max(pData(cds_subset)$Pseudotime))
 
 
  if (length(trajectory_states) != 2)
    stop("calILRs can only work for two Lineages")
  if(!all(trajectory_states %in% pData(cds_subset)[, "Lineage"]))
    stop("state(s) in trajectory_states are not included in 'Lineage'")
  
  if(verbose)
    message(paste("the pseudotime range for the calculation of ILRs:", overlap_rng[1], overlap_rng[2], sep = ' '))
  
  cds_branchA <- cds_subset[, pData(cds_subset)[, "Lineage"] ==
                              trajectory_states[1]]
  cds_branchB <- cds_subset[, pData(cds_subset)[, "Lineage"] ==
                              trajectory_states[2]]
  
  #use the genSmoothCurves to generate smooth curves for calculating the ABC scores:
  formula_all_variables <- all.vars(as.formula(trend_formula))
  
  t_rng <- range(pData(cds_branchA)$Pseudotime)
  str_new_cds_branchA <- data.frame(Pseudotime = seq(overlap_rng[1], overlap_rng[2],
                                                     length.out = num), Lineage = as.factor(trajectory_states[1]))
  colnames(str_new_cds_branchA)[2] <- formula_all_variables[2] #interaction term can be terms rather than Lineage
  if (verbose)
    print(paste("Check the whether or not Pseudotime scaled from 0 to 100: ",
                sort(pData(cds_branchA)$Pseudotime)))
  str_new_cds_branchB <- data.frame(Pseudotime = seq(overlap_rng[1], overlap_rng[2],
                                                     length.out = num), Lineage = as.factor(trajectory_states[2]))
  colnames(str_new_cds_branchB)[2] <- formula_all_variables[2]
  
  if (verbose)
    print(paste("Check the whether or not Pseudotime scaled from 0 to 100: ",
                sort(pData(cds_branchB)$Pseudotime)))
  
  str_branchAB_expression_curve_matrix <- genSmoothCurves(cds_subset, cores=cores, trend_formula = trend_formula,
                                                          relative_expr = relative_expr, new_data = rbind(str_new_cds_branchA, str_new_cds_branchB))
  
  str_branchA_expression_curve_matrix <- str_branchAB_expression_curve_matrix[, 1:num]
  str_branchB_expression_curve_matrix <- str_branchAB_expression_curve_matrix[, (num + 1):(2 * num)]
  
  ABCs_res <- str_branchA_expression_curve_matrix - str_branchB_expression_curve_matrix
  
  ABCs_res <- apply(ABCs_res, 1, function(x, num, ABC_method) {
    avg_delta_x <- (x[1:(num - 1)] + x[2:(num)])/2
    step <- (100/(num - 1))
    
    if (ABC_method == "integral") {
      res <- round(sum(avg_delta_x * step), 3)
    }
    else {
      stop('Current the ABC method only supports integral')
    }
    # else if (ABC_method == "global_normalization") {
    #   max <- max(max(predictBranchOri), max(x))
    #   res <- round(sum(avg_delta_x/max * step), 3)
    # }
    # else if (ABC_method == "local_normalization") {
    #   pair_wise_max <- apply(data.frame(x = x, y = predictBranchOri),
    #                          1, max)
    #   res <- round(sum((((predictBranchOri - x)/pair_wise_max)[1:(num -
    #                                                                 1)] + ((predictBranchOri - x)/pair_wise_max)[2:(num)])/2 *
    #                      step), 3)
    # }
    # else if (ABC_method == "four_values") {
    #   ori_ABCs <- round(sum((x[1:(num - 1)] + x[2:(num)])/2 *
    #                           step), 3)
    #   other_ABCs <- round(sum((predictBranchOri[1:(num -
    #                                                  1)] + predictBranchOri[2:(num)])/2 * step),
    #                       3)
    #   ori_ABCs_H <- round(sum(avg_delta_x[avg_delta_x >
    #                                         0] * step), 3)
    #   other_ABCs_H <- round(sum(avg_delta_x[avg_delta_x <
    #                                           0] * step), 3)
    #   res <- c(ori_ABCs = ori_ABCs, other_ABCs = other_ABCs,
    #            ori_ABCs_H = ori_ABCs_H, other_ABCs_H = other_ABCs_H)
    # }
    # else if (ABC_method == "ILRs") {
    #   str_logfc_df <- log2((predictBranchOri + 1)/(x +
    #                                                  1))
    #   res <- sum(str_logfc_df)
    # }
    return(res)}, num = num, ABC_method = ABC_method
  )
  
  ABCs_res <- merge(ABCs_res, fData(cds), by = "row.names")
  row.names(ABCs_res) <- ABCs_res[, 1]
  ABCs_res[, 1] <- NULL
  colnames(ABCs_res)[1] <- "ABCs"
  
  return(ABCs_res)
}

#' Calculate the Instant Log Ratio between two branching lineages
#' 
#' This function is used to calculate the Instant Log Ratio between two branching lineages which can be used to prepare the heatmap demonstrating the lineage gene expression divergence hirearchy. If "stretch" is specifified, each  
#' lineage will be firstly stretched into maturation level from 0-100. Since the results when we use "stretching" are always better and 
#' IRLs for non-stretched spline curves are often mismatched, we may only turn down "non-stretch" functionality in future versions. Then, we fit two separate nature spline curves for each 
#' individual linages. The log-ratios of the value on each spline curve corresponding to each lineages are calculated, which can be  
#' used as a measure for the magnitude of divergence between two branching lineages. 
#'
#' @param cds CellDataSet for the experiment
#' @param trend_formula trend_formula a formula string specifying the full model in differential expression tests (i.e. likelihood ratio tests) for each gene/feature.
#' @param trajectory_states states corresponding to two branches 
#' @param relative_expr A logic flag to determine whether or not the relative expressed should be used when we fitting the spline curves 
#' @param stretch a logic flag to determine whether or not each lineage should be stretched
#' @param cores Number of cores when fitting the spline curves
#' @param ILRs_limit the minimum Instant Log Ratio used to make the heatmap plot
#' @param label_by_short_name label the rows of the returned matrix by gene_short_name (TRUE) or feature id (FALSE)
#' @param useVST A logic flag to determine whether or not the Variance Stablization Transformation should be used to stablize the gene expression.
#' When VST is used, the difference between two lineages are used instead of the log-ratio.
#' @param round_exprs A logic flag to determine whether or not the expression value should be rounded into integer
#' @param output_type A character either of "all" or "after_bifurcation". If "after_bifurcation" is used, only the time points after the bifurcation point will be selected
#' @param lineage_labels the name for each lineage, for example, AT1 or AT2  
#' @param file the name for storing the data. Since the calculation of the Instant Log Ratio is very time consuming, so by default the result will be stored
#' @param return_all A logic flag to determine whether or not all the results from the analysis should be returned, this includes 
#' a dataframe for the log fold change, normalized log fold change, raw divergence, normalized divergence, fitting curves for each lineage 
#' @param verbose Whether or not detailed running information should be returned 
#' @param ... Additional arguments passed to buildLineageBranchCellDataSet
#' @return a ggplot2 plot object
#' @import ggplot2
#' @importFrom reshape2 melt
#' @export 
#' 
calILRs <- function (cds, 
          trend_formula = "~sm.ns(Pseudotime, df = 3)*Lineage",
          trajectory_states = c(2, 3), 
          relative_expr = TRUE, 
          stretch = T, 
          cores = 1, 
          ILRs_limit = 3, 
          label_by_short_name = TRUE,
          useVST = FALSE, 
          round_exprs = FALSE, 
          output_type = "all", 
          lineage_labels = NULL, 
          file = NULL, 
          return_all = F, 
          verbose = FALSE, 
          ...){
  
 
    cds_subset <- buildLineageBranchCellDataSet(cds = cds, #lineage_states = trajectory_states,
                                                lineage_labels = lineage_labels, stretch = stretch)
    overlap_rng <- c(0, max(pData(cds_subset)$Pseudotime))
  

  if (length(trajectory_states) != 2)
    stop("calILRs can only work for two Lineages")
  if(!all(pData(cds_subset)[pData(cds_subset)$State %in% trajectory_states, "Lineage"] %in% pData(cds_subset)[, "Lineage"]))
      stop("state(s) in trajectory_states are not included in 'Lineage'")
    
  if(verbose)
    message(paste("the pseudotime range for the calculation of ILRs:", overlap_rng[1], overlap_rng[2], sep = ' '))
  
  if(!is.null(lineage_labels)){
    trajectory_states <- lineage_labels
  }
  cds_branchA <- cds_subset[, pData(cds_subset)[, "Lineage"] ==
                              trajectory_states[1]]
  cds_branchB <- cds_subset[, pData(cds_subset)[, "Lineage"] ==
                              trajectory_states[2]]
  
  formula_all_variables <- all.vars(as.formula(trend_formula))
  
  if(!all(formula_all_variables %in% colnames(pData(cds_subset))))
    stop('All the variables in the model formula has to be included in the pData columns (excepting Lineage)')
  
  t_rng <- range(pData(cds_branchA)$Pseudotime)
  str_new_cds_branchA <- data.frame(Pseudotime = seq(overlap_rng[1], overlap_rng[2],
                                                     length.out = 100), Lineage = as.factor(trajectory_states[1]))
  if (verbose)
    message(paste("Check the whether or not Pseudotime scaled from 0 to 100: ",
                  sort(pData(cds_branchA)$Pseudotime)))
  colnames(str_new_cds_branchA)[2] <- formula_all_variables[2] #interaction term can be terms rather than Lineage
  
  str_new_cds_branchB <- data.frame(Pseudotime = seq(overlap_rng[1], overlap_rng[2],
                                                     length.out = 100), Lineage = as.factor(trajectory_states[2]))
  if (verbose)
    message(paste("Check the whether or not Pseudotime scaled from 0 to 100: ",
                  sort(pData(cds_branchB)$Pseudotime)))
  
  colnames(str_new_cds_branchB)[2] <- formula_all_variables[2] #interaction term can be terms rather than Lineage
  
  str_branchAB_expression_curve_matrix <- genSmoothCurves(cds_subset, cores=cores, trend_formula = trend_formula,
                                                          relative_expr = relative_expr, new_data = rbind(str_new_cds_branchA, str_new_cds_branchB))
  
  str_branchA_expression_curve_matrix <- str_branchAB_expression_curve_matrix[, 1:nrow(str_new_cds_branchA)]
  str_branchB_expression_curve_matrix <- str_branchAB_expression_curve_matrix[, 
                                                                              (nrow(str_new_cds_branchA) + 1):(nrow(str_new_cds_branchA) + nrow(str_new_cds_branchB))]
  
  if (useVST) {
    str_branchA_expression_curve_matrix <- vstExprs(cds,
                                                    expr_matrix = str_branchA_expression_curve_matrix,
                                                    round_vals = round_exprs)
    str_branchB_expression_curve_matrix <- vstExprs(cds,
                                                    expr_matrix = str_branchB_expression_curve_matrix,
                                                    round_vals = round_exprs)
    str_logfc_df <- str_branchA_expression_curve_matrix -
      str_branchB_expression_curve_matrix
  }
  else {
    str_logfc_df <- log2((str_branchA_expression_curve_matrix +
                            1)/(str_branchB_expression_curve_matrix + 1))
  }
  if (label_by_short_name) {
    row.names(str_logfc_df) <- fData(cds[, ])$gene_short_name
  }
  str_logfc_df[which(str_logfc_df <= -ILRs_limit, arr.ind = T)] <- -ILRs_limit
  str_logfc_df[which(str_logfc_df >= ILRs_limit, arr.ind = T)] <- ILRs_limit
  if (output_type == "after_bifurcation") {
    t_bifurcation_ori <- min(pData(cds[, c(which(pData(cds)$State ==
                                                   trajectory_states[1]), which(pData(cds)$State == trajectory_states[2]))])$Pseudotime)
    t_bifurcation <- pData(cds_subset[, pData(cds)$Pseudotime ==
                                        t_bifurcation_ori])$Pseudotime
    
    if (stretch)
      bif_index <- as.integer(pData(cds_subset[, pData(cds)$Pseudotime ==
                                                 t_bifurcation])$Pseudotime)
    else {
      bif_index <- as.integer(min(t_bifurcation/(max(pData(cds_branchA)$Pseudotime)/100),
                                  t_bifurcation/(max(pData(cds_branchB)$Pseudotime)/100)))
    }
    str_logfc_df[, bif_index:100] <- str_logfc_df
  }
  if (!is.null(file))
    save(str_logfc_df, str_branchA_expression_curve_matrix, str_branchB_expression_curve_matrix, file = file)
  
  if(return_all) {
    rMax <- function(df) {apply(df, 1, function(x) if(all(is.na(x))) NA else max(abs(x), na.rm = T))} #calculate row max
    
    str_raw_div_df <- str_branchA_expression_curve_matrix - str_branchB_expression_curve_matrix
    str_norm_div_df <- str_raw_div_df / rMax(str_raw_div_df) #calculate normalized divergence
    
    log_str_raw_div_df <- log2((str_branchA_expression_curve_matrix + .1)/(str_branchB_expression_curve_matrix + .1))
    norm_str_logfc_df <- str_logfc_df / rMax(log_str_raw_div_df) #calculate normalized divergence
    
    return(list(str_logfc_df = str_logfc_df, norm_str_logfc_df = norm_str_logfc_df,
                str_norm_div_df = str_norm_div_df, str_raw_div_df = str_raw_div_df, str_branchA_expression_curve_matrix = str_branchA_expression_curve_matrix, 
                str_branchB_expression_curve_matrix = str_branchB_expression_curve_matrix))
  }
  else
    return(str_logfc_df)
}

#' Detect the maturation time point where the gene expression starts to diverge 
#' 
#' This function is used to determine the bifurcation point for the gene expression between two distinct biological processes.
#' For processes we can not distinguish between lineages (or phenotype groups, like knockout VS un-knockout), this function will 
#' only detect bifurcation points after the inferenced bifurcatioin from the PQ-tree. The 
#'
#' @param str_log_df the ILRs dataframe calculated from calILRs function. If this data.frame is provided, all the following parameters are ignored. Note that we need to only use the ILRs after the bifurcation point if we duplicated the progenitor cell state.
#' @param ILRs_threshold the ILR value used to determine the earliest divergence time point
#' @param detect_all a logic flag to determine whether or not genes without ILRs pass the threshold will still report a bifurcation point
#' @param cds CellDataSet for the experiment
#' @param Lineage The column in pData used for calculating the ILRs (If not equal to "Lineage", a warning will report)
#' @param branch_point the id for the branch point choosed to select cell branches and progenitors 
#' @param lineage_states The states for two branching lineages
#' @param stretch a logic flag to determine whether or not each lineage should be stretched
#' @param cores Number of cores when fitting the spline curves
#' @param trend_formula the model formula to be used for fitting the expression trend over pseudotime
#' @param ILRs_limit the minimum Instant Log Ratio used to make the heatmap plot
#' @param relative_expr A logic flag to determine whether or not the relative expressed should be used when we fitting the spline curves 
#' @param label_by_short_name label the rows of the returned matrix by gene_short_name (TRUE) or feature id (FALSE)
#' @param useVST A logic flag to determine whether or not the Variance Stablization Transformation should be used to stablize the gene expression.
#' When VST is used, the difference between two lineages are used instead of the log-ratio.
#' @param round_exprs A logic flag to determine whether or not the expression value should be rounded into integer
#' @param output_type A character either of "all" or "after_bifurcation". If "after_bifurcation" is used, only the time points after the bifurcation point will be selected. Note that, if Lineage is set to "Lineage", we will only use "after_bifurcation" since we duplicated the progenitor cells and the bifurcation should only happen after the largest mature level from the progenitor cells
#' @param return_cross_point A logic flag to determine whether or not only return the cross point 
#' @param file the name for storing the data. Since the calculation of the Instant Log Ratio is very time consuming, so by default the result will be stored
#' @param verbose Whether to report verbose output
#' @param ... Additional arguments passed to calILRs
#' @return a vector containing the time for the bifurcation point with gene names for each value
#' @importFrom reshape2 melt
#' 
detectBifurcationPoint <- function(str_log_df = NULL, 
                                   ILRs_threshold = 0.1, 
                                   detect_all = T,
                                   cds = cds,
                                   Lineage = 'Lineage',
                                   branch_point=NULL,
                                   lineage_states = c(2, 3),
                                   stretch = T,
                                   cores = detectCores(),
                                   trend_formula = "~sm.ns(Pseudotime, df = 3)",
                                   ILRs_limit = 3,
                                   relative_expr = TRUE,
                                   label_by_short_name = TRUE,
                                   useVST = FALSE,
                                   round_exprs = FALSE,
                                   output_type = 'all', #'after_bifurcation
                                   return_cross_point = T, 
                                   file = "bifurcation_heatmap", verbose = FALSE, ...) {
  if(is.null(str_log_df)) {
    if(Lineage == 'Lineage') output_type = 'after_bifurcation'
    
    str_log_df <- calILRs(cds = cds,
                          Lineage,
                          lineage_states,
                          branch_point=branch_point,
                          stretch,
                          cores,
                          trend_formula,
                          ILRs_limit,
                          relative_expr,
                          label_by_short_name,
                          useVST,
                          round_exprs,
                          output_type = output_type,
                          file, verbose, ...)
  }
  # rMax <- function(df) {apply(df, 1, function(x) if(all(is.na(x))) NA else max(abs(x), na.rm = T))}
  
  else {
    bifurcation_time <- apply(str_log_df, 1, function(x) {
      # deriv <- diff(x) the ILRs are smooth, so use min is fine
      index <- Inf
      
      #new algorithm to bifurcation time point:
      if(all(is.na(x))) {
        return(NA)
      }
      
      max_ind <- which(abs(x) == max(abs(x)))
      
      # return(max(max_ind))
      if(length(max_ind) > 1) {
        max_ind <- min(max_ind)
        warning('multiple maximal time points detected ', max_ind)
      }
      
      #detect the cross point
      inflection_point_tmp <- which(x[1:(length(x) - 1)] * x[2:length(x)] <= 0)
      
      if(all(max_ind <= inflection_point_tmp)) return(NA) #remove all zero values and genes monotonically goes down
      
      inflection_point <- max(inflection_point_tmp[inflection_point_tmp < max_ind])
      
      if(return_cross_point == T) {
        return(inflection_point * sign(sum(x)))
      }
      
      else if (return_cross_point == F & !is.null(ILRs_threshold) ) { 
        rev_x <- rev(x[(inflection_point):max_ind])
        if(any(which(abs(rev_x) >= ILRs_threshold))){
          index_tmp <- max(which(abs(rev_x) > ILRs_threshold))
          index <- (max_ind - index_tmp + 1 ) * sign(sum(rev_x))
        }
        else if(detect_all & all(!is.na(rev_x))) {
          index_tmp <-  max(which(abs(rev_x) == max(abs(rev_x)))) #the earliest time point when the bifurcation is largest 
          index <- (max_ind - index_tmp + 1 ) * sign(sum(rev_x)) 
        }
        }
      index
    })
  }
  
  # print(bifurcation_time)
  # str_norm_div_df
  
  names(bifurcation_time) <- row.names(str_log_df)
  
  return(bifurcation_time)
}

#' Peform the beam analysis test
#'
#' Perform BEAM analysis
#'
#' @param cds a CellDataSet object upon which to perform this operation
#' @param fullModelFormulaStr a formula string specifying the full model in differential expression tests (i.e. likelihood ratio tests) for each gene/feature.
#' @param reducedModelFormulaStr a formula string specifying the reduced model in differential expression tests (i.e. likelihood ratio tests) for each gene/feature.
#' @param lineage_states ids for the immediate branch lineage which obtained from lineage construction based on MST
#' @param branch_point the id for the branch point choosed to select cell branches and progenitors 
#' @param relative_expr a logic flag to determine whether or not the relative gene expression should be used
#' @param stretch  a logic flag to determine whether or not each lineage should be stretched
#' @param q_thrsld The threshold for the qval 
#' @param lineage_labels the name for each lineage, for example, "AT1" or "AT2"  
#' @param verbose Whether to generate verbose output
#' @param cores the number of cores to be used while testing each gene for differential expression
#' @param ... additional arguments to be passed to differentialGeneTest
#' @return a data frame containing the p values and q-values from the likelihood ratio tests on the parallel arrays of models, as well as the time point where the gene starts to bifurcate between two branches
#' @export
#'
BEAM <- function(cds, fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3)*Lineage", 
					reducedModelFormulaStr = "~sm.ns(Pseudotime, df = 3)", 
					lineage_states = NULL,
					branch_point=1,
					relative_expr = TRUE, 
					stretch = TRUE, 
					q_thrsld = 0.05, 
					lineage_labels = NULL, 
					verbose = FALSE,
					cores = 1, 
					...) {

	branchTest_res <- branchTest(cds, fullModelFormulaStr = fullModelFormulaStr,
	                       reducedModelFormulaStr = reducedModelFormulaStr, 
	                       lineage_states = lineage_states, 
	                       branch_point=branch_point,
	                       relative_expr = relative_expr,
	                       stretch = stretch,
	                       cores = cores, 
	                       lineage_labels = lineage_labels, ...)
	cmbn_df <- branchTest_res[, 1:4] 

  #make a newCellDataSet object with the smoothed data? 
	if(verbose)
   message('pass branchTest')

	# if(draw_branched_kinetics) {
	# 	plot_genes_branched_pseudotime(cds, 
	# 		lineage_states = lineage_states, 
	# 		lineage_labels = lineage_labels,
	# 		stretch = TRUE, 
	# 		min_expr = NULL, 
	# 		cell_size = 0.75,
	# 		nrow = NULL, 
	# 		ncol = 1, 
	# 		panel_order = NULL, 
	# 		color_by = "State",
	# 		cell_color_by = "State",
	# 		trajectory_color_by = "State", 
	# 		trend_formula = fullModelFormulaStr, 
	# 		reducedModelFormulaStr = reducedModelFormulaStr, 
	# 		label_by_short_name = TRUE,
	# 		weighted = TRUE, 
	# 		add_ABC = FALSE, 
	# 		add_pval = FALSE,
	# 		normalize = TRUE,
	# 		bifurcation_time = NULL, 
	# 		#gene_pairs = NULL,
	# 	...)
	# }

	# if(draw_branched_heatmap) {
	# 	plot_genes_branched_heatmap(cds_subset, 
	# 	  num_clusters = 6,
	# 	  ABC_df = NULL, 
	# 	  branchTest_df = NULL, 
	# 	  lineage_labels = lineage_labels, 
	# 	  stretch = T, 
	# 	  scaling = T,
	# 	  norm_method = c("vstExprs", "log"), 
	# 	  use_fitting_curves = T, 
	# 	  dist_method = NULL, 
	# 	  hclust_method = "ward", 
	# 	  heatmap_height = 3, 
	# 	  heatmap_width = 4,
	# 	  ABC_lowest_thrsd = 0, 
	# 	  ABC_highest_thrsd = 2,
	# 	  qval_lowest_thrsd = 1, 
	# 	  qval_highest_thrsd = 5,
	# 	  hmcols = NULL, 
	# 	  Cell_type_color = c('#979797', '#F05662', '#7990C8'), 
	# 	  trend_formula = '~sm.ns(Pseudotime, df=3) * Lineage',
	# 	  pseudo_cnt = 0, 
	# 	  add_annotation_row = NULL,
	# 	  add_annotation_col = NULL,
	# 	  show_rownames = F, 
	# 	  cores = cores,
	# 	  use_gene_short_name = F,
	# 	  file_name = 'branched_heatmap.pdf')
	# }

# 
#   ILRs_res <- calILRs(cds = cds, 
#                       branch_point=branch_point,
#   			              trajectory_states = lineage_states, 
#   			              lineage_labels = lineage_labels, 
#   			              stretch = stretch, 
#   			              cores = cores, 
#   			              trend_formula = fullModelFormulaStr,
#   			              ILRs_limit = 3, 
#   			              relative_expr = relative_expr, 
#   			              return_all = T,
#   			              ...)
# 
#   if(verbose)
#    message('pass calILRs')
#   
#   BifurcationTimePoint_res <- detectBifurcationPoint(str_log_df = ILRs_res$str_norm_div_df,
#                                                      branch_point=branch_point,
#                                                      lineage_states = lineage_states, 
#                                                      stretch = stretch, 
#                                                      cores = cores, 
#                                                      trend_formula = fullModelFormulaStr, 
#                                                      relative_expr = relative_expr, 
#   	                                                 ...)
#       
  #if(verbose)
  # message('pass detectBifurcationPoint')
  # print('pass detectBifurcationPoint')
  
  #cmbn_df <- cbind(cmbn_df, data.frame(Bifurcation_time_point = BifurcationTimePoint_res))

	fd <- fData(cds)[row.names(cmbn_df),]

	#combined dataframe: 
	cmbn_df <- cbind(cmbn_df, fd)

  if(verbose)
   message('return results')

	return(cmbn_df)
}

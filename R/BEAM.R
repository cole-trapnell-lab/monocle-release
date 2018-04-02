#' Build a CellDataSet that splits cells among two branches
#' 
#' Analyzing branches with \code{\link{BEAM}()} requires fitting two models to
#' the expression data for each gene. The full model assigns each cell to one of
#' the two outcomes of the branch, and the reduced model excludes this
#' assignment. \code{buildBranchBranchCellDataSet()} takes a CellDataSet object
#' and returns a version where the cells are assigned to one of two branches.
#' The branch for each cell is encoded in a new column, "Branch", in the pData
#' table in the returned CellDataSet.
#' 
#' @param cds CellDataSet for the experiment
#' @param progenitor_method The method to use for dealing with the cells prior to the branch
#' @param branch_point The ID of the branch point to analyze. Can only be used
#'   when \code{\link{reduceDimension}()} is called with \code{reduction_method
#'   = "DDRTree"}.
#' @param branch_states The states for two branching branches
#' @param branch_labels The names for each branching branch
#' @import methods
#' @importFrom Biobase pData<- exprs
#' @importFrom stats setNames
#' @importFrom igraph V degree shortest_paths bfs
#' @param stretch A logical flag to determine whether or not the pseudotime trajectory for each branch should be stretched to the same range or not
#' @return a CellDataSet with the duplicated cells and stretched branches
#' @export
buildBranchCellDataSet <- function(cds,
                                   progenitor_method = c('sequential_split', 'duplicate'), 
                                   branch_states = NULL, 
                                   branch_point = 1,
                                   branch_labels = NULL, 
                                   stretch = TRUE)
{
  # TODO: check that branches are on the same paths
  if(is.null(pData(cds)$State) | is.null(pData(cds)$Pseudotime)) 
    stop('Please first order the cells in pseudotime using orderCells()')
  if(is.null(branch_point) & is.null(branch_states)) 
    stop('Please either specify the branch_point or branch_states to select subset of cells')
  #if(ncol(cds@reducedDimS) != ncol(cds))
  #  stop('You probably used clusterCells function which should be used together with buildBranchCellDataSet, try re-run reduceDimension without clustering cells again')
  
  if (!is.null(branch_labels) & !is.null(branch_states)) {
    if(length(branch_labels) != length(branch_states))
      stop("length of branch_labels doesn't match with that of branch_states")
    branch_map <- setNames(branch_labels, as.character(branch_states))
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
  if (is.null(branch_states) == FALSE){
    
    # If the user didn't specify a branch point,
    # let's walk back from the branch states
    for (leaf_state in branch_states){
      
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
        ancestor_cells_for_branch <- row.names(closest_vertex)[which(V(pr_graph_cell_proj_mst)[closest_vertex]$name %in% path_to_ancestor)]
      }else if (cds@dim_reduce_type == "ICA"){
        ancestor_cells_for_branch <- path_to_ancestor
      }
      ancestor_cells_for_branch <- intersect(ancestor_cells_for_branch, colnames(cds))
      paths_to_root[[as.character(leaf_state)]] <- ancestor_cells_for_branch
    }
  }else{
    if(cds@dim_reduce_type == "DDRTree")
      pr_graph_cell_proj_mst <- minSpanningTree(cds)
    else
      pr_graph_cell_proj_mst <- cds@auxOrderingData$ICA$cell_ordering_tree
    
    mst_branch_nodes <- cds@auxOrderingData[[cds@dim_reduce_type]]$branch_points
    branch_cell <- mst_branch_nodes[branch_point]
    mst_no_branch_point <- pr_graph_cell_proj_mst - V(pr_graph_cell_proj_mst)[branch_cell]
    
    path_to_ancestor <- shortest_paths(pr_graph_cell_proj_mst, branch_cell, root_cell)
    path_to_ancestor <- names(unlist(path_to_ancestor$vpath))
    
    #post_branch_cells <- c()
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
        
        closest_vertex <- cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex
        #branch_state <- unique(pData(cds)[backbone_nei, ]$State)[1]
        
        path_to_root <- intersect(path_to_root, colnames(cds))
        paths_to_root[[backbone_nei]] <- path_to_root
        #post_branch_cells <- c(post_branch_cells, backbone_nei)
      }
    }
  }
  all_cells_in_subset <- c()
  
  if (is.null(branch_labels) == FALSE){
    if (length(branch_labels) != 2)
      stop("Error: branch_labels must have exactly two entries")
    names(paths_to_root) <- branch_labels
  }
  
  for (path_to_ancestor in paths_to_root){
    if (length(path_to_ancestor) == 0){
      stop("Error: common ancestors between selected State values on path to root State")
    }
    all_cells_in_subset <- c(all_cells_in_subset, path_to_ancestor)
  }
  all_cells_in_subset <- unique(all_cells_in_subset)
  
  common_ancestor_cells <- intersect(paths_to_root[[1]], paths_to_root[[2]])
  # if (length(paths_to_root) > 2){
  #   for (i in seq(3,length(paths_to_root))){
  #     common_ancestor_cells <- intersect(intersect(paths_to_root[i], paths_to_root[i-1]), common_ancestor_cells)
  #   }
  # }
  
  #when n-center used, this creates problems
  cds <- cds[, row.names(pData(cds[,all_cells_in_subset]))] #or just union(ancestor_cells, branch_cells)
  
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
        branch_cells <- setdiff(path_to_ancestor, common_ancestor_cells)
        pData[branch_cells,]$Pseudotime <- ((pData[branch_cells,]$Pseudotime - branch_pseudotime) * path_scaling_factor + branch_pseudotime)
      }
    }
    #pData[common_ancestor_cells,]$Pseudotime <- pData[common_ancestor_cells,]$Pseudotime / max_pseudotime
    
    pData$Pseudotime <- 100 * pData$Pseudotime / max_pseudotime
  }
  pData$original_cell_id <- row.names(pData)
  
  pData$original_cell_id <- row.names(pData)
  
  if(length(paths_to_root) != 2)
    stop('more than 2 branch states are used!')
  
  pData[common_ancestor_cells, "Branch"] <- names(paths_to_root)[1] #set progenitors to the branch 1
  
  progenitor_pseudotime_order <- order(pData[common_ancestor_cells, 'Pseudotime'])
  
  if (progenitor_method == 'duplicate') {
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
      weight_vec <- c(weight_vec, rep(1, length(common_ancestor_cells)))
      weight_vec_block <- rep(1, length(common_ancestor_cells))
      
      #pData <- rbind(pData, pData[common_ancestor_cells, ])
      new_pData_block <- ancestor_pData_block
      # new_pData_block$Lineage <- lineage_states[i]
      # new_pData_block$State <- lineage_states[i]
      
      row.names(new_pData_block) <- paste('duplicate', i, 1:length(common_ancestor_cells), sep = '_')
      
      pData_lineage_cells <- pData[setdiff(paths_to_root[[i]], common_ancestor_cells),]
      # pData_lineage_cells$Lineage <- lineage_states[i]
      # pData_lineage_cells$State <- lineage_states[i]
      
      weight_vec_block <- c(weight_vec_block, rep(1, nrow(pData_lineage_cells)))
      
      weight_vec <- c(weight_vec, weight_vec_block)
      
      new_pData_block <- rbind(new_pData_block, pData_lineage_cells)
      new_pData_block$Branch <- names(paths_to_root)[i]
      pData_blocks[[i]] <- new_pData_block
    }
    pData <- do.call(rbind, pData_blocks)
    exprs_data <- do.call(cbind, expr_blocks)
  }
  else if(progenitor_method == 'sequential_split') {
    pData$Branch <- names(paths_to_root)[1]
    
    branchA <- progenitor_pseudotime_order[seq(1, length(common_ancestor_cells), by = 2)]
    pData[common_ancestor_cells[branchA], 'Branch'] <- names(paths_to_root)[1]
    branchB <- progenitor_pseudotime_order[seq(2, length(common_ancestor_cells), by = 2)]
    pData[common_ancestor_cells[branchB], 'Branch'] <- names(paths_to_root)[2]   
    
    # Duplicate the root cell to make sure both regression start at pseudotime zero:
    zero_pseudotime_root_cell <- common_ancestor_cells[progenitor_pseudotime_order[1]]
    exprs_data <- cbind(exprs(cds), 'duplicate_root' = exprs(cds)[, zero_pseudotime_root_cell])
    pData <- rbind(pData, pData[zero_pseudotime_root_cell, ])
    row.names(pData)[nrow(pData)] <- 'duplicate_root'
    pData[nrow(pData), 'Branch'] <- names(paths_to_root)[2]
    
    weight_vec <- rep(1, nrow(pData))
    
    for (i in 1:length(paths_to_root)){
      path_to_ancestor <- paths_to_root[[i]]
      branch_cells <- setdiff(path_to_ancestor, common_ancestor_cells)
      pData[branch_cells,]$Branch <- names(paths_to_root)[i]
    }
  }
  
  pData$Branch <- as.factor(pData$Branch)
  
  pData$State <- factor(pData$State)
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

#' Test for branch-dependent expression
#'
#' Testing for branch-dependent expression with \code{\link{BEAM}()} first
#' involves constructing a CellDataSet that assigns each cell to a branch, and 
#' then performing a likelihood ratio test to see if the branch assignments 
#' significantly improves the fit over a null model that does not split the cells.
#' \code{branchTest()} implements these two steps.
#'
#' @param cds a CellDataSet object upon which to perform this operation
#' @param fullModelFormulaStr a formula string specifying the full model in differential expression tests (i.e. likelihood ratio tests) for each gene/feature.
#' @param reducedModelFormulaStr a formula string specifying the reduced model in differential expression tests (i.e. likelihood ratio tests) for each gene/feature.
#' @param branch_states  states corresponding to two branches 
#' @param branch_point The ID of the branch point to analyze. Can only be used when reduceDimension is called with method = "DDRTree".
#' @param relative_expr a logic flag to determine whether or not the relative gene expression should be used
#' @param cores the number of cores to be used while testing each gene for differential expression
#' @param branch_labels the name for each branch, for example, AT1 or AT2  
#' @param verbose Whether to show VGAM errors and warnings. Only valid for cores = 1. 
#' @param ... Additional arguments passed to differentialGeneTest
#' @import methods
#' @importFrom stats as.formula
#' @return a data frame containing the p values and q-values from the likelihood ratio tests on the parallel arrays of models.
#' @export
branchTest <- function(cds, fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3)*Branch",
                       reducedModelFormulaStr = "~sm.ns(Pseudotime, df = 3)", 
                       branch_states = NULL, 
                       branch_point=1,
                       relative_expr = TRUE,
                       cores = 1, 
                       branch_labels = NULL, 
                       verbose = FALSE,
                       ...) {
  
  if("Branch" %in% all.vars(terms(as.formula(fullModelFormulaStr)))) {
    cds_subset <- buildBranchCellDataSet(cds = cds, 
                                                branch_states = branch_states,
                                                branch_point=branch_point,
                                                branch_labels = branch_labels,
                                                ...)
  }
  else
    cds_subset <- cds
  
  branchTest_res <- differentialGeneTest(cds_subset, 
                                         fullModelFormulaStr = fullModelFormulaStr, 
                                         reducedModelFormulaStr = reducedModelFormulaStr, 
                                         cores = cores, 
                                         relative_expr = relative_expr, 
                                         verbose=verbose)
  
  return(branchTest_res)
}


#' Compute the area between curves (ABC) for branch-dependent genes
#' 
#' This function is used to calculate the ABC score based on the the nature spline curves fitted for each branch. ABC score is used to 
#' quantify the total magnitude of divergence between two branchs. By default, the ABC score is the area between two fitted spline curves. 
#' The ABC score can be used to rank gene divergence. When coupled with p-val calculated from the branchTest, it can be used to identify
#' potential major regulators for branch bifurcation. 
#'
#' @param cds a CellDataSet object upon which to perform this operation
#' @param trend_formula a formula string specifying the full model in differential expression tests (i.e. likelihood ratio tests) for each gene/feature.
#' @param branch_point the point where two branches diverge
#' @param trajectory_states States corresponding to two branches 
#' @param relative_expr a logic flag to determine whether or not the relative gene expression should be used
#' @param stretch a logic flag to determine whether or not each branch should be stretched
#' @param cores the number of cores to be used while testing each gene for differential expression
#' @param verbose a logic flag to determine whether or not we should output detailed running information 
#' @param min_expr the lower limit for the expressed gene
#' @param integer_expression the logic flag to determine whether or not the integer numbers are used for calculating the ABCs. Default is False. 
#' @param num number of points on the fitted branch trajectories used for calculating the ABCs. Default is 5000. 
#' @param branch_labels the name for each branch, for example, AT1 or AT2  
#' @param ... Additional arguments passed to buildBranchCellDataSet
#' @import methods
#' @importFrom Biobase pData fData
#' @return a data frame containing the ABCs (Area under curves) score as the first column and other meta information from fData
#' @export 
calABCs <- function(cds,
                    trend_formula = "~sm.ns(Pseudotime, df = 3)*Branch",
                    branch_point = 1,
                    trajectory_states = NULL,
                    relative_expr = TRUE, 
                    stretch = TRUE, 
                    cores = 1, 
                    verbose = F,
                    min_expr = 0.5, 
                    integer_expression = FALSE, 
                    num = 5000, 
                    branch_labels = NULL,
                    ...){
  ABC_method = "integral"
  if(!is.null(trajectory_states)){
    if (length(trajectory_states) != 2)
      stop("Sorry, this function only supports the calculation of ABCs between TWO branch trajectories")
  }
  
    cds_subset <- buildBranchCellDataSet(cds = cds, 
                                                progenitor_method = 'duplicate',
                                                branch_point = branch_point, 
                                                branch_states = trajectory_states,
                                                branch_labels = branch_labels, stretch = stretch, ...)
    overlap_rng <- c(0, max(pData(cds_subset)$Pseudotime))
 
 
  # if (length(trajectory_states) != 2)
  #   stop("calILRs can only work for two branches")
  # if(!all(trajectory_states %in% pData(cds_subset)[, "Branch"]))
  #   stop("state(s) in trajectory_states are not included in 'Branch'")
  # 
  trajectory_states <- unique(pData(cds_subset)[, "Branch"])

  if(verbose)
    message(paste("the pseudotime range for the calculation of ILRs:", overlap_rng[1], overlap_rng[2], sep = ' '))
  
  cds_branchA <- cds_subset[, pData(cds_subset)[, "Branch"] ==
                              trajectory_states[1]]
  cds_branchB <- cds_subset[, pData(cds_subset)[, "Branch"] ==
                              trajectory_states[2]]
  
  #use the genSmoothCurves to generate smooth curves for calculating the ABC scores:
  formula_all_variables <- all.vars(as.formula(trend_formula))
  
  t_rng <- range(pData(cds_branchA)$Pseudotime)
  str_new_cds_branchA <- data.frame(Pseudotime = seq(overlap_rng[1], overlap_rng[2],
                                                     length.out = num), Branch = as.factor(as.character(trajectory_states[1])))
  colnames(str_new_cds_branchA)[2] <- formula_all_variables[2] #interaction term can be terms rather than Branch
  if (verbose)
    print(paste("Check the whether or not Pseudotime scaled from 0 to 100: ",
                sort(pData(cds_branchA)$Pseudotime)))
  str_new_cds_branchB <- data.frame(Pseudotime = seq(overlap_rng[1], overlap_rng[2],
                                                     length.out = num), Branch = as.factor(as.character(trajectory_states[2])))
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

#' Calculate the Instantaneous Log Ratio between two branches
#' 
#' This function is used to calculate the Instant Log Ratio between two branches which can be used to prepare the heatmap demonstrating the branch gene expression divergence hirearchy. If "stretch" is specifified, each  
#' branch will be firstly stretched into maturation level from 0-100. Since the results when we use "stretching" are always better and 
#' IRLs for non-stretched spline curves are often mismatched, we may only turn down "non-stretch" functionality in future versions. Then, we fit two separate nature spline curves for each 
#' individual linages. The log-ratios of the value on each spline curve corresponding to each branch are calculated, which can be  
#' used as a measure for the magnitude of divergence between two branching branchs. 
#'
#' @param cds CellDataSet for the experiment
#' @param trend_formula trend_formula a formula string specifying the full model in differential expression tests (i.e. likelihood ratio tests) for each gene/feature.
#' @param branch_point the point where two branches diverge
#' @param trajectory_states states corresponding to two branches 
#' @param relative_expr A logic flag to determine whether or not the relative expressed should be used when we fitting the spline curves 
#' @param stretch a logic flag to determine whether or not each branch should be stretched
#' @param cores Number of cores when fitting the spline curves
#' @param ILRs_limit the minimum Instant Log Ratio used to make the heatmap plot
#' @param label_by_short_name label the rows of the returned matrix by gene_short_name (TRUE) or feature id (FALSE)
#' @param useVST A logic flag to determine whether or not the Variance Stablization Transformation should be used to stablize the gene expression.
#' When VST is used, the difference between two branchs are used instead of the log-ratio.
#' @param round_exprs A logic flag to determine whether or not the expression value should be rounded into integer
#' @param output_type A character either of "all" or "after_bifurcation". If "after_bifurcation" is used, only the time points after the bifurcation point will be selected
#' @param branch_labels the name for each branch, for example, AT1 or AT2  
#' @param file the name for storing the data. Since the calculation of the Instant Log Ratio is very time consuming, so by default the result will be stored
#' @param return_all A logic flag to determine whether or not all the results from the analysis should be returned, this includes 
#' a dataframe for the log fold change, normalized log fold change, raw divergence, normalized divergence, fitting curves for each branch 
#' @param verbose Whether or not detailed running information should be returned 
#' @param ... Additional arguments passed to buildBranchCellDataSet
#' @return a ggplot2 plot object
#' @import ggplot2
#' @import methods
#' @importFrom Biobase pData fData
#' @importFrom reshape2 melt
#' @export 
calILRs <- function (cds, 
          trend_formula = "~sm.ns(Pseudotime, df = 3)*Branch",
          branch_point = 1,
          trajectory_states = NULL, 
          relative_expr = TRUE, 
          stretch = TRUE, 
          cores = 1, 
          ILRs_limit = 3, 
          label_by_short_name = TRUE,
          useVST = FALSE, 
          round_exprs = FALSE, 
          output_type = "all", 
          branch_labels = NULL, 
          file = NULL, 
          return_all = F, 
          verbose = FALSE, 
          ...){
    if(!is.null(trajectory_states)){
      if (length(trajectory_states) != 2)
        stop("Sorry, this function only supports the calculation of ILRs between TWO branch trajectories")
    }

    cds_subset <- buildBranchCellDataSet(cds = cds, branch_states = trajectory_states,
                                                branch_point = branch_point,
                                                progenitor_method = 'duplicate',
                                                branch_labels = branch_labels, stretch = stretch, ...)
    overlap_rng <- c(0, max(pData(cds_subset)$Pseudotime))
  

  # if (length(trajectory_states) != 2)
  #   stop("calILRs can only work for two Branches")
  # if(!all(pData(cds_subset)[pData(cds_subset)$State %in% trajectory_states, "Branch"] %in% pData(cds_subset)[, "Branch"]))
  #     stop("state(s) in trajectory_states are not included in 'Branch'")
    
  # if(verbose)
  #   message(paste("the pseudotime range for the calculation of ILRs:", overlap_rng[1], overlap_rng[2], sep = ' '))
  if(is.null(trajectory_states))
    trajectory_states <- unique(pData(cds_subset)[, "Branch"])

  if(!is.null(branch_labels)){
    trajectory_states <- branch_labels
  }
  else if(is.null(trajectory_states)){
    trajectory_states_tmp <- as.character(trajectory_states)
    branch_stats <- table(pData(cds_subset)[row.names(subset(pData(cds), as.character(State) == trajectory_states_tmp[1])), 'Branch'])
    trajectory_states[1] <- names(which.max(branch_stats))
    branch_stats <- table(pData(cds_subset)[row.names(subset(pData(cds), as.character(State) == trajectory_states_tmp[2])), 'Branch'])
    trajectory_states[2] <- names(which.max(branch_stats))
  }
      
  cds_branchA <- cds_subset[, pData(cds_subset)[, "Branch"] ==
                              trajectory_states[1]]
  cds_branchB <- cds_subset[, pData(cds_subset)[, "Branch"] ==
                              trajectory_states[2]]
  
  formula_all_variables <- all.vars(as.formula(trend_formula))
  
  if(!all(formula_all_variables %in% colnames(pData(cds_subset))))
    stop('All the variables in the model formula has to be included in the pData columns (excepting Branch)')
  
  t_rng <- range(pData(cds_branchA)$Pseudotime)
  str_new_cds_branchA <- data.frame(Pseudotime = seq(overlap_rng[1], overlap_rng[2],
                                                     length.out = 100), Branch = as.factor(as.character(trajectory_states[1])))
  if (verbose)
    message(paste("Check the whether or not Pseudotime scaled from 0 to 100: ",
                  sort(pData(cds_branchA)$Pseudotime)))
  colnames(str_new_cds_branchA)[2] <- formula_all_variables[2] #interaction term can be terms rather than Branch
  
  str_new_cds_branchB <- data.frame(Pseudotime = seq(overlap_rng[1], overlap_rng[2],
                                                     length.out = 100), Branch = as.factor(as.character(trajectory_states[2])))
  if (verbose)
    message(paste("Check the whether or not Pseudotime scaled from 0 to 100: ",
                  sort(pData(cds_branchB)$Pseudotime)))
  
  colnames(str_new_cds_branchB)[2] <- formula_all_variables[2] #interaction term can be terms rather than Branch
  
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

#' Calculate divergence times for branch-dependent genes
#' 
#' Branch-dependent genes may diverge at different points in pseudotime. \code{detectBifurcationPoint()}
#' calculates these times. Although the branch times will be shaped by and distributed
#' around the branch point in the trajectory, upstream regulators tend to branch
#' earlier in pseudotime than their targets. 
#'
#' @param str_log_df the ILRs dataframe calculated from calILRs function. If this data.frame is provided, all the following parameters are ignored. Note that we need to only use the ILRs after the bifurcation point if we duplicated the progenitor cell state.
#' @param ILRs_threshold the ILR value used to determine the earliest divergence time point
#' @param detect_all a logic flag to determine whether or not genes without ILRs pass the threshold will still report a bifurcation point
#' @param cds CellDataSet for the experiment
#' @param Branch The column in pData used for calculating the ILRs (If not equal to "Branch", a warning will report)
#' @param branch_point The ID of the branch point to analyze. Can only be used when reduceDimension is called with method = "DDRTree".
#' @param branch_states The states for two branching branchs
#' @param stretch a logic flag to determine whether or not each branch should be stretched
#' @param cores Number of cores when fitting the spline curves
#' @param trend_formula the model formula to be used for fitting the expression trend over pseudotime
#' @param ILRs_limit the minimum Instant Log Ratio used to make the heatmap plot
#' @param relative_expr A logic flag to determine whether or not the relative expressed should be used when we fitting the spline curves 
#' @param label_by_short_name label the rows of the returned matrix by gene_short_name (TRUE) or feature id (FALSE)
#' @param useVST A logic flag to determine whether or not the Variance Stablization Transformation should be used to stablize the gene expression.
#' When VST is used, the difference between two branchs are used instead of the log-ratio.
#' @param round_exprs A logic flag to determine whether or not the expression value should be rounded into integer
#' @param output_type A character either of "all" or "after_bifurcation". If "after_bifurcation" is used, only the time points after the bifurcation point will be selected. Note that, if Branch is set to "Branch", we will only use "after_bifurcation" since we duplicated the progenitor cells and the bifurcation should only happen after the largest mature level from the progenitor cells
#' @param return_cross_point A logic flag to determine whether or not only return the cross point 
#' @param file the name for storing the data. Since the calculation of the Instant Log Ratio is very time consuming, so by default the result will be stored
#' @param verbose Whether to report verbose output
#' @param ... Additional arguments passed to calILRs
#' @return a vector containing the time for the bifurcation point with gene names for each value
#' @import methods
#' @importFrom reshape2 melt
#' @importFrom parallel detectCores
#' @export 
detectBifurcationPoint <- function(str_log_df = NULL, 
                                   ILRs_threshold = 0.1, 
                                   detect_all = T,
                                   cds = cds,
                                   Branch = 'Branch',
                                   branch_point=NULL,
                                   branch_states = c(2, 3),
                                   stretch = T,
                                   cores = 1,
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
    if(Branch == 'Branch') output_type = 'after_bifurcation'
    
    str_log_df <- calILRs(cds = cds,
                          Branch,
                          branch_states,
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

#' Branched expression analysis modeling (BEAM).
#'
#' Identify genes with branch-dependent expression.
#' Branches in single-cell trajectories are generated by cell fate decisions 
#' in development and also arise when analyzing genetic, chemical, or environmental
#' perturbations. Branch expression analysis modeling is a statistical approach
#' for finding genes that are regulated in a manner that depends on the branch. 
#' Consider a progenitor cell that generates two distinct cell types. A single-cell
#' trajectory that includes progenitor cells and both differentiated cell types
#' will capture the "decision" as a branch point, with progenitors upstream of the branch
#' and the differentiated cells positioned along distinct branches. These branches
#' will be characterized by distinct gene expression programs. BEAM aims to find
#' all genes that differ between the branches. Such "branch-dependent" genes
#' can help identify the mechanism by which the fate decision is made.  
#' \code{BEAM()} Takes a CellDataSet and either a specified branch point, or a pair of 
#' trajectory outcomes (as States). If a branch point is provided, the function 
#' returns a dataframe of test results for dependence on that branch. If a pair
#' of outcomes is provided, it returns test results for the branch that unifies
#' those outcomes into a common path to the trajectory's root state.
#' \code{BEAM()} compares two models with a likelihood ratio test for branch-dependent
#' expression. The full model is the product of smooth Pseudotime and the Branch a cell is assigned to.
#' The reduced model just includes Pseudotime. You can modify these to include 
#' arbitrary additional effects in the full or both models.
#'
#' @param cds a CellDataSet object upon which to perform this operation
#' @param fullModelFormulaStr a formula string specifying the full model in differential expression tests (i.e. likelihood ratio tests) for each gene/feature.
#' @param reducedModelFormulaStr a formula string specifying the reduced model in differential expression tests (i.e. likelihood ratio tests) for each gene/feature.
#' @param branch_states ids for the immediate branch branch which obtained from branch construction based on MST
#' @param branch_point The ID of the branch point to analyze. Can only be used when reduceDimension is called with method = "DDRTree".
#' @param relative_expr a logic flag to determine whether or not the relative gene expression should be used
#' @param branch_labels the name for each branch, for example, "AT1" or "AT2"  
#' @param verbose Whether to generate verbose output
#' @param cores the number of cores to be used while testing each gene for differential expression
#' @param ... additional arguments to be passed to differentialGeneTest
#' @return a data frame containing the p values and q-values from the BEAM test, with one row per gene.
#' @import methods
#' @importFrom Biobase fData
#' @export
BEAM <- function(cds, fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3)*Branch", 
					reducedModelFormulaStr = "~sm.ns(Pseudotime, df = 3)", 
					branch_states = NULL,
					branch_point=1,
					relative_expr = TRUE, 
					branch_labels = NULL, 
					verbose = FALSE,
					cores = 1, 
					...) {

	branchTest_res <- branchTest(cds, fullModelFormulaStr = fullModelFormulaStr,
	                       reducedModelFormulaStr = reducedModelFormulaStr, 
	                       branch_states = branch_states, 
	                       branch_point=branch_point,
	                       relative_expr = relative_expr,
	                       cores = cores, 
	                       branch_labels = branch_labels, 
	                       verbose=verbose, 
	                       ...)
	cmbn_df <- branchTest_res[, 1:4] 

  #make a newCellDataSet object with the smoothed data? 
	if(verbose)
   message('pass branchTest')

	fd <- fData(cds)[row.names(cmbn_df),]

	#combined dataframe: 
	cmbn_df <- cbind(cmbn_df, fd)

  if(verbose)
   message('return results')

	return(cmbn_df)
}

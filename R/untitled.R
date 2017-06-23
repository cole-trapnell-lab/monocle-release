
plot_multiple_branches_pseudotime <- function(cds, 
                                           branches, 
                                           branches_name = NULL, 
                                           
                                           min_expr = NULL,                                     
                                           cell_size = 0.75,                                           
                                           norm_method = c("vstExprs", "log"), 
                                           nrow = NULL, 
                                           ncol = 1, 
                                           panel_order = NULL, 
                                           color_by = "Branch",

                                           trend_formula = '~sm.ns(Pseudotime, df=3)',
                                           label_by_short_name = TRUE,
                                           TPM = FALSE, 
                                           cores=1){


  if(TPM) {
    exprs(cds) <- esApply(cds, 2, function(x) x / sum(x) * 1e6)
  }

  if(!(all(branches %in% pData(cds)$State)) & length(branches) == 1){
    stop('This function only allows to make multiple branch plots where branches is included in the pData')
  }
  
  branch_label <- branches
  if(!is.null(branches_name)){
    if(length(branches) != length(branches_name)){
      stop('branches_name should have the same length as branches')
    }
    branch_label <- branches_name
  }
  
  #test whether or not the states passed to branches are true branches (not truncks) or there are terminal cells 
  g <- cds@minSpanningTree
  m <- NULL
  cds_exprs <- NULL 
  # branche_cell_num <- c()
  for(branch_in in branches) {
    branches_cells <- row.names(subset(pData(cds), State == branch_in))
    root_state <- subset(pData(cds), Pseudotime == 0)[, 'State']
    root_state_cells <- row.names(subset(pData(cds), State == root_state))
    
    if(cds@dim_reduce_type != 'ICA') {
      root_state_cells <- unique(paste('Y_', cds@auxOrderingData$DDRTree$pr_graph_cell_proj_closest_vertex[root_state_cells, ], sep = ''))
      branches_cells <- unique(paste('Y_', cds@auxOrderingData$DDRTree$pr_graph_cell_proj_closest_vertex[branches_cells, ], sep = ''))
    }
    root_cell <- root_state_cells[which(degree(g, v = root_state_cells) == 1)]
    tip_cell <- branches_cells[which(degree(g, v = branches_cells) == 1)]
    
    traverse_res <- traverseTree(g, root_cell, tip_cell)
    path_cells <- names(traverse_res$shortest_path[[1]])
    
    if(cds@dim_reduce_type != 'ICA') {
      pc_ind <- cds@auxOrderingData$DDRTree$pr_graph_cell_proj_closest_vertex 
      path_cells <- row.names(pc_ind)[paste('Y_', pc_ind[, 1], sep = '') %in% path_cells]
    }
    
    if(is.null(pData(cds)$no_expression)) {
      cds_subset <- cds[, path_cells]      
    } else {
      cds_subset <- cds[, path_cells %in% colnames(cds)[!pData(cds)$no_expression]]      
    }

    newdata <- data.frame(Pseudotime = pData(cds_subset)$Pseudotime, row.names = colnames(cds_subset))
    
    tmp <- t(esApply(cds_subset, 1, function(x) lowess(x[order(pData(cds_subset)$Pseudotime)])$y))
    # tmp <- t(esApply(cds_subset, 1, function(x) {
    #   x <- x[order(pData(cds_subset)$Pseudotime)]
    #   c(smooth::sma(x, order = 100, h = 1, silent="all")$fitted)}) #, x[length(x)]
    # )

    
    colnames(tmp) <- colnames(cds_subset)[order(pData(cds_subset)$Pseudotime)]
    # tmp <- genSmoothCurves(cds_subset, cores=cores, trend_formula = trend_formula,  
    #                        relative_expr = T, new_data = newdata)

    cds_exprs_tmp <- reshape2::melt(log2(tmp + 1))
    cds_exprs_tmp <- reshape2::melt(tmp)
    colnames(cds_exprs_tmp) <- c("f_id", "Cell", "expression")
    cds_exprs_tmp$Branch <- branch_label[which(branches == branch_in)]

    if(is.null(cds_exprs))
      cds_exprs <- cds_exprs_tmp
    else
      cds_exprs <- rbind(cds_exprs, cds_exprs_tmp)

    if(is.null(m))
      m <- tmp
    else
      m <- cbind(m, tmp)
  }
  
  #remove genes with no expression in any condition
  m=m[!apply(m,1,sum)==0,]
  
  norm_method <- match.arg(norm_method)
  
  # FIXME: this needs to check that vst values can even be computed. (They can only be if we're using NB as the expressionFamily)
  if(norm_method == 'vstExprs' && is.null(cds@dispFitInfo[["blind"]]$disp_func) == FALSE) {
    m = vstExprs(cds, expr_matrix=m)
  }     
  else if(norm_method == 'log') {
    m = log10(m+pseudocount)
  }
  
  if (is.null(min_expr)) {
    min_expr <- cds@lowerDetectionLimit
  }

  cds_pData <- pData(cds)

  cds_fData <- fData(cds)
  cds_exprs <- merge(cds_exprs, cds_fData, by.x = "f_id", by.y = "row.names")
  cds_exprs <- merge(cds_exprs, cds_pData, by.x = "Cell", by.y = "row.names")

  cds_exprs <- plyr::ddply(cds_exprs, .(Branch), mutate, Pseudotime = (Pseudotime - min(Pseudotime)) * 100 / (max(Pseudotime) - min(Pseudotime)) )

  # if (integer_expression) {
  #   cds_exprs$adjusted_expression <- round(cds_exprs$expression)
  # }
  # else {
  #   cds_exprs$adjusted_expression <- log10(cds_exprs$expression)
  # }
  if (label_by_short_name == TRUE) {
    if (is.null(cds_exprs$gene_short_name) == FALSE) {
      cds_exprs$feature_label <- as.character(cds_exprs$gene_short_name)
      cds_exprs$feature_label[is.na(cds_exprs$feature_label)] <- cds_exprs$f_id
    }
    else {
      cds_exprs$feature_label <- cds_exprs$f_id
    }
  }
  else {
    cds_exprs$feature_label <- cds_exprs$f_id
  }
  cds_exprs$feature_label <- as.factor(cds_exprs$feature_label)
  # trend_formula <- paste("adjusted_expression", trend_formula,
  #     sep = "")
  cds_exprs$Branch <- as.factor(cds_exprs$Branch) 
  
  # new_data <- data.frame(Pseudotime = pData(cds_subset)$Pseudotime, Branch = pData(cds_subset)$Branch)
  
  # full_model_expectation <- genSmoothCurves(cds_subset, cores=1, trend_formula = trend_formula, 
  #                                           relative_expr = T, new_data = new_data)
  # colnames(full_model_expectation) <- colnames(cds_subset)
  
  # cds_exprs$full_model_expectation <- apply(cds_exprs,1, function(x) full_model_expectation[x[2], x[1]])
  # if(!is.null(reducedModelFormulaStr)){
  #   reduced_model_expectation <- genSmoothCurves(cds_subset, cores=1, trend_formula = reducedModelFormulaStr,
  #                                                relative_expr = T, new_data = new_data)
  #   colnames(reduced_model_expectation) <- colnames(cds_subset)
  #   cds_exprs$reduced_model_expectation <- apply(cds_exprs,1, function(x) reduced_model_expectation[x[2], x[1]])
  # }
  
  # # FIXME: If you want to show the bifurcation time for each gene, this function
  # # should just compute it. Passing it in as a dataframe is just too complicated
  # # and will be hard on the user. 
  # # if(!is.null(bifurcation_time)){
  # #     cds_exprs$bifurcation_time <- bifurcation_time[as.vector(cds_exprs$gene_short_name)]
  # # }
  # if (method == "loess")
  #   cds_exprs$expression <- cds_exprs$expression + cds@lowerDetectionLimit
  # if (label_by_short_name == TRUE) {
  #   if (is.null(cds_exprs$gene_short_name) == FALSE) {
  #     cds_exprs$feature_label <- as.character(cds_exprs$gene_short_name)
  #     cds_exprs$feature_label[is.na(cds_exprs$feature_label)] <- cds_exprs$f_id
  #   }
  #   else {
  #     cds_exprs$feature_label <- cds_exprs$f_id
  #   }
  # }
  # else {
  #   cds_exprs$feature_label <- cds_exprs$f_id
  # }
  # cds_exprs$feature_label <- factor(cds_exprs$feature_label)
  # if (is.null(panel_order) == FALSE) {
  #   cds_exprs$feature_label <- factor(cds_exprs$feature_label,
  #                                     levels = panel_order)
  # }
  # cds_exprs$expression[is.na(cds_exprs$expression)] <- min_expr
  # cds_exprs$expression[cds_exprs$expression < min_expr] <- min_expr
  # cds_exprs$full_model_expectation[is.na(cds_exprs$full_model_expectation)] <- min_expr
  # cds_exprs$full_model_expectation[cds_exprs$full_model_expectation < min_expr] <- min_expr
  
  # if(!is.null(reducedModelFormulaStr)){
  #   cds_exprs$reduced_model_expectation[is.na(cds_exprs$reduced_model_expectation)] <- min_expr
  #   cds_exprs$reduced_model_expectation[cds_exprs$reduced_model_expectation < min_expr] <- min_expr
  # }
  
  cds_exprs$State <- as.factor(cds_exprs$State)
  cds_exprs$Branch <- as.factor(cds_exprs$Branch)
  
  q <- ggplot(aes(Pseudotime, expression), data = cds_exprs)
  # if (!is.null(bifurcation_time)) {
  #   q <- q + geom_vline(aes(xintercept = bifurcation_time),
  #                       color = "black", linetype = "longdash")
  # }
  if (is.null(color_by) == FALSE) {
    q <- q + geom_line(aes_string(color = color_by), size = I(cell_size))
  }
  #if (is.null(reducedModelFormulaStr) == FALSE)
    q <- q + facet_wrap(~feature_label, nrow = nrow, ncol = ncol, scales = "free_y") #+ scale_y_log10() 
  #else q <- q + scale_y_log10() + facet_wrap(~feature_label,
  #                                           nrow = nrow, ncol = ncol, scales = "free_y")
  #if (method == "loess")
  #  q <- q + stat_smooth(aes(fill = Branch, color = Branch),
  #                       method = "loess")
  #else if (method == "fitting") {
  #  q <- q + geom_line(aes_string(x = "Pseudotime", y = "full_model_expectation",
  #                                linetype = "Branch"), data = cds_exprs) #+ scale_color_manual(name = "Type", values = c(colour_cell, colour), labels = c("Pre-branch", "AT1", "AT2", "AT1", "AT2")
  #}
  
  #if(!is.null(reducedModelFormulaStr)) {
  #  q <- q + geom_line(aes_string(x = "Pseudotime", y = "reduced_model_expectation"),
  #                     color = 'black', linetype = 2, data =  cds_exprs)   
  #}
  
  q <- q + ylab("Expression") + xlab("Pseudotime (stretched)")
  
  q <- q + monocle_theme_opts()
  q + expand_limits(y = min_expr)
}
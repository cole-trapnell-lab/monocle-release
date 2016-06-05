utils::globalVariables(c("Pseudotime", "value", "ids", "prin_graph_dim_1", "prin_graph_dim_2", "State", 
                         "value", "feature_label", "expectation", "colInd", "rowInd", "value", 
                         "source_prin_graph_dim_1", "source_prin_graph_dim_2"))

monocle_theme_opts <- function()
{
    theme(strip.background = element_rect(colour = 'white', fill = 'white')) +
    theme(panel.border = element_blank(), axis.line = element_line()) +
    theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) +
    theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) + 
    theme(panel.background = element_rect(fill='white'))
}

#' Plots the minimum spanning tree on cells.
#'
#' @param cds CellDataSet for the experiment
#' @param x the column of reducedDimS(cds) to plot on the horizontal axis
#' @param y the column of reducedDimS(cds) to plot on the vertical axis
#' @param color_by the cell attribute (e.g. the column of pData(cds)) to map to each cell's color
#' @param show_tree whether to show the links between cells connected in the minimum spanning tree
#' @param show_backbone whether to show the diameter path of the MST used to order the cells
#' @param backbone_color the color used to render the backbone.
#' @param markers a gene name or gene id to use for setting the size of each cell in the plot
#' @param show_cell_names draw the name of each cell in the plot
#' @param cell_name_size the size of cell name labels
#' @return a ggplot2 plot object
#' @importFrom grid unit
#' @import ggplot2
#' @importFrom reshape2 melt
#' @export
#' @examples
#' \dontrun{
#' data(HSMM)
#' plot_spanning_tree(HSMM)
#' plot_spanning_tree(HSMM, color_by="Pseudotime", show_backbone=FALSE)
#' plot_spanning_tree(HSMM, markers="MYH3")
#' }
plot_spanning_tree <- function(cds, 
                               x=1, 
                               y=2, 
                               color_by="State", 
                               show_tree=TRUE, 
                               show_backbone=TRUE, 
                               backbone_color="black", 
                               markers=NULL, 
                               show_cell_names=FALSE, 
                               cell_size=1.5,
                               cell_link_size=0.75,
                               cell_name_size=2,
                               show_branch_points=TRUE){
  gene_short_name <- NULL
  sample_name <- NULL
  
  #TODO: need to validate cds as ready for this plot (need mst, pseudotime, etc)
  lib_info_with_pseudo <- pData(cds)
  
  if (is.null(cds@dim_reduce_type)){
    stop("Error: dimensionality not yet reduced. Please call reduceDimension() before calling this function.")
  }
  
  if (cds@dim_reduce_type == "ICA"){
    reduced_dim_coords <- reducedDimS(cds)
  }else if (cds@dim_reduce_type == "DDRTree"){
    reduced_dim_coords <- reducedDimK(cds)
  }else {
    stop("Error: unrecognized dimensionality reduction method.")
  }
  
  if (is.null(reduced_dim_coords)){
    stop("You must first call reduceDimension() before using this function")
  }
  
  ica_space_df <- data.frame(Matrix::t(reduced_dim_coords[c(x,y),]))
  colnames(ica_space_df) <- c("prin_graph_dim_1", "prin_graph_dim_2")
  
  ica_space_df$sample_name <- row.names(ica_space_df)
  #ica_space_with_state_df <- merge(ica_space_df, lib_info_with_pseudo, by.x="sample_name", by.y="row.names")
  #print(ica_space_with_state_df)
  dp_mst <- minSpanningTree(cds)
  
  if (is.null(dp_mst)){
    stop("You must first call orderCells() before using this function")
  }
  
  edge_list <- as.data.frame(get.edgelist(dp_mst))
  colnames(edge_list) <- c("source", "target")
  
  edge_df <- merge(ica_space_df, edge_list, by.x="sample_name", by.y="source", all=TRUE)
  #edge_df <- ica_space_df
  edge_df <- plyr::rename(edge_df, c("prin_graph_dim_1"="source_prin_graph_dim_1", "prin_graph_dim_2"="source_prin_graph_dim_2"))
  edge_df <- merge(edge_df, ica_space_df[,c("sample_name", "prin_graph_dim_1", "prin_graph_dim_2")], by.x="target", by.y="sample_name", all=TRUE)
  edge_df <- plyr::rename(edge_df, c("prin_graph_dim_1"="target_prin_graph_dim_1", "prin_graph_dim_2"="target_prin_graph_dim_2"))
  
  S_matrix <- reducedDimS(cds)
  data_df <- data.frame(t(S_matrix[c(x,y),]))
  colnames(data_df) <- c("data_dim_1", "data_dim_2")
  data_df$sample_name <- row.names(data_df)
  data_df <- merge(data_df, lib_info_with_pseudo, by.x="sample_name", by.y="row.names")
  
  markers_exprs <- NULL
  if (is.null(markers) == FALSE){
    markers_fData <- subset(fData(cds), gene_short_name %in% markers)
    if (nrow(markers_fData) >= 1){
      markers_exprs <- reshape2::melt(as.matrix(exprs(cds[row.names(markers_fData),])))
      markers_exprs <- merge(markers_exprs, markers_fData, by.x = "Var1", by.y="row.names")
      #print (head( markers_exprs[is.na(markers_exprs$gene_short_name) == FALSE,]))
      markers_exprs$feature_label <- as.character(markers_exprs$gene_short_name)
      markers_exprs$feature_label[is.na(markers_exprs$feature_label)] <- markers_exprs$Var1
    }
  }
  if (is.null(markers_exprs) == FALSE && nrow(markers_exprs) > 0){
    data_df <- merge(data_df, markers_exprs, by.x="sample_name", by.y="Var2")
    #print (head(edge_df))
    g <- ggplot(data=data_df, aes(x=data_dim_1, y=data_dim_2, size=log10(value + 0.1))) + facet_wrap(~feature_label)
  }else{
    g <- ggplot(data=data_df, aes(x=data_dim_1, y=data_dim_2)) 
  }
  if (show_tree){
    g <- g + geom_segment(aes_string(x="source_prin_graph_dim_1", y="source_prin_graph_dim_2", xend="target_prin_graph_dim_1", yend="target_prin_graph_dim_2"), size=.3, linetype="solid", na.rm=TRUE, data=edge_df)
  }
  
  # FIXME: setting size here overrides the marker expression funtionality. 
  # Don't do it!
  if (is.null(markers_exprs) == FALSE && nrow(markers_exprs) > 0){
    g <- g + geom_point(aes_string(color = color_by), na.rm = TRUE)
  }else {
    g <- g + geom_point(aes_string(color = color_by), size=I(cell_size), na.rm = TRUE)
  }
  
  
  if (show_branch_points & cds@auxOrderingData[[cds@dim_reduce_type]] == 'DDRTree'){
    mst_branch_nodes <- cds@auxOrderingData[[cds@dim_reduce_type]]$branch_points
    branch_point_df <- subset(edge_df, sample_name %in% mst_branch_nodes)[,c("sample_name", "source_prin_graph_dim_1", "source_prin_graph_dim_2")]
    branch_point_df$branch_point_idx <- match(branch_point_df$sample_name, mst_branch_nodes)
    branch_point_df <- branch_point_df[!duplicated(branch_point_df$branch_point_idx), ]
    
    g <- g + geom_point(aes_string(x="source_prin_graph_dim_1", y="source_prin_graph_dim_2"), 
                        size=5, na.rm=TRUE, data=branch_point_df) +
      geom_text(aes_string(x="source_prin_graph_dim_1", y="source_prin_graph_dim_2", label="branch_point_idx"), 
                size=4, color="white", na.rm=TRUE, data=branch_point_df)
  }
  if (show_cell_names){
    g <- g +geom_text(aes(label=sample_name), size=cell_name_size)
  }
  g <- g + 
    #scale_color_brewer(palette="Set1") +
    theme(panel.border = element_blank(), axis.line = element_line()) +
    theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) +
    theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) + 
    ylab("Component 1") + xlab("Component 2") +
    theme(legend.position="top", legend.key.height=grid::unit(0.35, "in")) +
    #guides(color = guide_legend(label.position = "top")) +
    theme(legend.key = element_blank()) +
    theme(panel.background = element_rect(fill='white'))
  g
}



#' Plots expression for one or more genes as a jittered, grouped points
#'
#' @param cds_subset CellDataSet for the experiment
#' @param grouping the cell attribute (e.g. the column of pData(cds)) to group cells by on the horizontal axis
#' @param min_expr the minimum (untransformed) expression level to use in plotted the genes.
#' @param cell_size the size (in points) of each cell used in the plot
#' @param nrow the number of rows used when laying out the panels for each gene's expression
#' @param ncol the number of columns used when laying out the panels for each gene's expression
#' @param panel_order the order in which genes should be layed out (left-to-right, top-to-bottom)
#' @param color_by the cell attribute (e.g. the column of pData(cds)) to be used to color each cell  
#' @param plot_trend whether to plot a trendline tracking the average expression across the horizontal axis.
#' @param label_by_short_name label figure panels by gene_short_name (TRUE) or feature id (FALSE)
#' @param relative_expr Whether to transform expression into relative values
#' @return a ggplot2 plot object
#' @import ggplot2
#' @importFrom reshape2 melt
#' @export
#' @examples
#' \dontrun{
#' data(HSMM)
#' MYOG_ID1 <- HSMM[row.names(subset(fData(HSMM), gene_short_name %in% c("MYOG", "ID1"))),]
#' plot_genes_jitter(MYOG_ID1, grouping="Media", ncol=2)
#' }
plot_genes_jitter <- function(cds_subset, grouping = "State", 
                              min_expr=NULL, cell_size=0.75, nrow=NULL, ncol=1, panel_order=NULL, 
                              color_by=NULL,
                              plot_trend=FALSE,
                              label_by_short_name=TRUE,
                              relative_expr=TRUE){
  
  if (cds_subset@expressionFamily@vfamily == "negbinomial" | cds_subset@expressionFamily@vfamily == "negbinomial.size"){
    integer_expression <- TRUE
  }else{
    integer_expression <- FALSE
    relative_expr <- TRUE
  }
  
  if (integer_expression)
  {
    cds_exprs <- exprs(cds_subset)
    if (relative_expr){
      if (is.null(sizeFactors(cds_subset)))
      {
        stop("Error: to call this function with relative_expr=TRUE, you must call estimateSizeFactors() first")
      }
      cds_exprs <- Matrix::t(Matrix::t(cds_exprs) / sizeFactors(cds_subset))
    }
    cds_exprs <- reshape2::melt(round(as.matrix(cds_exprs)))
  }else{
    cds_exprs <- exprs(cds_subset)
    cds_exprs <- reshape2::melt(as.matrix(cds_exprs))
  }
  if (is.null(min_expr)){
    min_expr <- cds_subset@lowerDetectionLimit
  }
  
  colnames(cds_exprs) <- c("f_id", "Cell", "expression")
  cds_exprs$expression[cds_exprs$expression < min_expr] <- min_expr
  cds_pData <- pData(cds_subset)
  cds_fData <- fData(cds_subset)
  
  cds_exprs <- merge(cds_exprs, cds_fData, by.x="f_id", by.y="row.names")
  cds_exprs <- merge(cds_exprs, cds_pData, by.x="Cell", by.y="row.names")
  
  cds_exprs$adjusted_expression <- log10(cds_exprs$expression)
  #cds_exprs$adjusted_expression <- log10(cds_exprs$adjusted_expression + abs(rnorm(nrow(cds_exprs), min_expr, sqrt(min_expr))))
  
  if (label_by_short_name == TRUE){
    if (is.null(cds_exprs$gene_short_name) == FALSE){
      cds_exprs$feature_label <- cds_exprs$gene_short_name
      cds_exprs$feature_label[is.na(cds_exprs$feature_label)]  <- cds_exprs$f_id
    }else{
      cds_exprs$feature_label <- cds_exprs$f_id
    }
  }else{
    cds_exprs$feature_label <- cds_exprs$f_id
  }
  
  #print (head(cds_exprs))
  
  if (is.null(panel_order) == FALSE)
  {
    cds_subset$feature_label <- factor(cds_subset$feature_label, levels=panel_order)
  }
  
  q <- ggplot(aes_string(x=grouping, y="expression"), data=cds_exprs) 
  
  if (is.null(color_by) == FALSE){
    q <- q + geom_jitter(aes_string(color=color_by), size=I(cell_size))
  }else{
    q <- q + geom_jitter(size=I(cell_size))
  }
  if (plot_trend == TRUE){
    q <- q + stat_summary(aes_string(color=color_by), fun.data = "mean_cl_boot", size=0.35)
    q <- q + stat_summary(aes_string(x=grouping, y="expression", color=color_by, group=color_by), fun.data = "mean_cl_boot", size=0.35, geom="line")
  }
  
  q <- q + scale_y_log10() + facet_wrap(~feature_label, nrow=nrow, ncol=ncol, scales="free_y")
  
  # Need this to gaurd against plotting failures caused by non-expressed genes
  if (min_expr < 1)
  {
    q <- q + expand_limits(y=c(min_expr, 1))
  }

  q <- q + ylab("Expression") + xlab(grouping)
  q <- q + monocle_theme_opts()
  q
}

#' Plots the number of cells expressing one or more genes as a barplot 
#'
#' @param cds_subset CellDataSet for the experiment
#' @param grouping the cell attribute (e.g. the column of pData(cds)) to group cells by on the horizontal axis
#' @param min_expr the minimum (untransformed) expression level to use in plotted the genes.
#' @param nrow the number of rows used when laying out the panels for each gene's expression
#' @param ncol the number of columns used when laying out the panels for each gene's expression
#' @param panel_order the order in which genes should be layed out (left-to-right, top-to-bottom)
#' @param plot_as_fraction whether to show the percent instead of the number of cells expressing each gene 
#' @param label_by_short_name label figure panels by gene_short_name (TRUE) or feature id (FALSE)
#' @param relative_expr Whether to transform expression into relative values
#' @return a ggplot2 plot object
#' @import ggplot2
#' @importFrom plyr ddply
#' @importFrom reshape2 melt
#' @export
#' @examples
#' \dontrun{
#' data(HSMM)
#' MYOG_ID1 <- HSMM[row.names(subset(fData(HSMM), gene_short_name %in% c("MYOG", "ID1"))),]
#' plot_genes_positive_cells(MYOG_ID1, grouping="Media", ncol=2)
#' }
plot_genes_positive_cells <- function(cds_subset, 
                                      grouping = "State", 
                                      min_expr=0.1, 
                                      nrow=NULL, 
                                      ncol=1, 
                                      panel_order=NULL, 
                                      plot_as_fraction=TRUE,
                                      label_by_short_name=TRUE,
                                      relative_expr=TRUE){
  
  percent <- NULL
  
  if (cds_subset@expressionFamily@vfamily == "negbinomial" | cds_subset@expressionFamily@vfamily == "negbinomial.size"){
    integer_expression <- TRUE
  }else{
    integer_expression <- FALSE
    relative_expr <- TRUE
  }
  
  if (integer_expression)
  {
    marker_exprs <- exprs(cds_subset)
    if (relative_expr){
      if (is.null(sizeFactors(cds_subset)))
      {
        stop("Error: to call this function with relative_expr=TRUE, you must call estimateSizeFactors() first")
      }
      marker_exprs <- Matrix::t(Matrix::t(marker_exprs) / sizeFactors(cds_subset))
    }
    marker_exprs_melted <- reshape2::melt(round(as.matrix(marker_exprs)))
  }else{
    marker_exprs_melted <- reshape2::melt(exprs(marker_exprs))
  }
   
  colnames(marker_exprs_melted) <- c("f_id", "Cell", "expression")
  
  marker_exprs_melted <- merge(marker_exprs_melted, pData(cds_subset), by.x="Cell", by.y="row.names")
  marker_exprs_melted <- merge(marker_exprs_melted, fData(cds_subset), by.x="f_id", by.y="row.names")
  
  if (label_by_short_name == TRUE){
    if (is.null(marker_exprs_melted$gene_short_name) == FALSE){
      marker_exprs_melted$feature_label <- marker_exprs_melted$gene_short_name
      marker_exprs_melted$feature_label[is.na(marker_exprs_melted$feature_label)]  <- marker_exprs_melted$f_id
    }else{
      marker_exprs_melted$feature_label <- marker_exprs_melted$f_id
    }    
  }else{
    marker_exprs_melted$feature_label <- marker_exprs_melted$f_id
  }
  
  if (is.null(panel_order) == FALSE)
  {
    marker_exprs_melted$feature_label <- factor(marker_exprs_melted$feature_label, levels=panel_order)
  }
  
  print (head(marker_exprs_melted))
  
  marker_counts <- plyr::ddply(marker_exprs_melted, c("feature_label", grouping), function(x) { 
    data.frame(target=sum(x$expression > min_expr), 
               target_fraction=sum(x$expression > min_expr)/nrow(x)) } )
  
  print (head(marker_counts))
  if (plot_as_fraction){
    marker_counts$target_fraction <- marker_counts$target_fraction * 100
    qp <- ggplot(aes_string(x=grouping, y="target_fraction", fill=grouping), data=marker_counts) +
      ylab("Cells (percent)") + scale_y_continuous(limits=c(0,100)) 
  }else{
    qp <- ggplot(aes_string(x=grouping, y="target", fill=grouping), data=marker_counts) +
      ylab("Cells")
  }
  
  qp <- qp + facet_wrap(~feature_label, nrow=nrow, ncol=ncol, scales="free_y")
  qp <-  qp + geom_bar(stat="identity") + monocle_theme_opts()

  return(qp)
}


#' Plots expression for one or more genes as a function of pseudotime
#'
#' @param cds_subset CellDataSet for the experiment
#' @param min_expr the minimum (untransformed) expression level to use in plotted the genes.
#' @param cell_size the size (in points) of each cell used in the plot
#' @param nrow the number of rows used when laying out the panels for each gene's expression
#' @param ncol the number of columns used when laying out the panels for each gene's expression
#' @param panel_order the order in which genes should be layed out (left-to-right, top-to-bottom)
#' @param color_by the cell attribute (e.g. the column of pData(cds)) to be used to color each cell 
#' @param trend_formula the model formula to be used for fitting the expression trend over pseudotime 
#' @param label_by_short_name label figure panels by gene_short_name (TRUE) or feature id (FALSE)
#' @param relative_expr Whether to transform expression into relative values
#' @return a ggplot2 plot object
#' @import ggplot2
#' @importFrom plyr ddply .
#' @importFrom reshape2 melt
#' @export
#' @examples
#' \dontrun{
#' data(HSMM)
#' my_genes <- row.names(subset(fData(HSMM), gene_short_name %in% c("CDK1", "MEF2C", "MYH3"))) 
#' cds_subset <- HSMM[my_genes,]
#' plot_genes_in_pseudotime(cds_subset, color_by="Time")
#' }
plot_genes_in_pseudotime <-function(cds_subset, 
                                    min_expr=NULL, 
                                    cell_size=0.75, 
                                    nrow=NULL, 
                                    ncol=1, 
                                    panel_order=NULL, 
                                    color_by="State",
                                    trend_formula="~ sm.ns(Pseudotime, df=3)",
                                    label_by_short_name=TRUE,
                                    relative_expr=TRUE,
                                    method = "VGAM"){
  
    if (cds_subset@expressionFamily@vfamily == "negbinomial" | cds_subset@expressionFamily@vfamily == "negbinomial.size") {
        integer_expression <- TRUE
    }
    else {
        integer_expression <- FALSE
        relative_expr <- TRUE
    }
    if (integer_expression) {
        cds_exprs <- exprs(cds_subset)
        if (relative_expr) {
            if (is.null(sizeFactors(cds_subset))) {
                stop("Error: to call this function with relative_expr=TRUE, you must call estimateSizeFactors() first")
            }
            cds_exprs <- Matrix::t(Matrix::t(cds_exprs)/sizeFactors(cds_subset))
        }
        cds_exprs <- reshape2::melt(round(as.matrix(cds_exprs)))
    }
    else {
        cds_exprs <- reshape2::melt(as.matrix(exprs(cds_subset)))
    }
    if (is.null(min_expr)) {
        min_expr <- cds_subset@lowerDetectionLimit
    }
    colnames(cds_exprs) <- c("f_id", "Cell", "expression")
    cds_pData <- pData(cds_subset)
    cds_fData <- fData(cds_subset)
    cds_exprs <- merge(cds_exprs, cds_fData, by.x = "f_id", by.y = "row.names")
    cds_exprs <- merge(cds_exprs, cds_pData, by.x = "Cell", by.y = "row.names")
    if (integer_expression) {
        cds_exprs$adjusted_expression <- cds_exprs$expression
    }
    else {
        cds_exprs$adjusted_expression <- log10(cds_exprs$expression)
    }
    # trend_formula <- paste("adjusted_expression", trend_formula,
    #     sep = "")
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
    cds_exprs$feature_label <- factor(cds_exprs$feature_label)

    new_data <- data.frame(Pseudotime = pData(cds_subset)$Pseudotime)
    model_expectation <- genSmoothCurves(cds_subset, cores=1, trend_formula = trend_formula,
                        relative_expr = T, pseudocount = 0, new_data = new_data, weights = pData(cds_subset)$weight)
    colnames(model_expectation) <- colnames(cds_subset)

    cds_exprs$expectation <- apply(cds_exprs,1, function(x) model_expectation[x[2], x[1]])

    cds_exprs$expression[cds_exprs$expression < min_expr] <- min_expr
    cds_exprs$expectation[cds_exprs$expectation < min_expr] <- min_expr
    if (is.null(panel_order) == FALSE) {
        cds_subset$feature_label <- factor(cds_subset$feature_label,
            levels = panel_order)
    }
    q <- ggplot(aes(Pseudotime, expression), data = cds_exprs)
    if (is.null(color_by) == FALSE) {
        q <- q + geom_point(aes_string(color = color_by), size = I(cell_size))
    }
    else {
        q <- q + geom_point(size = I(cell_size))
    }
    if (method == "loess") {
        q <- q + stat_smooth(aes(group = 1), color = I("black"),
            method = "loess", se = F)
    }
    else if (method == "VGAM") {
        q <- q + geom_line(aes(x = Pseudotime, y = expectation),
            data = cds_exprs)
    }
    q <- q + scale_y_log10() + facet_wrap(~feature_label, nrow = nrow,
        ncol = ncol, scales = "free_y")
    if (min_expr < 1) {
        q <- q + expand_limits(y = c(min_expr, 1))
    }
    if (relative_expr) {
        q <- q + ylab("Relative Expression")
    }
    else {
        q <- q + ylab("Absolute Expression")
    }
    q <- q + xlab("Pseudo-time")
    q <- q + monocle_theme_opts()
    q
}

#' Plots the minimum spanning tree on cells.
#'
#' @param cds CellDataSet for the experiment
#' @param clustering a clustering object produced by clusterCells
#' @param drawSummary whether to draw the summary line for each cluster
#' @param sumFun whether the function used to generate the summary for each cluster
#' @param ncol number of columns used to layout the faceted cluster panels
#' @param nrow number of columns used to layout the faceted cluster panels
#' @param row_samples how many genes to randomly select from the data
#' @param callout_ids a vector of gene names or gene ids to manually render as part of the plot
#' @return a ggplot2 plot object
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom stringr str_join
#' @export
#' @examples
#' \dontrun{
#' full_model_fits <- fitModel(HSMM_filtered[sample(nrow(fData(HSMM_filtered)), 100),],  modelFormulaStr="~VGAM::bs(Pseudotime)")
#' expression_curve_matrix <- responseMatrix(full_model_fits)
#' clusters <- clusterGenes(expression_curve_matrix, k=4)
#' plot_clusters(HSMM_filtered[ordering_genes,], clusters)
#' }
plot_clusters<-function(cds, 
                        clustering,
                        drawSummary=TRUE, 
                        sumFun=mean_cl_boot,
                        ncol=NULL, 
                        nrow=NULL, 
                        row_samples=NULL, 
                        callout_ids=NULL,
                        conf_int=0.68){
  m <- as.data.frame(clustering$exprs)
  m$ids <- rownames(clustering$exprs)
  if (is.null(clustering$labels) == FALSE)
  {
    m$cluster = factor(clustering$labels[clustering$clustering], levels = levels(clustering$labels))
  }else{
    m$cluster <- factor(clustering$clustering)
  }
  
  cluster_sizes <- as.data.frame(table(m$cluster))    
  
  cluster_sizes$Freq <- paste("(", cluster_sizes$Freq, ")")   
  facet_labels <- str_join(cluster_sizes$Var1, cluster_sizes$Freq, sep=" ") #update the function
  
  facet_wrap_labeller <- function(gg.plot,labels=NULL) {
    #works with R 3.0.1 and ggplot2 0.9.3.1
    require(gridExtra)
    
    g <- ggplotGrob(gg.plot)
    gg <- g$grobs      
    strips <- grep("strip_t", names(gg))
    
    for(ii in seq_along(labels))  {
      modgrob <- getGrob(gg[[strips[ii]]], "strip.text", 
                         grep=TRUE, global=TRUE)
      gg[[strips[ii]]]$children[[modgrob$name]] <- editGrob(modgrob,label=labels[ii])
    }
    
    g$grobs <- gg
    class(g) = c("arrange", "ggplot",class(g)) 
    g
  }
  
  
  m.melt <- melt(m, id.vars = c("ids", "cluster"))
  
  m.melt <- merge(m.melt, pData(cds), by.x="variable", by.y="row.names")
  
  
  if (is.null(row_samples) == FALSE){
    m.melt <- m.melt[sample(nrow(m.melt), row_samples),]
  }
  
  c <- ggplot(m.melt) + facet_wrap("cluster", ncol=ncol, nrow=nrow, scales="free_y")
  #c <- c + stat_density2d(aes(x = Pseudotime, y = value), geom="polygon", fill="white", color="black", size=I(0.1)) + facet_wrap("cluster", ncol=ncol, nrow=nrow)
    
  if (drawSummary) {
    c <- c + stat_summary(aes(x = Pseudotime, y = value, group = 1),
                          fun.data = sumFun, color = "red",
                          alpha = 0.2, size = 0.5, geom = "smooth")
  }
  
  #cluster_medians <- subset(m.melt, ids %in% clustering$medoids)
  
  #c <- c + geom_line()
  #c <- c + geom_line(aes(x=Pseudotime, y=value), data=cluster_medians, color=I("red"))
  c <- c + scale_color_hue(l = 50, h.start = 200) + theme(axis.text.x = element_text(angle = 0, 
                                                                                     hjust = 0)) + xlab("Pseudo-time") + ylab("Expression")
  c <- c + theme(strip.background = element_rect(colour = 'white', fill = 'white')) + 
    theme(panel.border = element_blank()) +
    theme(legend.position="none") +
    theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) +
    theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank())
  
  #   if (draw_cluster_size){
  #     cluster_sizes <- as.data.frame(table(m$cluster))
  #     colnames(cluster_sizes) <- c("cluster", "Freq")
  #     cluster_sizes <- cbind (cluster_sizes, Pseudotime = cluster_label_text_x, value = cluster_label_text_y)
  #     c <- c + geom_text(aes(x=Pseudotime, y=value, label=Freq), data=cluster_sizes, size=cluster_label_text_size)
  #   }
  
  if (is.null(callout_ids) == FALSE)
  {
    callout_melt <- subset(m.melt, ids %in% callout_ids)
    c <- c + geom_line(aes(x=Pseudotime, y=value), data=callout_melt, color=I("steelblue"))
  }
  c <- c + monocle_theme_opts()
  #c <- facet_wrap_labeller(c, facet_labels)
  c
}

# #' Plots a pseudotime-ordered, row-centered heatmap
# #' @export 
# plot_genes_heatmap <- function(cds, 
#                                rescaling='row', 
#                                clustering='row', 
#                                labCol=FALSE, 
#                                labRow=TRUE, 
#                                logMode=TRUE, 
#                                pseudocount=0.1, 
#                                use_vst=TRUE,
#                                border=FALSE, 
#                                heatscale=c(low='steelblue',mid='white',high='tomato'), 
#                                heatMidpoint=0,
#                                method="none",
#                                scaleMax=2, 
#                                scaleMin=-2, 
#                                relative_expr=TRUE, 
#                                ...){
  
#   ## the function can be be viewed as a two step process
#   ## 1. using the rehape package and other funcs the data is clustered, scaled, and reshaped
#   ## using simple options or by a user supplied function
#   ## 2. with the now resahped data the plot, the chosen labels and plot style are built
#   FM <- exprs(cds)
  
#   if (cds@expressionFamily@vfamily == "negbinomial"){
#     integer_expression <- TRUE
#   }else{
#     integer_expression <- FALSE
#     relative_expr <- TRUE
#   }
  
#   if (integer_expression)
#   {
#     if (relative_expr){
#       if (is.null(sizeFactors(cds)))
#       {
#         stop("Error: you must call estimateSizeFactors() first")
#       }
#       FM <- Matrix::t(Matrix::t(FM) / sizeFactors(cds))
#     }
#     FM <- round(FM)
#   }
  
#   m=FM
  
#   if (is.null(fData(cds)$gene_short_name) == FALSE){
#     feature_labels <- fData(cds)$gene_short_name
#     feature_labels[is.na(feature_labels)]  <- fData(cds)$f_id
#     row.names(m) <- feature_labels
#   }
  
#   #remove genes with no expression in any condition
#   m=m[!apply(m,1,sum)==0,]
  
#   if (use_vst && is.null(cds@dispFitInfo[["blind"]]$disp_func) == FALSE){
#     m = vstExprs(cds, expr_matrix=m)
#   }else if(logMode){
#     m = log10(m+pseudocount)
#   }
  
#   #remove genes with no sd
#   #m=m[!apply(m,1,sd)==0,]

#   ## you can either scale by row or column not both! 
#   ## if you wish to scale by both or use a different scale method then simply supply a scale
#   ## function instead NB scale is a base funct
  
#   ## I have supplied the default cluster and euclidean distance (JSdist) - and chose to cluster after scaling
#   ## if you want a different distance/cluster method-- or to cluster and then scale
#   ## then you can supply a custom function 
  
#   if(!is.function(method)){
#     method = function(mat){as.dist((1 - cor(Matrix::t(mat)))/2)}	
#   }
  
#   ## this is just reshaping into a ggplot format matrix and making a ggplot layer
  
#   if(is.function(rescaling))
#   { 
#     m=rescaling(m)
#   } else {
#     if(rescaling=='column'){
#       m=m[!apply(m,2,sd)==0,]
#       m=scale(m, center=TRUE)
#       m[is.nan(m)] = 0
#       m[m>scaleMax] = scaleMax
#       m[m<scaleMin] = scaleMin
#     }
#     if(rescaling=='row'){ 
#       m=m[!apply(m,1,sd)==0,]
#       m=Matrix::t(scale(Matrix::t(m),center=TRUE))
#       m[is.nan(m)] = 0
#       m[m>scaleMax] = scaleMax
#       m[m<scaleMin] = scaleMin
#     }
#   }
  
#   # If we aren't going to re-ordering the columns, order them by Pseudotime
#   if (clustering %in% c("row", "none"))
#     m = m[,row.names(pData(cds)[order(-pData(cds)$Pseudotime),])]
  
#   if(clustering=='row')
#     m=m[hclust(method(m))$order, ]
#   if(clustering=='column')  
#     m=m[,hclust(method(Matrix::t(m)))$order]
#   if(clustering=='both')
#     m=m[hclust(method(m))$order ,hclust(method(Matrix::t(m)))$order]
  
  
#   rows=dim(m)[1]
#   cols=dim(m)[2]
  
  
  
#   # if(logMode) {
#   #   melt.m=cbind(rowInd=rep(1:rows, times=cols), colInd=rep(1:cols, each=rows), reshape2::melt( log10(m+pseudocount)))
#   # }else{
#   #   melt.m=cbind(rowInd=rep(1:rows, times=cols), colInd=rep(1:cols, each=rows), reshape2::melt(m))
#   # }
  
#   melt.m=cbind(rowInd=rep(1:rows, times=cols), colInd=rep(1:cols, each=rows), reshape2::melt(m))
  
#   g=ggplot(data=melt.m)
  
#   ## add the heat tiles with or without a white border for clarity
  
#   if(border==TRUE)
#     g2=g+geom_raster(aes(x=colInd,y=rowInd, fill=value),colour='grey')
#   if(border==FALSE)
#     g2=g+geom_raster(aes(x=colInd,y=rowInd,ymax=rowInd, fill=value))
  
#   ## add axis labels either supplied or from the colnames rownames of the matrix
  
#   if(labCol==TRUE) 
#   {
#     g2=g2+scale_x_continuous(breaks=(1:cols)-0.5, labels=colnames(m))
#   }
#   if(labCol==FALSE) 
#   {
#     g2=g2+scale_x_continuous(breaks=(1:cols)-0.5, labels=rep('',cols))
#   }
  
  
#   if(labRow==TRUE) 
#   {
#     g2=g2+scale_y_continuous(breaks=(1:rows)-0.5, labels=rownames(m))	
#   }
#   if(labRow==FALSE)
#   { 
#     g2=g2+scale_y_continuous(breaks=(1:rows)-0.5, labels=rep('',rows))	
#   }
  
#   # Get rid of the ticks, they get way too dense with lots of rows
#   g2 <- g2 + theme(axis.ticks = element_blank()) 
  
#   ## get rid of grey panel background and gridlines
  
#   g2=g2+theme(panel.grid.minor=element_line(colour=NA), panel.grid.major=element_line(colour=NA),
#               panel.background=element_rect(fill=NA, colour=NA))
  
#   ##adjust x-axis labels
#   g2=g2+theme(axis.text.x=element_text(angle=-90, hjust=0))
  
#   #write(paste(c("Length of heatscale is :", length(heatscale))), stderr())
  
#   if(is.function(rescaling))
#   {
    
#   }else{ 
#     if(rescaling=='row' || rescaling == 'column'){
#       legendTitle <- "Relative\nexpression"
#     }else{
#       if (logMode)
#       {
#         legendTitle <- bquote(paste(log[10]," FPKM + ",.(pseudocount),sep=""))
#         #legendTitle <- paste(expression(plain(log)[10])," FPKM + ",pseudocount,sep="")
#       } else {
#         legendTitle <- "FPKM"
#       }
#     }
#   }
  
#   if (length(heatscale) == 2){
#     g2 <- g2 + scale_fill_gradient(low=heatscale[1], high=heatscale[2], name=legendTitle)
#   } else if (length(heatscale) == 3) {
#     if (is.null(heatMidpoint))
#     {
#       heatMidpoint = (max(m) + min(m)) / 2.0
#       #write(heatMidpoint, stderr())
#     }
#     g2 <- g2 + theme(panel.border = element_blank())
#     g2 <- g2 + scale_fill_gradient2(low=heatscale[1], mid=heatscale[2], high=heatscale[3], midpoint=heatMidpoint, name=legendTitle)
#   }else {
#     g2 <- g2 + scale_fill_gradientn(colours=heatscale, name=legendTitle)
#   }
  
#   #g2<-g2+scale_x_discrete("",breaks=tracking_ids,labels=gene_short_names)
  
#   g2 <- g2 + theme(axis.title.x=element_blank(), axis.title.y=element_blank())
  
#   ## finally add the fill colour ramp of your choice (default is blue to red)-- and return
#   return (g2)
# }


plot_genes_heatmap <- function(...){
  .Deprecated("plot_pseudotime_heatmap")
  plot_pseudotime_heatmap(...)
}

#' Plots a pseudotime-ordered, row-centered heatmap
#' 
#' @param cds_subset CellDataSet for the experiment (normally only the branching genes detected with branchTest)
#' @param cluster_rows Whether to cluster the rows of the heatmap.
#' @param hclust_method The method used by pheatmap to perform hirearchical clustering of the rows. 
#' @param num_clusters Number of clusters for the heatmap of branch genes
#' @param hmcols The color scheme for drawing the heatmap.
#' @param add_annotation_row Additional annotations to show for each row in the heatmap. Must be a dataframe with one row for each row in the fData table of cds_subset, with matching IDs.
#' @param add_annotation_col Additional annotations to show for each column in the heatmap. Must be a dataframe with one row for each cell in the pData table of cds_subset, with matching IDs.
#' @param show_rownames Whether to show the names for each row in the table.
#' @param use_gene_short_name Whether to use the short names for each row. If FALSE, uses row IDs from the fData table.
#' @param scale_max The maximum value (in standard deviations) to show in the heatmap. Values larger than this are set to the max.
#' @param scale_min The minimum value (in standard deviations) to show in the heatmap. Values smaller than this are set to the min.
#' @param norm_method Determines how to transform expression values prior to rendering
#' @param trend_formula A formula string specifying the model used in fitting the spline curve for each gene/feature.
#' @param return_heatmap Whether to return the pheatmap object to the user. 
#' @param cores Number of cores to use when smoothing the expression curves shown in the heatmap.
#' @return A list of heatmap_matrix (expression matrix for the branch committment), ph (pheatmap heatmap object),
#' annotation_row (annotation data.frame for the row), annotation_col (annotation data.frame for the column). 
#' @import pheatmap
#' @export
#'

plot_pseudotime_heatmap <- function(cds_subset, 
                                    
                                    cluster_rows = TRUE,
                                    hclust_method = "ward.D2", 
                                    num_clusters = 6,
                                    
                                    hmcols = NULL, 
                                    
                                    add_annotation_row = NULL,
                                    add_annotation_col = NULL,
                                    show_rownames = FALSE, 
                                    use_gene_short_name = TRUE,
                                    
                                    norm_method = c("vstExprs", "log"), 
                                    scale_max=3, 
                                    scale_min=-3, 
                                    
                                    trend_formula = '~sm.ns(Pseudotime, df=3)',
                                    
                                    return_heatmap=FALSE,
                                    cores=1){
  
  newdata <- data.frame(Pseudotime = seq(0, max(pData(cds_subset)$Pseudotime),length.out = 100)) 
  
  m <- genSmoothCurves(cds_subset, cores=cores, trend_formula = trend_formula,  
                       relative_expr = T, new_data = newdata)
  

  #remove genes with no expression in any condition
  m=m[!apply(m,1,sum)==0,]
  
  norm_method <- match.arg(norm_method)
  
  # FIXME: this needs to check that vst values can even be computed. (They can only be if we're using NB as the expressionFamily)
  if(norm_method == 'vstExprs' && is.null(cds_subset@dispFitInfo[["blind"]]$disp_func) == FALSE) {
    m = vstExprs(cds_subset, expr_matrix=m)
  }     
  else if(norm_method == 'log') {
    m = log10(m+pseudocount)
  }
  
  # Row-center the data.
  m=m[!apply(m,1,sd)==0,]
  m=Matrix::t(scale(Matrix::t(m),center=TRUE))
  m=m[is.na(row.names(m)) == FALSE,]
  m[is.nan(m)] = 0
  m[m>scale_max] = scale_max
  m[m<scale_min] = scale_min

  heatmap_matrix <- m
  
  row_dist <- as.dist((1 - cor(Matrix::t(heatmap_matrix)))/2)
  row_dist[is.na(row_dist)] <- 1
  
  bks <- seq(-3.1,3.1, by=0.1)
  if(is.null(hmcols)) {
    hmcols <- blue2green2red(length(bks) - 1)
  }
  
  ph <- pheatmap(heatmap_matrix, 
                 useRaster = T,
                 cluster_cols=FALSE, 
                 cluster_rows=cluster_rows, 
                 show_rownames=F, 
                 show_colnames=F, 
                 clustering_distance_rows=row_dist,
                 clustering_method = hclust_method,
                 cutree_rows=num_clusters,
                 silent=TRUE,
                 filename=NA,
                 breaks=bks,
                 color=hmcols)

  annotation_row <- data.frame(Cluster=factor(cutree(ph$tree_row, num_clusters)))
 
  if(!is.null(add_annotation_row)) {
    old_colnames_length <- ncol(annotation_row)
    annotation_row <- cbind(annotation_row, add_annotation_row[row.names(annotation_row), ])  
    colnames(annotation_row)[(old_colnames_length+1):ncol(annotation_row)] <- colnames(add_annotation_row)
    # annotation_row$bif_time <- add_annotation_row[as.character(fData(absolute_cds[row.names(annotation_row), ])$gene_short_name), 1]
  }
  
  
  if (use_gene_short_name == TRUE) {
    if (is.null(fData(cds_subset)$gene_short_name) == FALSE) {
      feature_label <- as.character(fData(cds_subset)[row.names(heatmap_matrix), 'gene_short_name'])
      feature_label[is.na(feature_label)] <- row.names(heatmap_matrix)
      
      row_ann_labels <- as.character(fData(cds_subset)[row.names(annotation_row), 'gene_short_name'])
      row_ann_labels[is.na(row_ann_labels)] <- row.names(annotation_row)
    }
    else {
      feature_label <- row.names(heatmap_matrix)
      row_ann_labels <- row.names(annotation_row)
    }
  }
  else {
    feature_label <- row.names(heatmap_matrix)
    row_ann_labels <- row.names(annotation_row)
  }
  
  row.names(heatmap_matrix) <- feature_label
  row.names(annotation_row) <- row_ann_labels
  
  
  colnames(heatmap_matrix) <- c(1:ncol(heatmap_matrix))
  
  ph_res <- pheatmap(heatmap_matrix[, ], #ph$tree_row$order
                     useRaster = T,
                     cluster_cols = FALSE, 
                     cluster_rows = cluster_rows, 
                     show_rownames=show_rownames, 
                     show_colnames=F, 
                     #scale="row",
                     clustering_distance_rows=row_dist, #row_dist
                     clustering_method = hclust_method, #ward.D2
                     cutree_rows=num_clusters,
                     # cutree_cols = 2,
                     annotation_row=annotation_row,
                     treeheight_row = 20, 
                     breaks=bks,
                     fontsize = 6,
                     color=hmcols, 
                     silent=TRUE,
                     filename=NA
  )
  
  grid::grid.rect(gp=grid::gpar("fill", col=NA))
  grid::grid.draw(ph_res$gtable)
  if (return_heatmap){
    return(ph_res)
  }
}

#' Plot the branch genes in pseduotime with separate lineage curves.
#' 
#' This plotting function is used to make the branching plots for a lineage dependent gene goes through the progenitor state
#' and bifurcating into two distinct lineages (Similar to the pitch-fork bifurcation in dynamic systems). In order to make the  
#' bifurcation plot, we first duplicated the progenitor states and by default stretch each lineage into maturation level 0-100.  
#' Then we fit two nature spline curves for each lineages using VGAM package.  
#'
#' @param cds CellDataSet for the experiment
#' @param lineage_states The states for two branching lineages
#' @param lineage_labels The names for each branching lineage
#' @param method The method to draw the curve for the gene expression branching pattern, either loess ('loess') or VGLM fitting ('fitting') 
#' @param stretch A logic flag to determine whether or not the pseudotime trajectory for each lineage should be stretched to the same range or not 
#' @param min_expr the minimum (untransformed) expression level to use in plotted the genes.
#' @param cell_size the size (in points) of each cell used in the plot
#' @param nrow number of columns used to layout the faceted cluster panels
#' @param ncol number of columns used to layout the faceted cluster panels
#' @param panel_order the order in which genes should be layed out (left-to-right, top-to-bottom)
#' @param color_by the cell attribute (e.g. the column of pData(cds)) to be used to color each cell 
#' @param trend_formula the model formula to be used for fitting the expression trend over pseudotime
#' @param label_by_short_name label figure panels by gene_short_name (TRUE) or feature id (FALSE)
#' @param weighted  A logic flag to determine whether or not we should use the navie logLikelihood weight scheme for the duplicated progenitor cells
#' @param add_ABC A logic flag to determine whether or not we should add the ABC score for each gene 
#' @param relative_expr A logic flag to determine whether or not we should use the relative expression values
#' @return a ggplot2 plot object
#' @import ggplot2
#' @importFrom plyr ddply
#' @importFrom reshape2 melt
#' @export 
#' @examples
#' \dontrun{
#' full_model_fits <- fitModel(HSMM_filtered[sample(nrow(fData(HSMM_filtered)), 100),],  modelFormulaStr="~VGAM::bs(Pseudotime)")
#' expression_curve_matrix <- responseMatrix(full_model_fits)
#' clusters <- clusterGenes(expression_curve_matrix, k=4)
#' plot_clusters(HSMM_filtered[ordering_genes,], clusters)
#' }
#' 

plot_genes_branched_pseudotime <- function (cds, 
                                            lineage_states = c(2, 3), 
                                            branch_point=NULL,
                                            lineage_labels = NULL,
                                            method = "fitting", 
                                            stretch = TRUE, 
                                            min_expr = NULL, 
                                            cell_size = 0.75,
                                            nrow = NULL, 
                                            ncol = 1, 
                                            panel_order = NULL, 
                                            color_by = "State",
                                            trajectory_linetype_by = "Lineage", 
                                            trend_formula = "~ sm.ns(Pseudotime, df=3) * Lineage", 
                                            reducedModelFormulaStr = NULL, 
                                            label_by_short_name = TRUE,
                                            weighted = TRUE, 
                                            add_ABC = FALSE, 
                                            add_pval = FALSE,
                                            normalize = TRUE,
                                            bifurcation_time = NULL, 
                                            #gene_pairs = NULL,
    ...)
{
    if (add_ABC) {
        ABCs_df <- calABCs(cds, 
                           trend_formula = trend_formula,
                           lineage_states = lineage_states, 
                           branch_point=branch_point,
                           stretch = stretch,
                           weighted = weighted, 
                           min_expr = min_expr, 
                           lineage_labels = lineage_labels,
            ...)
        fData(cds)[, "ABCs"] <- ABCs_df$ABCs
    }
    if (add_pval) {
        pval_df <- branchTest(cds, 
                              lineage_states=lineage_states,
                              branch_point=branch_point,
                              fullModelFormulaStr = trend_formula,
            reducedModelFormulaStr = "~ sm.ns(Pseudotime, df=3)", ...)
        fData(cds)[, "pval"] <- pval_df[row.names(cds), 'pval']
    }
    if("Lineage" %in% all.vars(terms(as.formula(trend_formula)))) { #only when Lineage is in the model formula we will duplicate the "progenitor" cells
        cds_subset <- buildLineageBranchCellDataSet(cds = cds, 
                                                    lineage_states = lineage_states, 
                                                    branch_point=branch_point,
                                                    lineage_labels = lineage_labels, 
                                                    method = method, 
                                                    stretch = stretch,
                                                    weighted = weighted, 
                                                    ...)
    }
    else {
        cds_subset <- cds
        pData(cds_subset)$Lineage <- pData(cds_subset)$State
    }
    if (cds_subset@expressionFamily@vfamily %in% c("zanegbinomialff",
        "negbinomial", "poissonff", "quasipoissonff")) {
        integer_expression <- TRUE
    }
    else {
        integer_expression <- FALSE
    }
    if (integer_expression) {
        CM <- exprs(cds_subset)
        if (normalize)
            CM <- Matrix::t(Matrix::t(CM)/sizeFactors(cds_subset))
        cds_exprs <- reshape2::melt(round(as.matrix(CM)))
    }
    else {
        cds_exprs <- reshape2::melt(exprs(cds_subset))
    }
    if (is.null(min_expr)) {
        min_expr <- cds_subset@lowerDetectionLimit
    }
    colnames(cds_exprs) <- c("f_id", "Cell", "expression")
    cds_pData <- pData(cds_subset)
    if (add_ABC)
        cds_fData <- ABCs_df
    else cds_fData <- fData(cds_subset)
    cds_exprs <- merge(cds_exprs, cds_fData, by.x = "f_id", by.y = "row.names")
    cds_exprs <- merge(cds_exprs, cds_pData, by.x = "Cell", by.y = "row.names")
    if (integer_expression) {
        cds_exprs$adjusted_expression <- round(cds_exprs$expression)
    }
    else {
        cds_exprs$adjusted_expression <- log10(cds_exprs$expression)
    }
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
    cds_exprs$Lineage <- as.factor(cds_exprs$Lineage) 

    new_data <- data.frame(Pseudotime = pData(cds_subset)$Pseudotime, Lineage = pData(cds_subset)$Lineage)
 
    full_model_expectation <- genSmoothCurves(cds_subset, cores=1, trend_formula = trend_formula, 
                        relative_expr = T, pseudocount = 0, new_data = new_data, weights = pData(cds_subset)$weight)
    colnames(full_model_expectation) <- colnames(cds_subset)
    
    cds_exprs$full_model_expectation <- apply(cds_exprs,1, function(x) full_model_expectation[x[2], x[1]])
    if(!is.null(reducedModelFormulaStr)){
        reduced_model_expectation <- genSmoothCurves(cds_subset, cores=1, trend_formula = reducedModelFormulaStr,
                            relative_expr = T, pseudocount = 0, new_data = new_data, weights = pData(cds_subset)$weight)
        colnames(reduced_model_expectation) <- colnames(cds_subset)
        cds_exprs$reduced_model_expectation <- apply(cds_exprs,1, function(x) reduced_model_expectation[x[2], x[1]])
    }

    if(!is.null(bifurcation_time)){
        cds_exprs$bifurcation_time <- bifurcation_time[as.vector(cds_exprs$gene_short_name)]
    }
    if (method == "loess")
        cds_exprs$expression <- cds_exprs$expression + cds@lowerDetectionLimit
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
    cds_exprs$feature_label <- factor(cds_exprs$feature_label)
    if (is.null(panel_order) == FALSE) {
        cds_exprs$feature_label <- factor(cds_exprs$feature_label,
            levels = panel_order)
    }
    cds_exprs$expression[is.na(cds_exprs$expression)] <- min_expr
    cds_exprs$expression[cds_exprs$expression < min_expr] <- min_expr
    cds_exprs$full_model_expectation[is.na(cds_exprs$full_model_expectation)] <- min_expr
    cds_exprs$full_model_expectation[cds_exprs$full_model_expectation < min_expr] <- min_expr
    
    if(!is.null(reducedModelFormulaStr)){
        cds_exprs$reduced_model_expectation[is.na(cds_exprs$reduced_model_expectation)] <- min_expr
        cds_exprs$reduced_model_expectation[cds_exprs$reduced_model_expectation < min_expr] <- min_expr
    }

    cds_exprs$State <- as.factor(cds_exprs$State)
    cds_exprs$Lineage <- as.factor(cds_exprs$Lineage)

    q <- ggplot(aes(Pseudotime, expression), data = cds_exprs)
    if (!is.null(bifurcation_time)) {
      q <- q + geom_vline(aes(xintercept = bifurcation_time),
                          color = "black", linetype = "longdash")
    }
    if (is.null(color_by) == FALSE) {
        q <- q + geom_point(aes_string(color = color_by), size = I(cell_size))
    }
    if (add_ABC)
        q <- q + scale_y_log10() + facet_wrap(~feature_label +
            ABCs, nrow = nrow, ncol = ncol, scales = "free_y")
    else if (add_pval)
        q <- q + scale_y_log10() + facet_wrap(~feature_label +
            pval, nrow = nrow, ncol = ncol, scales = "free_y")
    else q <- q + scale_y_log10() + facet_wrap(~feature_label,
        nrow = nrow, ncol = ncol, scales = "free_y")
    if (method == "loess")
        q <- q + stat_smooth(aes(fill = Lineage, color = Lineage),
            method = "loess")
    else if (method == "fitting") {
        q <- q + geom_line(aes_string(x = "Pseudotime", y = "full_model_expectation",
            linetype = trajectory_linetype_by), data = cds_exprs) #+ scale_color_manual(name = "Type", values = c(colour_cell, colour), labels = c("Progenitor", "AT1", "AT2", "AT1", "AT2")
    }

    if(!is.null(reducedModelFormulaStr)) {
        q <- q + geom_line(aes_string(x = "Pseudotime", y = "reduced_model_expectation"),
            color = 'black', linetype = 2, data =  cds_exprs)   
    }
    if (stretch)
        q <- q + ylab("Expression") + xlab("Pseudotime (stretched)")
    else q <- q + ylab("Expression") + xlab("Pseudotime")
    q <- q + monocle_theme_opts()
    q + expand_limits(y = min_expr)
}
#' Plot the branch genes in pseduotime with separate lineage curves 
#' @param cds CellDataSet for the experiment
#' @param rowgenes Gene ids or short names to be arrayed on the vertical axis.
#' @param colgenes Gene ids or short names to be arrayed on the horizontal axis
#' @param relative_expr Whether to transform expression into relative values
#' @param min_expr The minimum level of expression to show in the plot
#' @return a ggplot2 plot object
#' @import ggplot2
#' @importFrom reshape2 melt
#' @export 
plot_coexpression_matrix <- function(cds, 
                                     rowgenes, 
                                     colgenes, 
                                     relative_expr=TRUE, 
                                     min_expr=NULL, 
                                     cell_size=0.85, 
                                     label_by_short_name=TRUE,
                                     show_density=TRUE,
                                     round_expr=FALSE){
  
  gene_short_name <- NULL
  adjusted_expression.x <- NULL
  adjusted_expression.y <- NULL
  ..density.. <- NULL
  f_id <- NULL
  
  row_gene_ids <- row.names(subset(fData(cds), gene_short_name %in% rowgenes))
  row_gene_ids <- union(row_gene_ids, intersect(rowgenes, row.names(fData(cds))))
  
  col_gene_ids <- row.names(subset(fData(cds), gene_short_name %in% colgenes))
  col_gene_ids <- union(col_gene_ids, intersect(colgenes, row.names(fData(cds))))
  
  cds_subset <- cds[union(row_gene_ids, col_gene_ids),]
  
  if (cds_subset@expressionFamily@vfamily == "negbinomial" | cds_subset@expressionFamily@vfamily == "negbinomial.size"){
    integer_expression <- TRUE
  }else{
    integer_expression <- FALSE
    relative_expr <- TRUE
  }
  
  if (integer_expression)
  {
    cds_exprs <- exprs(cds_subset)
    if (relative_expr){
      if (is.null(sizeFactors(cds_subset)))
      {
        stop("Error: to call this function with relative_expr=TRUE, you must call estimateSizeFactors() first")
      }
      cds_exprs <- Matrix::t(Matrix::t(cds_exprs) / sizeFactors(cds_subset))
    }
    if (round_expr){
      cds_exprs <- reshape2::melt(round(as.matrix(cds_exprs)))
    } else {
      cds_exprs <- reshape2::melt(cds_exprs)
    }
      
  }else{
    cds_exprs <- reshape2::melt(exprs(cds_subset))
  }
  if (is.null(min_expr)){
    min_expr <- cds_subset@lowerDetectionLimit
  }
  
  colnames(cds_exprs) <- c("f_id", "Cell", "expression")
  cds_exprs$expression[cds_exprs$expression < min_expr] <- min_expr
  
  cds_pData <- pData(cds_subset)
  cds_fData <- fData(cds_subset)
  
  cds_exprs <- merge(cds_exprs, cds_fData, by.x="f_id", by.y="row.names")
  
  cds_exprs$adjusted_expression <- cds_exprs$expression

  #cds_exprs$adjusted_expression <- log10(cds_exprs$adjusted_expression + abs(rnorm(nrow(cds_exprs), min_expr, sqrt(min_expr))))
  
  if (label_by_short_name == TRUE){
    if (is.null(cds_exprs$gene_short_name) == FALSE){
      cds_exprs$feature_label <- as.character(cds_exprs$gene_short_name)
      cds_exprs$feature_label[is.na(cds_exprs$feature_label)]  <- cds_exprs$f_id
    }else{
      cds_exprs$feature_label <- cds_exprs$f_id
    }
  }else{
    cds_exprs$feature_label <- cds_exprs$f_id
  }
  
  cds_exprs$feature_label <- factor(cds_exprs$feature_label)
  
  row_cds_exprs <- subset(cds_exprs, f_id %in% row_gene_ids)
  col_cds_exprs <- subset(cds_exprs, f_id %in% col_gene_ids)
  
  joined_exprs <- merge(row_cds_exprs, col_cds_exprs, by="Cell")
  cds_exprs <- joined_exprs
  
  cds_exprs <- merge(cds_exprs, cds_pData, by.x="Cell", by.y="row.names")
  
  cds_exprs <- subset(cds_exprs, adjusted_expression.x > min_expr | adjusted_expression.y > min_expr)
  
  q <- ggplot(aes(adjusted_expression.x, adjusted_expression.y), data=cds_exprs, size=I(1))
  
  if (show_density){
    q <- q + stat_density2d(geom="raster", aes(fill = ..density..), contour = FALSE) + 
      scale_fill_gradient(low="white", high="red") 
  }

  q <- q + scale_x_log10() + scale_y_log10() + 
    geom_point(color=I("black"), size=I(cell_size * 1.50)) +
    geom_point(color=I("white"), size=I(cell_size)) +
    facet_grid(feature_label.x ~ feature_label.y, scales="free") 
    #scale_color_brewer(palette="Set1") +
  
  if (min_expr < 1)
  {
    q <- q + expand_limits(y=c(min_expr, 1), x=c(min_expr, 1))
  }
  
  #q <- q + monocle_theme_opts()
  
  q
}

#'  Create a heatmap to demonstrate the lineage divergence hierarchy for branching genes 
#'
#' @param cds CellDataSet for the experiment
#' @param ILRs_df Matrix of Instant Log Ratio for the branching genes calculated using the calILRs function 
#' @param ABC_df Matrix of Area Between Curves for the branching genes calculated using the calABCs function 
#' @param rownames_type The type of row names for the ILRs matrix (either 'gene_short_name' or 'ensemble_id')
#' @param ABC_type The type of genes based on ABC score to be selected (either 'positive', 'negative' or 'all')
#' @param dist_method The method to calculate distance for each gene used in the hirearchical clustering, any one of them: "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski". 
#' @param hclust_method The method to perform hirearchical clustering, any one of them: ward", "single", "complete", "average", "mcquitty", "median", "centroid"
#' @param ILRs_limit the minimum Instant Log Ratio used to make the heatmap plot
#' @param cluster_num The minimum level of expression to show in the plot
#' @return a heatmap plot from based on the pheatmap package
#' @import pheatmap
#' @export
#'
#'
plot_ILRs_heatmap <- function (cds, 
      ILRs_df, 
      ABC_df, 
      branching_pval_df, 
      rownames_type = c('gene_short_name', 'ensemble_id'), 
      ABC_type = c('positive', 'negative', 'all'),
      dist_method = "euclidean", 
      hclust_method = "ward", 
      ILRs_limit = 3, 
      cluster_num = 4, ...) {

    #clean the ILR dataset (avoid failed fittign genes or genes don't have names and cannot be identified)
    ILRs_df<- ILRs_df[!is.na(ILRs_df[, 1]), ]
    ILRs_df <- ILRs_df[row.names(ILRs_df) != "-", ]

    #limit ILR to certain range: 
    ILRs_df[which(ILRs_df <= -ILRs_limit)] <- -ILRs_limit
    ILRs_df[which(ILRs_df >= ILRs_limit)] <- ILRs_limit

    if(rownames_type == 'rownames_type')
      branch_gene_ABCs <- subset(ABC_df, gene_short_name %in% row.names(ILRs_df))
    else if(rownames_type == 'ensemble_id')
      branch_gene_ABCs <- ABC_df[row.names(ILRs_df), ]

    #select genes for certain lineage based on the ABC score
    if(ABC_type == 'positive')
      ILRs_df <- ILRs_df[unique(as.character(subset(branch_gene_ABCs, 
          ABCs > 0)[, "gene_short_name"])), ]
    else if(ABC_type == 'negative')
      ILRs_df <- ILRs_df[unique(as.character(subset(branch_gene_ABCs, 
        ABCs < 0)[, "gene_short_name"])), ]
    else if(ABC_type == 'all')
       ILRs_df <- ILRs_df

    ph <- pheatmap(ILRs_df, cluster_cols = FALSE, clustering_distance_rows = dist_method, 
        clustering_method = hclust_method, 
        ...)

    #create annotations for each gene
    annotation <- data.frame(class = as.factor(cutree(ph$tree_row, 
        cluster_num)), row.names = names(cutree(ph$tree_row, cluster_num)))

    #add also -log10(qval)
    gene_names <- row.names(ILRs_df[ph$tree_row$order, ]) 
    if(rownames_type == 'gene_short_name'){
      ensemble_names <- row.names(subset(fData(cds), gene_short_name %in% gene_names))
      #remove duplication: 
      ensemble_names <- ensemble_names[!duplicated(fData(cds[ensemble_names, ])$gene_short_name)]

      log_qval <- log10(branching_pval_df[ensemble_names, 'qval'])      
    }
    else{
      log_qval <- log10(branching_pval_df[gene_names, 'qval'])
    }
    annotation$log_qval <- -log_qval
    annotation[, "log10(abs(ABCs))"] <- log10(abs(ABC_df[row.names(annotation), 'ABCs']))

    #rotate the plot so that maturation level is on the x-axis
    pheatmap(t(ILRs_df[ph$tree_row$order, ]), 
      cluster_cols = T, 
      cluster_rows = F, 
      show_rownames = F, 
      show_colnames = F, 
      border_color = NA, 
      clustering_distance_cols = dist_method, 
      clustering_method = hclust_method, 
      annotation = annotation, 
      annotation_legend = T, ...)
}

#The following code is swipped from colorRamps package which is used to make the pallette
table.ramp <- function(n, mid = 0.5, sill = 0.5, base = 1, height = 1)
{
    x <- seq(0, 1, length.out = n)
    y <- rep(0, length(x))
    sill.min <- max(c(1, round((n - 1) * (mid - sill / 2)) + 1))
    sill.max <- min(c(n, round((n - 1) * (mid + sill / 2)) + 1))
    y[sill.min:sill.max] <- 1
    base.min <- round((n - 1) * (mid - base / 2)) + 1
    base.max <- round((n - 1) * (mid + base / 2)) + 1
    xi <- base.min:sill.min
    yi <- seq(0, 1, length.out = length(xi))
    i <- which(xi > 0 & xi <= n)
    y[xi[i]] <- yi[i]
    xi <- sill.max:base.max
    yi <- seq(1, 0, length.out = length(xi))
    i <- which(xi > 0 & xi <= n)
    y[xi[i]] <- yi[i]
    height * y
}

rgb.tables <- function(n,
red = c(0.75, 0.25, 1),
green = c(0.5, 0.25, 1),
blue = c(0.25, 0.25, 1))
{
    rr <- do.call("table.ramp", as.list(c(n, red)))
    gr <- do.call("table.ramp", as.list(c(n, green)))
    br <- do.call("table.ramp", as.list(c(n, blue)))
    rgb(rr, gr, br)
}

matlab.like <- function(n) rgb.tables(n)

matlab.like2 <- function(n)
rgb.tables(n,
red = c(0.8, 0.2, 1),
green = c(0.5, 0.4, 0.8),
blue = c(0.2, 0.2, 1))

blue2green2red <- matlab.like2

#'  Create a heatmap to demonstrate the bifurcation of gene expression along two lineages
#'
#' @param cds_subset CellDataSet for the experiment (normally only the branching genes detected with branchTest)
#' @param num_clusters Number of clusters for the heatmap of branch genes
#' @param ABC_df Matrix of Area Between Curves for the branching genes calculated using the calABCs function
#' @param branchTest_df Matrix of branchTest result (normally only on the branching genes)
#' @param lineage_labels The label for the lineages (first (second) lineage corresponding to state 2 (3))
#' @param norm_method Either vstExprs or log, determine which normalization approach will be performed for the fitted gene expression values
#' @param dist_method The method to calculate distance for each gene used in the hirearchical clustering, any one of them: "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski". (This option is not in function for now)
#' @param hclust_method The method to perform hirearchical clustering, any one of them: ward", "single", "complete", "average", "mcquitty", "median", "centroid"
#' @param heatmap_height Height of the saved heatmap
#' @param heatmap_width Width of the saved heatmap
#' @param ABC_lowest_thrsd The minimum log10(abs(ABCs)) used to annotate the heatmap
#' @param ABC_highest_thrsd The maximum log10(abs(ABCs)) used to annotate the heatmap
#' @param qval_lowest_thrsd The minimum log10(qval) used to annotate the heatmap
#' @param qval_highest_thrsd The maximum log10(qval) used to annotate the heatmap
#' @param hmcols The color scheme for drawing the heatmap
#' @param cores Number of cores to run this function
#' @param Cell_type_color The color for the progenitors and two lineages
#' @param return_all 
#' @return A list of heatmap_matrix (expression matrix for the lineage committment), ph (pheatmap heatmap object),
#' annotation_row (annotation data.frame for the row), annotation_col (annotation data.frame for the column). 
#' Note that, in order to draw the heatmap generate by this function you need to use grid.newpage(); grid.draw(res$ph$gt) (assuming "res" is the variable name for the result of this function)
#' @import pheatmap
#' @export
#'
plot_genes_branched_heatmap <- function(cds_subset, 
  num_clusters = 6,
  ABC_df = NULL, 
  branchTest_df = NULL, 
  lineage_states=c(2,3),
  branch_point=NULL,
  lineage_labels = c("Cell fate 1", "Cell fate 2"), 
  stretch = T, 
  scaling = T,
  norm_method = "vstExprs", 
  use_fitting_curves = T, 
  dist_method = NULL, 
  hclust_method = "ward.D2", 
  heatmap_height = 3, 
  heatmap_width = 4,
  ABC_lowest_thrsd = 0, 
  ABC_highest_thrsd = 2,
  qval_lowest_thrsd = 1, 
  qval_highest_thrsd = 5,
  hmcols = NULL, 
  Cell_type_color = c('#979797', '#F05662', '#7990C8'), 
  trend_formula = '~sm.ns(Pseudotime, df=3) * Lineage',
  pseudo_cnt = 0, 
  add_annotation_row = NULL,
  add_annotation_col = NULL,
  show_rownames = F, 
  cores = 1,
  use_gene_short_name = F,
  # file_name = 'branched_heatmap.pdf', 
  return_all=FALSE) {
    
    new_cds <- buildLineageBranchCellDataSet(cds_subset, 
                                             lineage_states=lineage_states, 
                                             branch_point=branch_point,
                                             stretch = stretch)
    new_cds@dispFitInfo <- cds_subset@dispFitInfo

    if(use_fitting_curves) {
        col_gap_ind <- 101
        # newdataA <- data.frame(Pseudotime = seq(0, 100, length.out = 100))
        # newdataB <- data.frame(Pseudotime = seq(0, 100, length.out = 100))

        newdataA <- data.frame(Pseudotime = seq(0, 100,
            length.out = 100), Lineage = as.factor(unique(as.character(pData(new_cds)$Lineage))[1]))   
        newdataB <- data.frame(Pseudotime = seq(0, 100,
            length.out = 100), Lineage = as.factor(unique(as.character(pData(new_cds)$Lineage))[2]))

        LineageAB_exprs <- genSmoothCurves(new_cds[, ], cores=cores, trend_formula = trend_formula,  
                    relative_expr = T, pseudocount = 0, new_data = rbind(newdataA, newdataB), weights = pData(new_cds)$weight)

        LineageA_exprs <- LineageAB_exprs[, 1:100]
        LineageB_exprs <- LineageAB_exprs[, 101:200]

        #common_ancestor_cells <- row.names(pData(new_cds)[duplicated(pData(new_cds)$original_cell_id),])
        common_ancestor_cells <- row.names(pData(new_cds)[pData(new_cds)$State == setdiff(pData(new_cds)$State, lineage_states),])
        
#         LineageA_exprs <- monocle::responseMatrix(monocle::fitModel(new_cds[, pData(new_cds)$Lineage == 2],  
#                                               modelFormulaStr="~sm.ns(Pseudotime, df=3)", cores=cores), newdataA)
#         LineageB_exprs <- monocle::responseMatrix(monocle::fitModel(new_cds[, pData(new_cds)$Lineage == 3],  
#                                               modelFormulaStr="~sm.ns(Pseudotime, df=3)", cores=cores), newdataB)
        LineageP_num <- (100 - floor(max(pData(new_cds)[common_ancestor_cells, 'Pseudotime'])))
        LineageA_num <- floor(max(pData(new_cds)[common_ancestor_cells, 'Pseudotime']))
        LineageB_num <- LineageA_num
    }
    else {
        LineageA_exprs <- exprs(new_cds[, pData(new_cds)$Lineage == as.numeric(unique(as.character(pData(new_cds)$Lineage))[1])])[, sort(pData(new_cds[, pData(new_cds)$Lineage == as.numeric(unique(as.character(pData(new_cds)$Lineage))[1])])$Pseudotime, index.return = T)$ix]
        LineageB_exprs <- exprs(new_cds[, pData(new_cds)$Lineage == as.numeric(unique(as.character(pData(new_cds)$Lineage))[2])])[, sort(pData(new_cds[, pData(new_cds)$Lineage == as.numeric(unique(as.character(pData(new_cds)$Lineage))[2])])$Pseudotime, index.return = T)$ix]

        col_gap_ind <- sum(pData(new_cds)$Lineage == as.numeric(unique(as.character(pData(new_cds)$Lineage))[1])) + 1

        newdataA <- data.frame(Pseudotime = sort(pData(new_cds[, pData(new_cds)$Lineage == as.numeric(unique(as.character(pData(new_cds)$Lineage))[1])])$Pseudotime))
        newdataB <- data.frame(Pseudotime = sort(pData(new_cds[, pData(new_cds)$Lineage == as.numeric(unique(as.character(pData(new_cds)$Lineage))[2])])$Pseudotime))

        #change to half of the number of the progenitor cell 
        #common_ancestor_cells <- row.names(pData(new_cds)[duplicated(pData(new_cds)$original_cell_id),])
        common_ancestor_cells <- row.names(pData(cds_subset)[pData(cds_subset)$State == setdiff(pData(cds_subset)$State, lineage_states),])
        
        #account for the changes in the buildLineageBranchCellDataSet
        if(length(common_ancestor_cells) %% 2 == 0) {
          LineageP_num <- length(common_ancestor_cells) / 2
          LineageA_num <- sum(pData(cds_subset)$State == as.numeric(unique(as.character(pData(new_cds)$Lineage))[1]))
          LineageB_num <- sum(pData(cds_subset)$State == as.numeric(unique(as.character(pData(new_cds)$Lineage))[2])) + 1
        }
        else {
          LineageP_num <- (length(common_ancestor_cells) + 1) / 2
          LineageA_num <- sum(pData(cds_subset)$State == as.numeric(unique(as.character(pData(new_cds)$Lineage))[1]))
          LineageB_num <- sum(pData(cds_subset)$State == as.numeric(unique(as.character(pData(new_cds)$Lineage))[2]))        
        }
    }

    if(norm_method == 'vstExprs') {
        LineageA_exprs <- vstExprs(new_cds, expr_matrix=LineageA_exprs)
        LineageB_exprs <- vstExprs(new_cds, expr_matrix=LineageB_exprs)
    }     
    else if(norm_method == 'log') {
        LineageA_exprs <- log10(LineageA_exprs + 1)
        LineageB_exprs <- log10(LineageB_exprs + 1)
    }

    heatmap_matrix <- cbind(LineageA_exprs[, (col_gap_ind - 1):1], LineageB_exprs)
    
    if(scaling) {
        heatmap_matrix <- Matrix::t(scale(Matrix::t(heatmap_matrix))) 
    }
    
    heatmap_matrix[heatmap_matrix > 3] <- 3
    heatmap_matrix[heatmap_matrix < -3] <- -3    

    heatmap_matrix_ori <- heatmap_matrix
    heatmap_matrix <- heatmap_matrix[is.finite(heatmap_matrix[, 1]) & is.finite(heatmap_matrix[, col_gap_ind]), ] #remove the NA fitting failure genes for each lineage 

    row_dist <- as.dist((1 - cor(Matrix::t(heatmap_matrix)))/2)
    row_dist[is.na(row_dist)] <- 1

    exp_rng <- range(heatmap_matrix) #bks is based on the expression range
    bks <- seq(exp_rng[1] - 0.1, exp_rng[2] + 0.1, by=0.1)
    if(is.null(hmcols)) {
        hmcols <- blue2green2red(length(bks) - 1)
    }
    
    # prin  t(hmcols)
    ph <- pheatmap(heatmap_matrix, 
             useRaster = T,
             cluster_cols=FALSE, 
             cluster_rows=TRUE, 
             show_rownames=F, 
             show_colnames=F, 
             #scale="row",
             clustering_distance_rows=row_dist,
             clustering_method = hclust_method,
             cutree_rows=num_clusters,
             silent=TRUE,
             filename=NA,
             breaks=bks,
             color=hmcols
             #color=hmcols#,
             # filename="expression_pseudotime_pheatmap.pdf",
             )
    #save(heatmap_matrix, row_dist, num_clusters, hmcols, ph, branchTest_df, qval_lowest_thrsd, lineage_labels, LineageA_num, LineageP_num, LineageB_num, file = 'heatmap_matrix')

    annotation_row <- data.frame(Cluster=factor(cutree(ph$tree_row, num_clusters)))
    if(!is.null(branchTest_df)) {
      annotation_row[, "-log 10(qval)"] <- - log10(branchTest_df[row.names(annotation_row), 'qval'])
      annotation_row[which(annotation_row[, "-log 10(qval)"] < qval_lowest_thrsd), "-log 10(qval)"] <- qval_lowest_thrsd
      annotation_row[which(annotation_row[, "-log 10(qval)"] > qval_highest_thrsd), "-log 10(qval)"] <- qval_highest_thrsd
    }
    if(!is.null(add_annotation_row)) {
        annotation_row <- cbind(annotation_row, add_annotation_row[row.names(annotation_row), ])  
        annotation_row <- add_annotation_row 
        # annotation_row$bif_time <- add_annotation_row[as.character(fData(absolute_cds[row.names(annotation_row), ])$gene_short_name), 1]
    }

    if(!is.null(ABC_df)) {
        annotation_row[, "log10(abs(ABCs))"] <- log10(abs(ABC_df[row.names(annotation_row), 'ABCs']))
        annotation_row[which(annotation_row[, "log10(abs(ABCs))"] < ABC_lowest_thrsd), "log10(abs(ABCs))"] <- ABC_lowest_thrsd
        annotation_row[which(annotation_row[, "log10(abs(ABCs))"] > ABC_highest_thrsd), "log10(abs(ABCs))"] <- ABC_highest_thrsd
    }
    
    colnames(heatmap_matrix) <- c(1:ncol(heatmap_matrix))
    annotation_col <- data.frame(row.names = c(1:ncol(heatmap_matrix)), "Cell Type" = c(rep(lineage_labels[1], LineageA_num),
                                                rep("Progenitor",  2 * LineageP_num),
                                                rep(lineage_labels[2], LineageB_num)))

    colnames(annotation_col) <- "Cell Type"  
    
    if(!is.null(add_annotation_col)) {
        annotation_col <- cbind(annotation_col, add_annotation_col[fData(cds[row.names(annotation_col), ])$gene_short_name, 1])  
    }

    Cluster_color <- brewer.pal(length(unique(annotation_row$Cluster)),"Set1")
    names(Cluster_color) <- 1:length(unique(annotation_row$Cluster))
    names(Cell_type_color) <- c("Progenitor", lineage_labels[1], lineage_labels[2])

    annotation_colors=list("Cell Type"=Cell_type_color, 
                                'Cluster' = Cluster_color)

    names(annotation_colors$`Cell Type`) = c('Progenitor', lineage_labels)

    if(use_gene_short_name == T) { #use_gene_short_name
        row.names(heatmap_matrix) <- fData(cds_subset)[row.names(heatmap_matrix), 'gene_short_name']
        row.names(annotation_row) <- fData(cds_subset)[row.names(annotation_row), 'gene_short_name']
    }

    #print(annotation_row)
    # pdf(paste(elife_directory, 'AT2_branch_gene_str_norm_div_df_heatmap_cole.pdf', sep = ''))#, height = 4, width = 3)
    #save(heatmap_matrix, hmcols, annotation_row, annotation_col, annotation_colors, row_dist, hclust_method, num_clusters, col_gap_ind, file = 'heatmap_matrix')

    #dev.off()
    #pdf(file_name, height = heatmap_height, width = heatmap_width)
    ph_res <- pheatmap(heatmap_matrix[, ], #ph$tree_row$order
             useRaster = T,
             cluster_cols=FALSE, 
             cluster_rows=TRUE, 
             show_rownames=show_rownames, 
             show_colnames=F, 
             #scale="row",
             clustering_distance_rows=row_dist, #row_dist
             clustering_method = hclust_method, #ward.D2
             cutree_rows=num_clusters,
             # cutree_cols = 2,
             annotation_row=annotation_row,
             annotation_col=annotation_col,
             annotation_colors=annotation_colors,
             gaps_col = col_gap_ind,
             treeheight_row = 1.5, 
             breaks=bks,
             fontsize = 6,
             color=hmcols, 
             silent=TRUE,
             # filename=NA
             # filename="expression_pseudotime_pheatmap2.pdf",
             )

    grid::grid.rect(gp=grid::gpar("fill"))
    grid::grid.draw(ph_res$gtable)
    if (return_all){
      return(list(LineageA_exprs = LineageA_exprs, LineageB_exprs = LineageB_exprs, heatmap_matrix = heatmap_matrix, 
        heatmap_matrix_ori = heatmap_matrix_ori, ph = ph, col_gap_ind = col_gap_ind, row_dist = row_dist, hmcols = hmcols, 
        annotation_colors = annotation_colors, annotation_row = annotation_row, annotation_col = annotation_col, 
        ph_res = ph_res))
    }
}

#' Plots genes by mean vs. dispersion, highlighting those selected for ordering
#' @export
plot_ordering_genes <- function(cds){
  disp_table <- dispersionTable(cds)

  ordering_genes <- row.names(subset(fData(cds), use_for_ordering))
  
  g <- qplot(mean_expression, dispersion_empirical, data=disp_table, log="xy", color=I("darkgrey")) + 
    geom_line(aes(y=dispersion_fit), color="red") 
  if (length(ordering_genes) > 0){
    g <- g + geom_point(aes(mean_expression, dispersion_empirical), 
                        data=disp_table[ordering_genes,], color="black")
  }
  g <- g + monocle_theme_opts()
  g
}

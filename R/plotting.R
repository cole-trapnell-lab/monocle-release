utils::globalVariables(c("Pseudotime", "value", "ids", "prin_graph_dim_1", "prin_graph_dim_2", "State", 
                         "value", "feature_label", "expectation", "colInd", "rowInd", "value", 
                         "source_prin_graph_dim_1", "source_prin_graph_dim_2"))

monocle_theme_opts <- function()
{
    theme(strip.background = element_rect(colour = 'white', fill = 'white')) +
    theme(panel.border = element_blank()) +
    theme(axis.line.x = element_line(size=0.25, color="black")) +
    theme(axis.line.y = element_line(size=0.25, color="black")) +
    theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) +
    theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) + 
    theme(panel.background = element_rect(fill='white')) +
    theme(legend.key=element_blank())
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
#' @param use_color_gradient Whether or not to use color gradient instead of cell size to show marker expression level 
#' @param markers_linear a boolean used to indicate whether you want to scale the markers logarithimically or linearly
#' @param show_cell_names draw the name of each cell in the plot
#' @param show_state_number show state number
#' @param cell_size The size of the point for each cell
#' @param cell_link_size The size of the line segments connecting cells (when used with ICA) or the principal graph (when used with DDRTree)
#' @param cell_name_size the size of cell name labels
#' @param state_number_size the size of the state number
#' @param show_branch_points Whether to show icons for each branch point (only available when reduceDimension was called with DDRTree)
#' @param theta How many degrees you want to rotate the trajectory
#' @param ... Additional arguments passed into scale_color_viridis function 
#' @return a ggplot2 plot object
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom igraph get.edgelist
#' @importFrom tibble rownames_to_column
#' @importFrom viridis scale_color_viridis
#' @importFrom dplyr left_join mutate n slice
#' @export
#' @examples
#' \dontrun{
#' lung <- load_lung()
#' plot_cell_trajectory(lung)
#' plot_cell_trajectory(lung, color_by="Pseudotime", show_backbone=FALSE)
#' plot_cell_trajectory(lung, markers="MYH3")
#' }
plot_cell_trajectory <- function(cds, 
                                 x=1, 
                                 y=2, 
                                 color_by="State", 
                                 show_tree=TRUE, 
                                 show_backbone=TRUE, 
                                 backbone_color="black", 
                                 markers=NULL, 
                                 use_color_gradient = FALSE,
                                 markers_linear = FALSE,
                                 show_cell_names=FALSE,
                                 show_state_number = FALSE,
                                 cell_size=1.5,
                                 cell_link_size=0.75,
                                 cell_name_size=2,
                                 state_number_size = 2.9,
                                 show_branch_points=TRUE,
                                 theta = 0,
                                 ...) {
  requireNamespace("igraph")
  gene_short_name <- NA
  sample_name <- NA
  sample_state <- pData(cds)$State
  data_dim_1 <- NA
  data_dim_2 <- NA
  
  #TODO: need to validate cds as ready for this plot (need mst, pseudotime, etc)
  lib_info_with_pseudo <- pData(cds)
  
  if (is.null(cds@dim_reduce_type)){
    stop("Error: dimensionality not yet reduced. Please call reduceDimension() before calling this function.")
  }
  
  if (cds@dim_reduce_type == "ICA"){
    reduced_dim_coords <- reducedDimS(cds)
  } else if (cds@dim_reduce_type %in% c("simplePPT", "DDRTree") ){
    reduced_dim_coords <- reducedDimK(cds)
  } else {
    stop("Error: unrecognized dimensionality reduction method.")
  }
  
  ica_space_df <- Matrix::t(reduced_dim_coords) %>%
    as.data.frame() %>%
    select_(prin_graph_dim_1 = x, prin_graph_dim_2 = y) %>%
    mutate(sample_name = rownames(.), sample_state = rownames(.))
  
  dp_mst <- minSpanningTree(cds)
  
  if (is.null(dp_mst)){
    stop("You must first call orderCells() before using this function")
  }
  
  edge_df <- dp_mst %>%
    igraph::as_data_frame() %>%
    select_(source = "from", target = "to") %>%
    left_join(ica_space_df %>% select_(source="sample_name", source_prin_graph_dim_1="prin_graph_dim_1", source_prin_graph_dim_2="prin_graph_dim_2"), by = "source") %>%
    left_join(ica_space_df %>% select_(target="sample_name", target_prin_graph_dim_1="prin_graph_dim_1", target_prin_graph_dim_2="prin_graph_dim_2"), by = "target")
  
  data_df <- t(monocle::reducedDimS(cds)) %>%
    as.data.frame() %>%
    select_(data_dim_1 = x, data_dim_2 = y) %>%
    rownames_to_column("sample_name") %>%
    mutate(sample_state) %>%
    left_join(lib_info_with_pseudo %>% rownames_to_column("sample_name"), by = "sample_name")
  
  return_rotation_mat <- function(theta) {
    theta <- theta / 180 * pi
    matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), nrow = 2)
  }
  rot_mat <- return_rotation_mat(theta)
  
  cn1 <- c("data_dim_1", "data_dim_2")
  cn2 <- c("source_prin_graph_dim_1", "source_prin_graph_dim_2")
  cn3 <- c("target_prin_graph_dim_1", "target_prin_graph_dim_2")
  data_df[, cn1] <- as.matrix(data_df[, cn1]) %*% t(rot_mat)
  edge_df[, cn2] <- as.matrix(edge_df[, cn2]) %*% t(rot_mat)
  edge_df[, cn3] <- as.matrix(edge_df[, cn3]) %*% t(rot_mat)
  
  markers_exprs <- NULL
  if (is.null(markers) == FALSE) {
    markers_fData <- subset(fData(cds), gene_short_name %in% markers)
    if (nrow(markers_fData) >= 1) {
      markers_exprs <- reshape2::melt(as.matrix(exprs(cds[row.names(markers_fData),])))
      colnames(markers_exprs)[1:2] <- c('feature_id','cell_id')
      markers_exprs <- merge(markers_exprs, markers_fData, by.x = "feature_id", by.y="row.names")
      #print (head( markers_exprs[is.na(markers_exprs$gene_short_name) == FALSE,]))
      markers_exprs$feature_label <- as.character(markers_exprs$gene_short_name)
      markers_exprs$feature_label[is.na(markers_exprs$feature_label)] <- markers_exprs$Var1
    }
  }
  if (is.null(markers_exprs) == FALSE && nrow(markers_exprs) > 0){
    data_df <- merge(data_df, markers_exprs, by.x="sample_name", by.y="cell_id")
    if(use_color_gradient) {
      if(markers_linear){
        g <- ggplot(data=data_df, aes(x=data_dim_1, y=data_dim_2)) + geom_point(aes(color= value), size=I(cell_size), na.rm = TRUE) + 
          scale_color_viridis(name = paste0("value"), ...) + facet_wrap(~feature_label)
        } else {
          g <- ggplot(data=data_df, aes(x=data_dim_1, y=data_dim_2)) + geom_point(aes(color=log10(value + 0.1)), size=I(cell_size), na.rm = TRUE) + 
              scale_color_viridis(name = paste0("log10(value + 0.1)"), ...) + facet_wrap(~feature_label)
        }
    } else {
      if(markers_linear){
        g <- ggplot(data=data_df, aes(x=data_dim_1, y=data_dim_2, size= (value * 0.1))) + facet_wrap(~feature_label)
      } else {
        g <- ggplot(data=data_df, aes(x=data_dim_1, y=data_dim_2, size=log10(value + 0.1))) + facet_wrap(~feature_label)
      }
    }
  } else {
    g <- ggplot(data=data_df, aes(x=data_dim_1, y=data_dim_2)) 
  }
  if (show_tree){
    g <- g + geom_segment(aes_string(x="source_prin_graph_dim_1", y="source_prin_graph_dim_2", xend="target_prin_graph_dim_1", yend="target_prin_graph_dim_2"), size=cell_link_size, linetype="solid", na.rm=TRUE, data=edge_df)
  }
  
  # FIXME: setting size here overrides the marker expression funtionality. 
  # Don't do it!
  if (is.null(markers_exprs) == FALSE && nrow(markers_exprs) > 0){
    if(use_color_gradient) {
      # g <- g + geom_point(aes_string(color = color_by), na.rm = TRUE)
    } else {
      g <- g + geom_point(aes_string(color = color_by), na.rm = TRUE)
    }
  }else {
    if(use_color_gradient) {
      # g <- g + geom_point(aes_string(color = color_by), na.rm = TRUE)
    } else {
      g <- g + geom_point(aes_string(color = color_by), size=I(cell_size), na.rm = TRUE)
    }
  }
  
  
  if (show_branch_points && cds@dim_reduce_type == 'DDRTree'){
    mst_branch_nodes <- cds@auxOrderingData[[cds@dim_reduce_type]]$branch_points
    branch_point_df <- ica_space_df %>%
      slice(match(mst_branch_nodes, sample_name)) %>%
      mutate(branch_point_idx = seq_len(n()))
    
    g <- g +
      geom_point(aes_string(x="prin_graph_dim_1", y="prin_graph_dim_2"),
                 size=5, na.rm=TRUE, branch_point_df) +
      geom_text(aes_string(x="prin_graph_dim_1", y="prin_graph_dim_2", label="branch_point_idx"),
                size=4, color="white", na.rm=TRUE, branch_point_df)
  }
  if (show_cell_names){
    g <- g + geom_text(aes(label=sample_name), size=cell_name_size)
  }
  if (show_state_number){
    g <- g + geom_text(aes(label = sample_state), size = state_number_size)
  }
  
  g <- g + 
    #scale_color_brewer(palette="Set1") +
    monocle_theme_opts() + 
    xlab(paste("Component", x)) + 
    ylab(paste("Component", y)) +
    theme(legend.position="top", legend.key.height=grid::unit(0.35, "in")) +
    #guides(color = guide_legend(label.position = "top")) +
    theme(legend.key = element_blank()) +
    theme(panel.background = element_rect(fill='white'))
  g
}

#' @rdname package-deprecated
#' @title Plots the minimum spanning tree on cells.
#' This function is deprecated.
#' @description This function arranges all of the cells in the cds in a tree and
#' predicts their location based on their pseudotime value
#' @param cds CellDataSet for the experiment
#' @param x the column of reducedDimS(cds) to plot on the horizontal axis
#' @param y the column of reducedDimS(cds) to plot on the vertical axis
#' @param color_by the cell attribute (e.g. the column of pData(cds)) to map to each cell's color
#' @param show_tree whether to show the links between cells connected in the minimum spanning tree
#' @param show_backbone whether to show the diameter path of the MST used to order the cells
#' @param backbone_color the color used to render the backbone.
#' @param markers a gene name or gene id to use for setting the size of each cell in the plot
#' @param show_cell_names draw the name of each cell in the plot
#' @param cell_size The size of the point for each cell
#' @param cell_link_size The size of the line segments connecting cells (when used with ICA) or the principal graph (when used with DDRTree)
#' @param cell_name_size the size of cell name labels
#' @param show_branch_points Whether to show icons for each branch point (only available when reduceDimension was called with DDRTree)
#' @return a ggplot2 plot object
#' @import ggplot2
#' @importFrom reshape2 melt
#' @export
#' @seealso plot_cell_trajectory
#' @examples
#' \dontrun{
#' library(HSMMSingleCell)
#' HSMM <- load_HSMM()
#' plot_cell_trajectory(HSMM)
#' plot_cell_trajectory(HSMM, color_by="Pseudotime", show_backbone=FALSE)
#' plot_cell_trajectory(HSMM, markers="MYH3")
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
  .Deprecated("plot_cell_trajectory") #include a package argument, too
  plot_cell_trajectory(cds=cds, 
                       x=x, 
                       y=y, 
                       color_by=color_by, 
                       show_tree=show_tree, 
                       show_backbone=show_backbone, 
                       backbone_color=backbone_color, 
                       markers=markers, 
                       show_cell_names=show_cell_names, 
                       cell_size=cell_size,
                       cell_link_size=cell_link_size,
                       cell_name_size=cell_name_size,
                       show_branch_points=show_branch_points)
}


#' @title Plots expression for one or more genes as a violin plot
#' 
#' @description Accepts a subset of a CellDataSet and an attribute to group cells by,
#' and produces one or more ggplot2 objects that plots the level of expression for
#' each group of cells. 
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
#' @param log_scale a boolean that determines whether or not to scale data logarithmically
#' @return a ggplot2 plot object
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom BiocGenerics sizeFactors
#' @export
#' @examples
#' \dontrun{
#' library(HSMMSingleCell)
#' HSMM <- load_HSMM()
#' my_genes <- HSMM[row.names(subset(fData(HSMM), gene_short_name %in% c("ACTA1", "ID1", "CCNB2"))),]
#' plot_genes_violin(my_genes, grouping="Hours", ncol=2, min_expr=0.1)
#' }
plot_genes_violin <- function (cds_subset, grouping = "State", min_expr = NULL, cell_size = 0.75, 
                              nrow = NULL, ncol = 1, panel_order = NULL, color_by = NULL, 
                              plot_trend = FALSE, label_by_short_name = TRUE, relative_expr = TRUE, 
                              log_scale = TRUE) 
{
  if (cds_subset@expressionFamily@vfamily %in% c("negbinomial", 
                                                 "negbinomial.size")) {
    integer_expression = TRUE
  }
  else {
    integer_expression = FALSE
    relative_expr = TRUE
  }
  if (integer_expression) {
    cds_exprs = exprs(cds_subset)
    if (relative_expr) {
      if (is.null(sizeFactors(cds_subset))) {
        stop("Error: to call this function with relative_expr=TRUE, you must call estimateSizeFactors() first")
      }
      cds_exprs = Matrix::t(Matrix::t(cds_exprs)/sizeFactors(cds_subset))
    }
    #cds_exprs = reshape2::melt(round(as.matrix(cds_exprs)))
    cds_exprs = reshape2::melt(as.matrix(cds_exprs))
  }
  else {
    cds_exprs = exprs(cds_subset)
    cds_exprs = reshape2::melt(as.matrix(cds_exprs))
  }
  if (is.null(min_expr)) {
    min_expr = cds_subset@lowerDetectionLimit
  }
  colnames(cds_exprs) = c("f_id", "Cell", "expression")
  cds_exprs$expression[cds_exprs$expression < min_expr] = min_expr
  cds_pData = pData(cds_subset)
  
  # 
  # # Custom bit for adding in a group for 
  # if(! is.null(show_combined)) {
  #   for(combine_gene in show_combined) {
  #     cds_pData_all = subset(cds_pData, gene == combine_gene)
  #     cds_pData_all[, grouping] = paste("All", combine_gene)
  #     cds_pData = rbind(cds_pData, cds_pData_all)
  #   }
  # }
  
  cds_fData = fData(cds_subset)
  cds_exprs = merge(cds_exprs, cds_fData, by.x = "f_id", by.y = "row.names")
  cds_exprs = merge(cds_exprs, cds_pData, by.x = "Cell", by.y = "row.names")
  cds_exprs$adjusted_expression = log10(cds_exprs$expression)
  
  
  
  
  if (label_by_short_name == TRUE) {
    if (is.null(cds_exprs$gene_short_name) == FALSE) {
      cds_exprs$feature_label = cds_exprs$gene_short_name
      cds_exprs$feature_label[is.na(cds_exprs$feature_label)] = cds_exprs$f_id
    }
    else {
      cds_exprs$feature_label = cds_exprs$f_id
    }
  }
  else {
    cds_exprs$feature_label = cds_exprs$f_id
  }
  if (is.null(panel_order) == FALSE) {
    cds_exprs$feature_label = factor(cds_exprs$feature_label, 
                                     levels = panel_order)
  }
  q = ggplot(aes_string(x = grouping, y = "expression"), data = cds_exprs)
  if (is.null(color_by) == FALSE) {
    q = q + geom_violin(aes_string(fill = color_by))
  }
  else {
    q = q + geom_violin()
  }
  if (plot_trend == TRUE) {
    q = q + stat_summary(fun.data = "mean_cl_boot", 
                         size = 0.2)
    q = q + stat_summary(aes_string(x = grouping, y = "expression", 
                                    group = color_by), fun.data = "mean_cl_boot", 
                         size = 0.2, geom = "line")
  }
  q = q + facet_wrap(~feature_label, nrow = nrow, 
                     ncol = ncol, scales = "free_y")
  if (min_expr < 1) {
     q = q + expand_limits(y = c(min_expr, 1))
  }
  
  
  q = q + ylab("Expression") + xlab(grouping)
  
  if (log_scale == TRUE){
    
    q = q + scale_y_log10()
  }
  q
}


#' Plots expression for one or more genes as a jittered, grouped points
#' 
#' @description Accepts a subset of a CellDataSet and an attribute to group cells by,
#' and produces one or more ggplot2 objects that plots the level of expression for
#' each group of cells. 
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
#' @importFrom BiocGenerics sizeFactors
#' @export
#' @examples
#' \dontrun{
#' library(HSMMSingleCell)
#' HSMM <- load_HSMM()
#' my_genes <- HSMM[row.names(subset(fData(HSMM), gene_short_name %in% c("MYOG", "ID1", "CCNB2"))),]
#' plot_genes_jitter(my_genes, grouping="Media", ncol=2)
#' }
plot_genes_jitter <- function(cds_subset, 
                              grouping = "State", 
                              min_expr=NULL, 
                              cell_size=0.75, 
                              nrow=NULL, 
                              ncol=1, 
                              panel_order=NULL, 
                              color_by=NULL,
                              plot_trend=FALSE,
                              label_by_short_name=TRUE,
                              relative_expr=TRUE){
  
  if (cds_subset@expressionFamily@vfamily %in% c("negbinomial", "negbinomial.size")){

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
    cds_exprs$feature_label <- factor(cds_exprs$feature_label, levels=panel_order)
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
  
  # Need this to guard against plotting failures caused by non-expressed genes
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
#'  @description Accetps a CellDataSet and a parameter,"grouping", used for dividing cells into groups.
#'  Returns one or more bar graphs (one graph for each gene in the CellDataSet).
#'  Each graph shows the percentage of cells that express a gene in the in the CellDataSet for
#'  each sub-group of cells created by "grouping".
#'  
#'  Let's say the CellDataSet passed in included genes A, B, and C and the "grouping parameter divided
#'  all of the cells into three groups called X, Y, and Z. Then three graphs would be produced called A,
#'  B, and C. In the A graph there would be three bars one for X, one for Y, and one for Z. So X bar in the
#'  A graph would show the percentage of cells in the X group that express gene A.
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
#' @param plot_limits A pair of number specifying the limits of the y axis. If NULL, scale to the range of the data.
#' @return a ggplot2 plot object
#' @import ggplot2
#' @importFrom plyr ddply
#' @importFrom reshape2 melt
#' @importFrom BiocGenerics sizeFactors
#' @export
#' @examples
#' \dontrun{
#' library(HSMMSingleCell)
#' HSMM <- load_HSMM()
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
                                      relative_expr=TRUE,
                                      plot_limits=c(0,100)){
  
  percent <- NULL

  if (cds_subset@expressionFamily@vfamily %in% c("negbinomial", "negbinomial.size")){
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
    marker_exprs_melted <- reshape2::melt(as.matrix(exprs(cds_subset)))
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

  marker_counts <- plyr::ddply(marker_exprs_melted, c("feature_label", grouping), function(x) { 
    data.frame(target=sum(x$expression > min_expr), 
               target_fraction=sum(x$expression > min_expr)/nrow(x)) } )
  
  #print (head(marker_counts))
  if (plot_as_fraction){
    marker_counts$target_fraction <- marker_counts$target_fraction * 100
    qp <- ggplot(aes_string(x=grouping, y="target_fraction", fill=grouping), data=marker_counts) +
      ylab("Cells (percent)")
    if (is.null(plot_limits) == FALSE)
      qp <- qp + scale_y_continuous(limits=plot_limits) 
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
#' @description Plots expression for one or more genes as a function of pseudotime.
#' Plotting allows you determine if the ordering produced by orderCells() is correct
#' and it does not need to be flipped using the "reverse" flag in orderCells
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
#' @param vertical_jitter A value passed to ggplot to jitter the points in the vertical dimension. Prevents overplotting, and is particularly helpful for rounded transcript count data.
#' @param horizontal_jitter A value passed to ggplot to jitter the points in the horizontal dimension. Prevents overplotting, and is particularly helpful for rounded transcript count data.
#' @return a ggplot2 plot object
#' @import ggplot2
#' @importFrom plyr ddply .
#' @importFrom reshape2 melt
#' @importFrom ggplot2 Position
#' @export
#' @examples
#' \dontrun{
#' library(HSMMSingleCell)
#' HSMM <- load_HSMM()
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
                                    vertical_jitter=NULL,
                                    horizontal_jitter=NULL){
    
  f_id <- NA
  Cell <- NA
    if (cds_subset@expressionFamily@vfamily %in% c("negbinomial", "negbinomial.size")) {
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
    #cds_exprs$f_id <- as.character(cds_exprs$f_id)
    #cds_exprs$Cell <- as.character(cds_exprs$Cell)
    
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
    cds_exprs$f_id <- as.character(cds_exprs$f_id)
    cds_exprs$feature_label <- factor(cds_exprs$feature_label)

    new_data <- data.frame(Pseudotime = pData(cds_subset)$Pseudotime)
    model_expectation <- genSmoothCurves(cds_subset, cores=1, trend_formula = trend_formula,
                        relative_expr = T, new_data = new_data)
    colnames(model_expectation) <- colnames(cds_subset)
    expectation <- ddply(cds_exprs, .(f_id, Cell), function(x) data.frame("expectation"=model_expectation[x$f_id, x$Cell]))
    cds_exprs <- merge(cds_exprs, expectation)
    #cds_exprs$expectation <- expectation#apply(cds_exprs,1, function(x) model_expectation[x$f_id, x$Cell])

    cds_exprs$expression[cds_exprs$expression < min_expr] <- min_expr
    cds_exprs$expectation[cds_exprs$expectation < min_expr] <- min_expr
    if (is.null(panel_order) == FALSE) {
      cds_exprs$feature_label <- factor(cds_exprs$feature_label,
            levels = panel_order)
    }
    q <- ggplot(aes(Pseudotime, expression), data = cds_exprs)
    if (is.null(color_by) == FALSE) {
        q <- q + geom_point(aes_string(color = color_by), size = I(cell_size), position=position_jitter(horizontal_jitter, vertical_jitter))
    }
    else {
        q <- q + geom_point(size = I(cell_size), position=position_jitter(horizontal_jitter, vertical_jitter))
    }

    q <- q + geom_line(aes(x = Pseudotime, y = expectation), data = cds_exprs)

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

#' Plots kinetic clusters of genes.
#'
#' @description returns a ggplot2 object showing the shapes of the
#' expression patterns followed by a set of pre-selected genes.
#' The topographic lines highlight the distributions of the kinetic patterns
#' relative to overall trend lines.
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
#' @importFrom stringr str_c
#' @importFrom ggplot2 Position
#' @import grid
#' @export
#' @examples
#' \dontrun{
#' full_model_fits <- fitModel(HSMM_filtered[sample(nrow(fData(HSMM_filtered)), 100),],  
#'    modelFormulaStr="~VGAM::bs(Pseudotime)")
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
                        callout_ids=NULL){
  .Deprecated("plot_genes_heatmap")
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
  facet_labels <- str_c(cluster_sizes$Var1, cluster_sizes$Freq, sep=" ") #update the function

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
# 
# #' Plots a pseudotime-ordered, row-centered heatmap
# #' @export 
# plot_genes_heatmap <- function(cds, 
#                                rescaling='row', 
#                                clustering='row', 
#                                labCol=FALSE, 
#                                labRow=TRUE, 
#                                logMode=TRUE, 
#                                use_vst=TRUE,
#                                border=FALSE, 
#                                heatscale=c(low='steelblue',mid='white',high='tomato'), 
#                                heatMidpoint=0,
#                                method="none",
#                                scaleMax=2, 
#                                scaleMin=-2, 
#                                relative_expr=TRUE, 
#                                ...){
#   
#   ## the function can be be viewed as a two step process
#   ## 1. using the rehape package and other funcs the data is clustered, scaled, and reshaped
#   ## using simple options or by a user supplied function
#   ## 2. with the now resahped data the plot, the chosen labels and plot style are built
#   FM <- exprs(cds)
#   
#   if (cds@expressionFamily@vfamily %in% c("negbinomial", "negbinomial.size")){
#     integer_expression <- TRUE
#   }else{
#     integer_expression <- FALSE
#     relative_expr <- TRUE
#   }
#   
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
#   
#   m=FM
#   
#   if (is.null(fData(cds)$gene_short_name) == FALSE){
#     feature_labels <- fData(cds)$gene_short_name
#     feature_labels[is.na(feature_labels)]  <- fData(cds)$f_id
#     row.names(m) <- feature_labels
#   }
#   
#   #remove genes with no expression in any condition
#   m=m[!apply(m,1,sum)==0,]
#   
#   if (use_vst && is.null(cds@dispFitInfo[["blind"]]$disp_func) == FALSE){
#     m = vstExprs(cds, expr_matrix=m)
#   }else if(logMode){
#     m = log10(m+pseudocount)
#   }
#   
#   #remove genes with no sd
#   #m=m[!apply(m,1,sd)==0,]
# 
#   ## you can either scale by row or column not both! 
#   ## if you wish to scale by both or use a different scale method then simply supply a scale
#   ## function instead NB scale is a base funct
#   
#   ## I have supplied the default cluster and euclidean distance (JSdist) - and chose to cluster after scaling
#   ## if you want a different distance/cluster method-- or to cluster and then scale
#   ## then you can supply a custom function 
#   
#   if(!is.function(method)){
#     method = function(mat){as.dist((1 - cor(Matrix::t(mat)))/2)}	
#   }
#   
#   ## this is just reshaping into a ggplot format matrix and making a ggplot layer
#   
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
#   
#   # If we aren't going to re-ordering the columns, order them by Pseudotime
#   if (clustering %in% c("row", "none"))
#     m = m[,row.names(pData(cds)[order(-pData(cds)$Pseudotime),])]
#   
#   if(clustering=='row')
#     m=m[hclust(method(m))$order, ]
#   if(clustering=='column')  
#     m=m[,hclust(method(Matrix::t(m)))$order]
#   if(clustering=='both')
#     m=m[hclust(method(m))$order ,hclust(method(Matrix::t(m)))$order]
#   
#   
#   rows=dim(m)[1]
#   cols=dim(m)[2]
#   
#   
#   
#   # if(logMode) {
#   #   melt.m=cbind(rowInd=rep(1:rows, times=cols), colInd=rep(1:cols, each=rows), reshape2::melt( log10(m+pseudocount)))
#   # }else{
#   #   melt.m=cbind(rowInd=rep(1:rows, times=cols), colInd=rep(1:cols, each=rows), reshape2::melt(m))
#   # }
#   
#   melt.m=cbind(rowInd=rep(1:rows, times=cols), colInd=rep(1:cols, each=rows), reshape2::melt(m))
#   
#   g=ggplot(data=melt.m)
#   
#   ## add the heat tiles with or without a white border for clarity
#   
#   if(border==TRUE)
#     g2=g+geom_raster(aes(x=colInd,y=rowInd, fill=value),colour='grey')
#   if(border==FALSE)
#     g2=g+geom_raster(aes(x=colInd,y=rowInd,ymax=rowInd, fill=value))
#   
#   ## add axis labels either supplied or from the colnames rownames of the matrix
#   
#   if(labCol==TRUE) 
#   {
#     g2=g2+scale_x_continuous(breaks=(1:cols)-0.5, labels=colnames(m))
#   }
#   if(labCol==FALSE) 
#   {
#     g2=g2+scale_x_continuous(breaks=(1:cols)-0.5, labels=rep('',cols))
#   }
#   
#   
#   if(labRow==TRUE) 
#   {
#     g2=g2+scale_y_continuous(breaks=(1:rows)-0.5, labels=rownames(m))	
#   }
#   if(labRow==FALSE)
#   { 
#     g2=g2+scale_y_continuous(breaks=(1:rows)-0.5, labels=rep('',rows))	
#   }
#   
#   # Get rid of the ticks, they get way too dense with lots of rows
#   g2 <- g2 + theme(axis.ticks = element_blank()) 
#   
#   ## get rid of grey panel background and gridlines
#   
#   g2=g2+theme(panel.grid.minor=element_line(colour=NA), panel.grid.major=element_line(colour=NA),
#               panel.background=element_rect(fill=NA, colour=NA))
#   
#   ##adjust x-axis labels
#   g2=g2+theme(axis.text.x=element_text(angle=-90, hjust=0))
#   
#   #write(paste(c("Length of heatscale is :", length(heatscale))), stderr())
#   
#   if(is.function(rescaling))
#   {
#     
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
#   
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
#   
#   #g2<-g2+scale_x_discrete("",breaks=tracking_ids,labels=gene_short_names)
#   
#   g2 <- g2 + theme(axis.title.x=element_blank(), axis.title.y=element_blank())
#   
#   ## finally add the fill colour ramp of your choice (default is blue to red)-- and return
#   return (g2)
# }


plot_genes_heatmap <- function(...){
  .Deprecated("plot_pseudotime_heatmap")
  plot_pseudotime_heatmap(...)
}

#' Plots a pseudotime-ordered, row-centered heatmap
#' 
#' @description The function plot_pseudotime_heatmap takes a CellDataSet object 
#' (usually containing a only subset of significant genes) and generates smooth expression 
#' curves much like plot_genes_in_pseudotime. 
#' Then, it clusters these genes and plots them using the pheatmap package. 
#' This allows you to visualize modules of genes that co-vary across pseudotime.
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
#' @importFrom stats sd as.dist cor cutree
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
                                    
                                    norm_method = c("log", "vstExprs"), 
                                    scale_max=3, 
                                    scale_min=-3, 
                                    
                                    trend_formula = '~sm.ns(Pseudotime, df=3)',
                                    
                                    return_heatmap=FALSE,
                                    cores=1){
  num_clusters <- min(num_clusters, nrow(cds_subset))
  pseudocount <- 1
  newdata <- data.frame(Pseudotime = seq(min(pData(cds_subset)$Pseudotime), max(pData(cds_subset)$Pseudotime),length.out = 100)) 
  
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
  
  if(is.null(hmcols)) {
    bks <- seq(-3.1,3.1, by = 0.1)
    hmcols <- blue2green2red(length(bks) - 1)
  }
  else {
    bks <- seq(-3.1,3.1, length.out = length(hmcols))
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
                 border_color = NA,
                 color=hmcols)

  if(cluster_rows) {
    annotation_row <- data.frame(Cluster=factor(cutree(ph$tree_row, num_clusters)))
  } else {
    annotation_row <- NULL
  }
  
  if(!is.null(add_annotation_row)) {
    old_colnames_length <- ncol(annotation_row)
    annotation_row <- cbind(annotation_row, add_annotation_row[row.names(annotation_row), ])  
    colnames(annotation_row)[(old_colnames_length+1):ncol(annotation_row)] <- colnames(add_annotation_row)
    # annotation_row$bif_time <- add_annotation_row[as.character(fData(absolute_cds[row.names(annotation_row), ])$gene_short_name), 1]
  }
  
  if(!is.null(add_annotation_col)) {
    if(nrow(add_annotation_col) != 100) {
      stop('add_annotation_col should have only 100 rows (check genSmoothCurves before you supply the annotation data)!')
    }
    annotation_col <- add_annotation_col
  } else {
    annotation_col <- NA
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
    if(!is.null(annotation_row))
      row_ann_labels <- row.names(annotation_row)
  }
  
  row.names(heatmap_matrix) <- feature_label
  if(!is.null(annotation_row))
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
                     annotation_col=annotation_col,
                     treeheight_row = 20, 
                     breaks=bks,
                     fontsize = 6,
                     color=hmcols, 
                     border_color = NA,
                     silent=TRUE,
                     filename=NA
  )
  
  grid::grid.rect(gp=grid::gpar("fill", col=NA))
  grid::grid.draw(ph_res$gtable)
  if (return_heatmap){
    return(ph_res)
  }
}


#' Plot the branch genes in pseduotime with separate branch curves.
#' 
#' @description Works similarly to plot_genes_in_psuedotime esceptit shows 
#' one kinetic trend for each lineage. 
#' 
#' @details This plotting function is used to make the branching plots for a branch dependent gene goes through the progenitor state
#' and bifurcating into two distinct branchs (Similar to the pitch-fork bifurcation in dynamic systems). In order to make the  
#' bifurcation plot, we first duplicated the progenitor states and by default stretch each branch into maturation level 0-100.  
#' Then we fit two nature spline curves for each branchs using VGAM package.  
#'
#' @param cds CellDataSet for the experiment
#' @param branch_states The states for two branching branchs
#' @param branch_point The ID of the branch point to analyze. Can only be used when reduceDimension is called with method = "DDRTree".
#' @param branch_labels The names for each branching branch
#' @param method The method to draw the curve for the gene expression branching pattern, either loess ('loess') or VGLM fitting ('fitting') 
#' @param min_expr The minimum (untransformed) expression level to use in plotted the genes.
#' @param cell_size The size (in points) of each cell used in the plot
#' @param nrow Number of columns used to layout the faceted cluster panels
#' @param ncol Number of columns used to layout the faceted cluster panels
#' @param panel_order The a character vector of gene short names (or IDs, if that's what you're using), specifying order in which genes should be layed out (left-to-right, top-to-bottom)
#' @param color_by The cell attribute (e.g. the column of pData(cds)) to be used to color each cell 
#' @param expression_curve_linetype_by The cell attribute (e.g. the column of pData(cds)) to be used for the linetype of each branch curve
#' @param trend_formula The model formula to be used for fitting the expression trend over pseudotime
#' @param reducedModelFormulaStr A formula specifying a null model. If used, the plot shows a p value from the likelihood ratio test that uses trend_formula as the full model
#' @param label_by_short_name Whether to label figure panels by gene_short_name (TRUE) or feature id (FALSE)
#' @param relative_expr Whether or not the plot should use relative expression values (only relevant for CellDataSets using transcript counts)
#' @param ... Additional arguments passed on to branchTest. Only used when reducedModelFormulaStr is not NULL.
#' @return a ggplot2 plot object
#' @import ggplot2
#' @importFrom plyr ddply
#' @importFrom reshape2 melt
#' @importFrom BiocGenerics sizeFactors
#' @export 
plot_genes_branched_pseudotime <- function (cds, 
                                            branch_states = NULL, 
                                            branch_point=1,
                                            branch_labels = NULL,
                                            method = "fitting", 
                                            min_expr = NULL, 
                                            cell_size = 0.75,
                                            nrow = NULL, 
                                            ncol = 1, 
                                            panel_order = NULL, 
                                            color_by = "State",
                                            expression_curve_linetype_by = "Branch", 
                                            trend_formula = "~ sm.ns(Pseudotime, df=3) * Branch", 
                                            reducedModelFormulaStr = NULL, 
                                            label_by_short_name = TRUE,
                                            relative_expr = TRUE,
                                            #gene_pairs = NULL,
                                            ...)
{
  Branch <- NA  
  if (is.null(reducedModelFormulaStr) == FALSE) {
    pval_df <- branchTest(cds, 
                          branch_states=branch_states,
                          branch_point=branch_point,
                          fullModelFormulaStr = trend_formula,
                          reducedModelFormulaStr = "~ sm.ns(Pseudotime, df=3)", 
                          ...)
    fData(cds)[, "pval"] <- pval_df[row.names(cds), 'pval']
  }
  if("Branch" %in% all.vars(terms(as.formula(trend_formula)))) { #only when Branch is in the model formula we will duplicate the "progenitor" cells
    cds_subset <- buildBranchCellDataSet(cds = cds, 
                                         branch_states = branch_states, 
                                         branch_point=branch_point,
                                         branch_labels = branch_labels, 
                                         progenitor_method = 'duplicate',
                                         ...)
  }
  else {
    cds_subset <- cds
    pData(cds_subset)$Branch <- pData(cds_subset)$State
  }
  if (cds_subset@expressionFamily@vfamily %in% c("negbinomial", "negbinomial.size")) {
    integer_expression <- TRUE
  }
  else {
    integer_expression <- FALSE
  }
  if (integer_expression) {
    CM <- exprs(cds_subset)
    if (relative_expr){
      if (is.null(sizeFactors(cds_subset))) {
        stop("Error: to call this function with relative_expr=TRUE, you must call estimateSizeFactors() first")
      }
      CM <- Matrix::t(Matrix::t(CM)/sizeFactors(cds_subset))
    }
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
  
  cds_fData <- fData(cds_subset)
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
  cds_exprs$Branch <- as.factor(cds_exprs$Branch) 
  
  new_data <- data.frame(Pseudotime = pData(cds_subset)$Pseudotime, Branch = pData(cds_subset)$Branch)
  
  full_model_expectation <- genSmoothCurves(cds_subset, cores=1, trend_formula = trend_formula, 
                                            relative_expr = T, new_data = new_data)
  colnames(full_model_expectation) <- colnames(cds_subset)
  
  cds_exprs$full_model_expectation <- apply(cds_exprs,1, function(x) full_model_expectation[x[2], x[1]])
  if(!is.null(reducedModelFormulaStr)){
    reduced_model_expectation <- genSmoothCurves(cds_subset, cores=1, trend_formula = reducedModelFormulaStr,
                                                 relative_expr = T, new_data = new_data)
    colnames(reduced_model_expectation) <- colnames(cds_subset)
    cds_exprs$reduced_model_expectation <- apply(cds_exprs,1, function(x) reduced_model_expectation[x[2], x[1]])
  }
  
  # FIXME: If you want to show the bifurcation time for each gene, this function
  # should just compute it. Passing it in as a dataframe is just too complicated
  # and will be hard on the user. 
  # if(!is.null(bifurcation_time)){
  #     cds_exprs$bifurcation_time <- bifurcation_time[as.vector(cds_exprs$gene_short_name)]
  # }
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
  cds_exprs$Branch <- as.factor(cds_exprs$Branch)
  
  q <- ggplot(aes(Pseudotime, expression), data = cds_exprs)
  # if (!is.null(bifurcation_time)) {
  #   q <- q + geom_vline(aes(xintercept = bifurcation_time),
  #                       color = "black", linetype = "longdash")
  # }
  if (is.null(color_by) == FALSE) {
    q <- q + geom_point(aes_string(color = color_by), size = I(cell_size))
  }
  if (is.null(reducedModelFormulaStr) == FALSE)
    q <- q + scale_y_log10() + facet_wrap(~feature_label +
                                            pval, nrow = nrow, ncol = ncol, scales = "free_y")
  else q <- q + scale_y_log10() + facet_wrap(~feature_label,
                                             nrow = nrow, ncol = ncol, scales = "free_y")
  if (method == "loess")
    q <- q + stat_smooth(aes(fill = Branch, color = Branch),
                         method = "loess")
  else if (method == "fitting") {
    q <- q + geom_line(aes_string(x = "Pseudotime", y = "full_model_expectation",
                                  linetype = "Branch"), data = cds_exprs) #+ scale_color_manual(name = "Type", values = c(colour_cell, colour), labels = c("Pre-branch", "AT1", "AT2", "AT1", "AT2")
  }
  
  if(!is.null(reducedModelFormulaStr)) {
    q <- q + geom_line(aes_string(x = "Pseudotime", y = "reduced_model_expectation"),
                       color = 'black', linetype = 2, data =  cds_exprs)   
  }
  
  q <- q + ylab("Expression") + xlab("Pseudotime (stretched)")
  
  q <- q + monocle_theme_opts()
  q + expand_limits(y = min_expr)
}

#' Not sure we're ready to release this one quite yet:
#' Plot the branch genes in pseduotime with separate branch curves 
#' @param cds CellDataSet for the experiment
#' @param rowgenes Gene ids or short names to be arrayed on the vertical axis.
#' @param colgenes Gene ids or short names to be arrayed on the horizontal axis
#' @param relative_expr Whether to transform expression into relative values
#' @param min_expr The minimum level of expression to show in the plot
#' @param cell_size A number how large the cells should be in the plot
#' @param label_by_short_name a boolean that indicates whether cells should be labeled by their short name
#' @param show_density a boolean that indicates whether a 2D density estimation should be shown in the plot
#' @param round_expr a boolean that indicates whether cds_expr values should be rounded or not
#' @return a ggplot2 plot object
#' @import ggplot2
#' @importFrom reshape2 melt
plot_coexpression_matrix <- function(cds, 
                                     rowgenes, 
                                     colgenes, 
                                     relative_expr=TRUE, 
                                     min_expr=NULL, 
                                     cell_size=0.85, 
                                     label_by_short_name=TRUE,
                                     show_density=TRUE,
                                     round_expr=FALSE){
  
  gene_short_name <- NA
  f_id <- NA
  adjusted_expression.x <- NULL
  adjusted_expression.y <- NULL
  ..density.. <- NULL
  
  
  row_gene_ids <- row.names(subset(fData(cds), gene_short_name %in% rowgenes))
  row_gene_ids <- union(row_gene_ids, intersect(rowgenes, row.names(fData(cds))))
  
  col_gene_ids <- row.names(subset(fData(cds), gene_short_name %in% colgenes))
  col_gene_ids <- union(col_gene_ids, intersect(colgenes, row.names(fData(cds))))
  
  cds_subset <- cds[union(row_gene_ids, col_gene_ids),]
  
  if (cds_subset@expressionFamily@vfamily %in% c("negbinomial", "negbinomial.size")){
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
      cds_exprs <- reshape2::melt(as.matrix(cds_exprs))
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

#' @importFrom grDevices rgb
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

#'  Create a heatmap to demonstrate the bifurcation of gene expression along two branchs
#'  
#'  @description returns a heatmap that shows changes in both lineages at the same time. 
#'  It also requires that you choose a branch point to inspect. 
#'  Columns are points in pseudotime, rows are genes, and the beginning of pseudotime is in the middle of the heatmap. 
#'  As you read from the middle of the heatmap to the right, you are following one lineage through pseudotime. As you read left, the other. 
#'  The genes are clustered hierarchically, so you can visualize modules of genes that have similar lineage-dependent expression patterns.
#'
#' @param cds_subset CellDataSet for the experiment (normally only the branching genes detected with branchTest)
#' @param branch_point The ID of the branch point to visualize. Can only be used when reduceDimension is called with method = "DDRTree".
#' @param branch_states The two states to compare in the heatmap. Mutually exclusive with branch_point. 
#' @param branch_labels The labels for the branchs. 
#' @param cluster_rows Whether to cluster the rows of the heatmap.
#' @param hclust_method The method used by pheatmap to perform hirearchical clustering of the rows. 
#' @param num_clusters Number of clusters for the heatmap of branch genes
#' @param hmcols The color scheme for drawing the heatmap.
#' @param branch_colors The colors used in the annotation strip indicating the pre- and post-branch cells.
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
#' @param ... Additional arguments passed to buildBranchCellDataSet
#' @return A list of heatmap_matrix (expression matrix for the branch committment), ph (pheatmap heatmap object),
#' annotation_row (annotation data.frame for the row), annotation_col (annotation data.frame for the column). 
#' @import pheatmap
#' @importFrom stats sd as.dist cor cutree
#' @export
#'
plot_genes_branched_heatmap <- function(cds_subset, 
                                        
                                        branch_point=1,
                                        branch_states=NULL,
                                        branch_labels = c("Cell fate 1", "Cell fate 2"), 
                                        cluster_rows = TRUE,
                                        hclust_method = "ward.D2", 
                                        num_clusters = 6,
                                        hmcols = NULL, 
                                        branch_colors = c('#979797', '#F05662', '#7990C8'), 
                                        add_annotation_row = NULL,
                                        add_annotation_col = NULL,
                                        show_rownames = FALSE, 
                                        use_gene_short_name = TRUE,
                                        scale_max=3, 
                                        scale_min=-3, 
                                        norm_method = c("log", "vstExprs"), 
                                        
                                        trend_formula = '~sm.ns(Pseudotime, df=3) * Branch',
                                        
                                        return_heatmap=FALSE,
                                        cores = 1, ...) {
  
  cds <- NA
  new_cds <- buildBranchCellDataSet(cds_subset, 
                                    branch_states=branch_states, 
                                    branch_point=branch_point, 
                                    progenitor_method = 'duplicate',
                                    ...)
  
  new_cds@dispFitInfo <- cds_subset@dispFitInfo
  
  if(is.null(branch_states)) {
    progenitor_state <- subset(pData(cds_subset), Pseudotime == 0)[, 'State']
    branch_states <- setdiff(pData(cds_subset)$State, progenitor_state)
  }
  
  col_gap_ind <- 101
  # newdataA <- data.frame(Pseudotime = seq(0, 100, length.out = 100))
  # newdataB <- data.frame(Pseudotime = seq(0, 100, length.out = 100))
  
  newdataA <- data.frame(Pseudotime = seq(0, 100,
                                          length.out = 100), Branch = as.factor(unique(as.character(pData(new_cds)$Branch))[1]))   
  newdataB <- data.frame(Pseudotime = seq(0, 100,
                                          length.out = 100), Branch = as.factor(unique(as.character(pData(new_cds)$Branch))[2]))
  
  BranchAB_exprs <- genSmoothCurves(new_cds[, ], cores=cores, trend_formula = trend_formula,  
                                    relative_expr = T, new_data = rbind(newdataA, newdataB))
  
  BranchA_exprs <- BranchAB_exprs[, 1:100]
  BranchB_exprs <- BranchAB_exprs[, 101:200]
  
  #common_ancestor_cells <- row.names(pData(new_cds)[duplicated(pData(new_cds)$original_cell_id),])
  common_ancestor_cells <- row.names(pData(new_cds)[pData(new_cds)$State == setdiff(pData(new_cds)$State, branch_states),])
  BranchP_num <- (100 - floor(max(pData(new_cds)[common_ancestor_cells, 'Pseudotime'])))
  BranchA_num <- floor(max(pData(new_cds)[common_ancestor_cells, 'Pseudotime']))
  BranchB_num <- BranchA_num
  
  norm_method <- match.arg(norm_method)
  
  # FIXME: this needs to check that vst values can even be computed. (They can only be if we're using NB as the expressionFamily)
  if(norm_method == 'vstExprs') {
    BranchA_exprs <- vstExprs(new_cds, expr_matrix=BranchA_exprs)
    BranchB_exprs <- vstExprs(new_cds, expr_matrix=BranchB_exprs)
  }     
  else if(norm_method == 'log') {
    BranchA_exprs <- log10(BranchA_exprs + 1)
    BranchB_exprs <- log10(BranchB_exprs + 1)
  }
  
  heatmap_matrix <- cbind(BranchA_exprs[, (col_gap_ind - 1):1], BranchB_exprs)
  
  heatmap_matrix=heatmap_matrix[!apply(heatmap_matrix, 1, sd)==0,]
  heatmap_matrix=Matrix::t(scale(Matrix::t(heatmap_matrix),center=TRUE))
  heatmap_matrix=heatmap_matrix[is.na(row.names(heatmap_matrix)) == FALSE,]
  heatmap_matrix[is.nan(heatmap_matrix)] = 0
  heatmap_matrix[heatmap_matrix>scale_max] = scale_max
  heatmap_matrix[heatmap_matrix<scale_min] = scale_min
  
  heatmap_matrix_ori <- heatmap_matrix
  heatmap_matrix <- heatmap_matrix[is.finite(heatmap_matrix[, 1]) & is.finite(heatmap_matrix[, col_gap_ind]), ] #remove the NA fitting failure genes for each branch 
  
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
  #save(heatmap_matrix, row_dist, num_clusters, hmcols, ph, branchTest_df, qval_lowest_thrsd, branch_labels, BranchA_num, BranchP_num, BranchB_num, file = 'heatmap_matrix')
  
  annotation_row <- data.frame(Cluster=factor(cutree(ph$tree_row, num_clusters)))
  
  if(!is.null(add_annotation_row)) {
    annotation_row <- cbind(annotation_row, add_annotation_row[row.names(annotation_row), ])  
    # annotation_row$bif_time <- add_annotation_row[as.character(fData(absolute_cds[row.names(annotation_row), ])$gene_short_name), 1]
  }
  
  colnames(heatmap_matrix) <- c(1:ncol(heatmap_matrix))
  annotation_col <- data.frame(row.names = c(1:ncol(heatmap_matrix)), "Cell Type" = c(rep(branch_labels[1], BranchA_num),
                                                                                      rep("Pre-branch",  2 * BranchP_num),
                                                                                      rep(branch_labels[2], BranchB_num)))
  
  colnames(annotation_col) <- "Cell Type"  
  
  if(!is.null(add_annotation_col)) {
    annotation_col <- cbind(annotation_col, add_annotation_col[fData(cds[row.names(annotation_col), ])$gene_short_name, 1])  
  }
  
  names(branch_colors) <- c("Pre-branch", branch_labels[1], branch_labels[2])
  
  annotation_colors=list("Cell Type"=branch_colors)
  
  names(annotation_colors$`Cell Type`) = c('Pre-branch', branch_labels)
  
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
                     treeheight_row = 20, 
                     breaks=bks,
                     fontsize = 6,
                     color=hmcols, 
                     border_color = NA,
                     silent=TRUE)
  
  grid::grid.rect(gp=grid::gpar("fill", col=NA))
  grid::grid.draw(ph_res$gtable)
  if (return_heatmap){
    return(list(BranchA_exprs = BranchA_exprs, BranchB_exprs = BranchB_exprs, heatmap_matrix = heatmap_matrix, 
                heatmap_matrix_ori = heatmap_matrix_ori, ph = ph, col_gap_ind = col_gap_ind, row_dist = row_dist, hmcols = hmcols, 
                annotation_colors = annotation_colors, annotation_row = annotation_row, annotation_col = annotation_col, 
                ph_res = ph_res))
  }
}

#' Plots genes by mean vs. dispersion, highlighting those selected for ordering
#' 
#' Each gray point in the plot is a gene. The black dots are those that were included
#' in the last call to setOrderingFilter. The red curve shows the mean-variance 
#' model learning by estimateDispersions().
#' 
#' @param cds The CellDataSet to be used for the plot.
#' @export
plot_ordering_genes <- function(cds){
  if(class(cds)[1] != "CellDataSet") {
    stop("Error input object is not of type 'CellDataSet'")
  }
  disp_table <- dispersionTable(cds)
  use_for_ordering <- NA
  mean_expression <- NA
  dispersion_empirical <- NA
  dispersion_fit <- NA
  gene_id <- NA
  ordering_genes <- row.names(subset(fData(cds), use_for_ordering == TRUE))
  
  g <- qplot(mean_expression, dispersion_empirical, data=disp_table, log="xy", color=I("darkgrey")) + 
    geom_line(aes(y=dispersion_fit), color="red") 
  if (length(ordering_genes) > 0){
    g <- g + geom_point(aes(mean_expression, dispersion_empirical), 
                        data=subset(disp_table, gene_id %in% ordering_genes), color="black")
  }
  g <- g + monocle_theme_opts()
  g
}

#' Plots clusters of cells .
#'
#' @param cds CellDataSet for the experiment
#' @param x the column of reducedDimS(cds) to plot on the horizontal axis
#' @param y the column of reducedDimS(cds) to plot on the vertical axis
#' @param color_by the cell attribute (e.g. the column of pData(cds)) to map to each cell's color
#' @param markers a gene name or gene id to use for setting the size of each cell in the plot
#' @param show_cell_names draw the name of each cell in the plot
#' @param cell_size The size of the point for each cell
#' @param cell_name_size the size of cell name labels
#' @param ... additional arguments passed into the scale_color_viridis function
#' @return a ggplot2 plot object
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom viridis scale_color_viridis
#' @export
#' @examples
#' \dontrun{
#' library(HSMMSingleCell)
#' HSMM <- load_HSMM()
#' HSMM <- reduceD
#' plot_cell_clusters(HSMM)
#' plot_cell_clusters(HSMM, color_by="Pseudotime")
#' plot_cell_clusters(HSMM, markers="MYH3")
#' }
plot_cell_clusters <- function(cds, 
                               x=1, 
                               y=2, 
                               color_by="Cluster", 
                               markers=NULL, 
                               show_cell_names=FALSE, 
                               cell_size=1.5,
                               cell_name_size=2, 
                               ...){
  if (is.null(cds@reducedDimA) | length(pData(cds)$Cluster) == 0){
    stop("Error: Clustering is not performed yet. Please call clusterCells() before calling this function.")
  }

  gene_short_name <- NULL
  sample_name <- NULL
  data_dim_1 <- NULL
  data_dim_2 <- NULL

  #TODO: need to validate cds as ready for this plot (need mst, pseudotime, etc)
  lib_info <- pData(cds)
  
  tSNE_dim_coords <- reducedDimA(cds)
  data_df <- data.frame(t(tSNE_dim_coords[c(x,y),]))
  colnames(data_df) <- c("data_dim_1", "data_dim_2")
  data_df$sample_name <- colnames(cds)
  data_df <- merge(data_df, lib_info, by.x="sample_name", by.y="row.names")
  
  markers_exprs <- NULL
  if (is.null(markers) == FALSE){
    markers_fData <- subset(fData(cds), gene_short_name %in% markers)
    if (nrow(markers_fData) >= 1){
      cds_subset <- cds[row.names(markers_fData),]
      if (cds_subset@expressionFamily@vfamily %in% c("negbinomial", "negbinomial.size")) {
        integer_expression <- TRUE
      }
      else {
        integer_expression <- FALSE

      }
      if (integer_expression) {
        cds_exprs <- exprs(cds_subset)

        if (is.null(sizeFactors(cds_subset))) {
         stop("Error: to call this function with relative_expr=TRUE, you must call estimateSizeFactors() first")
        }
        cds_exprs <- Matrix::t(Matrix::t(cds_exprs)/sizeFactors(cds_subset))

        cds_exprs <- reshape2::melt(round(as.matrix(cds_exprs)))
      }
      else {
        cds_exprs <- reshape2::melt(as.matrix(exprs(cds_subset)))
      }
      markers_exprs <- cds_exprs
      #markers_exprs <- reshape2::melt(as.matrix(cds_exprs))
      colnames(markers_exprs)[1:2] <- c('feature_id','cell_id')
      markers_exprs <- merge(markers_exprs, markers_fData, by.x = "feature_id", by.y="row.names")
      #print (head( markers_exprs[is.na(markers_exprs$gene_short_name) == FALSE,]))
      markers_exprs$feature_label <- as.character(markers_exprs$gene_short_name)
      markers_exprs$feature_label[is.na(markers_exprs$feature_label)] <- markers_exprs$Var1
    }
  }
  if (is.null(markers_exprs) == FALSE && nrow(markers_exprs) > 0){
    data_df <- merge(data_df, markers_exprs, by.x="sample_name", by.y="cell_id")

    g <- ggplot(data=data_df, aes(x=data_dim_1, y=data_dim_2)) + facet_wrap(~feature_label) 
  }else{
    g <- ggplot(data=data_df, aes(x=data_dim_1, y=data_dim_2)) 
  }
  
  # FIXME: setting size here overrides the marker expression funtionality. 
  # Don't do it!
  if (is.null(markers_exprs) == FALSE && nrow(markers_exprs) > 0){
    g <- g + geom_point(aes(color=log10(value + 0.1)), size=I(cell_size), na.rm = TRUE) + 
      scale_color_viridis(name = paste0("log10(value + 0.1)"), ...)
  }else {
    g <- g + geom_point(aes_string(color = color_by), size=I(cell_size), na.rm = TRUE)
  }
  
  g <- g + 
    #scale_color_brewer(palette="Set1") +
    monocle_theme_opts() + 
    xlab(paste("Component", x)) + 
    ylab(paste("Component", y)) +
    theme(legend.position="top", legend.key.height=grid::unit(0.35, "in")) +
    #guides(color = guide_legend(label.position = "top")) +
    theme(legend.key = element_blank()) +
    theme(panel.background = element_rect(fill='white')) +
    theme(text = element_text(size = 15))
  g
}

#' Plots the decision map of density clusters .
#'
#' @param cds CellDataSet for the experiment after running clusterCells_Density_Peak
#' @param rho_threshold The threshold of local density (rho) used to select the density peaks for plotting 
#' @param delta_threshold The threshold of local distance (delta) used to select the density peaks for plotting 
#' @export
#' @examples
#' \dontrun{
#' library(HSMMSingleCell)
#' HSMM <- load_HSMM()
#' plot_rho_delta(HSMM)
#' }

plot_rho_delta <- function(cds, rho_threshold = NULL, delta_threshold = NULL){
    if(!is.null(cds@auxClusteringData[["tSNE"]]$densityPeak) 
    & !is.null(pData(cds)$Cluster)
    & !is.null(pData(cds)$peaks)
    & !is.null(pData(cds)$halo)
    & !is.null(pData(cds)$delta)
    & !is.null(pData(cds)$rho)) {
    rho <- NULL
    delta <- NULL

    # df <- data.frame(rho = as.numeric(levels(pData(cds)$rho))[pData(cds)$rho], 
    #   delta = as.numeric(levels(pData(cds)$delta))[pData(cds)$delta])
    if(!is.null(rho_threshold) & !is.null(delta_threshold)){
      peaks <- pData(cds)$rho >= rho_threshold & pData(cds)$delta >= delta_threshold
    }
    else
      peaks <- pData(cds)$peaks

    df <- data.frame(rho = pData(cds)$rho, delta = pData(cds)$delta, peaks = peaks)

    g <- qplot(rho, delta, data = df, alpha = I(0.5), color = peaks) +  monocle_theme_opts() + 
      theme(legend.position="top", legend.key.height=grid::unit(0.35, "in")) +
      scale_color_manual(values=c("grey","black")) + 
      theme(legend.key = element_blank()) +
      theme(panel.background = element_rect(fill='white'))
  }
  else {
    stop('Please run clusterCells_Density_Peak before using this plotting function')
  }
  g
}

#' Plots the percentage of variance explained by the each component based on PCA from the normalized expression
#' data using the same procedure used in reduceDimension function.
#'
#' @param cds CellDataSet for the experiment after running reduceDimension with reduction_method as tSNE 
#' @param max_components Maximum number of components shown in the scree plot (variance explained by each component)
#' @param norm_method Determines how to transform expression values prior to reducing dimensionality
#' @param residualModelFormulaStr A model formula specifying the effects to subtract from the data before clustering.
#' @param pseudo_expr amount to increase expression values before dimensionality reduction
#' @param return_all A logical argument to determine whether or not the variance of each component is returned
#' @param use_existing_pc_variance Whether to plot existing results for variance explained by each PC
#' @param verbose Whether to emit verbose output during dimensionality reduction
#' @param ... additional arguments to pass to the dimensionality reduction function
#' @export
#' @examples
#' \dontrun{
#' library(HSMMSingleCell)
#' HSMM <- load_HSMM()
#' plot_pc_variance_explained(HSMM)
#' }
plot_pc_variance_explained <- function(cds, 
                            max_components=100, 
                            # reduction_method=c("DDRTree", "ICA", 'tSNE'),
                            norm_method = c("log", "vstExprs", "none"), 
                            residualModelFormulaStr=NULL,
                            pseudo_expr=NULL, 
                            return_all = F, 
                            use_existing_pc_variance=FALSE,
                            verbose=FALSE, 
                            ...){
  set.seed(2016)
  if(!is.null(cds@auxClusteringData[["tSNE"]]$variance_explained) & use_existing_pc_variance == T){
    prop_varex <- cds@auxClusteringData[["tSNE"]]$variance_explained
  }
  else{
    FM <- normalize_expr_data(cds, norm_method, pseudo_expr)
    
    #FM <- FM[unlist(sparseApply(FM, 1, sd, convert_to_dense=TRUE)) > 0, ]
    xm <- Matrix::rowMeans(FM)
    xsd <- sqrt(Matrix::rowMeans((FM - xm)^2))
    FM <- FM[xsd > 0,]
    
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
    
    # FM <- convert2DRData(cds, norm_method = 'log') 
    # FM <- FM[rowSums(is.na(FM)) == 0, ]
    irlba_res <- prcomp_irlba(t(FM), n = min(max_components, min(dim(FM)) - 1),
                              center = TRUE, scale. = TRUE)
    prop_varex <- irlba_res$sdev^2 / sum(irlba_res$sdev^2)
    # 
    # cell_means <- Matrix::rowMeans(FM_t)
    # cell_vars <- Matrix::rowMeans((FM_t - cell_means)^2)
    # 
    # irlba_res <- irlba(FM,
    #                    nv= min(max_components, min(dim(FM)) - 1), #avoid calculating components in the tail
    #                    nu=0,
    #                    center=cell_means,
    #                    scale=sqrt(cell_vars),
    #                    right_only=TRUE)
    # prop_varex <- irlba_res$d / sum(irlba_res$d)
    # 
    # # pca_res <- prcomp(t(FM), center = T, scale = T)
    # # std_dev <- pca_res$sdev 
    # # pr_var <- std_dev^2
    # prop_varex <- pr_var/sum(pr_var)
  }
  
  p <- qplot(1:length(prop_varex), prop_varex, alpha = I(0.5)) +  monocle_theme_opts() + 
    theme(legend.position="top", legend.key.height=grid::unit(0.35, "in")) +
    theme(panel.background = element_rect(fill='white')) + xlab('components') + 
    ylab('Variance explained \n by each component')
  
  cds@auxClusteringData[["tSNE"]]$variance_explained <- prop_varex # update CDS slot for variance_explained 

  if(return_all) {
    return(list(variance_explained = prop_varex, p = p))
  }
  else
    return(p)  
}

#' @importFrom igraph shortest_paths degree shortest.paths
traverseTree <- function(g, starting_cell, end_cells){
  distance <- shortest.paths(g, v=starting_cell, to=end_cells)
  branchPoints <- which(degree(g) == 3)
  path <- shortest_paths(g, from = starting_cell, end_cells)
  
  return(list(shortest_path = path$vpath, distance = distance, branch_points = intersect(branchPoints, unlist(path$vpath))))
}

#' Plots the minimum spanning tree on cells.
#' @description Plots the minimum spanning tree on cells.
#' @param cds CellDataSet for the experiment
#' @param x the column of reducedDimS(cds) to plot on the horizontal axis
#' @param y the column of reducedDimS(cds) to plot on the vertical axis
#' @param root_states the state used to set as the root of the graph
#' @param color_by the cell attribute (e.g. the column of pData(cds)) to map to each cell's color
#' @param show_tree whether to show the links between cells connected in the minimum spanning tree
#' @param show_backbone whether to show the diameter path of the MST used to order the cells
#' @param backbone_color the color used to render the backbone.
#' @param markers a gene name or gene id to use for setting the size of each cell in the plot
#' @param show_cell_names draw the name of each cell in the plot
#' @param cell_size The size of the point for each cell
#' @param cell_link_size The size of the line segments connecting cells (when used with ICA) or the principal graph (when used with DDRTree)
#' @param cell_name_size the size of cell name labels
#' @param show_branch_points Whether to show icons for each branch point (only available when reduceDimension was called with DDRTree)
#' @param ... Additional arguments passed to the scale_color_viridis function
#' @return a ggplot2 plot object
#' @import ggplot2
#' @importFrom igraph V get.edgelist layout_as_tree
#' @importFrom reshape2 melt
#' @importFrom viridis scale_color_viridis
#' @export
#' @examples
#' \dontrun{
#' library(HSMMSingleCell)
#' HSMM <- load_HSMM()
#' plot_complex_cell_trajectory(HSMM)
#' plot_complex_cell_trajectory(HSMM, color_by="Pseudotime", show_backbone=FALSE)
#' plot_complex_cell_trajectory(HSMM, markers="MYH3")
#' }
plot_complex_cell_trajectory <- function(cds, 
                                         x=1, 
                                         y=2, 
                                         root_states = NULL,
                                         color_by="State", 
                                         show_tree=TRUE, 
                                         show_backbone=TRUE, 
                                         backbone_color="black", 
                                         markers=NULL, 
                                         show_cell_names=FALSE, 
                                         cell_size=1.5,
                                         cell_link_size=0.75,
                                         cell_name_size=2,
                                         show_branch_points=TRUE, 
                                         ...){
  gene_short_name <- NA
  sample_name <- NA
  data_dim_1 <- NA
  data_dim_2 <- NA
  
  #TODO: need to validate cds as ready for this plot (need mst, pseudotime, etc)
  lib_info_with_pseudo <- pData(cds)
  
  if (is.null(cds@dim_reduce_type)){
    stop("Error: dimensionality not yet reduced. Please call reduceDimension() before calling this function.")
  }
  
  if (cds@dim_reduce_type == "ICA"){
    reduced_dim_coords <- reducedDimS(cds)
  }else if (cds@dim_reduce_type %in% c("SimplePPT", "DDRTree", "SGL-tree") ){
    reduced_dim_coords <- reducedDimK(cds)
    closest_vertex <- cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex
  }else {
    stop("Error: unrecognized dimensionality reduction method.")
  }
  
  if (is.null(reduced_dim_coords)){
    stop("You must first call reduceDimension() before using this function")
  }
  
  dp_mst <- minSpanningTree(cds)
  
  
  if(is.null(root_states)) {
    if(is.null(lib_info_with_pseudo$Pseudotime)){
      root_cell <- row.names(lib_info_with_pseudo)[degree(dp_mst) == 1][1]
    }
    else
      root_cell <- row.names(subset(lib_info_with_pseudo, Pseudotime == 0))
    
    if(cds@dim_reduce_type != "ICA")
      root_cell <- V(dp_mst)$name[cds@auxOrderingData$DDRTree$pr_graph_cell_proj_closest_vertex[root_cell, ]] 
    
  }
  else {
    candidate_root_cells <- row.names(subset(pData(cds), State %in% root_states))
    if(cds@dim_reduce_type == "ICA") {
      root_cell <- candidate_root_cells[which(degree(dp_mst, candidate_root_cells) == 1)]
    } else {
      Y_candidate_root_cells <- V(dp_mst)$name[cds@auxOrderingData$DDRTree$pr_graph_cell_proj_closest_vertex[candidate_root_cells, ]] 
      root_cell <- Y_candidate_root_cells[which(degree(dp_mst, Y_candidate_root_cells) == 1)]
    }
    
  }
  
  # #root_cell <- cds@auxOrderingData[[cds@dim_reduce_type]]$root_cell
  # root_state <- pData(cds)[root_cell,]$State
  # #root_state <- V(pr_graph_cell_proj_mst)[root_cell,]$State
  
  # pr_graph_root <- subset(pData(cds), State == root_state)
  
  # closest_vertex <- cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex
  # root_cell_point_in_Y <- closest_vertex[row.names(pr_graph_root),]
  tree_coords <- layout_as_tree(dp_mst, root=root_cell)
  
  #ica_space_df <- data.frame(Matrix::t(reduced_dim_coords[c(x,y),]))
  ica_space_df <- data.frame(tree_coords)
  row.names(ica_space_df) <- colnames(reduced_dim_coords)
  colnames(ica_space_df) <- c("prin_graph_dim_1", "prin_graph_dim_2")
  
  ica_space_df$sample_name <- row.names(ica_space_df)
  #ica_space_with_state_df <- merge(ica_space_df, lib_info_with_pseudo, by.x="sample_name", by.y="row.names")
  #print(ica_space_with_state_df)
  
  
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
  
  #S_matrix <- reducedDimS(cds)
  #data_df <- data.frame(t(S_matrix[c(x,y),]))
  
  if(cds@dim_reduce_type == "ICA"){
    S_matrix <- tree_coords[,] #colnames(cds)
    
  } else if(cds@dim_reduce_type %in% c("DDRTree", "SimplePPT", "SGL-tree")){
    S_matrix <- tree_coords[closest_vertex,]
    closest_vertex <- cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex
  }
  
  data_df <- data.frame(S_matrix)
  row.names(data_df) <- colnames(reducedDimS(cds))
  colnames(data_df) <- c("data_dim_1", "data_dim_2")
  data_df$sample_name <- row.names(data_df)
  data_df <- merge(data_df, lib_info_with_pseudo, by.x="sample_name", by.y="row.names")
  
  markers_exprs <- NULL
  if (is.null(markers) == FALSE){
    markers_fData <- subset(fData(cds), gene_short_name %in% markers)
    if (nrow(markers_fData) >= 1){
      markers_exprs <- reshape2::melt(as.matrix(exprs(cds[row.names(markers_fData),])))
      colnames(markers_exprs)[1:2] <- c('feature_id','cell_id')
      markers_exprs <- merge(markers_exprs, markers_fData, by.x = "feature_id", by.y="row.names")
      #print (head( markers_exprs[is.na(markers_exprs$gene_short_name) == FALSE,]))
      markers_exprs$feature_label <- as.character(markers_exprs$gene_short_name)
      markers_exprs$feature_label[is.na(markers_exprs$feature_label)] <- markers_exprs$Var1
    }
  }
  if (is.null(markers_exprs) == FALSE && nrow(markers_exprs) > 0){
    data_df <- merge(data_df, markers_exprs, by.x="sample_name", by.y="cell_id")
    #print (head(edge_df))
    g <- ggplot(data=data_df, aes(x=data_dim_1, y=data_dim_2, I(cell_size))) + facet_wrap(~feature_label)
  }else{
    g <- ggplot(data=data_df, aes(x=data_dim_1, y=data_dim_2)) 
  }
  if (show_tree){
    g <- g + geom_segment(aes_string(x="source_prin_graph_dim_1", y="source_prin_graph_dim_2", xend="target_prin_graph_dim_1", yend="target_prin_graph_dim_2"), size=cell_link_size, linetype="solid", na.rm=TRUE, data=edge_df)
  }
  
  # FIXME: setting size here overrides the marker expression funtionality. 
  # Don't do it!
  if (is.null(markers_exprs) == FALSE && nrow(markers_exprs) > 0){
    if(class(data_df[, color_by]) == 'numeric') {
      g <- g + geom_jitter(aes_string(color = paste0("log10(", color_by, " + 0.1)")), size=I(cell_size), na.rm = TRUE, height=5) + 
                             scale_color_viridis(name = paste0("log10(", color_by, ")"), ...)
    } else {
      g <- g + geom_jitter(aes_string(color = color_by), size=I(cell_size), na.rm = TRUE, height=5) 
    }
  }else {
    if(class(data_df[, color_by]) == 'numeric') {
      g <- g + geom_jitter(aes_string(color = paste0("log10(", color_by, " + 0.1)")), size=I(cell_size), na.rm = TRUE, height=5) + 
        scale_color_viridis(name = paste0("log10(", color_by, " + 0.1)"), ...)
    } else {
      g <- g + geom_jitter(aes_string(color = color_by), size=I(cell_size), na.rm = TRUE, height=5)
    }
  }

  if (show_branch_points && cds@dim_reduce_type == 'DDRTree'){
    mst_branch_nodes <- cds@auxOrderingData[[cds@dim_reduce_type]]$branch_points
    branch_point_df <- subset(edge_df, sample_name %in% mst_branch_nodes)[,c("sample_name", "source_prin_graph_dim_1", "source_prin_graph_dim_2")]
    branch_point_df$branch_point_idx <- match(branch_point_df$sample_name, mst_branch_nodes)
    branch_point_df <- branch_point_df[!duplicated(branch_point_df$branch_point_idx), ]
    
    g <- g + geom_point(aes_string(x="source_prin_graph_dim_1", y="source_prin_graph_dim_2"), 
                        size=2 * cell_size, na.rm=TRUE, data=branch_point_df) +
      geom_text(aes_string(x="source_prin_graph_dim_1", y="source_prin_graph_dim_2", label="branch_point_idx"), 
                size=1.5 * cell_size, color="white", na.rm=TRUE, data=branch_point_df)
  }
  if (show_cell_names){
    g <- g +geom_text(aes(label=sample_name), size=cell_name_size)
  }
  g <- g + 
    #scale_color_brewer(palette="Set1") +
    theme(strip.background = element_rect(colour = 'white', fill = 'white')) +
    theme(panel.border = element_blank()) +
    # theme(axis.line.x = element_line(size=0.25, color="black")) +
    # theme(axis.line.y = element_line(size=0.25, color="black")) +
    theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) +
    theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) + 
    theme(panel.background = element_rect(fill='white')) +
    theme(legend.key=element_blank()) + 
    xlab('') + 
    ylab('') +
    theme(legend.position="top", legend.key.height=grid::unit(0.35, "in")) +
    #guides(color = guide_legend(label.position = "top")) +
    theme(legend.key = element_blank()) +
    theme(panel.background = element_rect(fill='white')) + 
    theme(line = element_blank(), 
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank()) 
  g
}

# Modified function: Plot heatmap of 3 branches with the same coloring. Each CDS subset has to have the same set of genes.
#'  Create a heatmap to demonstrate the bifurcation of gene expression along multiple branches
#'
#' @param cds CellDataSet for the experiment (normally only the branching genes detected with BEAM)
#' @param branches The terminal branches (states) on the developmental tree you want to investigate.
#' @param branches_name Name (for example, cell type) of branches you believe the cells on the branches are associated with. 
#' @param cluster_rows Whether to cluster the rows of the heatmap.
#' @param hclust_method The method used by pheatmap to perform hirearchical clustering of the rows. 
#' @param num_clusters Number of clusters for the heatmap of branch genes
#' @param hmcols The color scheme for drawing the heatmap.
#' @param add_annotation_row Additional annotations to show for each row in the heatmap. Must be a dataframe with one row for each row in the fData table of cds_subset, with matching IDs.
#' @param add_annotation_col Additional annotations to show for each column in the heatmap. Must be a dataframe with one row for each cell in the pData table of cds_subset, with matching IDs.
#' @param show_rownames Whether to show the names for each row in the table.
#' @param use_gene_short_name Whether to use the short names for each row. If FALSE, uses row IDs from the fData table.
#' @param norm_method Determines how to transform expression values prior to rendering
#' @param scale_max The maximum value (in standard deviations) to show in the heatmap. Values larger than this are set to the max.
#' @param scale_min The minimum value (in standard deviations) to show in the heatmap. Values smaller than this are set to the min.
#' @param trend_formula A formula string specifying the model used in fitting the spline curve for each gene/feature.
#' @param return_heatmap Whether to return the pheatmap object to the user. 
#' @param cores Number of cores to use when smoothing the expression curves shown in the heatmap.
#' @return A list of heatmap_matrix (expression matrix for the branch committment), ph (pheatmap heatmap object),
#' annotation_row (annotation data.frame for the row), annotation_col (annotation data.frame for the column). 
#' @import pheatmap
#' @export
#'
plot_multiple_branches_heatmap <- function(cds, 
                                           branches, 
                                           branches_name = NULL, 
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
  pseudocount <- 1
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
    
    cds_subset <- cds[, path_cells]
    
    newdata <- data.frame(Pseudotime = seq(0, max(pData(cds_subset)$Pseudotime),length.out = 100)) 
    
    tmp <- genSmoothCurves(cds_subset, cores=cores, trend_formula = trend_formula,  
                           relative_expr = T, new_data = newdata)
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
  
  if(is.null(hmcols)) {
    bks <- seq(-3.1,3.1, by = 0.1)
    hmcols <- blue2green2red(length(bks) - 1)
  }
  else {
    bks <- seq(-3.1,3.1, length.out = length(hmcols))
  } 
  
  ph <- pheatmap(heatmap_matrix, 
                 useRaster = T,
                 cluster_cols=FALSE, 
                 cluster_rows=T, 
                 show_rownames=F, 
                 show_colnames=F, 
                 clustering_distance_rows=row_dist,
                 clustering_method = hclust_method,
                 cutree_rows=num_clusters,
                 silent=TRUE,
                 filename=NA,
                 breaks=bks,
                 color=hmcols)
  
  annotation_col <- data.frame(Branch=factor(rep(rep(branch_label, each = 100))))
  annotation_row <- data.frame(Cluster=factor(cutree(ph$tree_row, num_clusters)))
  col_gaps_ind <- c(1:(length(branches) - 1)) * 100
  
  if(!is.null(add_annotation_row)) {
    old_colnames_length <- ncol(annotation_row)
    annotation_row <- cbind(annotation_row, add_annotation_row[row.names(annotation_row), ])  
    colnames(annotation_row)[(old_colnames_length+1):ncol(annotation_row)] <- colnames(add_annotation_row)
    # annotation_row$bif_time <- add_annotation_row[as.character(fData(absolute_cds[row.names(annotation_row), ])$gene_short_name), 1]
  }
  
  
  if (use_gene_short_name == TRUE) {
    if (is.null(fData(cds)$gene_short_name) == FALSE) {
      feature_label <- as.character(fData(cds)[row.names(heatmap_matrix), 'gene_short_name'])
      feature_label[is.na(feature_label)] <- row.names(heatmap_matrix)
      
      row_ann_labels <- as.character(fData(cds)[row.names(annotation_row), 'gene_short_name'])
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
  
  if(!(cluster_rows)) {
    annotation_row <- NA
  }
  
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
                     annotation_col=annotation_col,
                     gaps_col = col_gaps_ind,
                     treeheight_row = 20, 
                     breaks=bks,
                     fontsize = 12,
                     color=hmcols, 
                     silent=TRUE,
                     border_color = NA,
                     filename=NA
  )
  
  grid::grid.rect(gp=grid::gpar("fill", col=NA))
  grid::grid.draw(ph_res$gtable)
  if (return_heatmap){
    return(ph_res)
  }
}

#'  Create a kinetic curves to demonstrate the bifurcation of gene expression along multiple branches
#'
#' @param cds CellDataSet for the experiment (normally only the branching genes detected with BEAM)
#' @param branches The terminal branches (states) on the developmental tree you want to investigate. 
#' @param branches_name Name (for example, cell type) of branches you believe the cells on the branches are associated with. 
#' @param min_expr The minimum level of expression to show in the plot
#' @param cell_size A number how large the cells should be in the plot 
#' @param norm_method Determines how to transform expression values prior to rendering
#' @param nrow the number of rows used when laying out the panels for each gene's expression
#' @param ncol the number of columns used when laying out the panels for each gene's expression
#' @param panel_order the order in which genes should be layed out (left-to-right, top-to-bottom)
#' @param color_by the cell attribute (e.g. the column of pData(cds)) to be used to color each cell  
#' @param trend_formula the model formula to be used for fitting the expression trend over pseudotime 
#' @param label_by_short_name label figure panels by gene_short_name (TRUE) or feature id (FALSE)
#' @param TPM Whether to convert the expression value into TPM values. 
#' @param cores Number of cores to use when smoothing the expression curves shown in the heatmap.
#' @return a ggplot2 plot object
#' 
#' @importFrom Biobase esApply exprs<-
#' @importFrom stats lowess
#' 
#' @export
#'
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
        
        #if(is.null(pData(cds)$no_expression)) {
        cds_subset <- cds[, path_cells]      
        # } else {
        #     cds_subset <- cds[, path_cells %in% colnames(cds)[!pData(cds)$no_expression]]      
        # }
        
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

  

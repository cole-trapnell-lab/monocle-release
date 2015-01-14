utils::globalVariables(c("Pseudotime", "value", "ids", "ICA_dim_1", "ICA_dim_2", "State", 
                         "value", "feature_label", "expectation", "colInd", "rowInd", "value", 
                         "source_ICA_dim_1", "source_ICA_dim_2"))

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
                               cell_name_size=1){
  #TODO: need to validate cds as ready for this plot (need mst, pseudotime, etc)
  lib_info_with_pseudo <- pData(cds)

  
  #print (lib_info_with_pseudo)
  S_matrix <- reducedDimS(cds)

  if (is.null(S_matrix)){
    stop("You must first call reduceDimension() before using this function")
  }
  
  ica_space_df <- data.frame(t(S_matrix[c(x,y),]))
  colnames(ica_space_df) <- c("ICA_dim_1", "ICA_dim_2")

  ica_space_df$sample_name <- row.names(ica_space_df)
  ica_space_with_state_df <- merge(ica_space_df, lib_info_with_pseudo, by.x="sample_name", by.y="row.names")
  #print(ica_space_with_state_df)
  dp_mst <- minSpanningTree(cds)
 
  if (is.null(dp_mst)){
    stop("You must first call orderCells() before using this function")
  }
  
  
  edge_list <- as.data.frame(get.edgelist(dp_mst))
  colnames(edge_list) <- c("source", "target")

  edge_df <- merge(ica_space_with_state_df, edge_list, by.x="sample_name", by.y="source", all=TRUE)
  
  edge_df <- rename(edge_df, c("ICA_dim_1"="source_ICA_dim_1", "ICA_dim_2"="source_ICA_dim_2"))
  edge_df <- merge(edge_df, ica_space_with_state_df[,c("sample_name", "ICA_dim_1", "ICA_dim_2")], by.x="target", by.y="sample_name", all=TRUE)
  edge_df <- rename(edge_df, c("ICA_dim_1"="target_ICA_dim_1", "ICA_dim_2"="target_ICA_dim_2"))
  
  diam <- as.data.frame(as.vector(V(dp_mst)[get.diameter(dp_mst, weights=NA)]$name))
  colnames(diam) <- c("sample_name")
  diam <- plyr::arrange(merge(ica_space_with_state_df,diam, by.x="sample_name", by.y="sample_name"), Pseudotime)

  markers_exprs <- NULL
  if (is.null(markers) == FALSE){
    markers_fData <- subset(fData(cds), gene_short_name %in% markers)
    if (nrow(markers_fData) >= 1){
      markers_exprs <- reshape2::melt(exprs(cds[row.names(markers_fData),]))
      markers_exprs <- merge(markers_exprs, markers_fData, by.x = "Var1", by.y="row.names")
      #print (head( markers_exprs[is.na(markers_exprs$gene_short_name) == FALSE,]))
      markers_exprs$feature_label <- as.character(markers_exprs$gene_short_name)
      markers_exprs$feature_label[is.na(markers_exprs$feature_label)] <- markers_exprs$Var1
    }
  }
  if (is.null(markers_exprs) == FALSE && nrow(markers_exprs) > 0){
    edge_df <- merge(edge_df, markers_exprs, by.x="sample_name", by.y="Var2")
    #print (head(edge_df))
    g <- ggplot(data=edge_df, aes(x=source_ICA_dim_1, y=source_ICA_dim_2, size=log10(value + 0.1))) + facet_wrap(~feature_label)
  }else{
    g <- ggplot(data=edge_df, aes(x=source_ICA_dim_1, y=source_ICA_dim_2)) 
  }
  if (show_tree){
    g <- g + geom_segment(aes_string(xend="target_ICA_dim_1", yend="target_ICA_dim_2", color=color_by), size=.3, linetype="solid", na.rm=TRUE)
  }
  
  g <- g +geom_point(aes_string(color=color_by), na.rm=TRUE) 
  if (show_backbone){
    #print (diam)
    g <- g +geom_path(aes(x=ICA_dim_1, y=ICA_dim_2), color=I(backbone_color), size=0.75, data=diam, na.rm=TRUE) + 
    geom_point(aes_string(x="ICA_dim_1", y="ICA_dim_2", color=color_by), size=I(1.5), data=diam, na.rm=TRUE)
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
                              min_expr=0.1, cell_size=0.75, nrow=NULL, ncol=1, panel_order=NULL, 
                              color_by=NULL,
                              plot_trend=FALSE,
                              label_by_short_name=TRUE,
                              relative_expr=TRUE){
  
  if (cds_subset@expressionFamily@vfamily == "negbinomial"){
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
      cds_exprs <- t(t(cds_exprs) / sizeFactors(cds_subset))
    }
    cds_exprs <- reshape2::melt(round(cds_exprs))
  }else{
    cds_exprs <- exprs(cds_subset)
    cds_exprs <- reshape2::melt(cds_exprs)
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
#' @return a ggplot2 plot object
#' @import ggplot2
#' @importFrom plyr dplyr
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
  
  
  if (cds_subset@expressionFamily@vfamily == "negbinomial"){
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
      marker_exprs <- t(t(marker_exprs) / sizeFactors(cds_subset))
    }
    marker_exprs_melted <- reshape2::melt(round(marker_exprs))
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
    qp <- ggplot(aes_string(x=grouping, y="target_fraction", fill=grouping), data=marker_counts) +
      scale_y_continuous(labels=percent, limits=c(0,1)) 
  }else{
    qp <- ggplot(aes_string(x=grouping, y="target", fill=grouping), data=marker_counts)
  }
  
  qp <- qp + facet_wrap(~feature_label, nrow=nrow, ncol=ncol, scales="free_y")
  qp <-  qp + geom_bar(stat="identity") +
    ylab("Cells") + monocle_theme_opts()

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
#' @return a ggplot2 plot object
#' @import ggplot2
#' @importFrom plyr ddply
#' @importFrom reshape2 melt
#' @importFrom VGAM vgam predict
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
                                    trend_formula="adjusted_expression ~ sm.ns(Pseudotime, df=3)",
                                    label_by_short_name=TRUE,
                                    relative_expr=TRUE){
  
  if (cds_subset@expressionFamily@vfamily == "negbinomial"){
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
      cds_exprs <- t(t(cds_exprs) / sizeFactors(cds_subset))
    }
    cds_exprs <- reshape2::melt(round(cds_exprs))
  }else{
    cds_exprs <- reshape2::melt(exprs(cds_subset))
  }
  if (is.null(min_expr)){
    min_expr <- cds_subset@lowerDetectionLimit
  }
  colnames(cds_exprs) <- c("f_id", "Cell", "expression")
  
  cds_pData <- pData(cds_subset)
  cds_fData <- fData(cds_subset)
  
  cds_exprs <- merge(cds_exprs, cds_fData, by.x="f_id", by.y="row.names")
  cds_exprs <- merge(cds_exprs, cds_pData, by.x="Cell", by.y="row.names")
  
  if (integer_expression)
  {
    cds_exprs$adjusted_expression <- cds_exprs$expression
  }else{
    cds_exprs$adjusted_expression <- log10(cds_exprs$expression)
  }
  
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
   
  merged_df_with_vgam <- plyr::ddply(cds_exprs, .(feature_label), function(x) { 
    fit_res <- tryCatch({
      #Extra <- list(leftcensored = with(x, adjusted_fpkm <= min_fpkm), rightcencored = rep(FALSE, nrow(x)))
      vg <- suppressWarnings(VGAM::vgam(formula = as.formula(trend_formula), 
                                  family = cds_subset@expressionFamily, 
                                  data = x, maxit=30, checkwz=FALSE))
      if (integer_expression){
        res <- predict(vg, type="response")
        res[res < min_expr] <- min_expr
      }else{
        res <- 10^(predict(vg, type="response"))
        res[res < log10(min_expr)] <- log10(min_expr)
      }
      res
    }
    ,error = function(e) {
      print("Error!")
      print(e)
      res <- rep(NA, nrow(x))
      res
    }
    )
    
    expectation = fit_res

    data.frame(Pseudotime=x$Pseudotime, expectation=expectation)
  })
  
  if (is.null(panel_order) == FALSE)
  {
    cds_subset$feature_label <- factor(cds_subset$feature_label, levels=panel_order)
  }
  
  q <- ggplot(aes(Pseudotime, expression), data=cds_exprs) 
  
  if (is.null(color_by) == FALSE){
    q <- q + geom_point(aes_string(color=color_by), size=I(cell_size))
  }else{
    q <- q + geom_point(size=I(cell_size))
  }
    
  q <- q + geom_line(aes(Pseudotime, expectation), data=merged_df_with_vgam)
  
  #q <- q + geom_ribbon(aes(x=pseudo_time, ymin=conf_lo, ymax=conf_hi), alpha=I(0.15), data=merged_df_with_vgam)
  
  q <- q + scale_y_log10() + facet_wrap(~feature_label, nrow=nrow, ncol=ncol, scales="free_y")
  
  if (relative_expr){
    q <- q + ylab("Relative Expression")
  }else{
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
#' @export
#' @examples
#' \dontrun{
#' full_model_fits <- fitModel(HSMM_filtered[sample(nrow(fData(HSMM_filtered)), 100),],  modelFormulaStr="expression~VGAM::bs(Pseudotime)")
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
  m <- as.data.frame(clustering$exprs)
  m$ids <- rownames(clustering$exprs)
  if (is.null(clustering$labels) == FALSE)
  {
    m$cluster = factor(clustering$labels[clustering$clustering], levels = levels(clustering$labels))
  }else{
    m$cluster <- factor(clustering$clustering)
  }
  m.melt <- reshape2::melt(m, id.vars = c("ids", "cluster"))

  m.melt <- merge(m.melt, pData(cds), by.x="variable", by.y="row.names")
  
  if (is.null(row_samples) == FALSE){
    m.melt <- m.melt[sample(nrow(m.melt), row_samples),]
  }
  
  c <- ggplot(m.melt)
  c <- c + stat_density2d(aes(x = Pseudotime, y = value), geom="polygon", fill="white", color="black", size=I(0.1)) + facet_wrap("cluster", ncol=ncol, nrow=nrow)
  
  if (drawSummary) {
    c <- c + stat_summary(aes(x = Pseudotime, y = value, group = 1), 
                          fun.data = sumFun, color = "red", fill = "black", 
                          alpha = 0.2, size = 0.5, geom = "smooth")
  }
  
  
  c <- c + scale_color_hue(l = 50, h.start = 200) + theme(axis.text.x = element_text(angle = 0, 
                                                                                     hjust = 0)) + xlab("Pseudo-time") + ylab("Expression")
  c <- c + theme(strip.background = element_rect(colour = 'white', fill = 'white')) + 
    theme(panel.border = element_blank(), axis.line = element_line(size=0.2)) +
    theme(axis.ticks = element_line(size = 0.2)) +
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
  c
}

plot_genes_heatmap <- function(cds, 
                               rescaling='row', 
                               clustering='row', 
                               labCol=FALSE, 
                               labRow=TRUE, 
                               logMode=TRUE, 
                               pseudocount=0.1, 
                               border=FALSE, 
                               heatscale=c(low='steelblue',mid='white',high='tomato'), 
                               heatMidpoint=0,
                               method="none",
                               scaleMax=2, 
                               scaleMin=-2, 
                               relative_expr=TRUE, 
                               ...){
  
  
  ## the function can be be viewed as a two step process
  ## 1. using the rehape package and other funcs the data is clustered, scaled, and reshaped
  ## using simple options or by a user supplied function
  ## 2. with the now resahped data the plot, the chosen labels and plot style are built
  FM <- exprs(cds)
  
  if (cds@expressionFamily@vfamily == "negbinomial"){
    integer_expression <- TRUE
  }else{
    integer_expression <- FALSE
    relative_expr <- TRUE
  }
  
  if (integer_expression)
  {
    if (relative_expr){
      if (is.null(sizeFactors(cds)))
      {
        stop("Error: you must call estimateSizeFactors() first")
      }
      FM <- t(t(FM) / sizeFactors(cds_subset))
    }
    FM <- round(FM)
  }
  
  m=FM
  
  if (is.null(fData(cds)$gene_short_name) == FALSE){
    feature_labels <- fData(cds)$gene_short_name
    feature_labels[is.na(feature_labels)]  <- fData(cds)$f_id
    row.names(m) <- feature_labels
  }
  
  #remove genes with no expression in any condition
  m=m[!apply(m,1,sum)==0,]
  
  if(logMode) 
  {
    m = log10(m+pseudocount)
  }
  
  #remove genes with no sd
  #m=m[!apply(m,1,sd)==0,]

  ## you can either scale by row or column not both! 
  ## if you wish to scale by both or use a different scale method then simply supply a scale
  ## function instead NB scale is a base funct
  
  ## I have supplied the default cluster and euclidean distance (JSdist) - and chose to cluster after scaling
  ## if you want a different distance/cluster method-- or to cluster and then scale
  ## then you can supply a custom function 
  
  if(!is.function(method)){
    method = function(mat){as.dist((1 - cor(t(mat)))/2)}	
  }
  
  ## this is just reshaping into a ggplot format matrix and making a ggplot layer
  
  if(is.function(rescaling))
  { 
    m=rescaling(m)
  } else {
    if(rescaling=='column'){
      m=scale(m, center=TRUE)
      m[is.nan(m)] = 0
      m[m>scaleMax] = scaleMax
      m[m<scaleMin] = scaleMin
    }
    if(rescaling=='row'){ 
      m=t(scale(t(m),center=TRUE))
      m[is.nan(m)] = 0
      m[m>scaleMax] = scaleMax
      m[m<scaleMin] = scaleMin
    }
  }
  
  # If we aren't going to re-ordering the columns, order them by Pseudotime
  if (clustering %in% c("row", "none"))
    m = m[,as.character(plyr::arrange(pData(cds), Pseudotime)$Cell)]
  
  if(clustering=='row')
    m=m[hclust(method(m))$order, ]
  if(clustering=='column')  
    m=m[,hclust(method(t(m)))$order]
  if(clustering=='both')
    m=m[hclust(method(m))$order ,hclust(method(t(m)))$order]
  
  
  rows=dim(m)[1]
  cols=dim(m)[2]
  
  
  
  # if(logMode) {
  #   melt.m=cbind(rowInd=rep(1:rows, times=cols), colInd=rep(1:cols, each=rows), reshape2::melt( log10(m+pseudocount)))
  # }else{
  #   melt.m=cbind(rowInd=rep(1:rows, times=cols), colInd=rep(1:cols, each=rows), reshape2::melt(m))
  # }
  
  melt.m=cbind(rowInd=rep(1:rows, times=cols), colInd=rep(1:cols, each=rows), reshape2::melt(m))
  
  g=ggplot(data=melt.m)
  
  ## add the heat tiles with or without a white border for clarity
  
  if(border==TRUE)
    g2=g+geom_rect(aes(xmin=colInd-1,xmax=colInd,ymin=rowInd-1,ymax=rowInd, fill=value),colour='grey')
  if(border==FALSE)
    g2=g+geom_rect(aes(xmin=colInd-1,xmax=colInd,ymin=rowInd-1,ymax=rowInd, fill=value))
  
  ## add axis labels either supplied or from the colnames rownames of the matrix
  
  if(labCol==TRUE) 
  {
    g2=g2+scale_x_continuous(breaks=(1:cols)-0.5, labels=colnames(m))
  }
  if(labCol==FALSE) 
  {
    g2=g2+scale_x_continuous(breaks=(1:cols)-0.5, labels=rep('',cols))
  }
  
  
  if(labRow==TRUE) 
  {
    g2=g2+scale_y_continuous(breaks=(1:rows)-0.5, labels=rownames(m))	
  }
  if(labRow==FALSE)
  { 
    g2=g2+scale_y_continuous(breaks=(1:rows)-0.5, labels=rep('',rows))	
  }
  
  # Get rid of the ticks, they get way too dense with lots of rows
  g2 <- g2 + theme(axis.ticks = element_blank()) 
  
  ## get rid of grey panel background and gridlines
  
  g2=g2+theme(panel.grid.minor=element_line(colour=NA), panel.grid.major=element_line(colour=NA),
              panel.background=element_rect(fill=NA, colour=NA))
  
  ##adjust x-axis labels
  g2=g2+theme(axis.text.x=element_text(angle=-90, hjust=0))
  
  #write(paste(c("Length of heatscale is :", length(heatscale))), stderr())
  
  if(is.function(rescaling))
  {
    
  }else{ 
    if(rescaling=='row' || rescaling == 'column'){
      legendTitle <- "Relative\nexpression"
    }else{
      if (logMode)
      {
        legendTitle <- bquote(paste(log[10]," FPKM + ",.(pseudocount),sep=""))
        #legendTitle <- paste(expression(plain(log)[10])," FPKM + ",pseudocount,sep="")
      } else {
        legendTitle <- "FPKM"
      }
    }
  }
  
  if (length(heatscale) == 2){
    g2 <- g2 + scale_fill_gradient(low=heatscale[1], high=heatscale[2], name=legendTitle)
  } else if (length(heatscale) == 3) {
    if (is.null(heatMidpoint))
    {
      heatMidpoint = (max(m) + min(m)) / 2.0
      #write(heatMidpoint, stderr())
    }
    g2 <- g2 + theme(panel.border = element_blank())
    g2 <- g2 + scale_fill_gradient2(low=heatscale[1], mid=heatscale[2], high=heatscale[3], midpoint=heatMidpoint, name=legendTitle)
  }else {
    g2 <- g2 + scale_fill_gradientn(colours=heatscale, name=legendTitle)
  }
  
  #g2<-g2+scale_x_discrete("",breaks=tracking_ids,labels=gene_short_names)
  
  
  ## finally add the fill colour ramp of your choice (default is blue to red)-- and return
  return (g2)
}

#' Plot the branch genes in pseduotime with separate lineage curves 
#' TO DO: add the confidence interval, add p-val, fix error when we cannot generate fitting for the data
#' @param cds CellDataSet for the experiment
#' @param ncol number of columns used to layout the faceted cluster panels
#' @param nrow number of columns used to layout the faceted cluster panels
#' @param row_samples how many genes to randomly select from the data
#' @return a ggplot2 plot object
#' @import ggplot2
#' @importFrom plyr ddply
#' @importFrom reshape2 melt
#' @export 
#' @examples
#' \dontrun{
#' full_model_fits <- fitModel(HSMM_filtered[sample(nrow(fData(HSMM_filtered)), 100),],  modelFormulaStr="expression~VGAM::bs(Pseudotime)")
#' expression_curve_matrix <- responseMatrix(full_model_fits)
#' clusters <- clusterGenes(expression_curve_matrix, k=4)
#' plot_clusters(HSMM_filtered[ordering_genes,], clusters)
#' }
#' 
#FIXME: This function is hopelessly buggy and broken.  Needs a re-write and THOROUGH TESTING
plot_genes_branched_pseudotime <- function (cds, lineage_states = c(2, 3), method = 'fitting', stretch = FALSE, 
                               min_expr = NULL, cell_size = 0.75, nrow = NULL, ncol = 1, panel_order = NULL, color_by = "State", 
                               fullModelFormulaStr = "expression ~ s(Pseudotime, df=3) * State", label_by_short_name = TRUE) 
{
  if(is.null(pData(cds)$State) | is.null(pData(cds)$Pseudotime)) 
    stop('Please first order the cells in pseudotime using orderCells()')
  
  lineage_cells <- row.names(pData(cds[,pData(cds)$State %in% lineage_states]))

  curr_cell <- setdiff(pData(cds[,lineage_cells])$Parent, lineage_cells)
  ancestor_cells <- c()
  
  while (1) {
    ancestor_cells <- c(ancestor_cells, curr_cell)
    if (is.na(pData(cds[,curr_cell])$Parent))
      break
    curr_cell <- as.character(pData(cds[,curr_cell])$Parent)
  }
  
  cds <- cds[, row.names(pData(cds[,union(ancestor_cells, lineage_cells)]))]
  
  State <- pData(cds)$State 
  Pseudotime <- pData(cds)$Pseudotime 
  

  progenitor_ind <- which(State %in% lineage_states == FALSE)
  
  progenitor_states <- unique(State[progenitor_ind])
  
  pData <- pData(cds)
  exprs_data <- exprs(cds)
  
  weight_vec <- rep(1, nrow(pData(cds)))
  weight_vec[progenitor_ind] <- 0.5
  
  range_df <- plyr::ddply(pData(cds), .(State), function(x) { range(x$Pseudotime)}) #pseudotime range for each state
  row.names(range_df) <- as.character(range_df$State)
  colnames(range_df) <- c("State","min", "max")
  
  print (range_df)
  #longest_branch <- groups[which(range_df[, 2] == max(range_df[, 2]))] 
  longest_lineage_branch <- tail(plyr::arrange(range_df[as.character(lineage_states),], max)$State, n=1)

  #first stretch pseudotime and then duplicate 
  if(stretch) { #should we stretch each lineage's pseudotime to have the same length (assume the same maturation real time)
    short_branches <- setdiff(row.names(range_df), as.character(c(progenitor_states, longest_lineage_branch)))
    
    #stretch (condense) the pseudotime from 0 to 100 
    #longest_branch_multiplier <- 100 / (range_df[as.character(longest_lineage_branch), 'max'] - range_df[as.character(progenitor_state), 'min']) 
    longest_branch_range <- range_df[as.character(longest_lineage_branch), 'max'] - min(range_df[, 'min'])
    longest_branch_multiplier <- 100 / longest_branch_range
    
    T_0 <- (max(range_df[as.character(progenitor_states), 'max']) - min(range_df[as.character(progenitor_states), 'min']))  * longest_branch_multiplier
    print (T_0)
    short_branches_multipliers <- (100 - T_0) / (range_df[short_branches , 'max'] - max(range_df[as.character(progenitor_states), 'max']))
    
    heterochrony_constant <- c(longest_branch_multiplier, short_branches_multipliers)
    names(heterochrony_constant) <- c(longest_lineage_branch, short_branches)
    message('stretch ratio are: ')
    print(heterochrony_constant)
    
    pData[pData$State %in% c(progenitor_states, longest_lineage_branch), 'Pseudotime'] <- 
      (pData[pData$State %in% c(progenitor_states, longest_lineage_branch), 'Pseudotime'] - min(range_df[progenitor_states, 'min'])) * longest_branch_multiplier
    print ("&&&&")
    print (pData[pData$State %in% c(progenitor_states, longest_lineage_branch), 'Pseudotime'])
    for(i in 1:length(short_branches)) { #stretch short branches
      print ("***")
      print (i)
      print ("!!!")
      print (pData[pData$State  == short_branches[i], 'Pseudotime'])
      pData[pData$State  == short_branches[i], 'Pseudotime'] <- 
        (pData[pData$State == short_branches[i], 'Pseudotime'] - max(range_df[as.character(progenitor_states), 'max'])) * 
        short_branches_multipliers[i] + T_0
      print(pData[pData$State  == short_branches[i], 'Pseudotime'])
    }
  }
  
  pData$State[progenitor_ind] <- lineage_states[1] #set progenitors to the lineage 1
  for (i in 1:(length(lineage_states) - 1)) { #duplicate progenitors for multiple branches
    exprs_data <- cbind(exprs_data, exprs_data[, progenitor_ind])
    weight_vec <- c(weight_vec, rep(0.5, length( progenitor_ind)))
    
    colnames(exprs_data)[(ncol(exprs_data) - length(progenitor_ind) + 1):ncol(exprs_data)] <- 
      paste('duplicate', lineage_states[i], 1:length(progenitor_ind), sep = '_')
    pData <- rbind(pData, pData[progenitor_ind, ])
    
    pData$State[(length(pData$State) - length(progenitor_ind) + 1):length(pData$State)] <- lineage_states[i + 1]
    row.names(pData)[(length(pData$State) - length(progenitor_ind) + 1):length(pData$State)] <- 
      paste('duplicate', lineage_states[i], 1:length(progenitor_ind), sep = '_')
  }
  pData$State <- as.factor(pData$State)
  
  pData$weight <- weight_vec
  Size_Factor <- pData$Size_Factor
  
  fData <- fData(cds)
  
  cds_subset <- newCellDataSet(as.matrix(exprs_data),
                               phenoData = new("AnnotatedDataFrame", data = pData),
                               featureData = new("AnnotatedDataFrame", data = fData),
                               expressionFamily=negbinomial(),
                               lowerDetectionLimit=1)
  pData(cds_subset)$State <- as.factor(pData(cds_subset)$State)
  pData(cds_subset)$Size_Factor <- Size_Factor
  
  if (cds_subset@expressionFamily@vfamily %in% c("zanegbinomialff", 
                                                 "negbinomial", "poissonff", "quasipoissonff")) {
    integer_expression <- TRUE
  }
  else {
    integer_expression <- FALSE
  }
  
  if (integer_expression) {
    CM <- exprs(cds_subset)
    CM <- t(t(CM) / sizeFactors(cds_subset)) 
    cds_exprs <- reshape2::melt(round(CM))
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
  cds_exprs$feature_label <- factor(cds_exprs$feature_label)
  
  #weight_vec <- rep(1, nrow(pData(cds_subset)))
  #names(weight_vec) <- row.names(pData(cds_subset))

  #weight_vec[row.names(subset(pData(cds), State = progenitor_state))] <- 0.5
  #weight_vec[nrow(pData(cds_subset)):length(pData)] <- 0.5
  print (weight_vec)
  #fit by VGAM: 
  merged_df_with_vgam <- plyr::ddply(cds_exprs, .(feature_label), function(x) {
    fit_res <- tryCatch({
      #expression <- x$adjusted_expression
      
      #df <- data.frame(expression = round(x))
      #df <- cbind(df, pData)
      
      vg <- suppressWarnings(vglm(formula = as.formula(fullModelFormulaStr), 
                                  family = cds_subset@expressionFamily, data = x, 
                                  maxit = 300, checkwz = FALSE,
                                  weights=pData$weight))
      if (integer_expression) {
        res <- predict(vg, type = "response")
        res[res < min_expr] <- min_expr
      }
      else {
        res <- 10^(predict(vg, type = "response"))
        res[res < log10(min_expr)] <- log10(min_expr)
      }
      res
    }, error = function(e) {
      print("Error!")
      print(e)
      res <- rep(NA, nrow(x))
      res
    })
    
    data.frame(Pseudotime=x$Pseudotime, 
               expectation=fit_res,
               State=x$State)
  })
  
  if(method == 'loess')
    cds_exprs$expression <-  cds_exprs$expression + cds@lowerDetectionLimit
  
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
    cds_subset$feature_label <- factor(cds_subset$feature_label, 
                                       levels = panel_order)
  }
  
  #mlt_df <- data.frame(expectation = as.vector(res$expectation), Pseudotime = pData$Pseudotime, State = pData$State, 
  #                     feature_label = rep(as.character(fData$gene_short_name), each = nrow(res))) 
  #mlt_df$State <- as.factor(mlt_df$State)
  
  
  q <- ggplot(aes(Pseudotime, expression), data = cds_exprs)
  if (is.null(color_by) == FALSE) {
    q <- q + geom_point(aes_string(color = I(color_by)), size = I(cell_size))
  }
  # else {
  #     q <- q + geom_point(size = I(cell_size))
  # }
  q <- q + scale_y_log10() + facet_wrap(~feature_label, nrow = nrow, 
                                        ncol = ncol, scales = "free_y")
  
  if(method == 'loess')
    q  <- q + stat_smooth(aes(fill = State, color = State), method = 'loess')
  else if(method == 'fitting') {
    q <- q + geom_line(aes(x=Pseudotime, y=expectation, color = State), data = merged_df_with_vgam)
  }
  
  if(stretch)
    q <- q + ylab("Expression") + xlab("Maturation levels")
  else 
    q <- q + ylab("Expression") + xlab("Pseudotime")
  q <- q + monocle_theme_opts()
  q
}

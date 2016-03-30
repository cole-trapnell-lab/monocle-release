cth_classifier_cds <- function(cds_subset, cth, curr_node, frequency_thresh) {
  #curr_cell_vertex <-  V(cth@classificationTree)[curr_node]
  next_nodes <- c()
  print (unique(pData(cds_subset)$Cluster))
  for (child in V(cth@classificationTree) [ suppressWarnings(nei(curr_node, mode="out")) ]){
    
    child_cell_class_func <- V(cth@classificationTree) [ child ]$classify_func[[1]]
    #type_res <- sparseApply(exprs(cds_subset), 2, child_cell_class_func, convert_to_dense=FALSE)
    type_res <- child_cell_class_func(exprs(cds_subset))
    #print(type_res)
    type_res <- unlist(type_res)
    
    names(type_res) <- row.names(pData(cds_subset))

    if ((sum(type_res) / length(type_res)) > frequency_thresh){
      next_nodes <- c(next_nodes, V(cth@classificationTree) [ child ]$name)
    }
    print (paste(V(cth@classificationTree) [ child ]$name, ":", sum(type_res),  " of ", length(type_res) ))
  }
  print (next_nodes)
  if (length(next_nodes) == 1){
    CellType <- cth_classifier_cds(cds_subset, cth, next_nodes[1], frequency_thresh)
  }else if(length(next_nodes) == 0){
    if (curr_node == "root")
      CellType = "Unknown"
    else
      CellType = curr_node
  }else if(length(next_nodes) > 1){
    CellType = "Ambiguous"
  }else{
    CellType = "Unknown"
  }
  return (CellType)
}

classifyCellsHelperCds <- function(cds_subset, cth, frequency_thresh){
  CellType <- cth_classifier_cds(cds_subset, cth, "root", frequency_thresh)
}


cth_classifier_cell <- function(cell_name, cth, curr_node, gate_res) {
  next_nodes <- c()
  for (child in V(cth@classificationTree) [ suppressWarnings(nei(curr_node, mode="out")) ]){
    type_res <- gate_res[[V(cth@classificationTree) [ child ]$name]]
    #print (class(type_res[cell_name]))
    #print (cell_name)
    if (type_res[cell_name] == TRUE)
      next_nodes <- c(next_nodes, V(cth@classificationTree) [ child ]$name)
  }
  
  if (length(next_nodes) == 1){
    CellType <- cth_classifier_cell(cell_name, cth, next_nodes[1], gate_res)
  }else if(length(next_nodes) == 0){
    if (curr_node == "root")
      CellType = "Unknown"
    else
      CellType = curr_node
  }else if(length(next_nodes) > 1){
    CellType = "Ambiguous"
  }else{
    CellType = "Unknown"
  }
  return (CellType)
}

classifyCellsHelperCell <- function(cds, cth){
  #next_node_list <- rep(list(), ncol(cds)) 
  
  gate_res <- list()
  for (v in V(cth@classificationTree)){
    cell_class_func <- V(cth@classificationTree) [ v ]$classify_func[[1]]
    type_res <- cell_class_func(exprs(cds))
    gate_res[[ V(cth@classificationTree) [ v ]$name]] <- type_res
  }
  cds_pdata <- dplyr::group_by_(dplyr::select_(add_rownames(pData(cds)), "rowname"), "rowname") 
  class_df <- as.data.frame(cds_pdata %>% do(CellType = cth_classifier_cell(.$rowname, cth, "root", gate_res)))
  CellType <- factor(unlist(class_df$CellType))
  names(CellType) <- class_df$rowname
  #CellType <- cth_classifier_cell(cds_subset, cth, "root", gate_res)
  return(CellType)
}

#' Classify cells according to a set of markers
#' @importFrom dplyr add_rownames select_ do group_by_ inner_join
#' @export 
classifyCells <- function(cds, cth, frequency_thresh, cores=1, ...) {
  progress_opts <- options()$dplyr.show_progress
  options(dplyr.show_progress = T)
  if (length(list(...)) > 0){
    cds_pdata <- dplyr::group_by_(dplyr::select_(add_rownames(pData(cds)), "rowname", ...), ...) 
    class_df <- as.data.frame(cds_pdata %>% do(CellType = classifyCellsHelperCds(cds[,.$rowname], cth, frequency_thresh)))
    class_df$CellType <- factor(unlist(class_df$CellType))
  }else{
    type_res <- classifyCellsHelperCell(cds, cth)
    class_df <- data.frame(rowname = names(type_res), CellType = type_res)
    class_df$CellType <- factor(class_df$CellType)
  }
  
  
  
  options(dplyr.show_progress = progress_opts)
  
  pData(cds) <- pData(cds)[!(names(pData(cds)) %in% "CellType")]
  
  #pData(cds)$cell_type <- cds_types
  
  
  pData(cds) <- as.data.frame(suppressMessages(inner_join(add_rownames(pData(cds)), class_df)))
  
  row.names(pData(cds)) <- pData(cds)$rowname
  pData(cds) <- pData(cds)[,-1]
  cds
}   

#' Calculate each gene's specificity for each cell type 
#' @export
calculateMarkerSpecificity <- function(cds, cth){
  markerSpecificityHelper <- function(cds, cth){
    averageExpression <- Matrix::rowMeans(exprs(cds))
    averageExpression <- unlist(averageExpression)
    averageExpression[is.na(averageExpression)] <- 0
    #names(averageExpression) <- row.names(fData(cds))
    return (data.frame(gene_id = row.names(fData(cds)), expr_val=averageExpression))
  }
  progress_opts <- options()$dplyr.show_progress
  options(dplyr.show_progress = T)
  
  cds <- cds[,row.names(subset(pData(cds), CellType %in% c("Unknown", "Ambiguous") == FALSE))]
  cds_pdata <- dplyr::group_by_(dplyr::select_(add_rownames(pData(cds)), "rowname", "CellType"), "CellType") 
  class_df <- as.data.frame(cds_pdata %>% do(markerSpecificityHelper(cds[,.$rowname], cth)))
  class_df <- dcast(class_df, CellType ~ gene_id, value.var = "expr_val")
  row.names(class_df) <- class_df$CellType
  class_df <- class_df[,-1]
  class_df <- t(as.matrix(class_df))
  
  marker_specificities <- lapply(1:ncol(class_df), function(cell_type_i){
    perfect_specificity <- rep(0.0, ncol(class_df))
    perfect_specificity[cell_type_i] <- 1.0
    apply(class_df, 1, function(x) { 
      if (sum(x) > 0) 1 - JSdistVec(makeprobsvec(x), perfect_specificity)
      else 0
    })
  })
  marker_specificities <- t(do.call(rbind, marker_specificities))
  colnames(marker_specificities) <- colnames(class_df)
  marker_specificities <- melt(marker_specificities)
  colnames(marker_specificities) <- c("gene_id", "CellType", "specificity")
  return (marker_specificities)
}

#' Select the most cell type specific markers
#' @importFrom dplyr top_n
#' @export
selectTopMarkers <- function(marker_specificities, num_markers = 10){
  as.data.frame(marker_specificities %>%
    group_by_("CellType") %>%
    top_n(n = num_markers, wt = specificity))
}

#' Test genes for cell type-dependent expression
#' 
#' @param cds A CellDataSet object containing cells to classify
#' @param cth The CellTypeHierarchy object to use for classification
#' @param residualModelFormulaStr A model formula string specify effects you
#' want to exclude when testing for cell type dependent expression
#' @param cores The number of cores to use when testing
#' @return A table of differential expression test results
#' @importFrom stringr str_replace_all
#' @importFrom dplyr sample_n
#' @export 
markerDiffTable <- function (cds, cth, residualModelFormulaStr="~1", balanced=TRUE, verbose=FALSE, cores=1) {
  
  if (verbose)
    message("Classifying cells according to markers")
  cds <- classifyCells(cds, cth, 0.05, cores=cores)
  cds <- cds[,pData(cds)$CellType %in% c("Unknown", "Ambiguous") == FALSE]
  
  if (balanced){
    cell_type_counts <- table(pData(cds)$CellType)
    cell_type_counts <- cell_type_counts[cell_type_counts > 0]
    least_frequent_type <- which(cell_type_counts == min(cell_type_counts))
    least_frequent_type <- names(cell_type_counts)[least_frequent_type]
    n_cells <- cell_type_counts[least_frequent_type]
    
    message(paste("Least frequent cell type is '", least_frequent_type, "', randomly selecting ", n_cells, " cells for marker identification test", sep=""))
    selected_cells <- c()

    for (cell_type in names(cell_type_counts)){
      cell_type_sample <- sample_n(add_rownames(subset(pData(cds), CellType == cell_type)), n_cells)$rowname
      selected_cells <- c(selected_cells, cell_type_sample)
    }
    
    cds <- cds[,selected_cells]
    # if(is.null(max_cells) == FALSE && max_cells < ncol(cds)) 
    #   selected_cells <- sample(ncol(cds), max_cells)
    # else
    #   selected_cells <- colnames(cds)
  }
  
  fullModelFormulaStr <- paste("CellType")
  fullModelFormulaStr <- paste("~", fullModelFormulaStr,sep = "")
  if (residualModelFormulaStr != "~1"){
    residual_terms <- str_replace_all(residualModelFormulaStr, "~", "")
    fullModelFormulaStr <- paste(fullModelFormulaStr, residual_terms, sep = "+")
  }

  if (verbose)
    message("Testing for marker-dependent expression")
  
  marker_diff <- differentialGeneTest(cds, 
                                      fullModelFormulaStr=fullModelFormulaStr,
                                      reducedModelFormulaStr=residualModelFormulaStr,
                                      verbose=T,
                                      cores=cores)
  
  return(marker_diff)
}

#' Creates a new CellTypeHierarcy object
#' @return A new CellTypeHierarchy object
#' @export
newCellTypeHierarchy <- function()
{
  cth <- new( "CellTypeHierarchy",
              classificationTree = graph.empty())
  
  root_node_id <- "root"
  
  cth@classificationTree <- cth@classificationTree + vertex(root_node_id, classify_func=list(function(x) TRUE))
  #cth@classificationTree %>% add_vertices(1, name = root_node_id, "classify_func"=list(function(x) TRUE))
  return(cth)
}

#' Adds a new cell type to a CellTypeHierarchy
#' @param cth The CellTypeHierarchy object
#' @param cell_type_name The name of the new cell type. Can't already exist in cth
#' @param classify_func A function that returns true when a cell is of the new type
#' @param parent_cell_typ_name If this cell type is a subtype of another, provide its name here
#' @return The updated CellTypeHierarchy
#' @export
addCellType <- function(cth, cell_type_name, classify_func, parent_cell_type_name="root") 
{
  if (cell_type_name %in% V(cth@classificationTree)$name){
    stop(paste("Error: cell type",cell_type_name, "already exists."))
  }
  
  # TODO: verify that classify_func has the right signature/call semantics?
  cth@classificationTree <- cth@classificationTree + vertex(cell_type_name, classify_func=list(classify_func))
  
  cth@classificationTree <- cth@classificationTree + edge(parent_cell_type_name, cell_type_name)
  return (cth)
}
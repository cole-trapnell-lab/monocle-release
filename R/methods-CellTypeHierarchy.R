#' Classify cells according to a set of markers
#' @export 
classifyCells <- function(cds, cth, frequency_thresh, ...) {
  
  classifyCellsHelper <- function(cds_subset, cth){
    
    cth_classifier <- function(cth, curr_node) {
      #curr_cell_vertex <-  V(cth@classificationTree)[curr_node]
      next_nodes <- c()
      for (child in V(cth@classificationTree) [ suppressWarnings(nei(curr_node, mode="out")) ]){
        
        child_cell_class_func <- V(cth@classificationTree) [ child ]$classify_func[[1]]
        type_res <- smartEsApply(cds_subset, 2, child_cell_class_func, convert_to_dense=FALSE)
        #print(type_res)
        type_res <- unlist(type_res)
        
        #type_res <- as.matrix(do.call(rbind, type_res))
        names(type_res) <- row.names(pData(cds_subset))
        #print(head(type_res))
        #print ("----")
        if (sum(type_res) / length(type_res) > frequency_thresh){
          next_nodes <- c(next_nodes, V(cth@classificationTree) [ child ]$name)
        }
      }
      
      if (length(next_nodes) == 1){
        CellType <- cth_classifier(cth, next_nodes[1])
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
    
    CellType <- cth_classifier(cth, "root")
    
  }
  progress_opts <- options()$dplyr.show_progress
  options(dplyr.show_progress = F)
  if (length(list(...)) > 0){
    cds_pdata <- dplyr::group_by_(dplyr::select_(add_rownames(pData(cds)), "rowname", ...), ...) 
  }else{
    cds_pdata <- dplyr::group_by(dplyr::select(add_rownames(pData(cds)), rowname), rowname)
  }
  
  class_df <- as.data.frame(cds_pdata %>% do(CellType = classifyCellsHelper(cds[,.$rowname], cth)))
  class_df$CellType <- factor(unlist(class_df$CellType))
  options(dplyr.show_progress = progress_opts)
  
  pData(cds) <- pData(cds)[!(names(pData(cds)) %in% "CellType")]
  
  #pData(cds)$cell_type <- cds_types
  
  
  pData(cds) <- as.data.frame(suppressMessages(inner_join(add_rownames(pData(cds)), class_df)))
  
  row.names(pData(cds)) <- pData(cds)$rowname
  pData(cds) <- pData(cds)[,-1]
  cds
}   

#' Test genes for cell type-dependent expression
#' 
#' @param cds A CellDataSet object containing cells to classify
#' @param cth The CellTypeHierarchy object to use for classification
#' @param residualModelFormulaStr A model formula string specify effects you
#' want to exclude when testing for cell type dependent expression
#' @param cores The number of cores to use when testing
#' @return A table of differential expression test results
#' @export 
markerDiffTable <- function (cds, cth, residualModelFormulaStr="~1", max_cells=NULL, cores=1) {
  
  cds <- classifyCells(cds, cth, 0.05)
  cds <- cds[,pData(cds)$CellType %in% c("Unknown", "Ambiguous") == FALSE]
  
  fullModelFormulaStr <- paste("CellType")
  fullModelFormulaStr <- paste("~", fullModelFormulaStr,sep = "")
  if (residualModelFormulaStr != "~1"){
    residual_terms <- str_replace_all(residualModelFormulaStr, "~", "")
    fullModelFormulaStr <- paste(fullModelFormulaStr, residual_terms, sep = "+")
  }
  
  if(max_cells < ncol(cds)) 
    selected_cells <- sample(ncol(cds), max_cells)
  else
    selected_cells <- colnames(cds)
  
  marker_diff <- differentialGeneTest(cds[,selected_cells], 
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
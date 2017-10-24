
#' @import methods
#' @importFrom Biobase exprs pData
#' @importFrom igraph V
cth_classifier_cds <- function(cds_subset, cth, curr_node, frequency_thresh) {
  #curr_cell_vertex <-  V(cth@classificationTree)[curr_node]
  next_nodes <- c()
  #print (unique(pData(cds_subset)$Cluster))
  for (child in V(cth@classificationTree) [ suppressWarnings(nei(curr_node, mode="out")) ]){
    
    child_cell_class_func <- V(cth@classificationTree) [ child ]$classify_func[[1]]
    #type_res <- sparseApply(exprs(cds_subset), 2, child_cell_class_func, convert_to_dense=FALSE)
    type_res <- child_cell_class_func(exprs(cds_subset))
    #print(type_res)
    type_res <- unlist(type_res)
    
    names(type_res) <- row.names(pData(cds_subset))
    cell_type_name <- V(cth@classificationTree) [ child ]$name
    if (length(frequency_thresh) > 1)
      required_thresh <- frequency_thresh[cell_type_name]
    else
      required_thresh <- frequency_thresh
    if ((sum(type_res) / length(type_res)) > frequency_thresh){
      next_nodes <- c(next_nodes, cell_type_name)
    }
    #print (paste(V(cth@classificationTree) [ child ]$name, ":", sum(type_res),  " of ", length(type_res) ))
  }
  
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

#' @importFrom igraph V
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

#' @importFrom Biobase exprs pData
#' @importFrom igraph V
#' @importFrom dplyr %>%
classifyCellsHelperCell <- function(cds, cth){
  #next_node_list <- rep(list(), ncol(cds)) 
  
  gate_res <- list()
  for (v in V(cth@classificationTree)){
    cell_class_func <- V(cth@classificationTree) [ v ]$classify_func[[1]]
    type_res <- cell_class_func(exprs(cds))
    gate_res[[ V(cth@classificationTree) [ v ]$name]] <- type_res
  }
  cds_pdata <- dplyr::group_by_(dplyr::select_(rownames_to_column(pData(cds)), "rowname"), "rowname") 
  class_df <- as.data.frame(cds_pdata %>% do(CellType = cth_classifier_cell(.$rowname, cth, "root", gate_res)))
  CellType <- factor(unlist(class_df$CellType))
  names(CellType) <- class_df$rowname
  #CellType <- cth_classifier_cell(cds_subset, cth, "root", gate_res)
  return(CellType)
}

#' @title Classify cells according to a set of markers
#' 
#' @description Creates a CellTypeHierarchy object which can store
#' cell types with the addCellType() function. When classifyCells
#' is used with a CellDataSet and a CellTypeHierarchy cells in the 
#' CellDataSet can be classified as cell types found in the CellTypeHierarchy
#' 
#' @details CellTypeHierarchy objects are Monocle's mechanism for
#'   classifying cells into types based on known markers. To classify the cells
#'   in a CellDataSet object according to known markers, first construct a
#'   CellTypeHierachy with \code{newCellTypeHierarchy()} and 
#'   \code{addCellType()} and then provide both the \code{CellDataSet}
#'   and the \code{CellTypeHierachy} to \code{classifyCells()}. Each
#'   call to \code{addCellType()} registers a classification function
#'   that accepts the expression data from a CellDataSet object as input, and
#'   returns a boolean vector indicating whether each cell is of the given type.
#'   When you call \code{classifyCells()}, each cell will be checked against the classification functions in the
#'   \code{CellTypeHierachy}.  If you wish to make a cell type a subtype of
#'   another that's already been registered with a CellTypeHierarchy object,
#'   make that one the "parent" type with the \code{cell_type_name} argument. If
#'   you want two types to be mutually exclusive, make them "siblings" by giving
#'   them the same parent. The classifcation functions in a CellTypeHierarchy must take a single argument, a matrix of
#'   expression values, as input. Note that this matrix could either be a 
#'   \code{\link[Matrix]{sparseMatrix}} or a dense matrix. Explicitly casting the input to a dense
#'   matrix inside a classification function is likely to drastically slow down 
#'   classifyCells and other routines that use CellTypeHierarhcy objects.
#'   Successive calls to \code{addCellType} build up a tree of classification
#'   functions inside a CellTypeHierarchy. When two functions are siblings in 
#'   the tree, classifyCells expects that a cell will meet the classification
#'   criteria for at most one of them. For example, you might place 
#'   classification functions for T cells and B cells as siblings, because
#'   a cell cannot be both of these at the same time. When a cell meets the 
#'   criteria for more than one function, it will be tagged as "Ambiguous". If
#'   \code{classifyCells} reports a large number of ambiguous cells, consider
#'   adjusting your classification functions. For example, some cells are 
#'   defined by very high expression of a key gene that is expressed at lower
#'   levels in other cell types. Raising the threshold for this gene in a 
#'   classification could resolve the ambiguities. A classification function
#'   can also have child functions. You can use this to specify subtypes of 
#'   cells. For example, T cells express the gene CD3, and there are many
#'   subtypes. You can encode each subset by first adding a general T cell
#'   classification function that recognizes CD3, and then adding an additional
#'   function that recognizes CD4 (for CD4+ helper T cells), one for CD8 (to
#'   identify CD8+ cytotoxic T cells), and so on. \code{classifyCells} will
#'   aim to assign each cell to its most specific subtype in the "CellType" 
#'   column. By default, \code{classifyCells} applies the classification functions to
#'   individual cells, but you can also apply it to cells in a "grouped" mode to 
#'   impute the type of cells that are missing expression of your known markers.
#'   You can specify additional (quoted) grouping variables to \code{classifyCells}.
#'   The function will group the cells according to these factors, and then 
#'   classify the cells. It will compute the frequency of each cell type in each
#'   group, and if a cell type is present at the frquency specified in 
#'   \code{frequency_thresh}, all the cells in the group are classified as that 
#'   type. If group contains more one cell type at this frequency, all the cells
#'   are marked "Ambiguous". This allows you to impute cell type based on 
#'   unsupervised clustering results (e.g. with \code{\link{clusterCells}()}) or
#'   some other grouping criteria.
#' 
#' 
#' @return \code{newCellTypeHierarchy} and \code{addCellType} both return an 
#'   updated CellTypeHierarchy object. \code{classifyCells} returns an updated 
#'   \code{CellDataSet} with a new column, "CellType", in the pData table.
#'   
#' @importFrom igraph vertex graph.empty
#'   
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

#' Add a new cell type
#' @description adds a cell type to a pre-existing CellTypeHierarchy and produces a function that accepts
#' expression data from a CellDataSet. When the function is called on a CellDataSet a boolean vector is returned
#' that indicates whether each cell is or is not the cell type that was added by addCellType.
#' @param cth The CellTypeHierarchy object 
#' @param cell_type_name The name of the new cell type. Can't already exist in
#'   cth
#' @param classify_func A function that returns true when a cell is of the new
#'   type
#' @param parent_cell_type_name If this cell type is a subtype of another,
#'   provide its name here
#'   
#' @importFrom igraph V edge
#'   
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

#' @title Classify cells according to a set of markers
#' 
#' @description classifyCells accepts a cellDataSet and and a cellTypeHierarchy.
#' Each cell in the cellDataSet is checked against the functions in the cellTypeHierarchy
#' to determine each cell's type
#' 
#' @describeIn newCellTypeHierarchy Add a cell type to a CellTypeHierarchy
#' @param cds The CelllDataSet you want to classify
#' @param ... character strings that you wish to pass to dplyr's group_by_ routine
#' @param enrichment_thresh fraction to be multipled by each cell type percentage. Only used if frequency_thresh is NULL, both cannot be NULL
#' @param frequency_thresh If at least this fraction of group of cells meet a cell types marker criteria, impute them all to be of that type.  
#' @importFrom dplyr select_ do group_by_ inner_join %>%
#' @importFrom tibble rownames_to_column
#' @importFrom Biobase pData pData<-
#' @export 
#' @examples
#' \dontrun{
#' # Initialize a new CellTypeHierachy
#' 
#' # Register a set of classification functions. There are multiple types of T cells
#' # A cell cannot be both a B cell and a T cell, a T cell and a Monocyte, or
#' # a B cell and a Monocyte.
#' cth <- newCellTypeHierarchy()
#' 
#' cth <- addCellType(cth, "T cell", 
#'                    classify_func=function(x) {x["CD3D",] > 0})
#'                    
#' cth <- addCellType(cth, "CD4+ T cell", 
#'                    classify_func=function(x) {x["CD4",] > 0}, 
#'                    parent_cell_type_name = "T cell")
#'                    
#' cth <- addCellType(cth, "CD8+ T cell", 
#'                    classify_func=function(x) {
#'                      x["CD8A",] > 0 | x["CD8B",] > 0
#'                    }, 
#'                    parent_cell_type_name = "T cell")
#'                    
#' cth <- addCellType(cth, "B cell", 
#'                    classify_func=function(x) {x["MS4A1",] > 0})
#'                    
#' cth <- addCellType(cth, "Monocyte", 
#'                    classify_func=function(x) {x["CD14",] > 0})
#' 
#' # Classify each cell in the CellDataSet "mix" according to these types
#' mix <- classifyCells(mix, cth)
#'
#' # Group the cells by the pData table column "Cluster". Apply the classification
#' functions to the cells groupwise. If a group is at least 5% of a type, make
#' them all that type. If the group is 5% one type, and 5% a different, mutually
#' exclusive type, mark the whole cluster "Ambiguous"
#' mix <- classifyCells(mix, Cluster, 0.05)
#' }
#' 
classifyCells <- function(cds, cth, frequency_thresh=NULL, enrichment_thresh=NULL, ...) {
  progress_opts <- options()$dplyr.show_progress
  options(dplyr.show_progress = F)
  if (length(list(...)) > 0){
    if (is.null(enrichment_thresh) && is.null(frequency_thresh))
      stop("Error: to use classifyCells in grouped mode, you must also set frequency_thresh")
    
    cds <- classifyCells(cds, cth)
    if (is.null(frequency_thresh)){
      frequency_thresholds <- prop.table(table(pData(cds)$CellType))
      frequency_thresholds <- frequency_thresholds * enrichment_thresh
      frequency_thresholds <- unlist(lapply(frequency_thresholds, min, 1.0))
    }else
      frequency_thresholds <- frequency_thresh

    cds_pdata <- dplyr::group_by_(dplyr::select_(rownames_to_column(pData(cds)), "rowname", ...), ...) 
    class_df <- as.data.frame(cds_pdata %>% dplyr::do(CellType = classifyCellsHelperCds(cds[,.$rowname], cth, frequency_thresh)))
    class_df$CellType <-  as.character(unlist(class_df$CellType))
    #class_df$rowname <- as.character(class_df$rowname)
  }else{
    type_res <- classifyCellsHelperCell(cds, cth)
    class_df <- data.frame(rowname = names(type_res), CellType = type_res)
    class_df$CellType <- as.character(class_df$CellType)
    class_df$rowname <- as.character(class_df$rowname)
  }
  
  options(dplyr.show_progress = progress_opts)
  
  pData(cds) <- pData(cds)[!(names(pData(cds)) %in% "CellType")]
  
  #pData(cds)$cell_type <- cds_types
  
  
  pData(cds) <- as.data.frame(suppressMessages(inner_join(rownames_to_column(pData(cds)), class_df)))
  
  pData(cds)$CellType <- factor(pData(cds)$CellType)
  
  row.names(pData(cds)) <- pData(cds)$rowname
  pData(cds) <- pData(cds)[,-1]
  cds
}   

#' @describeIn newCellTypeHierarchy Calculate each gene's specificity for each cell type 
#' 
#' Computes the Jensen-Shannon distance between the distribution of a gene's 
#' expression across cells and a hypothetical gene that is perfectly restricted
#' to each cell type. The Jensen-Shannon distance is an information theoretic
#' metric between two probability distributions. It is a widely accepted measure
#' of cell-type specificity. For a complete description see Cabili \emph{et. al},
#' Genes & Development (2011). 
#' 
#' @param cth CellTypeHierarchy
#' @param remove_ambig a boolean that determines if ambiguous cells should be removed
#' @param remove_unknown a boolean that determines whether unknown cells should be removed
#' @return For a CellDataset with N genes, and a CellTypeHierarchy with k types,
#' returns a dataframe with N x k rows. Each row contains a gene and a specifity
#' score for one of the types.
#' @importFrom reshape2 dcast
#' @importFrom dplyr %>%
#' @importFrom Biobase exprs fData pData
#' @export
calculateMarkerSpecificity <- function(cds, cth, remove_ambig=TRUE, remove_unknown=TRUE){
  
  if(class(cds)[1] != "CellDataSet") {
    stop("Error cds is not of type 'CellDataSet'")
  }
  
  if(class(cth)[1] != "CellTypeHierarchy") {
    stop("Error cth is not of type 'CellTypeHierarchy'")
  }
  
  CellType <- NA
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
  cds_pdata <- dplyr::group_by_(dplyr::select_(rownames_to_column(pData(cds)), "rowname", "CellType"), "CellType") 
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
  marker_specificities$gene_id <- as.character(marker_specificities$gene_id)
  return (marker_specificities)
}

#' Select the most cell type specific markers
#' 
#' This is a handy wrapper function around dplyr's top_n function to extract
#' the most specific genes for each cell type. Convenient, for example, for 
#' selecting a balanced set of genes to be used in semi-supervised clustering 
#' or ordering.
#' 
#' @param marker_specificities The dataframe of specificity results produced by \code{\link{calculateMarkerSpecificity}()}
#' @param num_markers The number of markers that will be shown for each cell type
#' @return A data frame of specificity results
#' @importFrom dplyr top_n %>%
#' @export
selectTopMarkers <- function(marker_specificities, num_markers = 10){
  specificity <- NA
  as.data.frame(marker_specificities %>%
    group_by_("CellType") %>%
    top_n(n = num_markers, wt = specificity))
}

#' Test genes for cell type-dependent expression
#' 
#' @description takes a CellDataSet and a CellTypeHierarchy and classifies all cells into types passed
#' functions passed into the CellTypeHierarchy. The function will remove all "Unknown" and "Ambiguous" types
#' before identifying genes that are differentially expressed between types.
#' 
#' @param cds A CellDataSet object containing cells to classify
#' @param cth The CellTypeHierarchy object to use for classification
#' @param residualModelFormulaStr A model formula string specify effects you
#' want to exclude when testing for cell type dependent expression
#' @param balanced Whether to downsample the cells so that there's an equal number of each type prior to performing the test
#' @param verbose Whether to emit verbose output during the the search for cell-type dependent genes
#' @param cores The number of cores to use when testing
#' @param reclassify_cells a boolean that indicates whether or not the cds and cth should be run through classifyCells again
#' @param remove_ambig a boolean that indicates whether or not ambiguous cells should be removed the cds
#' @param remove_unknown a boolean that indicates whether or not unknown cells should be removed from the cds
#' @return A table of differential expression test results
#' @importFrom stringr str_replace_all
#' @importFrom dplyr sample_n
#' @importFrom Biobase pData pData<-
#' @export 
markerDiffTable <- function (cds, cth, residualModelFormulaStr="~1", balanced=FALSE, reclassify_cells=TRUE, remove_ambig=TRUE, remove_unknown=TRUE, verbose=FALSE, cores=1) {
  if(class(cds)[1] != "CellDataSet") {
    stop("Error cds is not of type 'CellDataSet'")
  }
  
  if(class(cth)[1] != "CellTypeHierarchy") {
    stop("Error cth is not of type 'CellTypeHierarchy'")
  }
  
  CellType <- NULL
  if (verbose)
    message("Classifying cells according to markers")
  if (reclassify_cells)
    cds <- classifyCells(cds, cth, 0.05)
  if (remove_ambig)
    cds <- cds[,pData(cds)$CellType %in% c("Ambiguous") == FALSE]
  if (remove_unknown)
    cds <- cds[,pData(cds)$CellType %in% c("Unknown") == FALSE]
  pData(cds)$CellType <- droplevels(pData(cds)$CellType)
  
  if (balanced){
    cell_type_counts <- table(pData(cds)$CellType)
    cell_type_counts <- cell_type_counts[cell_type_counts > 0]
    least_frequent_type <- which(cell_type_counts == min(cell_type_counts))
    least_frequent_type <- names(cell_type_counts)[least_frequent_type]
    n_cells <- cell_type_counts[least_frequent_type]
    
    message(paste("Least frequent cell type is '", least_frequent_type, "', randomly selecting ", n_cells, " cells for marker identification test", sep=""))
    selected_cells <- c()

    for (cell_type in names(cell_type_counts)){
      cell_type_sample <- sample_n(rownames_to_column(subset(pData(cds), CellType == cell_type)), n_cells)$rowname
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
    fullModelFormulaStr <- paste(fullModelFormulaStr, residual_terms, sep = " + ")
  }

  if (verbose)
    message("Testing for marker-dependent expression")
  
  marker_diff <- differentialGeneTest(cds, 
                                      fullModelFormulaStr=fullModelFormulaStr,
                                      reducedModelFormulaStr=residualModelFormulaStr,
                                      verbose=verbose,
                                      cores=cores)
  
  return(marker_diff)
}

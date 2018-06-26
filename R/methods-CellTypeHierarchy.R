
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

count_cell_conflicts <- function(cth, gate_res){
  conflicts = data.frame()
  for (i in seq(2, length(gate_res) - 1)){ # Skip the root
    for (j in seq(i + 1, length(gate_res))){
      if (i < j){
        num_conflicts = sum(gate_res[[i]] & gate_res[[j]])
        conflicts = rbind(conflicts, data.frame("cell_type_1"=names(gate_res)[i], "cell_type_2"=names(gate_res)[j], "conflicts"=num_conflicts))
        
      }
    }
  }
  return(conflicts)
}

# 
# cth_classifier_cell <- function(cell_name, cth, curr_node, gate_res, max_depth = NULL) {
#   next_nodes <- c()
#   for (child in V(cth@classificationTree) [ suppressWarnings(nei(curr_node, mode="out")) ]){
#     type_res <- gate_res[[V(cth@classificationTree) [ child ]$name]]
#     if (V(cth@classificationTree)[curr_node]$name != "root")
#     {
#       type_res <- gate_res[[ V(cth@classificationTree)[curr_node]$name]] & type_res 
#     }
#     #print (class(type_res[cell_name]))
#     #print (cell_name)
#     if (type_res[cell_name] == TRUE)
#       next_nodes <- c(next_nodes, V(cth@classificationTree) [ child ]$name)
#   }
#   
#   if (length(next_nodes) == 1){
#     if (is.null(max_depth)){
#       CellType <- cth_classifier_cell(cell_name, cth, next_nodes[1], gate_res)
#     }else if(max_depth > 1){
#       CellType <- cth_classifier_cell(cell_name, cth, next_nodes[1], gate_res, max_depth=max_depth - 1)
#     }else{
#       CellType = V(cth@classificationTree)[next_nodes[1]]$name
#     }
#     
#   }else if(length(next_nodes) == 0){
#     if (V(cth@classificationTree)[curr_node]$name == "root")
#       CellType = "Unknown"
#     else
#       CellType = V(cth@classificationTree)[curr_node]$name
#   }else if(length(next_nodes) > 1){
#     CellType = "Ambiguous"
#   }else{
#     CellType = "Unknown"
#   }
#   return (CellType)
# }
cth_classifier_cell <- function(cth, gate_res, curr_node=1, max_depth = NULL) {
  
  ambiguities = rep(FALSE, length(gate_res[[2]]))
  names(ambiguities) = row.names(gate_res[[2]])
  
  # First, mark off the ambiguous cells:
  for (v in V(cth@classificationTree)){
    if (length(V(cth@classificationTree) [suppressWarnings(nei(v, mode="out")) ]) > 1){
      type_res <- gate_res[ V(cth@classificationTree) [ suppressWarnings(nei(v, mode="out")) ]$name]
      mat <- do.call("cBind",type_res)
      #row.names(mat) = gate_res[[1]]
      
      tryCatch({ambiguities[Matrix::rowSums(mat) > 1] = TRUE}, error = function(e) {})
      
    }
  }
 
  
  assignments = rep("Unknown", length(gate_res[[2]]))
  names(assignments) = row.names(gate_res[[2]])
  
  fill_in_assignments <- function(curr_assignments, cth, v, gate_res, max_depth=NULL){
    for (child in V(cth@classificationTree) [ suppressWarnings(nei(v, mode="out")) ]){
      # type_res <- gate_res[[V(cth@classificationTree) [ child ]$name]]
      # if (V(cth@classificationTree)[v]$name != "root")
      # {
      #   type_res <- gate_res[[ V(cth@classificationTree)[v]$name]] & type_res 
      # }
      parents = V(cth@classificationTree)[igraph::all_simple_paths(cth@classificationTree, v, to = child, mode = "out")[[1]]]$name
      parents = setdiff(parents, "root")
      if (length(intersect(parents, names(gate_res))) > 0){
        type_res <- gate_res[parents]
        if (length(type_res) > 1){
          mat <- do.call("cBind",type_res)
          type_res = apply(mat, 1, function(x) { prod(x) })
        }else{
          type_res = type_res[[1]]
        }
        
        curr_assignments[which(type_res == TRUE)] = V(cth@classificationTree) [ child ]$name
        if(is.null(max_depth) == FALSE){
          if (max_depth > 1){
            curr_assignments = fill_in_assignments(curr_assignments, cth, child, gate_res, max_depth-1)
          }
        }else{
          curr_assignments = fill_in_assignments(curr_assignments, cth, child, gate_res)
        }
      }
    }
    return (curr_assignments)
  }
  
  assignments = fill_in_assignments(assignments, cth, curr_node, gate_res, max_depth)
  assignments[which(ambiguities == TRUE)] = "Ambiguous"
  return(assignments)
  # next_nodes <- c()
  # for (child in V(cth@classificationTree) [ suppressWarnings(nei(curr_node, mode="out")) ]){
  #   type_res <- gate_res[[V(cth@classificationTree) [ child ]$name]]
  #   if (V(cth@classificationTree)[curr_node]$name != "root")
  #   {
  #     type_res <- gate_res[[ V(cth@classificationTree)[curr_node]$name]] & type_res 
  #   }
  #   #print (class(type_res[cell_name]))
  #   #print (cell_name)
  #   if (type_res[cell_name] == TRUE)
  #     next_nodes <- c(next_nodes, V(cth@classificationTree) [ child ]$name)
  # }
  # 
  # if (length(next_nodes) == 1){
  #   if (is.null(max_depth)){
  #     CellType <- cth_classifier_cell(cell_name, cth, next_nodes[1], gate_res)
  #   }else if(max_depth > 1){
  #     CellType <- cth_classifier_cell(cell_name, cth, next_nodes[1], gate_res, max_depth=max_depth - 1)
  #   }else{
  #     CellType = V(cth@classificationTree)[next_nodes[1]]$name
  #   }
  #   
  # }else if(length(next_nodes) == 0){
  #   if (V(cth@classificationTree)[curr_node]$name == "root")
  #     CellType = "Unknown"
  #   else
  #     CellType = V(cth@classificationTree)[curr_node]$name
  # }else if(length(next_nodes) > 1){
  #   CellType = "Ambiguous"
  # }else{
  #   CellType = "Unknown"
  # }
  # return (CellType)
}


#' @importFrom Biobase exprs pData
#' @importFrom igraph V
#' @importFrom dplyr %>%
classifyCellsHelperCell <- function(cds, cth){
  #next_node_list <- rep(list(), ncol(cds)) 
  
  gate_res <- list()
  for (v in V(cth@classificationTree)){
    cell_class_func <- V(cth@classificationTree) [ v ]$classify_func[[1]]
  
    parent <- environment(cell_class_func)
    if (is.null(parent))
      parent <- emptyenv()
    e1 <- new.env(parent=parent)
    multiassign(names(pData(cds)), pData(cds), envir=e1)
    environment(cell_class_func) <- e1
    
    type_res <- cell_class_func(exprs(cds))
    if (length(type_res)!= ncol(cds)){
      message(paste("Error: classification function for", V(cth@classificationTree) [ v ]$name, "returned a malformed result."))
      stop()
    }
    type_res = as(as(type_res,"sparseVector"), "sparseMatrix")
    row.names(type_res) = row.names(pData(cds))
    colnames(type_res) =V(cth@classificationTree) [ v ]$name
    gate_res[[ V(cth@classificationTree) [ v ]$name]] <- type_res
  }
  
  CellType = cth_classifier_cell(cth, gate_res)
  conflicts = count_cell_conflicts(cth, gate_res)
  cds@auxClusteringData[["class_conflicts"]] = conflicts
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
  
  cth@classificationTree <- cth@classificationTree + vertex(root_node_id, classify_func=list(function(x) {rep(TRUE, ncol(x))}))
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


#' @title classifyCellsHelperCellGlmNet
#' @description Description of classifyCellsHelperCellGlmNet -- to be added 
#' @param cds CellDataSet containing cells that will be clustered
#' @param cth CellTypeHierarchy that dictates cell types present and requirements of cell to be considered  a certain cell type
#' @importFrom glmnet cv.glmnet
#' @export
classifyCellsHelperCellGlmNet <- function(cds, cth){
  
  cds_sub = cds[,pData(cds)$CellType %in% c("Unknown", "Ambiguous") == FALSE]
  cds_sub = cds_sub[fData(cds_sub)$use_for_ordering,]
  test_cells = row.names(pData(cds_sub))[sample(ncol(cds_sub), round(0.2 * ncol(cds_sub)))]
  training_cells = setdiff(row.names(pData(cds_sub)), test_cells)
  cds_sub = cds_sub[,training_cells]
  
  y = droplevels(pData(cds_sub)$CellType)
  x = t(t(exprs(cds_sub)) / pData(cds_sub)$Size_Factor)
  
  x@x = log(x@x[x@x > 0] + 1)
  x = t(x)
  cvfit = cv.glmnet(x, y, family = "multinomial")
  
  #selected_lambda = cvfit$lambda.min
  feature_genes = coef(cvfit, s = "lambda.min")
  
  test_cds = cds[fData(cds)$use_for_ordering, test_cells]
  test_x = t(t(exprs(test_cds)) / pData(test_cds)$Size_Factor)
  test_x@x = log(test_x@x[test_x@x > 0] + 1)
  test_x = t(test_x)
  predictions = predict(cvfit, newx = test_x, s = "lambda.min", type = "class")
  
  #fData(cds_sub)[colnames(x)[which(feature_genes[["Myoblast"]][,1] != 0) + 1],]$gene_short_name
  return(cds)
}

# Train a multinomial classifier at each internal node of a CellTypeHierarchy
#' @importFrom igraph V
#' @import doParallel
cth_train_glmnet <- function(cds, cth, curr_node, gate_res, rank_prob_ratio = 2, min_observations = 8, max_training_samples = 10000, cores = 1) {
  message(paste("Classifying children of", V(cth@classificationTree)[curr_node]$name))
  child_cell_types = V(cth@classificationTree) [ suppressWarnings(nei(curr_node, mode="out")) ]$name
  
  cds_pdata <- dplyr::group_by_(dplyr::select_(rownames_to_column(pData(cds)), "rowname"), "rowname") 
  ctf_cell_type = cth_classifier_cell(cth, gate_res, curr_node, max_depth=1)
  conflicts = count_cell_conflicts(cth, gate_res)
  cds@auxClusteringData[["class_conflicts"]] = conflicts
  outgroup_samples =  ctf_cell_type %in% child_cell_types == FALSE & 
                      gate_res[[V(cth@classificationTree)[curr_node]$name]][,1] == FALSE &
                      ctf_cell_type != "Ambiguous" 
  if (sum(outgroup_samples) > 0){
    outgroup_samples = ctf_cell_type[outgroup_samples]
  }else{
    outgroup_samples = c()
  }
  
  ctf_cell_type = ctf_cell_type[ctf_cell_type %in% child_cell_types]
  #ctf_cell_type[ctf_cell_type %in% child_cell_types == FALSE ]
  if (length(ctf_cell_type) > max_training_samples){
    message("Downsampling training data.")
    obs_counts = table(ctf_cell_type)
    #obs_prob = 1/obs_counts
    target_obs_per_cell_type =  ceiling(max_training_samples / length(obs_counts))
    training_sample = c()
    for(i in names(obs_counts)){
      num_obs_for_type_i = min(target_obs_per_cell_type, obs_counts[i])
      obs_for_type_i = sample(which(ctf_cell_type == i), num_obs_for_type_i)
      training_sample = append(training_sample, obs_for_type_i)
    }
    training_sample = ctf_cell_type[training_sample]
  }else{
    training_sample = ctf_cell_type
  }
  
  if (length(outgroup_samples) > 0){
    num_outgroup = min(max_training_samples, length(outgroup_samples))
    outgroup = sample(outgroup_samples, num_outgroup)
    training_sample = append(training_sample, outgroup)
  }
  
  training_sample = factor(training_sample)
  training_sample = droplevels(training_sample)
  obs_counts = table(training_sample)
  
  gm_mean = function(x, na.rm=TRUE){
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
  
  obs_weights = gm_mean(obs_counts) / obs_counts
  print(obs_counts)
  excluded_cell_types = names(which(obs_counts < min_observations))
  training_sample = training_sample[training_sample %in% excluded_cell_types == FALSE]
  
  cds_sub = cds[,names(training_sample)]
  
  
  # chi_sq_tests = smartEsApply(cds_sub, 1, function(x) { 
  #   expr_x = as.factor(x > 0)
  #   if(length(levels(expr_x)) < 2){
  #     return (1); 
  #   }
  #   return(suppressWarnings(chisq.test(expr_x, ctf_cell_type)$p.value));
  #   },
  #   convert_to_dense=FALSE)
  # chi_sq_tests = unlist(chi_sq_tests)
  # chi_sq_tests = p.adjust(chi_sq_tests)
  # candidate_model_genes = sort(names(chi_sq_tests[chi_sq_tests < 0.05]))
  
  y = droplevels(training_sample)
  
  candidate_model_genes = c()
  for (cell_type in levels(y)){
    genes_in_cell_type = names(which(Matrix::rowSums(exprs(cds_sub[,y == cell_type]) > 0) > 0.01 * sum(y == cell_type)))
    candidate_model_genes = append(candidate_model_genes, genes_in_cell_type)
  }
  candidate_model_genes = unique(candidate_model_genes)
  #candidate_model_genes = names(which(Matrix::rowSums(exprs(cds_sub) > 0) > 0))
  
  cds_sub = cds_sub[candidate_model_genes,]
  
  
  #obs_weights = obs_weights*obs_weights
  
  #x = t(t(exprs(cds_sub)) / pData(cds_sub)$Size_Factor)
  
  #x@x = log(x@x + 1)
  #x = x > 0
  x = t(exprs(cds_sub))
  
  
  #require(doMC)
  #registerDoMC(cores=4)
  predictions = tryCatch({
    if (cores > 1){
      registerDoParallel(cores=cores)
      cvfit = cv.glmnet(x, y, 
                        #weights=pData(cds_sub)$Size_Factor, 
                        weights=obs_weights[y],
                        #offset=matrix(rep(log(pData(cds_sub)$Size_Factor), length(levels(y))),ncol=length(levels(y))),
                        alpha=1, family = "multinomial", type.multinomial = "grouped", 
                        type.measure="class",
                        type.logistic = "modified.Newton",
                        lambda.min.ratio=0.001,
                        standardize=FALSE,
                        parallel=TRUE,
                        thresh=1e-6,
                        nfolds=3,
                        nlambda=20)
    }else{
      cvfit = cv.glmnet(x, y, 
                        #weights=pData(cds_sub)$Size_Factor, 
                        weights=obs_weights[y],
                        #offset=matrix(rep(log(pData(cds_sub)$Size_Factor), length(levels(y))),ncol=length(levels(y))),
                        alpha=1, family = "multinomial", type.multinomial = "grouped", 
                        type.logistic = "modified.Newton",
                        type.measure="class",
                        lambda.min.ratio=0.001,
                        standardize=FALSE,
                        parallel=FALSE,
                        thresh=1e-6,
                        nfolds=3,
                        nlambda=20)
    }
    feature_genes = coef(cvfit, s = "lambda.min")
    #sort(names(which(Matrix::rowSums(abs(feature_genes$`T cell`)) > 0)))
    message("Model training finished.")
    
    #x = t(t(exprs(cds[candidate_model_genes,])) / pData(cds)$Size_Factor)
    x = t(exprs(cds[candidate_model_genes,]))
    
    message("Making predictions...")
    prediction_probs = as.data.frame(predict(cvfit, 
                                             newx = x,
                                             #newoffset=matrix(rep(log(pData(cds)$Size_Factor), length(levels(y))),ncol=length(levels(y))),
                                             s = "lambda.min", type = "response"))
    prediction_probs = melt(rownames_to_column(prediction_probs))
    prediction_probs = prediction_probs %>% dplyr::group_by(rowname) %>% dplyr::mutate(assignment_prob = value / max(value))
    prediction_probs = prediction_probs %>% dplyr::arrange(rowname, desc(assignment_prob)) %>% dplyr::top_n(2) %>% dplyr::mutate(odds_ratio = value / min(value))
    
    assignments = prediction_probs %>% dplyr::filter(odds_ratio > rank_prob_ratio)
    
    random_guess_thresh = 1.0 / length(levels(y))
    #min_guess_prob = 0.8
    assignments = assignments %>% dplyr::filter(value > random_guess_thresh)
    
    assignments = dplyr::left_join(data.frame(rowname=row.names(pData(cds))), assignments)
    #assignments = assignments, class_df)
    #assignments$CurrLevelCellType = unlist(assignments$CurrLevelCellType)
    assignments$variable = as.character(unlist(assignments$variable))
    
    assignments$variable = str_replace_all(assignments$variable, "\\.1", "")
    
    #assignments = assignments %>% dplyr::mutate(FinalCellType = dplyr::if_else(dplyr::if_else(variable %in% row.names(type_distances) == TRUE && CellType %in% row.names(type_distances) == TRUE,
    #                                                     is.finite(type_distances[CellType, variable]), CellType, variable), variable))
    
    predictions = dcast(assignments, rowname ~ variable)
    row.names(predictions) = predictions$rowname
    
    if (ncol(predictions) > 2){
      predictions = predictions[,setdiff(colnames(predictions), "NA")]
      predictions = predictions[,-1]
      predictions = predictions[rownames(pData(cds)),]
      predictions = as.matrix(predictions)
      predictions[is.na(predictions)] = FALSE
      predictions[predictions != 0] = TRUE
      cell_type_names = colnames(predictions)
      
      predictions = split(predictions, rep(1:ncol(predictions), each = nrow(predictions)))
      excluded_cell_types = append(excluded_cell_types, setdiff(levels(y), cell_type_names))
      if (length(excluded_cell_types > 0)){
        excluded_cell_type_predictions = matrix(FALSE, nrow=nrow(pData(cds)), ncol=length(excluded_cell_types), dimnames=list(row.names(pData(cds)), excluded_cell_types))
        excluded_cell_type_predictions = split(excluded_cell_type_predictions, rep(1:ncol(excluded_cell_type_predictions), each = nrow(excluded_cell_type_predictions)))
        predictions = append(predictions, excluded_cell_type_predictions)
        names(predictions) = c(cell_type_names, excluded_cell_types)
      }else{
        names(predictions) = cell_type_names
      }
    }else{
      cell_type_names = levels(y)
      predictions = matrix(FALSE, nrow=nrow(pData(cds)), ncol=length(cell_type_names), dimnames=list(row.names(pData(cds)), cell_type_names))
      predictions = split(predictions, rep(1:ncol(predictions), each = nrow(predictions)))
      names(predictions) = cell_type_names
    }
    predictions
  }, error = function(e) { 
    print (e)
    cell_type_names = levels(y)
    predictions = matrix(FALSE, nrow=nrow(pData(cds)), ncol=length(cell_type_names), dimnames=list(row.names(pData(cds)), cell_type_names))
    predictions = split(predictions, rep(1:ncol(predictions), each = nrow(predictions)))
    names(predictions) = cell_type_names
    predictions
  })
  
  #predictions = data.frame(row.names = assignments$rowname, CellType = assignments$variable)
  #colnames(predictions) = "CellType"
  #row.names(predictions) = rownames(x)
  #predictions = model.matrix(~.+0,predictions)

  for (i in 1:length(predictions)){
    p = as(as(predictions[[i]], "sparseVector"), "sparseMatrix")
    row.names(p) = row.names(pData(cds))
    predictions[[i]] = p
  }
  
  return(predictions)
}


#' @title classifyCellsGlmNet
#' @description 
#' Description of classifyCellsGlmNet function -- To be added 
#' 
#' @param cds CellDataSet containing cells that will be clustered
#' @param cth CellTypeHierarchy that dictates cell types present and requirements of cell to be considered  a certain cell type
#' @param rank_prob_ratio The probability ratio of the rank 
#' @param min_observations Minimun of the observation 
#' @param max_training_samples Maximum training samples 
#' @param cores Number of cores computer should use to execute function
#' @export
classifyCellsGlmNet <- function(cds, cth, rank_prob_ratio = 2, min_observations = 8,  max_training_samples = 10000, cores=1){
  
  gate_res <- list()
  for (v in V(cth@classificationTree)){
    cell_class_func <- V(cth@classificationTree) [ v ]$classify_func[[1]]
    
    parent <- environment(cell_class_func)
    if (is.null(parent))
      parent <- emptyenv()
    e1 <- new.env(parent=parent)
    multiassign(names(pData(cds)), pData(cds), envir=e1)
    environment(cell_class_func) <- e1
    
    type_res <- cell_class_func(exprs(cds))
    if (length(type_res)!= ncol(cds)){
      message(paste("Error: classification function for", V(cth@classificationTree) [ v ]$name, "returned a malformed result."))
      stop()
    }
    type_res = as(as(type_res,"sparseVector"), "sparseMatrix")
    row.names(type_res) = row.names(pData(cds))
    colnames(type_res) =V(cth@classificationTree) [ v ]$name
    gate_res[[ V(cth@classificationTree) [ v ]$name]] <- type_res
  }
  
  imputed_gate_res <- list()
  for (v in V(cth@classificationTree)){
    child_cell_types = V(cth@classificationTree) [ suppressWarnings(nei(v, mode="out")) ]$name
    if (length(child_cell_types) > 0){
      new_gate_res = cth_train_glmnet(cds, cth, v, gate_res, 
                                      rank_prob_ratio=rank_prob_ratio, 
                                      min_observations=min_observations, 
                                      max_training_samples=max_training_samples, 
                                      cores=cores)
      imputed_gate_res = append(imputed_gate_res, new_gate_res)
    }
  }
  
  cds_pdata <- dplyr::group_by_(dplyr::select_(rownames_to_column(pData(cds)), "rowname"), "rowname") 
  CellType = cth_classifier_cell(cth, imputed_gate_res)
  MarkerCellType = cth_classifier_cell(cth, gate_res)

  type_distances = igraph::distances(cth@classificationTree, mode = "out")
  cells_to_assign = names(CellType)[which(CellType %in% c("Unknown", "Ambiguous") == FALSE & MarkerCellType %in% row.names(type_distances))]
  cells_to_assign = names(CellType[cells_to_assign])[which(is.finite(type_distances[cbind(
                                                                         CellType[cells_to_assign], 
                                                                         MarkerCellType[cells_to_assign])]))]
  if (length(cells_to_assign > 0)){
    CellType[cells_to_assign] = MarkerCellType[cells_to_assign]
  }
  
  cells_to_assign = names(CellType)[which(CellType == "Unknown")]
  if (length(cells_to_assign > 0)){
    CellType[cells_to_assign] = MarkerCellType[cells_to_assign]
  }
  
  #CellType <- cth_classifier_cell(cds_subset, cth, "root", gate_res)
  return(CellType)
  
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
classifyCells <- function(cds, cth, method=c("glmnet", "markers-only"), rank_prob_ratio = 2, min_observations=8, max_training_samples=10000, cores=1) {
  if(is.null(method)) {
    method = "glmnet"
  }
  progress_opts <- options()$dplyr.show_progress
  options(dplyr.show_progress = F)
  if (method == "glmnet"){
    type_res <- classifyCellsGlmNet(cds, cth, rank_prob_ratio, min_observations, max_training_samples, cores)
  } else if (method == "markers-only"){
      type_res <- classifyCellsHelperCell(cds, cth)
  } else {
      stop("method should either be set to \"glmnet\" or \"markers-only\"")
  }
 
  class_df <- data.frame(rowname = names(type_res), CellType = type_res)
  class_df$CellType <- as.character(class_df$CellType)
  class_df$rowname <- as.character(class_df$rowname)
  
  options(dplyr.show_progress = progress_opts)
  
  pData(cds) <- pData(cds)[!(names(pData(cds)) %in% "CellType")]
  
  pData(cds) <- as.data.frame(suppressMessages(inner_join(rownames_to_column(pData(cds)), class_df)))
  
  pData(cds)$CellType <- factor(pData(cds)$CellType)
  
  row.names(pData(cds)) <- pData(cds)$rowname
  pData(cds) <- pData(cds)[,-1]
  cds
}   

#' @import methods
#' @importFrom Biobase exprs pData
#' @importFrom igraph V
cth_classifier_cluster_cds <- function(cds_subset, cth, curr_node, frequency_thresh) {
  #curr_cell_vertex <-  V(cth@classificationTree)[curr_node]
  next_nodes <- c()
  print (paste ("Cluster", unique(pData(cds_subset)$Cluster)))
  for (child in V(cth@classificationTree) [ suppressWarnings(nei(curr_node, mode="out")) ]){
    
    #child_cell_class_func <- V(cth@classificationTree) [ child ]$classify_func[[1]]
    #type_res <- sparseApply(exprs(cds_subset), 2, child_cell_class_func, convert_to_dense=FALSE)
    #type_res <- child_cell_class_func(exprs(cds_subset))
    type_res <- pData(cds_subset)$CellType %in% V(cth@classificationTree)[subcomponent(cth@classificationTree, child, "out")]$name
    
    #type_res <- unlist(type_res)
    if (sum(type_res) > 0){
      names(type_res) <- row.names(pData(cds_subset))
      cell_type_name <- V(cth@classificationTree) [ child ]$name
      if (length(frequency_thresh) > 1)
        required_thresh <- frequency_thresh[cell_type_name]
      else
        required_thresh <- frequency_thresh
      fraction_present = sum(type_res) / length(type_res)
      print (paste("       ", cell_type_name, fraction_present))
      if (fraction_present > frequency_thresh){
        next_nodes <- c(next_nodes, cell_type_name)
      }
    }
    #print (paste(V(cth@classificationTree) [ child ]$name, ":", sum(type_res),  " of ", length(type_res) ))
  }
  
  if (length(next_nodes) == 1){
    CellType <- cth_classifier_cluster_cds(cds_subset, cth, next_nodes[1], frequency_thresh)
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

classifyClustersHelperCds <- function(cds_subset, cth, frequency_thresh){
  CellType <- cth_classifier_cluster_cds(cds_subset, cth, "root", frequency_thresh)
}


#' @title classifyCellClusters
#' @param cds CellDataSet containing cells that will be clustered
#' @param cth CellTypeHierarchy that dictates cell types present and requirements of cell to be considered  a certain cell type
#' @param frequency_thresh When a CellTypeHierarchy is provided, cluster cells will impute cell types in clusters that are composed of at least this much of exactly one cell type.
#' @param grouping_var variable that controls how cells are grouped
#' @description Impute cell type classifications using clustering results
#' @export
classifyCellClusters <- function(cds, cth, frequency_thresh=0.80, grouping_var = "Cluster") {
  progress_opts <- options()$dplyr.show_progress
  options(dplyr.show_progress = F)
  
  if (is.null(frequency_thresh))
    stop("Error: to use classifyCells in grouped mode, you must also set frequency_thresh")
  
  cds_pdata <- dplyr::group_by_(dplyr::select_(rownames_to_column(pData(cds)), "rowname", grouping_var), grouping_var) 
  class_df <- as.data.frame(cds_pdata %>% dplyr::do(CellType = classifyClustersHelperCds(cds[,.$rowname], cth, frequency_thresh)))
  class_df$CellType <-  as.character(unlist(class_df$CellType))
  #class_df$rowname <- as.character(class_df$rowname)
  
  options(dplyr.show_progress = progress_opts)
  
  old_cell_types = as.character(pData(cds)$CellType)
  
  pData(cds) <- pData(cds)[!(names(pData(cds)) %in% "CellType")]
  
  #pData(cds)$cell_type <- cds_types
  
  pData(cds) <- as.data.frame(suppressMessages(inner_join(rownames_to_column(pData(cds)), class_df)))
  
  row.names(pData(cds)) <- pData(cds)$rowname
  pData(cds) <- pData(cds)[,-1]
  
  pData(cds)$CellType[pData(cds)$CellType == "Unknown"] = old_cell_types[pData(cds)$CellType == "Unknown"]
  
  pData(cds)$CellType <- factor(pData(cds)$CellType)
  
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
    cds <- classifyCells(cds, cth)
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

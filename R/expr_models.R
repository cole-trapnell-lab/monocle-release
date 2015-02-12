#' Helper function for parallel VGAM fitting
#' 
#' @param relative_expr Whether to transform expression into relative values
fit_model_helper <- function(x, 
                             modelFormulaStr, 
                             expressionFamily, 
                             relative_expr, 
                             disp_func=NULL, 
                             pseudocount=0){
  modelFormulaStr <- paste("f_expression", modelFormulaStr, sep="")
  
  orig_x <- x
  x <- x + pseudocount
  
  if (expressionFamily@vfamily == "negbinomial"){
    if (relative_expr)
    {
      x <- x / Size_Factor
    }
    f_expression <- round(x)
    if (is.null(disp_func) == FALSE){
      disp_guess <- calulate_NB_dispersion_hint(disp_func, round(x_orig))
      if (is.null(disp_guess) == FALSE && disp_guess > 0 && is.na(disp_guess) == FALSE   ) {
        # FIXME: In theory, we could lose some user-provided parameters here
        # e.g. if users supply zero=NULL or something.    
        expressionFamily <- negbinomial(isize=1/disp_guess)
      }
    }
  }else if (expressionFamily@vfamily %in% c("gaussianff", "uninormal")){
    f_expression <- x
  }else{
    f_expression <- log10(x)
  }
  
  tryCatch({
    FM_fit <-  suppressWarnings(VGAM::vglm(as.formula(modelFormulaStr), family=expressionFamily))
    FM_fit
  }, 
  #warning = function(w) { FM_fit },
  error = function(e) { print (e); NULL }
  )
}


#' Fits a model for each gene in a CellDataSet object.
#' @param cds the CellDataSet upon which to perform this operation
#' @param modelFormulaStr a formula string specifying the model to fit for the genes.
#' @param cores the number of processor cores to be used during fitting.
#' @return a list of VGAM model objects
#' @export
#' @details
#' 
#' This function fits a Tobit-family vector generalized additive model (VGAM) from the VGAM package for each gene in a CellDataSet. 
#' The default formula string speficies that the (log transformed) expression values follow a Tobit distribution with upper and lower bounds
#' specificed by max_expr and min_expr, respectively. By default, expression levels are modeled as smooth functions of the Pseudotime value of each 
#' cell. That is, expression is a function of progress through the biological process.  More complicated formulae can be provided to account for
#' additional covariates (e.g. day collected, genotype of cells, media conditions, etc).
fitModel <- function(cds,
                     modelFormulaStr="~sm.ns(Pseudotime, df=3)",
                     relative_expr=TRUE,
                     pseudocount=0,
                     cores=1){
  if (cores > 1){
    f<-mcesApply(cds, 1, fit_model_helper, required_packages=c("BiocGenerics", "VGAM", "plyr"), cores=cores, 
                 modelFormulaStr=modelFormulaStr, 
                 expressionFamily=cds@expressionFamily,
                 relative_expr=relative_expr,
                 disp_func=cds@dispFitInfo[["blind"]]$disp_func,
                 pseudocount=pseudocount)
    f
  }else{
    f<-esApply(cds,1,fit_model_helper, 
               modelFormulaStr=modelFormulaStr, 
               expressionFamily=cds@expressionFamily,
               relative_expr=relative_expr,
               disp_func=cds@dispFitInfo[["blind"]]$disp_func,
               pseudocount=pseudocount)
    f
  }
}

#' Response values
#' 
#' Generates a matrix of response values for a set of fitted models
#' @param models a list of models, e.g. as returned by fitModels()
#' @return a matrix where each row is a vector of response values for a particular feature's model, and columns are cells.
#' @export
responseMatrix <- function(models){
  res_list <- lapply(models, function(x) { 
    if (is.null(x)) { NA } else { 
      if (x@family@vfamily  == "negbinomial"){
        VGAM::predict(x, type="response") 
      }else if (x@family@vfamily %in% c("gaussianff", "uninormal")){
        VGAM::predict(x, type="response")
      }else{
        10^predict(x, type="response") 
      }
    } 
  } )
  res_list_lengths <- lapply(res_list[is.na(res_list) == FALSE], length)
  
  stopifnot(length(unique(res_list_lengths)) == 1)
  num_na_fits <- length(res_list[is.na(res_list)])
  if (num_na_fits > 0){
    na_matrix<- matrix(rep(rep(NA, res_list_lengths[[1]]), num_na_fits),nrow=num_na_fits) 
    row.names(na_matrix) <- names(res_list[is.na(res_list)])
    
    non_na_matrix <- t(do.call(cbind, lapply(res_list[is.na(res_list) == FALSE], unlist)))
    row.names(non_na_matrix) <- names(res_list[is.na(res_list) == FALSE])
    res_matrix <- rbind(non_na_matrix, na_matrix)
    res_matrix <- res_matrix[names(res_list),]
  }else{
    res_matrix <- t(do.call(cbind, lapply(res_list, unlist)))
    row.names(res_matrix) <- names(res_list[is.na(res_list) == FALSE])
  }
  
  res_matrix
}


## This function was swiped from DESeq (Anders and Huber) and modified for our purposes
parametricDispersionFit <- function( means, disps )
{
  coefs <- c( .1, 1 )
  iter <- 0
  while(TRUE) {
    residuals <- disps / ( coefs[1] + coefs[2] / means )
    good <- which( (residuals > 1e-4) & (residuals < 15) )
    fit <- glm( disps[good] ~ I(1/means[good]),
                family=Gamma(link="identity"), start=coefs )
    oldcoefs <- coefs
    coefs <- coefficients(fit)
    if( !all( coefs > 0 ) ){
      #print(data.frame(means,disps))
      stop( "Parametric dispersion fit failed. Try a local fit and/or a pooled estimation. (See '?estimateDispersions')" )
      
      
    }
    if( sum( log( coefs / oldcoefs )^2 ) < 1e-6 )
      break
    iter <- iter + 1
    if( iter > 10 ) {
      warning( "Dispersion fit did not converge." )
      break }
  }
  
  names( coefs ) <- c( "asymptDisp", "extraPois" )
  ans <- function( q )
    coefs[1] + coefs[2] / q
  #ans
  coefs
}

## This function was swiped from DESeq (Anders and Huber) and modified for our purposes
#' @export
vstExprs <- function(cds, dispModelName="blind", expr_matrix=NULL, round_vals=TRUE ) {
  fitInfo <- cds@dispFitInfo[[dispModelName]]
  if (is.null(fitInfo)){
    stop(paste("Error: No dispersion model named '",dispModelName,"'. You must call estimateSizeFactors(...,dispModelName='",dispModelName,
               "') before calling this function.", sep=""))
  }
  
  coefs <- attr( fitInfo$disp_func, "coefficients" )
  if (is.null(expr_matrix)){
    ncounts <- exprs(cds)
    ncounts <- t(t(ncounts) / sizeFactors(cds))
    if (round_vals)
      ncounts <- round(ncounts)
  }else{
    ncounts <- expr_matrix
  }
  vst <- function( q )
    log( (1 + coefs["extraPois"] + 2 * coefs["asymptDisp"] * q +
            2 * sqrt( coefs["asymptDisp"] * q * ( 1 + coefs["extraPois"] + coefs["asymptDisp"] * q ) ) )
         / ( 4 * coefs["asymptDisp"] ) ) / log(2)
  vst( ncounts )
}


#' Helper function for parallel dispersion modeling
#' 
#' @importFrom dplyr distinct
disp_calc_helper <- function(x, modelFormulaStr, expressionFamily){
  disp <- NULL
  modelFormulaStr <- paste("f_expression", modelFormulaStr, sep="")
  f_expression <- x

  if (expressionFamily@vfamily == "negbinomial"){
    x <- x / Size_Factor
    f_expression <- round(x)
  }else if (expressionFamily@vfamily %in% c("gaussianff", "uninormal")){
    f_expression <- x
  }else{
    f_expression <- log10(x)
  }
  
  disp_vals <- tryCatch({
    fitted_model <-  suppressWarnings(VGAM::vgam(as.formula(modelFormulaStr), family=expressionFamily))
    disp_vals <- as.data.frame(VGAM::predict(fitted_model))
    colnames(disp_vals) <- c("mu", "disp")
    disp_vals$disp <- signif(disp_vals$disp)
    disp_vals <- dplyr::distinct(disp_vals, disp)
    
    disp_vals$mu <- exp(disp_vals$mu)
    disp_vals$disp <- 1.0/exp(disp_vals$disp)
    df_res <- data.frame(mu=disp_vals$mu, disp=disp_vals$disp)
  },
  #warning = function(w) { print (w) },
  error = function(e) { print (e); NULL }
  )
  disp_vals
}

estimateDispersionsForCellDataSet <- function(cds, modelFormulaStr, relative_expr, cores)
{
  
  if (cores > 1){
    disp_table<-mcesApply(cds, 1, disp_calc_helper, c("BiocGenerics", "Biobase", "VGAM", "dplyr"), cores=cores, 
                          modelFormulaStr=modelFormulaStr, 
                          expressionFamily=cds@expressionFamily)
  }else{
    disp_table<-esApply(cds,1,disp_calc_helper, 
                        modelFormulaStr=modelFormulaStr, 
                        expressionFamily=cds@expressionFamily)
  }
  #print (disp_table)
  disp_table <- do.call(rbind.data.frame, disp_table)
  
  coefs <- parametricDispersionFit(disp_table$mu, disp_table$disp)
  
  names( coefs ) <- c( "asymptDisp", "extraPois" )
  ans <- function( q )
    coefs[1] + coefs[2] / q
  attr( ans, "coefficients" ) <- coefs
  
  res <- list(disp_table = disp_table, disp_func = ans)
  return(res)
}

buildLineageBranchCellDataSet <- function(cds, lineage_states = c(2, 3), lineage_labels=NULL, method = 'fitting', stretch = FALSE)
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

    short_branches_multipliers <- (100 - T_0) / (range_df[short_branches , 'max'] - max(range_df[as.character(progenitor_states), 'max']))
    
    heterochrony_constant <- c(longest_branch_multiplier, short_branches_multipliers)
    names(heterochrony_constant) <- c(longest_lineage_branch, short_branches)

    pData[pData$State %in% c(progenitor_states, longest_lineage_branch), 'Pseudotime'] <- 
      (pData[pData$State %in% c(progenitor_states, longest_lineage_branch), 'Pseudotime'] - min(range_df[progenitor_states, 'min'])) * longest_branch_multiplier

    
    for(i in 1:length(short_branches)) { #stretch short branches
      pData[pData$State  == short_branches[i], 'Pseudotime'] <- 
        (pData[pData$State == short_branches[i], 'Pseudotime'] - max(range_df[as.character(progenitor_states), 'max'])) * 
        short_branches_multipliers[i] + T_0
    }
  }
  pData$original_cell_id <- pData$Cell
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
  pData$Lineage <- as.factor(pData$State)
  pData$State <- factor(pData(cds)[as.character(pData$original_cell_id),]$State, 
                            levels =levels(cds$State))
  pData$weight <- weight_vec
  Size_Factor <- pData$Size_Factor
  
  fData <- fData(cds)
  
  cds_subset <- newCellDataSet(as.matrix(exprs_data),
                               phenoData = new("AnnotatedDataFrame", data = pData),
                               featureData = new("AnnotatedDataFrame", data = fData),
                               expressionFamily=cds@expressionFamily,
                               lowerDetectionLimit=cds@lowerDetectionLimit)
  pData(cds_subset)$State <- as.factor(pData(cds_subset)$State)
  pData(cds_subset)$Size_Factor <- Size_Factor
  return (cds_subset)
}

calulate_NB_dispersion_hint <- function(disp_func, f_expression)
{
  expr_median <- median(f_expression[f_expression > 0])
  if (is.null(expr_median) == FALSE) {
    disp_guess_fit <- disp_func(expr_median)
    
    # For NB: Var(Y)=mu*(1+mu/k)
    f_expression_var <- var(f_expression)
    f_expression_mean <- mean(f_expression)
    
    disp_guess_meth_moments <- f_expression_var - f_expression_mean
    disp_guess_meth_moments <- f_expression_var / (f_expression_mean^2)
    
    return (max(disp_guess_fit, disp_guess_meth_moments))
  }
  return (NULL)
}

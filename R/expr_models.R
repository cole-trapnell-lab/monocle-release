#' Helper function for parallel VGAM fitting
#' 
#' @importFrom VGAM vgam
fit_model_helper <- function(x, modelFormulaStr, expressionFamily, relative_expr){
  if (expressionFamily@vfamily == "negbinomial"){
    if (relative_expr)
    {
      x <- x / Size_Factor
    }
    expression <- round(x)
  }else if (expressionFamily@vfamily %in% c("gaussianff", "uninormal")){
    expression <- x
  }else{
    expression <- log10(x)
  }
  
  tryCatch({
    FM_fit <-  suppressWarnings(VGAM::vgam(as.formula(modelFormulaStr), family=expressionFamily))
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
                     modelFormulaStr="expression~sm.ns(Pseudotime, df=3)",
                     relative_expr=TRUE,
                     cores=1){
  if (cores > 1){
    f<-mcesApply(cds, 1, fit_model_helper, required_packages=c("BiocGenerics", "VGAM", "plyr"), cores=cores, 
                 modelFormulaStr=modelFormulaStr, 
                 expressionFamily=cds@expressionFamily,
                 relative_expr=relative_expr)
    f
  }else{
    f<-esApply(cds,1,fit_model_helper, 
               modelFormulaStr=modelFormulaStr, 
               expressionFamily=cds@expressionFamily,
               relative_expr=relative_expr)
    f
  }
}

#' Response values
#' 
#' Generates a matrix of response values for a set of fitted models
#' @param models a list of models, e.g. as returned by fitModels()
#' @return a matrix where each row is a vector of response values for a particular feature's model, and columns are cells.
#' @importFrom VGAM predict
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
#' @importFrom VGAM predict vgam
disp_calc_helper <- function(x, modelFormulaStr, expressionFamily){
  if (expressionFamily@vfamily == "negbinomial"){
    x <- x / Size_Factor
    expression <- round(x)
  }else if (expressionFamily@vfamily %in% c("gaussianff", "uninormal")){
    expression <- x
  }else{
    expression <- log10(x)
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
    disp_table<-mcesApply(cds, 1, disp_calc_helper, c("BiocGenerics", "VGAM", "dplyr"), cores=cores, 
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


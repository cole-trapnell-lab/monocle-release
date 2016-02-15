# Helper function for parallel differential expression testing
diff_test_helper <- function(x, 
                             fullModelFormulaStr, 
                             reducedModelFormulaStr, 
                             expressionFamily, 
                             relative_expr,
                             weights,
                             disp_func=NULL,
                             exprs_thrsld_percentage = NULL, 
                             verbose=FALSE){
  if(is.null(exprs_thrsld_percentage) == FALSE) {
    if((sum(x > 0) / length(x)) < exprs_thrsld_percentage) {
      test_res <- data.frame(status = "NO_TEST", family=NA, pval=1.0, qval=1.0)
      return(test_res)
    }
  }
  
  reducedModelFormulaStr <- paste("f_expression", reducedModelFormulaStr, sep="")
  fullModelFormulaStr <- paste("f_expression", fullModelFormulaStr, sep="")
  
  x_orig <- x
  
  disp_guess <- 0
  
  if (expressionFamily@vfamily == "negbinomial"){
    if (relative_expr == TRUE)
    {
      x <- x / Size_Factor
    }
    f_expression <- round(x)
    if (is.null(disp_func) == FALSE){
      disp_guess <- calulate_NB_dispersion_hint(disp_func, round(x_orig))
      if (is.null(disp_guess) == FALSE && disp_guess > 0 && is.na(disp_guess) == FALSE  ) {
        # FIXME: In theory, we could lose some user-provided parameters here
        # e.g. if users supply zero=NULL or something. 
        expressionFamily <- negbinomial(isize=1/disp_guess)
      }
    }
  }else if (expressionFamily@vfamily %in% c("gaussianff", "uninormal")){
    f_expression <- x
  }else if (expressionFamily@vfamily %in% c("binomialff")){
    f_expression <- x
    f_expression[f_expression > 1] <- 1
  }else{
    f_expression <- log10(x)
  }
  
  test_res <- tryCatch({
    if (verbose){
      full_model_fit <- VGAM::vglm(as.formula(fullModelFormulaStr), family=expressionFamily, checkwz=TRUE)
      reduced_model_fit <- VGAM::vglm(as.formula(reducedModelFormulaStr), family=expressionFamily, checkwz=TRUE)                         
    }else{
      full_model_fit <- suppressWarnings(VGAM::vglm(as.formula(fullModelFormulaStr), family=expressionFamily))
      reduced_model_fit <- suppressWarnings(VGAM::vglm(as.formula(reducedModelFormulaStr), family=expressionFamily))                    
    }
    #print(full_model_fit)
    #print(coef(reduced_model_fit))
    compareModels(list(full_model_fit), list(reduced_model_fit))
  }, 
  #warning = function(w) { FM_fit },
  error = function(e) { 
    #print (e);
    # If we threw an exception, re-try with a simpler model.  Which one depends on
    # what the user has specified for expression family
    #print(disp_guess)
    backup_expression_family <- NULL
    if (expressionFamily@vfamily == "negbinomial"){
        disp_guess <- calulate_NB_dispersion_hint(disp_func, round(x_orig), expr_selection_func = max)
        backup_expression_family <- negbinomial(isize=1/disp_guess, parallel=TRUE, zero=NULL)
    }else if (expressionFamily@vfamily %in% c("gaussianff", "uninormal")){
      backup_expression_family <- NULL
    }else if (expressionFamily@vfamily %in% c("binomialff")){
      backup_expression_family <- NULL
    }else{
      backup_expression_family <- NULL
    }
    if (FALSE){ #is.null(backup_expression_family) == FALSE #remove the backup fitting
      test_res <- tryCatch({
      if (verbose){
        full_model_fit <- VGAM::vglm(as.formula(fullModelFormulaStr), family=backup_expression_family, epsilon=1e-1, checkwz=FALSE)
        reduced_model_fit <- VGAM::vglm(as.formula(reducedModelFormulaStr), family=backup_expression_family, epsilon=1e-1,checkwz=FALSE)                         
      }else{
        full_model_fit <- suppressWarnings(VGAM::vglm(as.formula(fullModelFormulaStr), family=backup_expression_family, epsilon=1e-1,checkwz=FALSE))
        reduced_model_fit <- suppressWarnings(VGAM::vglm(as.formula(reducedModelFormulaStr), family=backup_expression_family, epsilon=1e-1, checkwz=FALSE))                    
      }
      #print(full_model_fit)
      #print(coef(reduced_model_fit))
      compareModels(list(full_model_fit), list(reduced_model_fit))
      }, 
      #warning = function(w) { FM_fit },
      error = function(e) { 
        #print (e);
        data.frame(status = "FAIL", family=NA, pval=1.0, qval=1.0)
      })
      #print(test_res)
      test_res
    } else {
      data.frame(status = "FAIL", family=NA, pval=1.0, qval=1.0)
    }
    #data.frame(status = "FAIL", pval=1.0) 
  }
  )
  test_res
}

#' Compare model fits
#' 
#' Performs likelihood ratio tests on nested vector generalized additive models 
#' @param full_models a list of models, e.g. as returned by fitModels(), forming the numerators of the L.R.Ts.
#' @param reduced_models a list of models, e.g. as returned by fitModels(), forming the denominators of the L.R.Ts.
#' @return a data frame containing the p values and q-values from the likelihood ratio tests on the parallel arrays of models.
#' @export
compareModels <- function(full_models, reduced_models){
  stopifnot(length(full_models) == length(reduced_models))
  test_res <- mapply(function(x,y) { 
    if (is.null(x) == FALSE && is.null(y) == FALSE) {
      lrt <- VGAM::lrtest(x,y) 
      pval=lrt@Body["Pr(>Chisq)"][2,]
      data.frame(status = "OK", family=x@family@vfamily, pval=pval)
    } else { data.frame(status = "FAIL", family=NA, pval=1.0) } 
  } , full_models, reduced_models, SIMPLIFY=FALSE, USE.NAMES=TRUE)
  
  test_res <- do.call(rbind.data.frame, test_res)
  test_res$qval <- p.adjust(test_res$pval, method="BH")
  test_res
}

#' Tests each gene for differential expression as a function of progress through a biological process, or according to other covariates as specified. 
#' @param cds a CellDataSet object upon which to perform this operation
#' @param fullModelFormulaStr a formula string specifying the full model in differential expression tests (i.e. likelihood ratio tests) for each gene/feature.
#' @param reducedModelFormulaStr a formula string specifying the reduced model in differential expression tests (i.e. likelihood ratio tests) for each gene/feature.
#' @param cores the number of cores to be used while testing each gene for differential expression.
#' @param relative_expr Whether to transform expression into relative values.
#' @param exprs_thrsld_percentage Skip testing for genes expressed in fewer than this percentage of cells. Default is NULL, which excludes no genes from testing.
#' @param verbose Whether to show VGAM errors and warnings. Only valid for cores = 1. 
#' @return a data frame containing the p values and q-values from the likelihood ratio tests on the parallel arrays of models.
#' @export
differentialGeneTest <- function(cds, 
                                 fullModelFormulaStr="~sm.ns(Pseudotime, df=3)",
                                 reducedModelFormulaStr="~1", 
                                 cores=1, 
                                 relative_expr=TRUE,
                                 exprs_thrsld_percentage = NULL, 
                                 verbose=FALSE){
  if (relative_expr && cds@expressionFamily@vfamily == "negbinomial"){
    if (is.null(sizeFactors(cds))){
      stop("Error: to call this function with relative_expr==TRUE, you must first call estimateSizeFactors() on the CellDataSet.")
    }
  }
  
  if (cores > 1){
    diff_test_res<-mcesApply(cds, 1, diff_test_helper, 
                             c("BiocGenerics", "VGAM", "Matrix"), 
                             cores=cores, 
                             fullModelFormulaStr=fullModelFormulaStr,
                             reducedModelFormulaStr=reducedModelFormulaStr,
                             expressionFamily=cds@expressionFamily,
                             relative_expr=relative_expr,
                             disp_func=cds@dispFitInfo[["blind"]]$disp_func,
                             exprs_thrsld_percentage = exprs_thrsld_percentage, 
                             verbose=verbose)
    diff_test_res
  }else{
    diff_test_res<-smartEsApply(cds,1,diff_test_helper, 
                                fullModelFormulaStr=fullModelFormulaStr,
                                reducedModelFormulaStr=reducedModelFormulaStr, 
                                expressionFamily=cds@expressionFamily, 
                                relative_expr=relative_expr,
                                disp_func=cds@dispFitInfo[["blind"]]$disp_func,
                                exprs_thrsld_percentage = exprs_thrsld_percentage, 
                                verbose=verbose)
    diff_test_res
  }
  
  diff_test_res <- do.call(rbind.data.frame, diff_test_res)

  diff_test_res$qval <- 1
  diff_test_res$qval[which(diff_test_res$status == 'OK')] <- p.adjust(subset(diff_test_res, status == 'OK')[, 'pval'], method="BH")
  
  diff_test_res <- merge(diff_test_res, fData(cds), by="row.names")
  row.names(diff_test_res) <- diff_test_res[, 1] #remove the first column and set the row names to the first column
  diff_test_res[, 1] <- NULL 

  diff_test_res
}

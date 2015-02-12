#' Helper function for parallel differential expression testing
#' 
#' @param relative_expr Whether to transform expression into relative values
diff_test_helper <- function(x, 
                             fullModelFormulaStr, 
                             reducedModelFormulaStr, 
                             expressionFamily, 
                             relative_expr,
                             weights,
                             disp_func=NULL,
                             pseudocount=0){
  reducedModelFormulaStr <- paste("f_expression", reducedModelFormulaStr, sep="")
  fullModelFormulaStr <- paste("f_expression", fullModelFormulaStr, sep="")
  
  x_orig <- x
  x <- x + pseudocount
  
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
  }else{
    f_expression <- log10(x)
  }
  
  test_res <- tryCatch({
    full_model_fit <- VGAM::vglm(as.formula(fullModelFormulaStr), family=expressionFamily, weights=weights)
    reduced_model_fit <- suppressWarnings(VGAM::vglm(as.formula(reducedModelFormulaStr), family=expressionFamily, weights=weights))
    compareModels(list(full_model_fit), list(reduced_model_fit))
  }, 
  #warning = function(w) { FM_fit },
  error = function(e) { 
    print (e); 
    NULL
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
      lrt <- lrtest(x,y) 
      pval=lrt@Body["Pr(>Chisq)"][2,]
      data.frame(status = "OK", pval=pval)
    } else { data.frame(status = "FAIL", pval=1.0) } 
  } , full_models, reduced_models, SIMPLIFY=FALSE, USE.NAMES=TRUE)
  
  test_res <- do.call(rbind.data.frame, test_res)
  test_res$qval <- p.adjust(test_res$pval, method="BH")
  test_res
}

#' Tests each gene for differential expression as a function of progress through a biological process, or according to other covariates as specified. 
#' @param cds a CellDataSet object upon which to perform this operation
#' @param fullModelFormulaStr a formula string specifying the full model in differential expression tests (i.e. likelihood ratio tests) for each gene/feature.
#' @param reducedModelFormulaStr a formula string specifying the reduced model in differential expression tests (i.e. likelihood ratio tests) for each gene/feature.
#' @param cores the number of cores to be used while testing each gene for differential expression
#' @param relative_expr Whether to transform expression into relative values
#' @return a data frame containing the p values and q-values from the likelihood ratio tests on the parallel arrays of models.
#' @export
differentialGeneTest <- function(cds, 
                                 fullModelFormulaStr="~sm.ns(Pseudotime, df=3)",
                                 reducedModelFormulaStr="~1", 
                                 cores=1, 
                                 relative_expr=TRUE,
                                 weights=NULL,
                                 pseudocount=0){
  if (relative_expr && cds@expressionFamily@vfamily == "negbinomial"){
    if (is.null(sizeFactors(cds))){
      stop("Error: to call this function with relative_expr==TRUE, you must first call estimateSizeFactors() on the CellDataSet.")
    }
  }
  
  if (cores > 1){
    diff_test_res<-mcesApply(cds, 1, diff_test_helper, 
                             c("BiocGenerics", "VGAM"), 
                             cores=cores, 
                             fullModelFormulaStr=fullModelFormulaStr,
                             reducedModelFormulaStr=reducedModelFormulaStr,
                             expressionFamily=cds@expressionFamily,
                             relative_expr=relative_expr,
                             weights=weights,
                             disp_func=cds@dispFitInfo[["blind"]]$disp_func,
                             pseudocount=pseudocount)
    diff_test_res
  }else{
    diff_test_res<-esApply(cds,1,diff_test_helper, 
                           fullModelFormulaStr=fullModelFormulaStr,
                           reducedModelFormulaStr=reducedModelFormulaStr, 
                           expressionFamily=cds@expressionFamily, 
                           relative_expr=relative_expr,
                           weights=weights,
                           disp_func=cds@dispFitInfo[["blind"]]$disp_func,
                           pseudocount=pseudocount)
    diff_test_res
  }
  
  diff_test_res <- do.call(rbind.data.frame, diff_test_res)
  diff_test_res$qval <- p.adjust(diff_test_res$pval, method="BH")
  
  diff_test_res <- merge(diff_test_res, fData(cds), by="row.names")
  
  diff_test_res
}


#' find the branch genes
#' @param cds a CellDataSet object upon which to perform this operation
#' @param fullModelFormulaStr a formula string specifying the full model in differential expression tests (i.e. likelihood ratio tests) for each gene/feature.
#' @param reducedModelFormulaStr a formula string specifying the reduced model in differential expression tests (i.e. likelihood ratio tests) for each gene/feature.
#' @param cores the number of cores to be used while testing each gene for differential expression
#' @param lineage_states ids for the immediate branch lineage which obtained from lineage construction based on MST
#' @return a data frame containing the p values and q-values from the likelihood ratio tests on the parallel arrays of models.
#' @export
branchTest <- function(cds, fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3)*Lineage",
                       reducedModelFormulaStr = "~sm.ns(Pseudotime, df = 3)", 
                       lineage_states = c(2, 3), 
                       relative_expr = TRUE,
                       stretch = TRUE,
                       pseudocount=0,
                       cores = 1) {
  
  cds_subset <- buildLineageBranchCellDataSet(cds, lineage_states, lineage_labels, method, stretch)
  
  branchTest_res <- differentialGeneTest(cds_subset, 
                                         fullModelFormulaStr = fullModelFormulaStr, 
                                         reducedModelFormulaStr = reducedModelFormulaStr, 
                                         cores = cores, 
                                         relative_expr = relative_expr, 
                                         weights = pData(cds_subset)$weight,
                                         pseudocount = pseudocount)
  
  return(branchTest_res)
}




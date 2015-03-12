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
  }else if (expressionFamily@vfamily %in% c("binomialff")){
    f_expression <- x
    f_expression[f_expression > 1] <- 1
  }else{
    f_expression <- log10(x)
  }
  
  test_res <- tryCatch({
    #print (f_expression)
    full_model_fit <- suppressWarnings(VGAM::vglm(as.formula(fullModelFormulaStr), family=expressionFamily, weights=weights))
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
      lrt <- VGAM::lrtest(x,y) 
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
  row.names(diff_test_res) <- diff_test_res[, 1] #remove the first column and set the row names to the first column
  diff_test_res[, 1] <- NULL 

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
                       cores = 1, 
                       weighted = TRUE, 
                       lineage_labels = NULL, ...) {

  cds_subset <- buildLineageBranchCellDataSet(cds, lineage_states, lineage_labels, stretch, weighted, ...)
  cds_subset <- estimateSizeFactors(cds_subset) 

  branchTest_res <- differentialGeneTest(cds_subset, 
                                         fullModelFormulaStr = fullModelFormulaStr, 
                                         reducedModelFormulaStr = reducedModelFormulaStr, 
                                         cores = cores, 
                                         relative_expr = relative_expr, 
                                         weights = pData(cds_subset)$weight,
                                         pseudocount = pseudocount)
  
  return(branchTest_res)
}

#to do: 
# 1. think about how calculating ABCs for multiple lineages (use a common reference as implemented?)
# 2. how to store ABCs? 
# 3. how to use fit_model_helper directly for calculating the model fitting?
#' calculate the area between TWO fitted lineage trajectories  
#' @param cds a CellDataSet object upon which to perform this operation
#' @param fullModelFormulaStr a formula string specifying the full model in differential expression tests (i.e. likelihood ratio tests) for each gene/feature.
#' @param reducedModelFormulaStr a formula string specifying the reduced model in differential expression tests (i.e. likelihood ratio tests) for each gene/feature.
#' @param branchTest a logic flag to determine whether or not to perform the branchTest inside the function. Because of the long computations, we recommend to first perform the branchTest and then calculate the ABCs for the genes of interests. Otherwise the ABCs will be appended to the last column of branchTest results. 
#' @param cores the number of cores to be used while testing each gene for differential expression
#' @param lineage_states ids for the immediate branch lineage which obtained from lineage construction based on MST (only two lineages are allowed for this function)
#' @param min_expr the lower limit for the expressed gene
#' @param integer_expression the logic flag to determine whether or not the integer numbers are used for calculating the ABCs. Default is False. 
#' @param num number of points on the fitted lineage trajectories used for calculating the ABCs. Default is 5000. 
#' @return a data frame containing the p values and q-values from the likelihood ratio tests on the parallel arrays of models.
#' @export
calABCs <- function(cds, fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3)*Lineage",
                       reducedModelFormulaStr = "~sm.ns(Pseudotime, df = 3)", 
                       method = 'fitting', 
                       ABC_method = 'integral', 
                       points_num = 1000, 
                       fc_limit = 3, 
                       branchTest = FALSE, 
                       lineage_states = c(2, 3), 
                       relative_expr = TRUE,
                       stretch = TRUE,
                       pseudocount=0,
                       cores = 1, 
                       weighted = TRUE, 
                       min_expr = 0.5, 
                       integer_expression = FALSE, 
                       num = 5000, lineage_labels = NULL, ...) {
  if(branchTest)
    branchTest_res <- branchTest(cds, fullModelFormulaStr = fullModelFormulaStr,
                       reducedModelFormulaStr = reducedModelFormulaStr, 
                       lineage_states = lineage_states, 
                       relative_expr = relative_expr,
                       stretch = stretch,
                       pseudocount = pseudocount,
                       cores = cores, 
                       weighted = weighted, lineage_labels = lineage_labels, ...)
  if(length(lineage_states) != 2)
    stop('Sorry, this function only supports the calculation of ABCs between TWO lineage trajectories')

  new_cds <- buildLineageBranchCellDataSet(cds, lineage_states, lineage_labels, stretch, weighted, ...)
  new_cds <- estimateSizeFactors(new_cds)

  #parallelize the calculation of ABCs 
  res_ABC <- mcesApply(new_cds, 1, function(x, modelFormulaStr, ABC_method, expressionFamily, relative_expr, disp_func, pseudocount, num, lineage_states, lineage_labels, fc_limit, points_num, ...) {
    fit_res <- tryCatch({
      
      #how to pass the enviroment to fit_model_helper function?
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
          if (is.null(disp_guess) == FALSE ) {
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
      
      ##add loess curve fitting here 
      FM_fit <- tryCatch({
        FM_fit <-  suppressWarnings(VGAM::vglm(as.formula(modelFormulaStr), family=expressionFamily, maxit = 30, checkwz = FALSE))
      }, 
      #warning = function(w) { FM_fit },
      error = function(e) { print (e); NULL }
      )
  
    #calculate the area between curve: 
    #return the matrix for the difference 
    modelFormulaStr <- paste("f_expression", modelFormulaStr, sep="")

    if(!is.null(lineage_labels))
      lineage <- lineage_labels
    else 
      lineage <- lineage_states

    predictBranchOri <- VGAM::predict(FM_fit, newdata = data.frame(Pseudotime = seq(0, 100, length.out = num), Lineage = as.factor(rep(lineage[1], num))), type = 'response')
    predictBranchOthers <- lapply(lineage[2:length(lineage)], function(x) {
        VGAM::predict(FM_fit, newdata = data.frame(Pseudotime = seq(0, 100, length.out = num), Lineage = as.factor(rep(x, num))), type = 'response')        
    })
    ABCs <- lapply(predictBranchOthers, function(x){
        avg_delta_x <- ((predictBranchOri - x)[1:(num - 1)] + (predictBranchOri - x)[2:(num)]) / 2
        step <- (100 / (num - 1))

        if(ABC_method == 'integral'){
          res <- round(sum( avg_delta_x * step), 3)
        }
        else if(ABC_method == 'global_normalization'){
          max <- max(max(predictBranchOri), max(x))
          res <- round(sum( avg_delta_x / max * step), 3)
        }
        else if(ABC_method == 'local_normalization'){
          pair_wise_max <- apply(data.frame(x = x, y = predictBranchOri), 1, max)
          res <- round(sum( (((predictBranchOri - x) / pair_wise_max)[1:(num - 1)] + ((predictBranchOri - x) / pair_wise_max)[2:(num)]) / 2 * step), 3)
        }
        else if(ABC_method == 'four_values'){ #check this 
          ori_ABCs <- round(sum( (x[1:(num - 1)] + x[2:(num)]) / 2 * step), 3)
          other_ABCs <- round(sum( (predictBranchOri[1:(num - 1)] + predictBranchOri[2:(num)]) / 2 * step), 3)
          ori_ABCs_H <- round(sum( avg_delta_x[avg_delta_x > 0] * step), 3)
          other_ABCs_H <- round(sum( avg_delta_x[avg_delta_x < 0] * step), 3)
          res <- c(ori_ABCs = ori_ABCs, other_ABCs = other_ABCs, ori_ABCs_H = ori_ABCs_H, other_ABCs_H = other_ABCs_H)
        }
        else if(ABC_method == 'difference') {#copy from the plot_heatmap function 
          str_logfc_df <- log2((predictBranchOri + 1) / (x + 1))

          str_logfc_df[which(str_logfc_df <= -fc_limit)] <- -fc_limit
          str_logfc_df[which(str_logfc_df >= fc_limit)] <- fc_limit

          res <- str_logfc_df[c(seq(1, num, length.out = points_num - 1), num)] #only return 100 points for clustering 
        }
        return(res)
      })
    ABCs
    # list(res = res, ABC = unlist(ABCs)) predictBranchOri - x
  }, error = function(e) {
    print("Error!")
    print(e)
    res <- rep(NA, length(x))
    ABCs <- rep(NA, length(lineage_states) - 1)
    # list(res = res, ABC = ABCs)
    ABCs
  })

  return(fit_res)
}, required_packages = c("BiocGenerics", "VGAM", "base", "stats"), cores = cores, 
    fullModelFormulaStr, ABC_method, new_cds@expressionFamily, relative_expr, new_cds@dispFitInfo[['blind']]$disp_func, 
    pseudocount, num, lineage_states, lineage_labels, fc_limit, points_num)


  if(ABC_method %in% c('integral', 'global_normalization', 'local_normalization')) {
    ABCs_res <- do.call(rbind.data.frame, lapply(res_ABC, function(x) x[[1]]))
    row.names(ABCs_res) <- names(res_ABC)
    
    if(branchTest) {
      branchTest_res[, 'ABCs'] <- ABCs_res[row.names(branchTest_res)]
      return(branchTest_res)
    }

    ABCs_res <- merge(ABCs_res, fData(cds), by="row.names")
    row.names(ABCs_res) <- ABCs_res[, 1] #remove the first column and set the row names to the first column
    ABCs_res[, 1] <- NULL 
    colnames(ABCs_res)[1] <- 'ABCs'
  }
  else if(ABC_method %in% c('four_values', 'difference')) {
    ABCs_res <- do.call(rbind, lapply(res_ABC,  unlist))
  }
  
  return(ABCs_res)
}

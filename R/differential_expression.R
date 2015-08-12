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
                             pseudocount=0,
                             verbose=FALSE){
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
    if (verbose){
      full_model_fit <- VGAM::vglm(as.formula(fullModelFormulaStr), family=expressionFamily, weights=weights)
      reduced_model_fit <- VGAM::vglm(as.formula(reducedModelFormulaStr), family=expressionFamily, weights=weights)                         
    }else{
      full_model_fit <- suppressWarnings(VGAM::vglm(as.formula(fullModelFormulaStr), family=expressionFamily, weights=weights))
      reduced_model_fit <- suppressWarnings(VGAM::vglm(as.formula(reducedModelFormulaStr), family=expressionFamily, weights=weights))                    
    }
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
#' @param verbose Whether to show VGAM errors and warnings. Only valid for cores = 1. 
#' @return a data frame containing the p values and q-values from the likelihood ratio tests on the parallel arrays of models.
#' @export
differentialGeneTest <- function(cds, 
                                 fullModelFormulaStr="~sm.ns(Pseudotime, df=3)",
                                 reducedModelFormulaStr="~1", 
                                 cores=1, 
                                 relative_expr=TRUE,
                                 weights=NULL,
                                 pseudocount=0,
                                 verbose=FALSE){
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
                             pseudocount=pseudocount,
                             verbose=verbose)
    diff_test_res
  }else{
    diff_test_res<-esApply(cds,1,diff_test_helper, 
                           fullModelFormulaStr=fullModelFormulaStr,
                           reducedModelFormulaStr=reducedModelFormulaStr, 
                           expressionFamily=cds@expressionFamily, 
                           relative_expr=relative_expr,
                           weights=weights,
                           disp_func=cds@dispFitInfo[["blind"]]$disp_func,
                           pseudocount=pseudocount,
                           verbose=verbose)
    diff_test_res
  }
  
  diff_test_res <- do.call(rbind.data.frame, diff_test_res)
  diff_test_res$qval <- p.adjust(diff_test_res$pval, method="BH")
  
  diff_test_res <- merge(diff_test_res, fData(cds), by="row.names")
  row.names(diff_test_res) <- diff_test_res[, 1] #remove the first column and set the row names to the first column
  diff_test_res[, 1] <- NULL 

  diff_test_res
}


#' Peform the branching test
#'
#' This function is used to perform the branching test to ask the question about whether or not the genes under tests are significant 
#' lineage dependent genes. If stretch equal to TRUE, each lineage is firstly stretched into maturation level 0-100 and the progenitor 
#' cells are duplicated and assigned to each lineage. This test can be used to detect lineage dependent genes. 
#'
#' @param cds a CellDataSet object upon which to perform this operation
#' @param fullModelFormulaStr a formula string specifying the full model in differential expression tests (i.e. likelihood ratio tests) for each gene/feature.
#' @param reducedModelFormulaStr a formula string specifying the reduced model in differential expression tests (i.e. likelihood ratio tests) for each gene/feature.
#' @param lineage_states ids for the immediate branch lineage which obtained from lineage construction based on MST
#' @param relative_expr a logic flag to determine whether or not the relative gene expression should be used
#' @param stretch  a logic flag to determine whether or not each lineage should be stretched
#' @param pseudocount pseudo count added before fitting the spline curves 
#' @param weighted  A logic flag to determine whether or not we should use the navie logLikelihood weight scheme for the duplicated progenitor cells
#' @param cores the number of cores to be used while testing each gene for differential expression
#' @param lineage_labelsthe name for each lineage, for example, AT1 or AT2  
#' @param gene_names gene names used to make the bifurcation plots for two genes. 
#' @return a data frame containing the p values and q-values from the likelihood ratio tests on the parallel arrays of models.
#' @export
#'
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
#' Calculate the area between TWO fitted lineage trajectories  
#' 
#' This function is used to calculate the ABC score based on the the nature spline curves fitted for each lineage. ABC score is used to 
#' quantify the magnitude of divergence between two lineages. By default, the ABC score is the area between two fitted spline curves. 
#' The ABC score can be used to rank gene divergence. When coupled with p-val calculated from the branchTest, it can be used to identify
#' potential major regulators for lineage bifurcation. 
#'
#' @param cds a CellDataSet object upon which to perform this operation
#' @param fullModelFormulaStr a formula string specifying the full model in differential expression tests (i.e. likelihood ratio tests) for each gene/feature.
#' @param reducedModelFormulaStr a formula string specifying the reduced model in differential expression tests (i.e. likelihood ratio tests) for each gene/feature.
#' @param ABC_method the method used to calculate the ABC scores. It can be either one of "integral", "global_normalization", "local_normalization", "four_values", "ILRs".
#' We use "integral" by default, which is defined as the area between two spline curves. "global_normalization" or "local_normalization" are similar measures between normalized by the global maximum or current maximum. 
#' "ILRs" is similar to calculate the Instant Log Ratio used in calILRs function  
#' @param branchTest a logic flag to determine whether or not to perform the branchTest inside the function. Because of the long computations, we recommend to first perform the branchTest and then calculate the ABCs for the genes of interests. Otherwise the ABCs will be appended to the last column of branchTest results. 
#' @param lineage_states ids for the immediate branch lineage which obtained from lineage construction based on MST (only two lineages are allowed for this function)
#' @param relative_expr a logic flag to determine whether or not the relative gene expression should be used
#' @param stretch a logic flag to determine whether or not each lineage should be stretched
#' @param pseudocount pseudo count added before fitting the spline curves 
#' @param cores the number of cores to be used while testing each gene for differential expression
#' @param weighted A logic flag to determine whether or not we should use the navie logLikelihood weight scheme for the duplicated progenitor cells
#' @param min_expr the lower limit for the expressed gene
#' @param integer_expression the logic flag to determine whether or not the integer numbers are used for calculating the ABCs. Default is False. 
#' @param num number of points on the fitted lineage trajectories used for calculating the ABCs. Default is 5000. 
#' @param lineage_labels the name for each lineage, for example, AT1 or AT2  
#' @return a data frame containing the p values and q-values from the likelihood ratio tests on the parallel arrays of models.
#' @export
calABCs <- function(cds, fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3)*Lineage",
                       reducedModelFormulaStr = "~sm.ns(Pseudotime, df = 3)", 
                       ABC_method = 'integral', 
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
  res_ABC <- mcesApply(new_cds, 1, function(x, modelFormulaStr, ABC_method, expressionFamily, relative_expr, disp_func, pseudocount, num, lineage_states, lineage_labels, ...) {
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
          disp_guess <- calulate_NB_dispersion_hint(disp_func, round(orig_x))
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

    #generalize the algoirthm to work for non-Lineage variables: 
    formula_all_variables <- all.vars(as.formula(fullModelFormulaStr))
    if(formula_all_variables[1] != 'Pseudotime')
      stop("Pseudotime is not in your model formula!")
    if(formula_all_variables[2] != 'Lineage'){
      warning("Warning: Lineage is not in your model formula!")      
      lineage <- pData(new_cds)[, formula_all_variables[2]]
    }

    newdata_df_ori <- data.frame(Pseudotime = seq(0, 100, length.out = num), formula_all_variables = as.factor(rep(lineage[1], num)))
    colnames(newdata_df)[2] <- formula_all_variables[2] #generalize the function to handle arbitrary variables
    predictBranchOri <- VGAM::predict(FM_fit, newdata = newdata_df, type = 'response')

    predictBranchOthers <- lapply(lineage[2:length(lineage)], function(x) {
        newdata_df_other <- data.frame(Pseudotime = seq(0, 100, length.out = num), formula_all_variables = as.factor(rep(x, num)))
        colnames(newdata_df)[2] <- formula_all_variables[2]

        VGAM::predict(FM_fit, newdata = newdata_df_other, type = 'response')        
    })
    ABCs <- lapply(predictBranchOthers, function(x){
        avg_delta_x <- ((predictBranchOri - x)[1:(num - 1)] + (predictBranchOri - x)[2:(num)]) / 2
        step <- (100 / (num - 1))

        if(ABC_method == 'integral'){
          res <- round(sum(avg_delta_x * step), 3)
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
        else if (ABC_method == "ILRs") {
                  str_logfc_df <- log2((predictBranchOri + 1)/(x + 
                    1))
                  res <- sum(str_logfc_df)
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
    pseudocount, num, lineage_states, lineage_labels, ...)

  if(ABC_method %in% c('integral', 'global_normalization', 'local_normalization', 'ILRs')) {
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
  else if(ABC_method %in% c('four_values')) {
    ABCs_res <- do.call(rbind, lapply(res_ABC,  unlist))
  }
  
  return(ABCs_res)
}

#' Calculate the Instant Log Ratio between two branching lineages
#' 
#' This function is used to calculate the Instant Log Ratio between two branching lineages which can be used to prepare the heatmap demonstrating the lineage gene expression divergence hirearchy. If "stretch" is specifified, each  
#' lineage will be firstly stretched into maturation level from 0-100. Since the results when we use "stretching" are always better and 
#' IRLs for non-stretched spline curves are often mismatched, we may only turn down "non-stretch" functionality in future versions. Then, we fit two separate nature spline curves for each 
#' individual linages. The log-ratios of the value on each spline curve corresponding to each lineages are calculated, which can be  
#' used as a measure for the magnitude of divergence between two branching lineages. 
#'
#' @param cds CellDataSet for the experiment
#' @param Lineage The column in pData used for calculating the ILRs (If not equal to "Lineage", a warning will report)
#' @param lineage_states The states for two branching lineages
#' @param cores Number of cores when fitting the spline curves
#' @param trend_formula the model formula to be used for fitting the expression trend over pseudotime
#' @param ILRs_limit the minimum Instant Log Ratio used to make the heatmap plot
#' @param relative_expr A logic flag to determine whether or not the relative expressed should be used when we fitting the spline curves 
#' @param weighted A logic flag to determine whether or not we should use the navie logLikelihood weight scheme for the duplicated progenitor cells
#' @param label_by_short_name label the rows of the returned matrix by gene_short_name (TRUE) or feature id (FALSE)
#' @param useVST A logic flag to determine whether or not the Variance Stablization Transformation should be used to stablize the gene expression.
#' When VST is used, the difference between two lineages are used instead of the log-ratio.
#' @param round_exprs A logic flag to determine whether or not the expression value should be rounded into integer
#' @param pseudocount pseudo count added before fitting the spline curves 
#' @param output_type A character either of "all" or "after_bifurcation". If "after_bifurcation" is used, only the time points after the bifurcation point will be selected
#' @param file the name for storing the data. Since the calculation of the Instant Log Ratio is very time consuming, so by default the result will be stored
#' @return a ggplot2 plot object
#' @import ggplot2
#' @importFrom reshape2 melt
#' @export 
#' 
calILRs <- function (cds = cds,
    Lineage = 'Lineage', 
    lineage_states = c(2, 3), 
    stretch = T, 
    cores = detectCores(), 
    trend_formula = "~sm.ns(Pseudotime, df = 3)", 
    ILRs_limit = 3, 
    relative_expr = TRUE, 
    weighted = FALSE, 
    label_by_short_name = TRUE, 
    useVST = FALSE, 
    round_exprs = FALSE, 
    pseudocount = 0, 
    output_type = c('all', 'after_bifurcation'), 
    file = "bifurcation_heatmap", verbose = FALSE, ...) {   

    cds_subset <- buildLineageBranchCellDataSet(cds = cds, lineage_states = lineage_states, 
        lineage_labels = NULL, method = "fitting", stretch = stretch, weighted = weighted, ...)

    #generate cds for branches 
    #we may also need to u
    if(Lineage != 'Lineage')
      warning("Warning: You didn't choose Lineage to calculate the ILRs")
    if(length(lineage_states) != 2)
      stop('calILRs can only work for two Lineages')

    cds_branchA <- cds_subset[, pData(cds_subset)[, Lineage] == 
        lineage_states[1]]
    cds_branchB <- cds_subset[, pData(cds_subset)[, Lineage] == 
        lineage_states[2]]

    #fit nature spline curve for each branch
    branchA_full_model_fits <- fitModel(cds_branchA, modelFormulaStr = trend_formula, 
        cores = cores, relative_expr = relative_expr, pseudocount = pseudocount)
    branchB_full_model_fits <- fitModel(cds_branchB, modelFormulaStr = trend_formula, 
        cores = cores, relative_expr = relative_expr, pseudocount = pseudocount)

    t_rng <- range(pData(cds_branchA)$Pseudotime)

    str_new_cds_branchA <- data.frame(Pseudotime = seq(0, max(pData(cds_branchA)$Pseudotime), 
        length.out = 100))
    if(verbose)
      print(paste("Check the whether or not Pseudotime scaled from 0 to 100: ", sort(pData(cds_branchA)$Pseudotime)))
    str_new_cds_branchB <- data.frame(Pseudotime = seq(0, max(pData(cds_branchB)$Pseudotime), 
        length.out = 100))
    if(verbose)
      print(paste("Check the whether or not Pseudotime scaled from 0 to 100: ", sort(pData(cds_branchA)$Pseudotime)))

    str_branchA_expression_curve_matrix <- responseMatrix(branchA_full_model_fits, 
        newdata = str_new_cds_branchA)
    str_branchB_expression_curve_matrix <- responseMatrix(branchB_full_model_fits, 
        newdata = str_new_cds_branchB)

    #VST for the fitted spline curves
    if (useVST) {
        str_branchA_expression_curve_matrix <- vstExprs(cds, 
            expr_matrix = str_branchA_expression_curve_matrix, 
            round_vals = round_exprs)
        str_branchB_expression_curve_matrix <- vstExprs(cds, 
            expr_matrix = str_branchB_expression_curve_matrix, 
            round_vals = round_exprs)

        #when VST is used, the difference between two lineages are defined as ILRs: 
        str_logfc_df <- str_branchA_expression_curve_matrix - 
            str_branchB_expression_curve_matrix
    }
    else {
        str_logfc_df <- log2((str_branchA_expression_curve_matrix + 1) / 
          (str_branchB_expression_curve_matrix + 1))
    }

    #should we label the ILRs matrix with gene short names? 
    if (label_by_short_name) {
        row.names(str_logfc_df) <- fData(cds[, ])$gene_short_name
    }

    #limit the range of ILRs
    str_logfc_df[which(str_logfc_df <= -ILRs_limit)] <- -ILRs_limit
    str_logfc_df[which(str_logfc_df >= ILRs_limit)] <- ILRs_limit

    if(output_type == 'after_bifurcation') {
      t_bifurcation_ori <- min(pData(cds[, c(which(pData(cds)$State == lineage_states[1]), #the pseudotime for the bifurcation point
        which(pData(cds)$State == lineage_states[2]))])$Pseudotime)
      t_bifurcation <- pData(cds_subset[, pData(cds)$Pseudotime == t_bifurcation_ori])$Pseudotime #corresponding stretched pseudotime

      if(stretch)
        bif_index <- as.integer(pData(cds_subset[,  pData(cds)$Pseudotime == t_bifurcation])$Pseudotime)
      else { #earliest bifurcation point on the original pseudotime scale (no stretching)
        bif_index <- as.integer(min(t_bifurcation / (max(pData(cds_branchA)$Pseudotime) / 100), 
                        t_bifurcation / (max(pData(cds_branchB)$Pseudotime) / 100))) 
      }
      #select only ILRs data points after the bifuration point
      
      str_logfc_df[, bif_index:100] <- str_logfc_df
    }

    if(!is.null(file)) #save the data file calculated since it will take a lot time to generate
      save(str_logfc_df, str_branchA_expression_curve_matrix, str_branchB_expression_curve_matrix, 
          file = file)

    return(str_logfc_df)
}

#' Detect the maturation time point where the gene expression starts to diverge 
#' 
#' This function is used to determine the bifurcation point for the gene expression between two distinct biological processes.
#' For processes we can not distinguish between lineages (or phenotype groups, like knockout VS un-knockout), this function will 
#' only detect bifurcation points after the inferenced bifurcatioin from the PQ-tree. The 
#'
#' @param str_log_df the ILRs dataframe calculated from calILRs function. If this data.frame is provided, all the following parameters are ignored. Note that we need to only use the ILRs after the bifurcation point if we duplicated the progenitor cell state.
#' @param div_threshold the ILR value used to determine the earliest divergence time point
#' @param detect_all a logic flag to determine whether or not genes without ILRs pass the threshold will still report a bifurcation point
#' @param cds CellDataSet for the experiment
#' @param Lineage The column in pData used for calculating the ILRs (If not equal to "Lineage", a warning will report)
#' @param lineage_states The states for two branching lineages
#' @param cores Number of cores when fitting the spline curves
#' @param trend_formula the model formula to be used for fitting the expression trend over pseudotime
#' @param ILRs_limit the minimum Instant Log Ratio used to make the heatmap plot
#' @param relative_expr A logic flag to determine whether or not the relative expressed should be used when we fitting the spline curves 
#' @param weighted A logic flag to determine whether or not we should use the navie logLikelihood weight scheme for the duplicated progenitor cells
#' @param label_by_short_name label the rows of the returned matrix by gene_short_name (TRUE) or feature id (FALSE)
#' @param useVST A logic flag to determine whether or not the Variance Stablization Transformation should be used to stablize the gene expression.
#' When VST is used, the difference between two lineages are used instead of the log-ratio.
#' @param round_exprs A logic flag to determine whether or not the expression value should be rounded into integer
#' @param pseudocount pseudo count added before fitting the spline curves 
#' @param output_type A character either of "all" or "after_bifurcation". If "after_bifurcation" is used, only the time points after the bifurcation point will be selected. Note that, if Lineage is set to "Lineage", we will only use "after_bifurcation" since we duplicated the progenitor cells and the bifurcation should only happen after the largest mature level from the progenitor cells
#' @param file the name for storing the data. Since the calculation of the Instant Log Ratio is very time consuming, so by default the result will be stored
#' @return a vector containing the time for the bifurcation point with gene names for each value
#' @importFrom reshape2 melt
#' @export 
#' 
detectBifurcationPoint <- function(str_log_df = NULL, str_norm_div_df = NULL, ILRs_threshold = 0.5, div_threshold = 0.2, detect_all = T,
cds = cds,
Lineage = 'Lineage',
lineage_states = c(2, 3),
stretch = T,
cores = detectCores(),
trend_formula = "~sm.ns(Pseudotime, df = 3)",
ILRs_limit = 3,
relative_expr = TRUE,
weighted = FALSE,
label_by_short_name = TRUE,
useVST = FALSE,
round_exprs = FALSE,
pseudocount = 0,
output_type = c('all', 'after_bifurcation'),
file = "bifurcation_heatmap", verbose = FALSE, ...) {
    if(is.null(str_log_df)) {
        if(Lineage == 'Lineage') output_type = 'after_bifurcation'
        
        str_log_df <- calILRs(cds = cds,
        Lineage,
        lineage_states,
        stretch,
        cores,
        trend_formula,
        ILRs_limit,
        relative_expr,
        weighted,
        label_by_short_name,
        useVST,
        round_exprs,
        pseudocount,
        output_type = output_type,
        file, verbose, ...)
    }
    
    if(!is.null(str_log_df) & is.null(str_norm_div_df)) {
        bifurcation_time <- apply(str_log_df, 1, function(x) {
            # deriv <- diff(x) the ILRs are smooth, so use min is fine
            index <- NA
            if(any(which(abs(x) > div_threshold)))
            index <- min(which(abs(x) > ILRs_threshold)) * sign(sum(x))
            else if(detect_all & all(!is.na(x))) index <-  min(which(abs(x) == max(abs(x))))
            index
        }
        ) #detect the earliest divergence point
    }
    else if(is.null(str_log_df) & !is.null(str_norm_div_df)) {
        bifurcation_time <- apply(str_norm_div_df, 1, function(x) {
            # deriv <- diff(x) the ILRs are smooth, so use min is fine
            index <- NA
            if(any(which(abs(x) > div_threshold)))
            index <- min(which(abs(x) > ILRs_threshold)) * sign(sum(x))
            else if(detect_all & all(!is.na(x))) index <-  min(which(abs(x) == max(abs(x))))
            index
        }
        ) #detect the earliest divergence point
    }
    else { #use both of div threshold / ILRs threshold
        str_log_df_list <- split(str_log_df, row.names(str_log_df))
        str_norm_div_df_list <- split(str_norm_div_df, row.names(str_norm_div_df))
        
        save(logic_tmp, str_norm_div_df_list,  str_log_df_list, file = 'logic_tmp')
        bifurcation_time <- mapply(function(x, y, ILRs_thresh = ILRs_threshold, div_thresh = div_threshold) {
            # deriv <- diff(x) the ILRs are smooth, so use min is fine
            index <- NA
            
            logic_tmp <- (abs(x) > ILRs_thresh) & (abs(y) > div_thresh)
            if(all(is.logical(logic_tmp)) & any(logic_tmp, na.rm = T))
            index <- min(which(abs(x) > ILRs_thresh & abs(y) > div_thresh)) * sign(sum(x))
            else if(detect_all & all(!is.na(x))) index <-  min(which(abs(x) == max(abs(x))))
            index
        }, str_log_df_list, str_norm_div_df_list
        ) #detect the earliest divergence point
    }
    # print(bifurcation_time)
    # str_norm_div_df
    
    names(bifurcation_time) <- row.names(str_log_df)
    
    return(bifurcation_time)
}

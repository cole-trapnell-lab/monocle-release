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
      full_model_fit <- VGAM::vglm(as.formula(fullModelFormulaStr), family=expressionFamily, weights=weights, checkwz=TRUE)
      reduced_model_fit <- VGAM::vglm(as.formula(reducedModelFormulaStr), family=expressionFamily, weights=weights, checkwz=TRUE)                         
    }else{
      full_model_fit <- suppressWarnings(VGAM::vglm(as.formula(fullModelFormulaStr), family=expressionFamily, weights=weights))
      reduced_model_fit <- suppressWarnings(VGAM::vglm(as.formula(reducedModelFormulaStr), family=expressionFamily, weights=weights))                    
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
	f_expression <- x
        disp_guess <- calulate_QP_dispersion_hint(disp_func, x_orig)
        backup_expression_family <- poissonff(dispersion=disp_guess)
    }else if (expressionFamily@vfamily %in% c("gaussianff", "uninormal")){
      backup_expression_family <- NULL
    }else if (expressionFamily@vfamily %in% c("binomialff")){
      backup_expression_family <- NULL
    }else{
      backup_expression_family <- NULL
    }
    if (is.null(backup_expression_family) == FALSE){
      test_res <- tryCatch({
      if (verbose){
        full_model_fit <- VGAM::vglm(as.formula(fullModelFormulaStr), family=backup_expression_family, weights=weights, checkwz=TRUE)
        reduced_model_fit <- VGAM::vglm(as.formula(reducedModelFormulaStr), family=backup_expression_family, weights=weights, checkwz=TRUE)                         
      }else{
        full_model_fit <- suppressWarnings(VGAM::vglm(as.formula(fullModelFormulaStr), family=backup_expression_family, weights=weights, checkwz=TRUE))
        reduced_model_fit <- suppressWarnings(VGAM::vglm(as.formula(reducedModelFormulaStr), family=backup_expression_family, weights=weights, checkwz=TRUE))                    
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

  if("Lineage" %in% all.vars(terms(as.formula(fullModelFormulaStr)))) {
     cds_subset <- buildLineageBranchCellDataSet(cds = cds, lineage_states = lineage_states,
     lineage_labels = lineage_labels, method = method, stretch = stretch,
     weighted = weighted, ...)
  }
  else
    cds_subset <- cds

  branchTest_res <- differentialGeneTest(cds_subset, 
                                         fullModelFormulaStr = fullModelFormulaStr, 
                                         reducedModelFormulaStr = reducedModelFormulaStr, 
                                         cores = cores, 
                                         relative_expr = relative_expr, 
                                         weights = pData(cds_subset)$weight,
                                         pseudocount = pseudocount)
  
  return(branchTest_res)
}

#add genSmoothCurves function: 
genSmoothCurves <- function(cds, cores = 1, fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3)", weights = NULL, 
                        relative_expr = T, pseudocount = 0, new_data) { 
    
    expressionFamily <- cds@expressionFamily

    if(cores > 1) {
        expression_curve_matrix <- mcesApply(cds, 1, function(x, fullModelFormulaStr, expressionFamily, relative_expr, pseudocount, new_data){
            environment(fit_model_helper) <- environment()
            environment(responseMatrix) <- environment()
            model_fits <- fit_model_helper(x, modelFormulaStr = fullModelFormulaStr, expressionFamily = expressionFamily, weights = weights, 
                                       relative_expr = relative_expr, pseudocount = pseudocount, disp_func = cds@dispFitInfo[['blind']]$disp_func)
            if(is.null(model_fits))
                expression_curve <- matrix(rep(NA, length(x)), nrow = 1)
            else
                expression_curve <- responseMatrix(list(model_fits), newdata = new_data)

            }, required_packages=c("BiocGenerics", "VGAM", "plyr"), cores=cores, 
            fullModelFormulaStr = fullModelFormulaStr, expressionFamily = expressionFamily, relative_expr = relative_expr, pseudocount = pseudocount, new_data = new_data
            )
    }
    else {
        expression_curve_matrix <- esApply(cds, 1, function(x, fullModelFormulaStr, expressionFamily, relative_expr, pseudocount, new_data = new_data){
            environment(fit_model_helper) <- environment()
            environment(responseMatrix) <- environment()
            model_fits <- fit_model_helper(x, modelFormulaStr = fullModelFormulaStr, expressionFamily = expressionFamily, weights = weights, 
                                       relative_expr = relative_expr, pseudocount = pseudocount, disp_func = cds@dispFitInfo[['blind']]$disp_func)
            if(is.null(model_fits))
                expression_curve <- matrix(rep(NA, nrow(new_data)), nrow = 1)
            else
                expression_curve <- responseMatrix(list(model_fits), new_data)

            }, 
            fullModelFormulaStr = fullModelFormulaStr, expressionFamily = expressionFamily, relative_expr = relative_expr, pseudocount = pseudocount, new_data = new_data
            )
    }

    t(expression_curve_matrix)
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
calABCs <- function(cds, trajectory_type = "Lineage", 
  trajectory_states = c(2, 3),
  branchTest = FALSE, 
  relative_expr = TRUE, 
  fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3)*Lineage",
  stretch = TRUE, 
  pseudocount = 0, 
  cores = 1, 
  weighted = TRUE, 
  verbose = F,
  min_expr = 0.5, 
  integer_expression = FALSE, 
  num = 5000, 
  lineage_labels = NULL,
  ABC_method = c('integral'), #'global_normalization', 'local_normalization', 'four_values', 'ILRs'),
  ...){
    if (length(trajectory_states) != 2)
    stop("Sorry, this function only supports the calculation of ABCs between TWO lineage trajectories")
    
    if(trajectory_type == "Lineage") {
        cds_subset <- buildLineageBranchCellDataSet(cds = cds, #lineage_states = trajectory_states,
          lineage_labels = lineage_labels, method = method, stretch = stretch,
          weighted = weighted, ...)
        overlap_rng <- c(0, max(pData(cds_subset)$Pseudotime))
    }
    else{
        pd <- pData(cds)
        range_df <- ddply(pd, .(get(trajectory_type)), function(x) {
            range(x$Pseudotime)
        })
        
        overlap_rng <- c(max(range_df[, 2]), min(range_df[, 3])) #calculate the overlapping pseudotime range
        # cds <- cds[, row.names(subset(pd, Pseudotime > overlap_rng[1] & Pseudotime < overlap_rng[2]))]
        cds_subset <- cds
    }
    if (trajectory_type != "Lineage")
        warning("Warning: You didn't choose Lineage to calculate the ILRs")
    if (length(trajectory_states) != 2)
        stop("calILRs can only work for two Lineages")
    if(!all(trajectory_states %in% pData(cds_subset)[, trajectory_type]))
        stop("state(s) in trajectory_states are not included in trajectory_type")
    
    if(verbose)
      message(paste("the pseudotime range for the calculation of ILRs:", overlap_rng[1], overlap_rng[2], sep = ' '))
    
    cds_branchA <- cds_subset[, pData(cds_subset)[, trajectory_type] ==
      trajectory_states[1]]
    cds_branchB <- cds_subset[, pData(cds_subset)[, trajectory_type] ==
      trajectory_states[2]]
    
    #use the genSmoothCurves to generate smooth curves for calculating the ABC scores:
    formula_all_variables <- all.vars(as.formula(fullModelFormulaStr))
    
    t_rng <- range(pData(cds_branchA)$Pseudotime)
    str_new_cds_branchA <- data.frame(Pseudotime = seq(overlap_rng[1], overlap_rng[2],
      length.out = num), Lineage = as.factor(trajectory_states[1]))
    colnames(str_new_cds_branchA)[2] <- formula_all_variables[2] #interaction term can be terms rather than Lineage
    if (verbose)
        print(paste("Check the whether or not Pseudotime scaled from 0 to 100: ",
          sort(pData(cds_branchA)$Pseudotime)))
    str_new_cds_branchB <- data.frame(Pseudotime = seq(overlap_rng[1], overlap_rng[2],
      length.out = num), Lineage = as.factor(trajectory_states[2]))
    colnames(str_new_cds_branchB)[2] <- formula_all_variables[2]
    
    if (verbose)
      print(paste("Check the whether or not Pseudotime scaled from 0 to 100: ",
        sort(pData(cds_branchB)$Pseudotime)))
    
    str_branchAB_expression_curve_matrix <- genSmoothCurves(cds_subset, cores=cores, fullModelFormulaStr = fullModelFormulaStr,  weights = pData(cds_subset)$weight,
                    relative_expr = relative_expr, pseudocount = pseudocount, new_data = rbind(str_new_cds_branchA, str_new_cds_branchB))
    
    str_branchA_expression_curve_matrix <- str_branchAB_expression_curve_matrix[, 1:num]
    str_branchB_expression_curve_matrix <- str_branchAB_expression_curve_matrix[, (num + 1):(2 * num)]
    
    ABCs_res <- str_branchA_expression_curve_matrix - str_branchB_expression_curve_matrix
    
    ABCs_res <- apply(ABCs_res, 1, function(x, num, ABC_method) {
        avg_delta_x <- (x[1:(num - 1)] + x[2:(num)])/2
        step <- (100/(num - 1))
        
        if (ABC_method == "integral") {
            res <- round(sum(avg_delta_x * step), 3)
        }
        else if (ABC_method == "global_normalization") {
            max <- max(max(predictBranchOri), max(x))
            res <- round(sum(avg_delta_x/max * step), 3)
        }
        else if (ABC_method == "local_normalization") {
            pair_wise_max <- apply(data.frame(x = x, y = predictBranchOri),
            1, max)
            res <- round(sum((((predictBranchOri - x)/pair_wise_max)[1:(num -
            1)] + ((predictBranchOri - x)/pair_wise_max)[2:(num)])/2 *
            step), 3)
        }
        else if (ABC_method == "four_values") {
            ori_ABCs <- round(sum((x[1:(num - 1)] + x[2:(num)])/2 *
            step), 3)
            other_ABCs <- round(sum((predictBranchOri[1:(num -
            1)] + predictBranchOri[2:(num)])/2 * step),
            3)
            ori_ABCs_H <- round(sum(avg_delta_x[avg_delta_x >
            0] * step), 3)
            other_ABCs_H <- round(sum(avg_delta_x[avg_delta_x <
            0] * step), 3)
            res <- c(ori_ABCs = ori_ABCs, other_ABCs = other_ABCs,
            ori_ABCs_H = ori_ABCs_H, other_ABCs_H = other_ABCs_H)
        }
        else if (ABC_method == "ILRs") {
            str_logfc_df <- log2((predictBranchOri + 1)/(x +
            1))
            res <- sum(str_logfc_df)
        }
        return(res)}, num = num, ABC_method = ABC_method
    )
    
    ABCs_res <- merge(ABCs_res, fData(cds), by = "row.names")
    row.names(ABCs_res) <- ABCs_res[, 1]
    ABCs_res[, 1] <- NULL
    colnames(ABCs_res)[1] <- "ABCs"
    
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
#' @param fullModelFormulaStr the model formula to be used for fitting the expression trend over pseudotime
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
  trajectory_type = "Lineage", 
  trajectory_states = c(2, 3), 
  lineage_labels = NULL, 
  stretch = T, 
  cores = detectCores(), 
  fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3)*Lineage",
  ILRs_limit = 3, 
  relative_expr = TRUE, 
  weighted = FALSE, 
  label_by_short_name = TRUE,
  useVST = FALSE, 
  round_exprs = FALSE, 
  pseudocount = 0, 
  output_type = c("all", "after_bifurcation"), 
  file = "bifurcation_heatmap", 
  return_all = F, 
  verbose = FALSE, ...){
    
    if(trajectory_type == "Lineage") {
        cds_subset <- buildLineageBranchCellDataSet(cds = cds, #lineage_states = trajectory_states,
        lineage_labels = lineage_labels, method = method, stretch = stretch,
        weighted = weighted, ...)
        overlap_rng <- c(0, max(pData(cds_subset)$Pseudotime))
    }
    else{
        pd <- pData(cds)
        range_df <- ddply(pd, .(get(trajectory_type)), function(x) {
            range(x$Pseudotime)
        })
        
        overlap_rng <- c(max(range_df[, 2]), min(range_df[, 3])) #calculate the overlapping pseudotime range
        # cds <- cds[, row.names(subset(pd, Pseudotime > overlap_rng[1] & Pseudotime < overlap_rng[2]))]
        cds_subset <- cds
    }
    if (trajectory_type != "Lineage")
      warning("Warning: You didn't choose Lineage to calculate the ILRs")
    if (length(trajectory_states) != 2)
      stop("calILRs can only work for two Lineages")
    if(!all(trajectory_states %in% pData(cds_subset)[, trajectory_type]))
      stop("state(s) in trajectory_states are not included in trajectory_type")
    
    if(verbose)
      message(paste("the pseudotime range for the calculation of ILRs:", overlap_rng[1], overlap_rng[2], sep = ' '))
    
    cds_branchA <- cds_subset[, pData(cds_subset)[, trajectory_type] ==
      trajectory_states[1]]
    cds_branchB <- cds_subset[, pData(cds_subset)[, trajectory_type] ==
      trajectory_states[2]]
    
    formula_all_variables <- all.vars(as.formula(fullModelFormulaStr))
    
    t_rng <- range(pData(cds_branchA)$Pseudotime)
    str_new_cds_branchA <- data.frame(Pseudotime = seq(overlap_rng[1], overlap_rng[2],
                            length.out = 100), Lineage = as.factor(trajectory_states[1]))
    if (verbose)
      message(paste("Check the whether or not Pseudotime scaled from 0 to 100: ",
        sort(pData(cds_branchA)$Pseudotime)))
    colnames(str_new_cds_branchA)[2] <- formula_all_variables[2] #interaction term can be terms rather than Lineage
    
    str_new_cds_branchB <- data.frame(Pseudotime = seq(overlap_rng[1], overlap_rng[2],
                            length.out = 100), Lineage = as.factor(trajectory_states[2]))
    if (verbose)
     message(paste("Check the whether or not Pseudotime scaled from 0 to 100: ",
        sort(pData(cds_branchB)$Pseudotime)))
    
    colnames(str_new_cds_branchB)[2] <- formula_all_variables[2] #interaction term can be terms rather than Lineage
    
    str_branchAB_expression_curve_matrix <- genSmoothCurves(cds_subset, cores=cores, fullModelFormulaStr = fullModelFormulaStr,
                      relative_expr = relative_expr, pseudocount = pseudocount, new_data = rbind(str_new_cds_branchA, str_new_cds_branchB))
    
    str_branchA_expression_curve_matrix <- str_branchAB_expression_curve_matrix[, 1:nrow(str_new_cds_branchA)]
    str_branchB_expression_curve_matrix <- str_branchAB_expression_curve_matrix[, 
                                                  (nrow(str_new_cds_branchA) + 1):(nrow(str_new_cds_branchA) + nrow(str_new_cds_branchB))]
    
    if (useVST) {
        str_branchA_expression_curve_matrix <- vstExprs(cds,
          expr_matrix = str_branchA_expression_curve_matrix,
          round_vals = round_exprs)
        str_branchB_expression_curve_matrix <- vstExprs(cds,
          expr_matrix = str_branchB_expression_curve_matrix,
          round_vals = round_exprs)
        str_logfc_df <- str_branchA_expression_curve_matrix -
        str_branchB_expression_curve_matrix
    }
    else {
        str_logfc_df <- log2((str_branchA_expression_curve_matrix +
          1)/(str_branchB_expression_curve_matrix + 1))
    }
    if (label_by_short_name) {
        str_logfc_df <- t(str_logfc_df)
        row.names(str_logfc_df) <- fData(cds[, ])$gene_short_name
    }
    str_logfc_df[which(str_logfc_df <= -ILRs_limit, arr.ind = T)] <- -ILRs_limit
    str_logfc_df[which(str_logfc_df >= ILRs_limit, arr.ind = T)] <- ILRs_limit
    if (output_type == "after_bifurcation") {
        t_bifurcation_ori <- min(pData(cds[, c(which(pData(cds)$State ==
            lineage_states[1]), which(pData(cds)$State == lineage_states[2]))])$Pseudotime)
        t_bifurcation <- pData(cds_subset[, pData(cds)$Pseudotime ==
        t_bifurcation_ori])$Pseudotime
        
        if (stretch)
          bif_index <- as.integer(pData(cds_subset[, pData(cds)$Pseudotime ==
                t_bifurcation])$Pseudotime)
        else {
            bif_index <- as.integer(min(t_bifurcation/(max(pData(cds_branchA)$Pseudotime)/100),
            t_bifurcation/(max(pData(cds_branchB)$Pseudotime)/100)))
        }
        str_logfc_df[, bif_index:100] <- str_logfc_df
    }
    if (!is.null(file))
      save(str_logfc_df, str_branchA_expression_curve_matrix, str_branchB_expression_curve_matrix, file = file)

    if(return_all) {
        rMax <- function(df) {apply(df, 1, function(x) if(all(is.na(x))) NA else max(abs(x), na.rm = T))} #calculate row max

        str_raw_div_df <- str_branchA_expression_curve_matrix - str_branchB_expression_curve_matrix
        str_norm_div_df <- str_raw_div_df / rMax(str_raw_div_df) #calculate normalized divergence

        log_str_raw_div_df <- log2((str_branchA_expression_curve_matrix + .1)/(str_branchB_expression_curve_matrix + .1))
        norm_str_logfc_df <- str_logfc_df / rMax(log_str_raw_div_df) #calculate normalized divergence

        return(list(str_logfc_df = str_logfc_df, norm_str_logfc_df = norm_str_logfc_df,
                str_norm_div_df = str_norm_div_df, str_raw_div_df = str_raw_div_df, str_branchA_expression_curve_matrix = str_branchA_expression_curve_matrix, 
                str_branchB_expression_curve_matrix = str_branchB_expression_curve_matrix))
    }
    else
        return(str_logfc_df)
}

#' Detect the maturation time point where the gene expression starts to diverge 
#' 
#' This function is used to determine the bifurcation point for the gene expression between two distinct biological processes.
#' For processes we can not distinguish between lineages (or phenotype groups, like knockout VS un-knockout), this function will 
#' only detect bifurcation points after the inferenced bifurcatioin from the PQ-tree. The 
#'
#' @param str_log_df the ILRs dataframe calculated from calILRs function. If this data.frame is provided, all the following parameters are ignored. Note that we need to only use the ILRs after the bifurcation point if we duplicated the progenitor cell state.
#' @param ILRs_threshold the ILR value used to determine the earliest divergence time point
#' @param detect_all a logic flag to determine whether or not genes without ILRs pass the threshold will still report a bifurcation point
#' @param cds CellDataSet for the experiment
#' @param Lineage The column in pData used for calculating the ILRs (If not equal to "Lineage", a warning will report)
#' @param lineage_states The states for two branching lineages
#' @param cores Number of cores when fitting the spline curves
#' @param fullModelFormulaStr the model formula to be used for fitting the expression trend over pseudotime
#' @param ILRs_limit the minimum Instant Log Ratio used to make the heatmap plot
#' @param relative_expr A logic flag to determine whether or not the relative expressed should be used when we fitting the spline curves 
#' @param weighted A logic flag to determine whether or not we should use the navie logLikelihood weight scheme for the duplicated progenitor cells
#' @param label_by_short_name label the rows of the returned matrix by gene_short_name (TRUE) or feature id (FALSE)
#' @param useVST A logic flag to determine whether or not the Variance Stablization Transformation should be used to stablize the gene expression.
#' When VST is used, the difference between two lineages are used instead of the log-ratio.
#' @param round_exprs A logic flag to determine whether or not the expression value should be rounded into integer
#' @param pseudocount pseudo count added before fitting the spline curves 
#' @param output_type A character either of "all" or "after_bifurcation". If "after_bifurcation" is used, only the time points after the bifurcation point will be selected. Note that, if Lineage is set to "Lineage", we will only use "after_bifurcation" since we duplicated the progenitor cells and the bifurcation should only happen after the largest mature level from the progenitor cells
#' @param return_cross_point A logic flag to determine whether or not only return the cross point 
#' @param file the name for storing the data. Since the calculation of the Instant Log Ratio is very time consuming, so by default the result will be stored
#' @return a vector containing the time for the bifurcation point with gene names for each value
#' @importFrom reshape2 melt
#' @export 
#' 
detectBifurcationPoint <- function(str_log_df = NULL, 
  ILRs_threshold = 0.1, 
  detect_all = T,
  cds = cds,
  Lineage = 'Lineage',
  lineage_states = c(2, 3),
  stretch = T,
  cores = detectCores(),
  fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3)",
  ILRs_limit = 3,
  relative_expr = TRUE,
  weighted = FALSE,
  label_by_short_name = TRUE,
  useVST = FALSE,
  round_exprs = FALSE,
  pseudocount = 0,
  output_type = c('all', 'after_bifurcation'),
  return_cross_point = T, 
  file = "bifurcation_heatmap", verbose = FALSE, ...) {
    if(is.null(str_log_df)) {
        if(Lineage == 'Lineage') output_type = 'after_bifurcation'
        
        str_log_df <- calILRs(cds = cds,
          Lineage,
          lineage_states,
          stretch,
          cores,
          fullModelFormulaStr,
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
    # rMax <- function(df) {apply(df, 1, function(x) if(all(is.na(x))) NA else max(abs(x), na.rm = T))}
    
    else {
        bifurcation_time <- apply(str_log_df, 1, function(x) {
            # deriv <- diff(x) the ILRs are smooth, so use min is fine
            index <- Inf
            
            #new algorithm to bifurcation time point:
            if(all(is.na(x))) {
                return(NA)
            }
            
            max_ind <- which(abs(x) == max(abs(x)))
            
            # return(max(max_ind))
            if(length(max_ind) > 1) {
                max_ind <- min(max_ind)
                warning('multiple maximal time points detected ', max_ind)
            }
            
            #detect the cross point
            inflection_point_tmp <- which(x[1:(length(x) - 1)] * x[2:length(x)] <= 0)
            
            if(all(max_ind <= inflection_point_tmp)) return(NA) #remove all zero values and genes monotonically goes down
            
            inflection_point <- max(inflection_point_tmp[inflection_point_tmp < max_ind])
            

            if(any(abs(x) > ILRs_threshold) & return_cross_point == T)
                return(inflection_point * sign(sum(x)))
            else
                return(Inf)

            rev_x <- rev(x[(inflection_point):max_ind])
            if(any(which(abs(rev_x) >= ILRs_threshold))){
                index_tmp <- max(which(abs(rev_x) > ILRs_threshold))
                index <- (max_ind - index_tmp + 1 ) * sign(sum(rev_x))
            }
            else if(detect_all & all(!is.na(rev_x))) {
                index <-  min(which(abs(rev_x) == max(abs(rev_x)))) * sign(sum(rev_x))
            }
            
            index
        })
    }
    
    # print(bifurcation_time)
    # str_norm_div_df
    
    names(bifurcation_time) <- row.names(str_log_df)
    
    return(bifurcation_time)
}


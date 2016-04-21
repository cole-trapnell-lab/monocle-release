#Helper function for parallel VGAM fitting
fit_model_helper <- function(x, 
                             modelFormulaStr, 
                             expressionFamily, 
                             relative_expr, 
                             disp_func=NULL, 
                             verbose=FALSE,
                             ...){
    modelFormulaStr <- paste("f_expression", modelFormulaStr,
        sep = "")
    orig_x <- x
    if (expressionFamily@vfamily %in% c("negbinomial", "negbinomial.size")) {
        if (relative_expr) {
            x <- x/Size_Factor
        }
        f_expression <- round(x)
        if (is.null(disp_func) == FALSE) {
            disp_guess <- calulate_NB_dispersion_hint(disp_func,
                round(orig_x))
            if (is.null(disp_guess) == FALSE && disp_guess >
                0 && is.na(disp_guess) == FALSE) {
                size_guess <- 1/disp_guess
                if (expressionFamily@vfamily == "negbinomial")
                  expressionFamily <- negbinomial(isize=1/disp_guess, ...)
                else
                  expressionFamily <- negbinomial.size(size=1/disp_guess, ...)
            }
        }
    }
    else if (expressionFamily@vfamily %in% c("gaussianff", "uninormal", "binomialff", "betabinomial", "betabinomialff")) {
        f_expression <- x
    }
    else {
        f_expression <- log10(x)
    }
    tryCatch({
        if (verbose) {
            FM_fit <- VGAM::vglm(as.formula(modelFormulaStr),
                family = expressionFamily, epsilon=1e-1)
        }
        else {
            FM_fit <- suppressWarnings(VGAM::vglm(as.formula(modelFormulaStr),
                family = expressionFamily, epsilon=1e-1))
        }
        FM_fit
    }, error = function(e) {
        print (e);
        # If we threw an exception, re-try with a simpler model.  Which one depends on
        # what the user has specified for expression family
        #print(disp_guess)
        backup_expression_family <- NULL
        if (expressionFamily@vfamily %in% c("negbinomial", "negbinomial.size")){
            disp_guess <- calulate_NB_dispersion_hint(disp_func, round(orig_x), expr_selection_func = max)
            backup_expression_family <- negbinomial()
        }else if (expressionFamily@vfamily %in% c("gaussianff", "uninormal")){
          backup_expression_family <- NULL
        }else if (expressionFamily@vfamily %in% c("binomialff", "betabinomial", "betabinomialff")){
          backup_expression_family <- NULL
        }else{
          backup_expression_family <- NULL
        }
        if (is.null(backup_expression_family) == FALSE){

          #FM_fit <- VGAM::vglm(as.formula(modelFormulaStr), family=backup_expression_family, trace=T, epsilon=1e-1, checkwz=F)
          test_res <- tryCatch({
          if (verbose){
            FM_fit <- VGAM::vglm(as.formula(modelFormulaStr), family=backup_expression_family, epsilon = 1e-1, checkwz= TRUE)
          }else{
            FM_fit <- suppressWarnings(VGAM::vglm(as.formula(modelFormulaStr), family=backup_expression_family, epsilon = 1e-1, checkwz = TRUE))
          }
            FM_fit
          }, 
          #warning = function(w) { FM_fit },
          error = function(e) { 
            #print (e);
            NULL
          })
          #print(test_res)
          test_res
        } else {
          #print(e); 
          NULL
        }
  })
}

#' Fits a model for each gene in a CellDataSet object.
#' @param cds the CellDataSet upon which to perform this operation
#' @param modelFormulaStr a formula string specifying the model to fit for the genes.
#' @param relative_expr Whether to fit a model to relative or absolute expression. Only meaningful for count-based expression data. If TRUE, counts are normalized by Size_Factor prior to fitting.
#' @param cores the number of processor cores to be used during fitting.
#' @return a list of VGAM model objects
#' @importFrom qlcMatrix rowMax
#' @export
#' @details
#' 
#' This function fits a vector generalized additive model (VGAM) from the VGAM package for each gene in a CellDataSet. 
#' By default, expression levels are modeled as smooth functions of the Pseudotime value of each 
#' cell. That is, expression is a function of progress through the biological process.  More complicated formulae can be provided to account for
#' additional covariates (e.g. day collected, genotype of cells, media conditions, etc).
fitModel <- function(cds,
                     modelFormulaStr="~sm.ns(Pseudotime, df=3)",
                     relative_expr=TRUE,
                     cores=1){
  if (cores > 1){
    f<-mcesApply(cds, 1, fit_model_helper, required_packages=c("BiocGenerics", "Biobase", "VGAM", "plyr", "Matrix"), cores=cores, 
                 modelFormulaStr=modelFormulaStr, 
                 expressionFamily=cds@expressionFamily,
                 relative_expr=relative_expr,
                 disp_func=cds@dispFitInfo[["blind"]]$disp_func)
    f
  }else{
    f<-smartEsApply(cds,1,fit_model_helper, 
                    convert_to_dense=TRUE,
               modelFormulaStr=modelFormulaStr, 
               expressionFamily=cds@expressionFamily,
               relative_expr=relative_expr,
               disp_func=cds@dispFitInfo[["blind"]]$disp_func)
    f
  }
}

#' Response values
#' 
#' Generates a matrix of response values for a set of fitted models
#' @param models a list of models, e.g. as returned by fitModels()
#' @param newdata a dataframe used to generate new data for interpolation of time points
#' @param response_type the response desired, as accepted by VGAM's predict function
#' @param cores number of cores used for calculation
#' @return a matrix where each row is a vector of response values for a particular feature's model, and columns are cells.
#' @export
responseMatrix <- function(models, newdata = NULL, response_type="response", cores = detectCores()) {
    res_list <- mclapply(models, function(x) {
      if (is.null(x)) { NA } else {
          if (x@family@vfamily %in% c("negbinomial", "negbinomial.size")) {
              predict(x, newdata = newdata, type = response_type)
          } else if (x@family@vfamily %in% c("gaussianff")) {
              predict(x, newdata = newdata, type = response_type)
          }
          else {
              10^predict(x, newdata = newdata, type = response_type)
          }
      }
    }, mc.cores = cores)

    res_list_lengths <- lapply(res_list[is.na(res_list) == FALSE],
        length)
    stopifnot(length(unique(res_list_lengths)) == 1)
    num_na_fits <- length(res_list[is.na(res_list)])
    if (num_na_fits > 0) {
        na_matrix <- matrix(rep(rep(NA, res_list_lengths[[1]]),
            num_na_fits), nrow = num_na_fits)
        row.names(na_matrix) <- names(res_list[is.na(res_list)])
        non_na_matrix <- Matrix::t(do.call(cbind, lapply(res_list[is.na(res_list) ==
            FALSE], unlist)))
        row.names(non_na_matrix) <- names(res_list[is.na(res_list) ==
            FALSE])
        res_matrix <- rbind(non_na_matrix, na_matrix)
        res_matrix <- res_matrix[names(res_list), ]
    }
    else {
        res_matrix <- Matrix::t(do.call(cbind, lapply(res_list, unlist)))
        row.names(res_matrix) <- names(res_list[is.na(res_list) ==
            FALSE])
    }
    res_matrix
}

#' Response values
#' 
#' Generates a matrix of response values for a set of fitted models
#' @param models a list of models, e.g. as returned by fitModels()
#' @param residual_type the response desired, as accepted by VGAM's predict function
#' @param cores number of cores used for calculation
#' @return a matrix where each row is a vector of response values for a particular feature's model, and columns are cells.
residualMatrix <- function(models,  residual_type="response", cores = detectCores()) {
  res_list <- mclapply(models, function(x) {
    if (is.null(x)) { NA } else {
        resid(x, type = residual_type)
    }
  }, mc.cores = cores)
  
  res_list_lengths <- lapply(res_list[is.na(res_list) == FALSE],
                             length)
  stopifnot(length(unique(res_list_lengths)) == 1)
  num_na_fits <- length(res_list[is.na(res_list)])
  if (num_na_fits > 0) {
    na_matrix <- matrix(rep(rep(NA, res_list_lengths[[1]]),
                            num_na_fits), nrow = num_na_fits)
    row.names(na_matrix) <- names(res_list[is.na(res_list)])
    non_na_matrix <- Matrix::t(do.call(cbind, lapply(res_list[is.na(res_list) ==
                                                                FALSE], unlist)))
    row.names(non_na_matrix) <- names(res_list[is.na(res_list) ==
                                                 FALSE])
    res_matrix <- rbind(non_na_matrix, na_matrix)
    res_matrix <- res_matrix[names(res_list), ]
  }
  else {
    res_matrix <- Matrix::t(do.call(cbind, lapply(res_list, unlist)))
    row.names(res_matrix) <- names(res_list[is.na(res_list) ==
                                              FALSE])
  }
  res_matrix
}


#' Fit smooth spline curves and return the response matrix
#'
#' This function will fit smooth spline curves for the gene expression dynamics along pseudotime in a gene-wise manner and return
#' the corresponding response matrix. This function is build on other functions (fit_models and responseMatrix) and used in calILRs and calABCs functions
#'
#' @param cds a CellDataSet object upon which to perform this operation
#' @param new_data a data.frame object including columns (for example, Pseudotime) with names specified in the model formula. The values in the data.frame should be consist with the corresponding values from cds object.
#' @param trend_formula a formula string specifying the model formula used in fitting the spline curve for each gene/feature.
#' @param relative_expr a logic flag to determine whether or not the relative gene expression should be used
#' @param response_type the response desired, as accepted by VGAM's predict function
#' @param cores the number of cores to be used while testing each gene for differential expression
#' @return a data frame containing the data for the fitted spline curves.
#' @export
#'
genSmoothCurves <- function(cds,  new_data, trend_formula = "~sm.ns(Pseudotime, df = 3)",
                        relative_expr = T, response_type="response", cores = 1) { 
  
  expressionFamily <- cds@expressionFamily
  
    if(cores > 1) {
      expression_curve_matrix <- mcesApply(cds, 1, function(x, trend_formula, expressionFamily, relative_expr, new_data, fit_model_helper, responseMatrix, 
                                                              calulate_NB_dispersion_hint, calulate_QP_dispersion_hint){
            environment(fit_model_helper) <- environment()
            environment(responseMatrix) <- environment()
            model_fits <- fit_model_helper(x, modelFormulaStr = trend_formula, expressionFamily = expressionFamily, 
                                       relative_expr = relative_expr, disp_func = cds@dispFitInfo[['blind']]$disp_func)
            if(is.null(model_fits))
                expression_curve <- as.data.frame(matrix(rep(NA, nrow(new_data)), nrow = 1))
            else
                expression_curve <- as.data.frame(responseMatrix(list(model_fits), newdata = new_data, response_type=response_type))
            colnames(expression_curve) <- row.names(new_data)
            expression_curve
            #return(expression_curve)
            }, required_packages=c("BiocGenerics", "Biobase", "VGAM", "plyr"), cores=cores, 
            trend_formula = trend_formula, expressionFamily = expressionFamily, relative_expr = relative_expr, new_data = new_data, 
            fit_model_helper = fit_model_helper, responseMatrix = responseMatrix, calulate_NB_dispersion_hint = calulate_NB_dispersion_hint,
            calulate_QP_dispersion_hint = calulate_QP_dispersion_hint
            )
        expression_curve_matrix <- as.matrix(do.call(rbind, expression_curve_matrix))
        return(expression_curve_matrix)
    }
    else {
        expression_curve_matrix <- smartEsApply(cds, 1, function(x, trend_formula, expressionFamily, relative_expr, new_data){
            environment(fit_model_helper) <- environment()
            environment(responseMatrix) <- environment()
            model_fits <- fit_model_helper(x, modelFormulaStr = trend_formula, expressionFamily = expressionFamily,
                                       relative_expr = relative_expr, disp_func = cds@dispFitInfo[['blind']]$disp_func)
            if(is.null(model_fits))
                expression_curve <- as.data.frame(matrix(rep(NA, nrow(new_data)), nrow = 1))
            else
                expression_curve <- as.data.frame(responseMatrix(list(model_fits), new_data, response_type=response_type))
            colnames(expression_curve) <- row.names(new_data)
            expression_curve
            }, 
            convert_to_dense=TRUE,
            trend_formula = trend_formula, expressionFamily = expressionFamily, relative_expr = relative_expr, new_data = new_data
            )
        expression_curve_matrix <- as.matrix(do.call(rbind, expression_curve_matrix))
        row.names(expression_curve_matrix) <- row.names(fData(cds))
        return(expression_curve_matrix)
      }

}
#' Fit smooth spline curves and return the residuals matrix
#'
#' This function will fit smooth spline curves for the gene expression dynamics along pseudotime in a gene-wise manner and return
#' the corresponding residuals matrix. This function is build on other functions (fit_models and residualsMatrix)
#'
#' @param cds a CellDataSet object upon which to perform this operation
#' @param trend_formula a formula string specifying the model formula used in fitting the spline curve for each gene/feature.
#' @param relative_expr a logic flag to determine whether or not the relative gene expression should be used
#' @param residual_type the response desired, as accepted by VGAM's predict function
#' @param cores the number of cores to be used while testing each gene for differential expression
#' @return a data frame containing the data for the fitted spline curves.
#'
genSmoothCurveResiduals <- function(cds, trend_formula = "~sm.ns(Pseudotime, df = 3)",
                            relative_expr = T,  residual_type="response", cores = 1) { 
  
  expressionFamily <- cds@expressionFamily
  
  if(cores > 1) {
    expression_curve_matrix <- mcesApply(cds, 1, function(x, trend_formula, expressionFamily, relative_expr, fit_model_helper, residualMatrix, 
                                                    calulate_NB_dispersion_hint, calulate_QP_dispersion_hint){
      environment(fit_model_helper) <- environment()
      environment(responseMatrix) <- environment()
      model_fits <- fit_model_helper(x, modelFormulaStr = trend_formula, expressionFamily = expressionFamily,
                                     relative_expr = relative_expr, disp_func = cds@dispFitInfo[['blind']]$disp_func)
      if(is.null(model_fits))
        expression_curve <- as.data.frame(matrix(rep(NA, nrow(pData(cds))), nrow = 1))
      else
        expression_curve <- as.data.frame(residualMatrix(list(model_fits), residual_type=residual_type))
      #colnames(expression_curve) <- row.names(pData(cds))
      expression_curve
      #return(expression_curve)
    }, required_packages=c("BiocGenerics", "Biobase", "VGAM", "plyr"), cores=cores, 
    trend_formula = trend_formula, expressionFamily = expressionFamily, relative_expr = relative_expr,
    fit_model_helper = fit_model_helper, residualMatrix = residualMatrix, calulate_NB_dispersion_hint = calulate_NB_dispersion_hint,
    calulate_QP_dispersion_hint = calulate_QP_dispersion_hint
    )
    expression_curve_matrix <- do.call(rbind, expression_curve_matrix)
    return(expression_curve_matrix)
  }
  else {
    expression_curve_matrix <- smartEsApply(cds, 1, function(x, trend_formula, expressionFamily, relative_expr){
      environment(fit_model_helper) <- environment()
      environment(residualMatrix) <- environment()
      model_fits <- fit_model_helper(x, modelFormulaStr = trend_formula, expressionFamily = expressionFamily, 
                                     relative_expr = relative_expr, disp_func = cds@dispFitInfo[['blind']]$disp_func)
      if(is.null(model_fits))
        expression_curve <-  as.data.frame(matrix(rep(NA, nrow(pData(cds))), nrow = 1))
      else
        expression_curve <-  as.data.frame(residualMatrix(list(model_fits), residual_type=residual_type))
      #colnames(expression_curve) <- row.names(pData(cds))
      expression_curve
    }, 
    convert_to_dense=TRUE,
    trend_formula = trend_formula, expressionFamily = expressionFamily, relative_expr = relative_expr
    )
    expression_curve_matrix <- do.call(rbind, expression_curve_matrix)
    row.names(expression_curve_matrix) <- row.names(fData(cds))
    return(expression_curve_matrix)
  }
  
}


## This function was swiped from DESeq (Anders and Huber) and modified for our purposes
parametricDispersionFit <- function( disp_table, initial_coefs=c(1e-6, 1) )
{
  coefs <- initial_coefs
  iter <- 0
  while(TRUE) {
    residuals <- disp_table$disp / ( coefs[1] + coefs[2] / disp_table$mu )
    good <- disp_table[which( (residuals > initial_coefs[1]) & (residuals < 10000) ),]
    #good <- disp_table
    fit <- glm( disp ~ I(1/mu), data=good,
                family=Gamma(link="identity"), start=coefs )
    oldcoefs <- coefs
    coefs <- coefficients(fit)
    if (coefs[1] < initial_coefs[1]){
      coefs[1] <- initial_coefs[1]
    }
    if (coefs[2] < 0){
      stop( "Parametric dispersion fit failed. Try a local fit and/or a pooled estimation. (See '?estimateDispersions')" )
    }
    #     if( !all( coefs > 0 ) ){
    #       #print(data.frame(means,disps))
    #       stop( "Parametric dispersion fit failed. Try a local fit and/or a pooled estimation. (See '?estimateDispersions')" )
    #     }
    if( sum( log( coefs / oldcoefs )^2 ) < initial_coefs[1] )
      break
    iter <- iter + 1
    #print(coefs)
    if( iter > 10 ) {
      warning( "Dispersion fit did not converge." )
      break 
    }
  }
  
  if( !all( coefs > 0 ) ){
    #print(data.frame(means,disps))
    stop( "Parametric dispersion fit failed. Try a local fit and/or a pooled estimation. (See '?estimateDispersions')" )
  }
  
  #names( coefs ) <- c( "asymptDisp", "extraPois" )
  #ans <- function( q )
  #  coefs[1] + coefs[2] / q
  #ans
  #coefs
  list(fit, coefs)
}

# parametricDispersionFit <- function( means, disps )
# {
#   coefs <- c( 1e-6, 1 )
#   iter <- 0
#  
#     residuals <- disps / ( coefs[1] + coefs[2] / means )
#     good <- which( (residuals > 1e-4) & (residuals < 10000) )
#     fit <- vglm( log(disps[good]) ~ log(means[good]), family=gaussianff())
#     oldcoefs <- coefs
#     coefs <- coefficients(fit)
# 
#     iter <- iter + 1
#     print(coefs)
#   names( coefs ) <- c( "asymptDisp", "extraPois" )
#   ans <- function( q )
#     exp(coefs[1] + coefs[2] * log(q))
#   #ans
#   coefs
# }


#' Return a variance-stabilized matrix of expression values
#'
#' This function was taken from the DESeq package (Anders and Huber) and modified 
#' to suit Monocle's needs
#'
#' @param cds A CellDataSet to use for variance stabilization.
#' @param dispModelName The name of the dispersion function to use for VST.
#' @param expr_matrix An matrix of values to transform. Must be normalized (e.g. by size factors) already. This function doesn't do this for you.
#' @param round_vals Whether to round expression values to the nearest integer before applying the transformation. 
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
    ncounts <- Matrix::t(Matrix::t(ncounts) / sizeFactors(cds))
    if (round_vals)
      ncounts <- round(ncounts)
  }else{
    ncounts <- expr_matrix
  }
  vst <- function( q )
    log( (1 + coefs["extraPois"] + 2 * coefs["asymptDisp"] * q +
            2 * sqrt( coefs["asymptDisp"] * q * ( 1 + coefs["extraPois"] + coefs["asymptDisp"] * q ) ) )
         / ( 4 * coefs["asymptDisp"] ) ) / log(2)
  ## NOTE: this converts to a dense matrix
  vst( ncounts )
}

disp_calc_helper_NB <- function(cds, expressionFamily, min_cells_detected){
  
  rounded <- round(exprs(cds))
  nzGenes <- Matrix::rowSums(rounded > cds@lowerDetectionLimit)
  nzGenes <- names(nzGenes[nzGenes > min_cells_detected])
  
  x <- t(t(rounded[nzGenes,]) / pData(cds[nzGenes,])$Size_Factor)
  
  xim <- mean(1/ pData(cds[nzGenes,])$Size_Factor)

  if (isSparseMatrix(exprs(cds))){
    f_expression_mean <- as(Matrix::rowMeans(x), "sparseVector")
  }else{
    f_expression_mean <- Matrix::rowMeans(x)
  }
  
  
    # For NB: Var(Y)=mu*(1+mu/k)
  f_expression_var <- Matrix::rowMeans((x - f_expression_mean)^2)

  disp_guess_meth_moments <- f_expression_var - xim * f_expression_mean
    
  disp_guess_meth_moments <- disp_guess_meth_moments / (f_expression_mean^2) #fix the calculation of k
  

  res <- data.frame(mu=as.vector(f_expression_mean), disp=as.vector(disp_guess_meth_moments))
  res[res$mu == 0]$mu = NA
  res[res$mu == 0]$disp = NA  
  res$disp[res$disp < 0] <- 0

  res <- cbind(gene_id=row.names(fData(cds[nzGenes,])), res)
  res
}

disp_calc_helper_BB <- function(cds, expressionFamily, min_cells_detected){
  
  nzGenes <- Matrix::rowSums(exprs(cds) > 0)
  nzGenes <- names(nzGenes[nzGenes > min_cells_detected])
  
  x <- t(t(exprs(cds[nzGenes,])) / pData(cds[nzGenes,])$Size_Factor)
  m1 <- Matrix::rowMeans(x)
  m2 <- Matrix::rowMeans(x^2)
  n <- qlcMatrix::rowMax(x)
  n[n < 1] <- 1
  #n <- max(1, n)
  
  denom <- n * ((m2/m1) - m1 - 1) + m1
  alpha_hat <- (n*m1 - m2) / denom
  beta_hat <- (n - m1)*(n - (m2/m1)) / denom
  
  mu = alpha_hat/(alpha_hat + beta_hat)
  rho = 1/(1 + alpha_hat + beta_hat)
  
  res <- data.frame(mu=mu, disp=rho)
  res[res$mu == 0]$mu = NA
  res[res$mu == 0]$disp = NA
  res <- cbind(gene_id=row.names(fData(cds[nzGenes,])), res)
  res <- subset(res, disp < 1e5)

  res
}


#' Helper function to estimate dispersions
#' @importFrom stringr str_split str_trim
estimateDispersionsForCellDataSet <- function(cds, modelFormulaStr, relative_expr, min_cells_detected, removeOutliers, cores)
{
  
  # if (cores > 1){
  #     disp_table<-mcesApply(cds, 1, disp_calc_helper, c("BiocGenerics", "Biobase", "VGAM", "dplyr", "Matrix"), cores=cores, 
  #                         modelFormulaStr=modelFormulaStr, 
  #                         expressionFamily=cds@expressionFamily)
  # }else{
  #     disp_table<-smartEsApply(cds,1,disp_calc_helper, 
  #                              convert_to_dense=TRUE,
  #                              modelFormulaStr=modelFormulaStr, 
  #                              expressionFamily=cds@expressionFamily)
  # }
  
  model_terms <- unlist(lapply(str_split(modelFormulaStr, "~|\\+|\\*"), str_trim))
  model_terms <- model_terms[model_terms != ""]
  progress_opts <- options()$dplyr.show_progress
  options(dplyr.show_progress = T)
  
  # FIXME: this needs refactoring, badly.
  if (cds@expressionFamily@vfamily %in% c("negbinomial", "negbinomial.size")){
    if (length(model_terms) > 1 || (length(model_terms) == 1 && model_terms[1] != "1")){
      cds_pdata <- dplyr::group_by_(dplyr::select_(add_rownames(pData(cds)), "rowname", .dots=model_terms), .dots=model_terms) 
      disp_table <- as.data.frame(cds_pdata %>% do(disp_calc_helper_NB(cds[,.$rowname], cds@expressionFamily, min_cells_detected)))
    }else{
      cds_pdata <- dplyr::group_by_(dplyr::select_(add_rownames(pData(cds)), "rowname")) 
      disp_table <- as.data.frame(cds_pdata %>% do(disp_calc_helper_NB(cds[,.$rowname], cds@expressionFamily, min_cells_detected)))
      #disp_table <- data.frame(rowname = names(type_res), CellType = type_res)
    }
    
    #message("fitting disersion curves")
    #print (disp_table)
    if(!is.list(disp_table))
      stop("Parametric dispersion fitting failed, please set a different lowerDetectionLimit")
    #disp_table <- do.call(rbind.data.frame, disp_table)
    disp_table <- subset(disp_table, is.na(mu) == FALSE)
    res <- parametricDispersionFit(disp_table)
    fit <- res[[1]]
    coefs <- res[[2]]
    #removeOutliers = TRUE
    if (removeOutliers){
      CD <- cooks.distance(fit)
      #cooksCutoff <- qf(.99, 2, ncol(cds) - 2)
      cooksCutoff <- 4/nrow(disp_table)
      #print (head(CD[CD > cooksCutoff]))
      #print (head(names(CD[CD > cooksCutoff])))
      message (paste("Removing", length(CD[CD > cooksCutoff]), "outliers"))
      outliers <- union(names(CD[CD > cooksCutoff]), setdiff(row.names(disp_table), names(CD))) 
      res <- parametricDispersionFit(disp_table[row.names(disp_table) %in% outliers == FALSE,])
      fit <- res[[1]]
      coefs <- res[[2]]
    }
    names( coefs ) <- c( "asymptDisp", "extraPois" )
    ans <- function( q )
      coefs[1] + coefs[2] / q
    attr( ans, "coefficients" ) <- coefs
    
  }else if (cds@expressionFamily@vfamily %in% c("betabinomial")){
    if (length(model_terms) > 1 || (length(model_terms) == 1 && model_terms[1] != "1")){
      cds_pdata <- dplyr::group_by_(dplyr::select_(add_rownames(pData(cds)), "rowname", .dots=model_terms), .dots=model_terms) 
      disp_table <- as.data.frame(cds_pdata %>% do(disp_calc_helper_BB(cds[,.$rowname], cds@expressionFamily, min_cells_detected)))
    }else{
      cds_pdata <- dplyr::group_by_(dplyr::select_(add_rownames(pData(cds)), "rowname")) 
      disp_table <- as.data.frame(cds_pdata %>% do(disp_calc_helper_BB(cds[,.$rowname], cds@expressionFamily, min_cells_detected)))
      #disp_table <- data.frame(rowname = names(type_res), CellType = type_res)
    }
    
    #message("fitting disersion curves")
    #print (disp_table)
    if(!is.list(disp_table))
      stop("Parametric dispersion fitting failed, please set a different lowerDetectionLimit")
    #disp_table <- do.call(rbind.data.frame, disp_table)
    disp_table <- subset(disp_table, is.na(mu) == FALSE)
    res <- parametricDispersionFit(disp_table, c(1e-10, 1))
    fit <- res[[1]]
    coefs <- res[[2]]
    #removeOutliers = TRUE
    if (removeOutliers){
      CD <- cooks.distance(fit)
      #cooksCutoff <- qf(.99, 2, ncol(cds) - 2)
      cooksCutoff <- 4/nrow(disp_table)
      #print (head(CD[CD > cooksCutoff]))
      #print (head(names(CD[CD > cooksCutoff])))
      message (paste("Removing", length(CD[CD > cooksCutoff]), "outliers"))
      outliers <- union(names(CD[CD > cooksCutoff]), setdiff(row.names(disp_table), names(CD))) 
      res <- parametricDispersionFit(disp_table[row.names(disp_table) %in% outliers == FALSE,], c(1e-10, 1))
      fit <- res[[1]]
      coefs <- res[[2]]
    }
    
    names( coefs ) <- c( "asymptDisp", "extraPois" )
    ans <- function( q )
      coefs[1] + coefs[2] / q
    attr( ans, "coefficients" ) <- coefs
  }
  
  res <- list(disp_table = disp_table, disp_func = ans)
  return(res)
}

calulate_NB_dispersion_hint <- function(disp_func, f_expression, expr_selection_func=mean)
{
  expr_hint <- expr_selection_func(f_expression)
  if (expr_hint > 0 && is.null(expr_hint) == FALSE) {
    disp_guess_fit <- disp_func(expr_hint)
    
    # For NB: Var(Y)=mu*(1+mu/k)
    f_expression_var <- var(f_expression)
    f_expression_mean <- mean(f_expression)
    
    disp_guess_meth_moments <- f_expression_var - f_expression_mean 
    disp_guess_meth_moments <- disp_guess_meth_moments / (f_expression_mean^2) #fix the calculation of k 
    
    #return (max(disp_guess_fit, disp_guess_meth_moments))
    return (disp_guess_fit)
  }
  return (NULL)
}

# note that quasipoisson expects a slightly different format for the 
# dispersion parameter, hence the differences in return value between
# this function and calulate_NB_dispersion_hint
calulate_QP_dispersion_hint <- function(disp_func, f_expression, expr_selection_func=mean)
{
  expr_hint <- expr_selection_func(f_expression)
  if (expr_hint > 0 && is.null(expr_hint) == FALSE) {
    disp_guess_fit <- disp_func(expr_hint)
    
    # For NB: Var(Y)=mu*(1+mu/k)
    f_expression_var <- var(f_expression)
    f_expression_mean <- mean(f_expression)
    
    disp_guess_meth_moments <- f_expression_var - f_expression_mean 
    disp_guess_meth_moments <- disp_guess_meth_moments / (f_expression_mean^2) #fix the calculation of k 
    
    return (1 + f_expression_mean * max(disp_guess_fit, disp_guess_meth_moments))
  }
  return (NULL)
}

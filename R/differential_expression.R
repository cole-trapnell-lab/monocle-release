#' @title Helper function for parallel differential expression testing
#' @param x test
#' @param fullModelFormulaStr a formula string specifying the full model in differential expression tests (i.e. likelihood ratio tests) for each gene/feature.
#' @param reducedModelFormulaStr a formula string specifying the reduced model in differential expression tests (i.e. likelihood ratio tests) for each gene/feature.
#' @param expressionFamily specifies the VGAM family function used for expression responses
#' @param relative_expr Whether to transform expression into relative values
#' @param weights test
#' @param disp_func test
#' @param verbose Whether to show VGAM errors and warnings. Only valid for cores = 1.
#' @name diff_test_helper
#' @description test
diff_test_helper <- function(x, 
                             fullModelFormulaStr, 
                             reducedModelFormulaStr, 
                             expressionFamily, 
                             relative_expr,
                             weights,
                             disp_func=NULL,
                             verbose=FALSE
                             ){ 
  
  reducedModelFormulaStr <- paste("f_expression", reducedModelFormulaStr, sep="")
  fullModelFormulaStr <- paste("f_expression", fullModelFormulaStr, sep="")
  
  x_orig <- x
  disp_guess <- 0
  
  if (expressionFamily@vfamily %in% c("negbinomial", "negbinomial.size")){
    if (relative_expr == TRUE)
    {
      x <- x / Size_Factor
    }
    f_expression <- round(x)
    if (is.null(disp_func) == FALSE){
      disp_guess <- calculate_NB_dispersion_hint(disp_func, round(x_orig))
      if (is.null(disp_guess) == FALSE && disp_guess > 0 && is.na(disp_guess) == FALSE  ) {
        # FIXME: In theory, we could lose some user-provided parameters here
        # e.g. if users supply zero=NULL or something. 
        if (expressionFamily@vfamily == "negbinomial")
          expressionFamily <- negbinomial(isize=1/disp_guess)
        else
          expressionFamily <- negbinomial.size(size=1/disp_guess)
      }
    }
  }else if (expressionFamily@vfamily %in% c("gaussianff", "uninormal")){
    f_expression <- x
  }else if (expressionFamily@vfamily %in% c("binomialff")){
    f_expression <- x
    #f_expression[f_expression > 1] <- 1
  }else{
    f_expression <- log10(x)
  }
  
  test_res <- tryCatch({
    if (expressionFamily@vfamily %in% c("binomialff")){
      if (verbose){
        full_model_fit <- VGAM::vglm(as.formula(fullModelFormulaStr), epsilon=1e-1, family=expressionFamily)
        reduced_model_fit <- VGAM::vglm(as.formula(reducedModelFormulaStr), epsilon=1e-1, family=expressionFamily)                         
      }else{
        full_model_fit <- suppressWarnings(VGAM::vglm(as.formula(fullModelFormulaStr), epsilon=1e-1, family=expressionFamily))
        reduced_model_fit <- suppressWarnings(VGAM::vglm(as.formula(reducedModelFormulaStr), epsilon=1e-1, family=expressionFamily))                    
      }
    }else{
      if (verbose){
        full_model_fit <- VGAM::vglm(as.formula(fullModelFormulaStr), epsilon=1e-1, family=expressionFamily)
        reduced_model_fit <- VGAM::vglm(as.formula(reducedModelFormulaStr), epsilon=1e-1, family=expressionFamily)                         
      }else{
        full_model_fit <- suppressWarnings(VGAM::vglm(as.formula(fullModelFormulaStr), epsilon=1e-1, family=expressionFamily))
        reduced_model_fit <- suppressWarnings(VGAM::vglm(as.formula(reducedModelFormulaStr), epsilon=1e-1, family=expressionFamily))                    
      }
    }

    #print(full_model_fit)
    #print(coef(reduced_model_fit))
    compareModels(list(full_model_fit), list(reduced_model_fit))
  }, 
  #warning = function(w) { FM_fit },
  error = function(e) { 
    if(verbose)
      print (e);
      data.frame(status = "FAIL", family=expressionFamily@vfamily, pval=1.0, qval=1.0)
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
#' @importFrom stats p.adjust
#' @export
compareModels <- function(full_models, reduced_models){
  stopifnot(length(full_models) == length(reduced_models))
  test_res <- mapply(function(x,y) { 
    if (is.null(x) == FALSE && is.null(y) == FALSE) {
      lrt <- VGAM::lrtest(x,y) 
      pval=lrt@Body["Pr(>Chisq)"][2,]
      family = x@family@vfamily
      if (length(family) > 1)
        family = family[1]
      data.frame(status = "OK", family=family, pval=pval)
    } else { data.frame(status = "FAIL", family=NA, pval=1.0) } 
  } , full_models, reduced_models, SIMPLIFY=FALSE, USE.NAMES=TRUE)
  
  test_res <- do.call(rbind.data.frame, test_res)
  test_res$qval <- p.adjust(test_res$pval, method="BH")
  test_res
}

#' Test genes for differential expression
#' 
#' Tests each gene for differential expression as a function of pseudotime 
#' or according to other covariates as specified. \code{differentialGeneTest} is
#' Monocle's main differential analysis routine. 
#' It accepts a CellDataSet and two model formulae as input, which specify generalized
#' lineage models as implemented by the \code{VGAM} package. 
#' 
#' @param cds a CellDataSet object upon which to perform this operation
#' @param fullModelFormulaStr a formula string specifying the full model in differential expression tests (i.e. likelihood ratio tests) for each gene/feature.
#' @param reducedModelFormulaStr a formula string specifying the reduced model in differential expression tests (i.e. likelihood ratio tests) for each gene/feature.
#' @param relative_expr Whether to transform expression into relative values.
#' @param cores the number of cores to be used while testing each gene for differential expression.
#' @param verbose Whether to show VGAM errors and warnings. Only valid for cores = 1. 
#' @return a data frame containing the p values and q-values from the likelihood ratio tests on the parallel arrays of models.
#' @importFrom Biobase fData
#' @importFrom stats p.adjust
#' @seealso \code{\link[VGAM]{vglm}}
#' @export
differentialGeneTest <- function(cds, 
                                 fullModelFormulaStr="~sm.ns(Pseudotime, df=3)",
                                 reducedModelFormulaStr="~1", 
                                 relative_expr=TRUE,
                                 cores=1, 
                                 verbose=FALSE
                                 ){
  status <- NA
  
  if(class(cds)[1] != "CellDataSet") {
    stop("Error cds is not of type 'CellDataSet'")
  }
  
  all_vars <- c(all.vars(formula(fullModelFormulaStr)), all.vars(formula(reducedModelFormulaStr)))
   
  pd <- pData(cds)
  
  for(i in all_vars) {
    x <- pd[, i]
    if(any((c(Inf, NaN, NA) %in% x))){
      stop("Error: Inf, NaN, or NA values were located in pData of cds in columns mentioned in model terms")
    }
  }
  
  
  if (relative_expr && cds@expressionFamily@vfamily %in% c("negbinomial", "negbinomial.size")){
    if (is.null(sizeFactors(cds)) || sum(is.na(sizeFactors(cds)))){
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
                             verbose=verbose
                       #       ,
                       # backup_method = backup_method, 
                       # use_epislon = use_epislon, 
                       # stepsize = stepsize
                             )
    diff_test_res
  }else{
    diff_test_res<-smartEsApply(cds,1,diff_test_helper, 
                                convert_to_dense=TRUE,
                                fullModelFormulaStr=fullModelFormulaStr,
                                reducedModelFormulaStr=reducedModelFormulaStr, 
                                expressionFamily=cds@expressionFamily, 
                                relative_expr=relative_expr,
                                disp_func=cds@dispFitInfo[["blind"]]$disp_func,
                                verbose=verbose
                       #          ,
                       # backup_method = backup_method, 
                       # use_epislon = use_epislon,
                       # stepsize = stepsize

                                )
    diff_test_res
  }
  
  diff_test_res <- do.call(rbind.data.frame, diff_test_res)

  diff_test_res$qval <- 1
  diff_test_res$qval[which(diff_test_res$status == 'OK')] <- p.adjust(subset(diff_test_res, status == 'OK')[, 'pval'], method="BH")
  
  diff_test_res <- merge(diff_test_res, fData(cds), by="row.names")
  row.names(diff_test_res) <- diff_test_res[, 1] #remove the first column and set the row names to the first column
  diff_test_res[, 1] <- NULL 

  diff_test_res[row.names(cds), ] # make sure gene name ordering in the DEG test result is the same as the CDS
}

#' Test genes for differential expression
#' 
#' Tests each gene for differential expression as a function of pseudotime 
#' or according to other covariates as specified. \code{differentialGeneTest} is
#' Monocle's main differential analysis routine. 
#' It accepts a CellDataSet and two model formulae as input, which specify generalized
#' lineage models as implemented by the \code{VGAM} package. 
#' 
#' @param cds a CellDataSet object upon which to perform this operation
#' @param relative_expr Whether to transform expression into relative values.
#' @param k Number of nearest neighbors used for building the kNN graph which is passed to knn2nb function during the Moran's I test procedure
#' @param cores the number of cores to be used while testing each gene for differential expression.
#' @param verbose Whether to show VGAM errors and warnings. Only valid for cores = 1. 
#' @return a data frame containing the p values and q-values from the Moran's I test on the parallel arrays of models.
#' @importFrom spdep knn2nb nb2listw moran.test 
#' @importFrom stats p.adjust 
#' @seealso \code{\link[spdep]{moran.test}}
#' @export
spatialDifferentialTest <- function(cds, 
                                 relative_expr=TRUE,
                                 k = 5, 
                                 cores=1, 
                                 verbose=FALSE) {

  # 1. first retrieve the association from each cell to any principal points. Then for any two connected principal points, 
  # find all cells associated with the two principal points,  build kNN graph, then pool all kNN and then remove redundant 
  # points and finally use this kNN graph to calculate a global Moran’s I and get the p-value

  # cds@auxOrderingData[["L1graph"]]$adj_mat # graph from UMAP 

  data <- t(reducedDimA(cds)) # cell coordinates on low dimensional 
  cell2pp_map <- cds@auxOrderingData[["L1graph"]]$pr_graph_cell_proj_closest_vertex # mapping from each cell to the principal points 
  principal_g <- cds@auxOrderingData[["L1graph"]]$W 
  
  # Use all pairings of i and j
  # principal_g[upper.tri(principal_g)] <- 0
  i_vec <- rep(seq_len(ncol(principal_g)), times = apply(principal_g, 1, function(x) sum(x > 0)))
  j_vec <- as.vector(unlist(apply(principal_g, 1, function(x) which(x > 0))))
  
  # calculate kNN in parallel for cells belong to each pair of connected principal points
  conn_kNN_graph <- mcmapply(i_vec[which(j_vec != 0)], j_vec[which(j_vec != 0)], 
                       FUN = function(i, j) {
                         # message('i, j are ', i, ' ', j)
                         cells_id_to_pp <- which(cell2pp_map %in% c(i, j))
                         data <- data[cells_id_to_pp, ]
                         res <- RANN::nn2(data, data, min(k + 1, length(cells_id_to_pp)), searchtype = "standard")[[1]][, ]
                         res <- matrix(cells_id_to_pp[res], ncol = k + 1, byrow = F)
                      },
                      mc.cores = cores
              ) 
  conn_kNN_graph <- do.call(rbind, conn_kNN_graph)
  conn_kNN_graph <- unique(conn_kNN_graph)[, ]
  
  res <- list(nn = conn_kNN_graph[, -1], np = nrow(conn_kNN_graph), 
              k = k, dimension = ncol(data), x = data[conn_kNN_graph[, 1], ])

  class(res) <- "knn"
  attr(res, "call") <- match.call()

  k1 <- knn2nb(res)
  lw <- nb2listw(k1, zero.policy=TRUE)
  
  moran_test_res <- mclapply(row.names(cds), FUN = function(x) {
    exprs_val <- exprs(cds)[x, ]
    
    if (cds@expressionFamily@vfamily %in% c("gaussianff", "uninormal", "binomialff")){
      exprs_val <- x 
    }else{
      if(relative_expr) {
        exprs_val <- log10(exprs_val / sizeFactors(cds) + 1)
      } else {
        exprs_val <- log10(exprs_val + 1)
      }
    }
    
    test_res <- tryCatch({
      mt <- moran.test(exprs_val[conn_kNN_graph[, 1]], lw, zero.policy=TRUE)
      data.frame(status = 'OK', pval = mt$p.value, statistics = mt$statistic)
    }, 
      error = function(e) {
        data.frame(status = 'FAIL', pval = NA, statistics = NA)
      })
  }, mc.cores = cores)
  
  moran_test_res <- do.call(rbind.data.frame, moran_test_res)
  row.names(moran_test_res) <- row.names(cds)
  
  moran_test_res <- merge(moran_test_res, fData(cds), by="row.names")
  row.names(moran_test_res) <- moran_test_res[, 1] #remove the first column and set the row names to the first column
  moran_test_res[, 1] <- NULL 
  moran_test_res$qval <- 1
  moran_test_res$qval[which(moran_test_res$status == 'OK')] <- p.adjust(subset(moran_test_res, status == 'OK')[, 'pval'], method="BH")
  
  moran_test_res[row.names(cds), ] # make sure gene name ordering in the DEG test result is the same as the CDS
}


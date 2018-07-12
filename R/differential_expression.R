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

#' Function to calculate the a neighbours list with spatial weights for the chosen coding scheme from a cell dataset object
calculateLW <- function(cds, verbose = FALSE, k = 25, return_sparse_matrix = FALSE) {
  # first retrieve the association from each cell to any principal points, then build kNN graph for all cells
  # remove edges that connected between groups that disconnected in the corresponding principal graph and
  # finally use this kNN graph to calculate a global Moran’s I and get the p-value
  
  if(verbose) {
    message("retrieve the matrices for Moran's test...")
  }
  
  if(cds@dim_reduce_type == 'L1graph') {
    cell_coords <- t(reducedDimA(cds)) # cell coordinates on low dimensional
    principal_g <- cds@auxOrderingData[["L1graph"]]$W
  } else if(cds@dim_reduce_type %in% c('DDRTree', 'SimplePPT')) {
    cell_coords <- t(reducedDimS(cds))
    principal_g <-  igraph::get.adjacency(cds@minSpanningTree)[1:ncol(reducedDimK(cds)), 1:ncol(reducedDimK(cds))]
  } else if(cds@dim_reduce_type %in% c('UMAP')) {
    cell_coords <- t(reducedDimS(cds))
    knn_res <- RANN::nn2(cell_coords, cell_coords, min(k + 1, nrow(cell_coords)), searchtype = "standard")[[1]]
  }
  
  exprs_mat <- exprs(cds)
  cell2pp_map <- cds@auxOrderingData[[cds@dim_reduce_type]]$pr_graph_cell_proj_closest_vertex # mapping from each cell to the principal points
  
  if(is.null(cell2pp_map)) {
    links <- monocle:::jaccard_coeff(knn_res[, -1], F)
    links <- links[links[, 1] > 0, ]
    relations <- as.data.frame(links)
    colnames(relations) <- c("from", "to", "weight")
    knn_res_graph <- igraph::graph.data.frame(relations, directed = T)
    
    if(return_sparse_matrix) {
      tmp <- get.adjacency(knn_res_graph)
      dimnames(tmp) <- list(colnames(cds), colnames(cds))
      
      return(tmp)
    }
    

    knn_list <- lapply(1:nrow(knn_res), function(x) knn_res[x, -1])

  } else {
    # This cds object might be a subset of the one on which ordering was performed,
    # so we may need to subset the nearest vertex and low-dim coordinate matrices:
    cell2pp_map <-  cell2pp_map[row.names(cell2pp_map) %in% row.names(pData(cds)),, drop=FALSE]
    cell2pp_map <- cell2pp_map[colnames(cds), ]
    
    # cds@auxOrderingData[["L1graph"]]$adj_mat # graph from UMAP
    
    if(verbose) {
      message("Identify connecting principal point pairs ...")
    }
    
    # an alternative approach to make the kNN graph based on the principal graph
    knn_res <- RANN::nn2(cell_coords, cell_coords, min(k + 1, nrow(cell_coords)), searchtype = "standard")[[1]]
    # kNN_res_pp_map <- matrix(cell2pp_map[knn_res], ncol = k + 1, byrow = F) # convert the matrix of knn graph from the cell IDs into a matrix of principal points IDs
    
    principal_g_tmp <- principal_g # kNN can be built within group of cells corresponding to each principal points
    diag(principal_g_tmp) <- 1 # so set diagnol as 1
    
    cell_membership <- as.factor(cell2pp_map)
    uniq_member <- sort(unique(cell_membership))
    
    membership_matrix <- sparse.model.matrix( ~ cell_membership + 0)
    colnames(membership_matrix) <- levels(uniq_member)
    
    # sparse matrix multiplication for calculating the feasible space
    feasible_space <- membership_matrix %*% tcrossprod(principal_g_tmp[as.numeric(levels(uniq_member)), as.numeric(levels(uniq_member))], membership_matrix)
    
    links <- monocle:::jaccard_coeff(knn_res[, -1], F)
    links <- links[links[, 1] > 0, ]
    relations <- as.data.frame(links)
    colnames(relations) <- c("from", "to", "weight")
    knn_res_graph <- igraph::graph.data.frame(relations, directed = T)
    
    # remove edges across cells belong to two disconnected principal points
    tmp <- get.adjacency(knn_res_graph) * feasible_space
    
    if(return_sparse_matrix) {
      dimnames(tmp) <- list(colnames(cds), colnames(cds))
      return(tmp)
    }
    
    knn_list <- slam::rowapply_simple_triplet_matrix(slam::as.simple_triplet_matrix(tmp), function(x) {
      res <- which(as.numeric(x) > 0)
      if(length(res) == 0)
        res <- 0L
      res
    })
  }
  # create the lw list for moran.test
  class(knn_list) <- "nb"
  attr(knn_list, "region.id") <- colnames(cds)
  attr(knn_list, "call") <- match.call()
  # attr(knn_list, "type") <- "queen"
  lw <- nb2listw(knn_list, zero.policy = TRUE)
  
  lw
}

#' Test genes for differential expression based on the low dimensional embedding and the principal graph
#'
#' Tests each gene for differential expression as a function of pseudotime
#' or according to other covariates as specified. \code{differentialGeneTest} is
#' Monocle's main differential analysis routine.
#' It accepts a CellDataSet and two model formulae as input, which specify generalized
#' lineage models as implemented by the \code{VGAM} package.
#'
#' @param cds a CellDataSet object upon which to perform this operation
#' @param landmark_num Number of landmark cells selected for performing aggregate Moran's I test, default is NULL (no landmark selection and all cells are used)
#' @param relative_expr Whether to transform expression into relative values.
#' @param k Number of nearest neighbors used for building the kNN graph which is passed to knn2nb function during the Moran's I (Geary's C) test procedure.
#' @param method a character string specifying the method (the default 'Moran_I' or 'Geary_C') for detecting significant genes showing correlated genes along the principal graph embedded in the low dimensional space.
#' @param alternative a character string specifying the alternative hypothesis, must be one of greater (default), less or two.sided.
#' @param cores the number of cores to be used while testing each gene for differential expression.
#' @param verbose Whether to show VGAM errors and warnings. Only valid for cores = 1.
#' @return a data frame containing the p values and q-values from the Moran's I test on the parallel arrays of models.
#' @importFrom spdep knn2nb nb2listw moran spweights.constants
#' @importFrom stats p.adjust
#' @importFrom pbmcapply pbmclapply
#' @seealso \code{\link[spdep]{moran.test}}
#' @export
principalGraphTest <- function(cds,
                               relative_expr=TRUE,
                               k = 25,
                               method = c('Moran_I'),
                               alternative = 'greater',
                               cores=1,
                               verbose=FALSE) {
  lw <- calculateLW(cds, verbose = verbose, k = k)
  
  if(verbose) {
    message("Performing Moran's test: ...")
  }
  exprs_mat <- exprs(cds)
  
  wc <- spweights.constants(lw, zero.policy = TRUE, adjust.n = TRUE)
  test_res <- pbmclapply(row.names(exprs_mat), FUN = function(x, alternative, method) {
    exprs_val <- exprs_mat[x, ]
    
    if (cds@expressionFamily@vfamily %in% c("gaussianff", "uninormal", "binomialff")){
      exprs_val <- exprs_val
    }else{
      if(relative_expr) {
        exprs_val <- log10(exprs_val / sizeFactors(cds) + 0.1)
      } else {
        exprs_val <- log10(exprs_val + 0.1)
      }
    }
    
    test_res <- tryCatch({
      if(method == "Moran_I") {
        mt <- my.moran.test(exprs_val, lw, wc, alternative = alternative)
        data.frame(status = 'OK', pval = mt$p.value, morans_test_statistic = mt$statistic, morans_I = mt$estimate[["Moran I statistic"]])
      } else if(method == 'Geary_C') {
        gt <- my.geary.test(exprs_val, lw, wc, alternative = alternative)
        data.frame(status = 'OK', pval = gt$p.value, geary_test_statistic = gt$statistic, geary_C = gt$estimate[["Geary C statistic"]])
      }
    },
    error = function(e) {
      data.frame(status = 'FAIL', pval = NA, morans_test_statistic = NA, morans_I = NA)
    })
  }, alternative = alternative, method = method, mc.cores = cores, ignore.interactive = TRUE)
  
  if(verbose) {
    message("returning results: ...")
  }
  
  test_res <- do.call(rbind.data.frame, test_res)
  row.names(test_res) <- row.names(cds)
  
  test_res <- merge(test_res, fData(cds), by="row.names")
  row.names(test_res) <- test_res[, 1] #remove the first column and set the row names to the first column
  test_res[, 1] <- NULL
  test_res$qval <- 1
  test_res$qval[which(test_res$status == 'OK')] <- p.adjust(subset(test_res, status == 'OK')[, 'pval'], method="BH")
  
  test_res[row.names(cds), ] # make sure gene name ordering in the DEG test result is the same as the CDS
}

my.moran.test <- function (x, listw, wc, alternative = "greater", randomisation = TRUE)
{
  zero.policy = TRUE
  adjust.n = TRUE
  na.action = na.fail
  drop.EI2 = FALSE
  xname <- deparse(substitute(x))
  wname <- deparse(substitute(listw))
  NAOK <- deparse(substitute(na.action)) == "na.pass"
  x <- na.action(x)
  na.act <- attr(x, "na.action")
  if (!is.null(na.act)) {
    subset <- !(1:length(listw$neighbours) %in% na.act)
    listw <- subset(listw, subset, zero.policy = zero.policy)
  }
  n <- length(listw$neighbours)
  if (n != length(x))
    stop("objects of different length")
  
  S02 <- wc$S0 * wc$S0
  res <- spdep::moran(x, listw, wc$n, wc$S0, zero.policy = zero.policy,
                      NAOK = NAOK)
  I <- res$I
  K <- res$K
  
  EI <- (-1)/wc$n1
  if (randomisation) {
    VI <- wc$n * (wc$S1 * (wc$nn - 3 * wc$n + 3) - wc$n *
                    wc$S2 + 3 * S02)
    tmp <- K * (wc$S1 * (wc$nn - wc$n) - 2 * wc$n * wc$S2 +
                  6 * S02)
    if (tmp > VI)
      warning("Kurtosis overflow,\ndistribution of variable does not meet test assumptions")
    VI <- (VI - tmp)/(wc$n1 * wc$n2 * wc$n3 * S02)
    if (!drop.EI2)
      VI <- (VI - EI^2)
    if (VI < 0)
      warning("Negative variance,\ndistribution of variable does not meet test assumptions")
  }
  else {
    VI <- (wc$nn * wc$S1 - wc$n * wc$S2 + 3 * S02)/(S02 *
                                                      (wc$nn - 1))
    if (!drop.EI2)
      VI <- (VI - EI^2)
    if (VI < 0)
      warning("Negative variance,\ndistribution of variable does not meet test assumptions")
  }
  ZI <- (I - EI)/sqrt(VI)
  statistic <- ZI
  names(statistic) <- "Moran I statistic standard deviate"
  if (alternative == "two.sided")
    PrI <- 2 * pnorm(abs(ZI), lower.tail = FALSE)
  else if (alternative == "greater")
    PrI <- pnorm(ZI, lower.tail = FALSE)
  else PrI <- pnorm(ZI)
  if (!is.finite(PrI) || PrI < 0 || PrI > 1)
    warning("Out-of-range p-value: reconsider test arguments")
  vec <- c(I, EI, VI)
  names(vec) <- c("Moran I statistic", "Expectation", "Variance")
  method <- paste("Moran I test under", ifelse(randomisation,
                                               "randomisation", "normality"))
  
  res <- list(statistic = statistic, p.value = PrI, estimate = vec)
  if (!is.null(na.act))
    attr(res, "na.action") <- na.act
  class(res) <- "htest"
  res
}

my.geary.test <- function (x, listw, wc, randomisation = TRUE, alternative = "greater")
{
  zero.policy = TRUE
  adjust.n = TRUE
  spChk = NULL
  
  alternative <- match.arg(alternative, c("less", "greater",
                                          "two.sided"))
  if (!inherits(listw, "listw"))
    stop(paste(deparse(substitute(listw)), "is not a listw object"))
  if (!is.numeric(x))
    stop(paste(deparse(substitute(x)), "is not a numeric vector"))
  if (any(is.na(x)))
    stop("NA in X")
  n <- length(listw$neighbours)
  if (n != length(x))
    stop("objects of different length")
  if (is.null(spChk))
    spChk <- get.spChkOption()
  if (spChk && !chkIDs(x, listw))
    stop("Check of data and weights ID integrity failed")
  
  S02 <- wc$S0 * wc$S0
  res <- spdep::geary(x, listw, wc$n, wc$n1, wc$S0, zero.policy)
  C <- res$C
  if (is.na(C))
    stop("NAs generated in geary - check zero.policy")
  K <- res$K
  EC <- 1
  if (randomisation) {
    VC <- (wc$n1 * wc$S1 * (wc$nn - 3 * n + 3 - K * wc$n1))
    VC <- VC - ((1/4) * (wc$n1 * wc$S2 * (wc$nn + 3 * n -
                                            6 - K * (wc$nn - n + 2))))
    VC <- VC + (S02 * (wc$nn - 3 - K * (wc$n1^2)))
    VC <- VC/(n * wc$n2 * wc$n3 * S02)
  }
  else {
    VC <- ((2 * wc$S1 + wc$S2) * wc$n1 - 4 * S02)/(2 * (n +
                                                          1) * S02)
  }
  ZC <- (EC - C)/sqrt(VC)
  statistic <- ZC
  names(statistic) <- "Geary C statistic standard deviate"
  PrC <- NA
  if (is.finite(ZC)) {
    if (alternative == "two.sided")
      PrC <- 2 * pnorm(abs(ZC), lower.tail = FALSE)
    else if (alternative == "greater")
      PrC <- pnorm(ZC, lower.tail = FALSE)
    else PrC <- pnorm(ZC)
    if (!is.finite(PrC) || PrC < 0 || PrC > 1)
      warning("Out-of-range p-value: reconsider test arguments")
  }
  vec <- c(C, EC, VC)
  names(vec) <- c("Geary C statistic", "Expectation", "Variance")
  method <- paste("Geary C test under", ifelse(randomisation,
                                               "randomisation", "normality"))
  data.name <- paste(deparse(substitute(x)), "\nweights:",
                     deparse(substitute(listw)), "\n")
  res <- list(statistic = statistic, p.value = PrC, estimate = vec,
              alternative = ifelse(alternative == "two.sided", alternative,
                                   paste("Expectation", alternative, "than statistic")),
              method = method, data.name = data.name)
  class(res) <- "htest"
  res
}

#' Find marker genes for each group of cells
#'
#' Tests each gene for differential expression as a function of pseudotime
#' or according to other covariates as specified. \code{differentialGeneTest} is
#' Monocle's main differential analysis routine.
#' It accepts a CellDataSet and two model formulae as input, which specify generalized
#' lineage models as implemented by the \code{VGAM} package.
#'
#' @seealso \code{\link{principalGraphTest}} principalGraphTest
#' @param cds a CellDataSet object upon which to perform this operation
#' @param spatial_res the result returned from spatialDifferentialTest
#' @param group_by a column in the pData specifying the groups for calculating the specifities. By default it is Cluster
#' @param qval_threshold The q-value threshold for genes to be selected
#' @param morans_I_threshold The lowest Morans' I threshold for selecting genes
#' @param lower_threshold The lowest gene expression threshold for genes to be considered as expressed
#' @param pseudocount Pseduo-count added to gene expression before calculating the log
#' @param top_n_by_group Select top_n_by_group from each group based on the specificity
#' @param block_size the number of genes to process per bloclk which dramatically reduces the memory footprints. Default to be 100.
#' @param verbose Whether to show VGAM errors and warnings. Only valid for cores = 1.
#' @param ... Additional arguments passed to function
#' @return a data frame containing the p values and q-values from the likelihood ratio tests on the parallel arrays of models.
#' @importFrom dplyr group_by summarize desc arrange top_n do
#' @importFrom reshape2 melt
#' @export
#'
find_cluster_markers <- function(cds,
                                 spatial_res,
                                 group_by = 'Cluster',
                                 qval_threshold = 0.05,
                                 morans_I_threshold = 0.25,
                                 lower_threshold = 0,
                                 pseudocount = 1,
                                 top_n_by_group = NULL,
                                 block_size = 100,
                                 verbose = FALSE,
                                 ...) {
  if(!(group_by %in% colnames(pData(cds)))) {
    stop('Please ensure group_by is included in the pData')
  }
  
  gene_ids <- row.names(subset(spatial_res, qval < qval_threshold & morans_I > morans_I_threshold))
  num_blocks = ceiling(length(gene_ids) / block_size)
  
  specificity_res <- NULL
  for(i in 1:num_blocks) {
    if (i < num_blocks){
      exprs_mat = as.matrix(cds@assayData$exprs[gene_ids[((((i-1) * block_size)+1):(i*block_size))], ])
    }else{
      exprs_mat = as.matrix(cds@assayData$exprs[gene_ids[((((i-1) * block_size)+1):(length(gene_ids)))], ])
    }
    
    exprs_mat <- melt(exprs_mat)
    colnames(exprs_mat) <- c('Gene', 'Cell', 'Expression')
    exprs_mat$Gene <- as.character(exprs_mat$Gene)
    exprs_mat$Group <- pData(cds)[exprs_mat$Cell, group_by]
    
    ExpVal <- exprs_mat %>% group_by(Group, Gene) %>% summarize(mean = log(mean(Expression) + pseudocount), percentage = sum(Expression > lower_threshold) / length(Expression))
    
    ExpVal <- merge(ExpVal, spatial_res, by.x = 'Gene', by = "row.names")
    ExpVal$Group <- ExpVal$Group
    
    FUN <- function(df) {
      class_df <- data.frame(Group = df$Group, percentage = df$percentage)
      uniq_group <- unique(df$Group)
      specificity <- rep(0, length(uniq_group))
      for(cell_type_i in 1:length(uniq_group)) {
        perfect_specificity <- rep(0.0, nrow(class_df))
        perfect_specificity[cell_type_i] <- 1.0
        
        if(sum(class_df$percentage) > 0) {
          specificity[cell_type_i] <- 1 - JSdistVec(makeprobsvec(class_df$percentage), perfect_specificity)
        } else {
          specificity[cell_type_i] <- 0
        }
      }
      specificity
    }
    
    tmp <- ExpVal %>% group_by(Gene) %>% do({
      tmp <- dplyr::as_data_frame(.)
      tmp$specificity = FUN(tmp)
      tmp
    })
    
    specificity_res <- rbind(tmp, specificity_res)
  }
  
  specificity_res <- specificity_res %>% arrange(Group, desc(specificity), desc(-qval), desc(morans_I))
  if(!is.null(top_n_by_group)) {
    specificity_res <- specificity_res %>% group_by(Group) %>% top_n(n = top_n_by_group, wt = specificity)
  }
  specificity_res
}
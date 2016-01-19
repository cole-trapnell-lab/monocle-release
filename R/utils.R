
#' Creates a new CellDateSet object.
#'
#' @param cellData expression data matrix for an experiment
#' @param phenoData data frame containing attributes of individual cells
#' @param featureData data frame containing attributes of features (e.g. genes)
#' @param lowerDetectionLimit the minimum expression level that consistitutes true expression
#' @param expressionFamily the VGAM family function to be used for expression response variables
#' @return a new CellDataSet object
#' @export
#' @examples
#' \dontrun{
#' sample_sheet_small <- read.delim("../data/sample_sheet_small.txt", row.names=1)
#' sample_sheet_small$Time <- as.factor(sample_sheet_small$Time)
#' gene_annotations_small <- read.delim("../data/gene_annotations_small.txt", row.names=1)
#' fpkm_matrix_small <- read.delim("../data/fpkm_matrix_small.txt")
#' pd <- new("AnnotatedDataFrame", data = sample_sheet_small)
#' fd <- new("AnnotatedDataFrame", data = gene_annotations_small)
#' HSMM <- new("CellDataSet", exprs = as.matrix(fpkm_matrix_small), phenoData = pd, featureData = fd)
#' }
newCellDataSet <- function( cellData, 
                            phenoData = NULL, 
                            featureData = NULL, 
                            lowerDetectionLimit = 0.1, 
                            expressionFamily=VGAM::tobit(Lower=log10(lowerDetectionLimit), lmu="identitylink"))
{
  #cellData <- as.matrix( cellData )
  
  
  sizeFactors <- rep( NA_real_, ncol(cellData) )
  
  
  if( is.null( phenoData ) )
    phenoData <- annotatedDataFrameFrom( cellData, byrow=FALSE )
  if( is.null( featureData ) ) 
    featureData <- annotatedDataFrameFrom( cellData, byrow=TRUE )
  
  phenoData$`Size_Factor` <- sizeFactors
  
  cds <- new( "CellDataSet",
              assayData = assayDataNew( "environment", exprs=cellData ),
              phenoData=phenoData, 
              featureData=featureData, 
              lowerDetectionLimit=lowerDetectionLimit,
              expressionFamily=expressionFamily,
              dispFitInfo = new.env( hash=TRUE ))
  
  validObject( cds )
  cds
}




#' Multicore esApply wrapper
#'
#' @importFrom parallel makeCluster stopCluster clusterCall parRapply parCapply
mcesApply <- function(X, MARGIN, FUN, required_packages, cores=1, ...) {
  parent <- environment(FUN)
  if (is.null(parent))
    parent <- emptyenv()
  e1 <- new.env(parent=parent)
  multiassign(names(pData(X)), pData(X), envir=e1)
  environment(FUN) <- e1
  cl <- parallel::makeCluster(cores)
  
  cleanup <- function(){
    parallel::stopCluster(cl)
  }
  on.exit(cleanup)
  
  if (is.null(required_packages) == FALSE){
    parallel::clusterCall(cl, function(pkgs) {
      for (req in pkgs) {
        library(req, character.only=TRUE)
      }
    }, required_packages)
  }
  
  if (MARGIN == 1){
    res <- parRapply(cl, exprs(X), FUN, ...)
  }else{
    res <- parCapply(cl, exprs(X), FUN, ...)
  }
  
  res
}

smartEsApply <- function(X, MARGIN, FUN, ...) {
  parent <- environment(FUN)
  if (is.null(parent))
    parent <- emptyenv()
  e1 <- new.env(parent=parent)
  multiassign(names(pData(X)), pData(X), envir=e1)
  environment(FUN) <- e1
  sp_X <- asSlamMatrix(exprs(X))
  if (MARGIN == 1){
    res <- rowapply_simple_triplet_matrix(sp_X, FUN, ...)
  }else{
    res <- colapply_simple_triplet_matrix(sp_X, FUN, ...)
  }
  
  res
}


####
#' Filter genes with extremely high or low negentropy
#'
#' @param cds a CellDataSet object upon which to perform this operation
#' @param lower_negentropy_bound the centile below which to exclude to genes 
#' @param upper_negentropy_bound the centile above which to exclude to genes
#' @param expression_lower_thresh the expression level below which to exclude genes used to determine negentropy
#' @param expression_upper_thresh the expression level above which to exclude genes used to determine negentropy
#' @return a vector of gene names
#' @export
#' @examples
#' \dontrun{
#' reasonableNegentropy <- selectNegentropyGenes(HSMM, "07%", "95%", 1, 100)
#' }
selectNegentropyGenes <- function(cds, lower_negentropy_bound="0%",
                                  upper_negentropy_bound="99%", 
                                  expression_lower_thresh=0.1,
                                  expression_upper_thresh=Inf){
  
  log_expression <- NULL
  
  FM <- exprs(cds)
  if (cds@expressionFamily@vfamily == "negbinomial")
  {
    expression_lower_thresh <- expression_lower_thresh / colSums(FM)
    expression_upper_thresh <- expression_upper_thresh / colSums(FM)
    FM <- Matrix::t(Matrix::t(FM)/colSums(FM))
  }
  
  negentropy_exp <- apply(FM,1,function(x) { 
    
    expression <- x[x > expression_lower_thresh]
    expression <- log2(x); 
    expression <- expression[is.finite(expression)]
    
    if (length(expression)){
      expression <- scale(expression)
      mean(-exp(-(expression^2)/2))^2
    }else{
      0
    }
  }
  
  )
  
  
  means <- apply(FM,1,function(x) { 
    expression <- x[x > expression_lower_thresh]
    expression <- log2(x); 
    expression[is.finite(expression) == FALSE] <- NA; 
    
    if (length(expression)){
      mean(expression, na.rm=T)
    }else{
      NA
    }
  }
  )
  ordering_df <- data.frame(log_expression = means, negentropy = negentropy_exp)
  
  negentropy <- NULL
  log_express <- NULL
  negentropy_residual <- NULL
  
  ordering_df <- subset(ordering_df, 
                        is.na(log_expression) == FALSE &
                          is.nan(log_expression) == FALSE &
                          is.na(negentropy) == FALSE &
                          is.nan(negentropy) == FALSE)
  negentropy_fit <- vglm(negentropy~sm.ns(log_expression, df=4),data=ordering_df, family=VGAM::gaussianff())
  ordering_df$negentropy_response <- predict(negentropy_fit, newdata=ordering_df, type="response")
  ordering_df$negentropy_residual <- ordering_df$negentropy - ordering_df$negentropy_response
  lower_negentropy_thresh <- quantile(ordering_df$negentropy_residual, probs=seq(0,1,0.01), na.rm=T)[lower_negentropy_bound]
  upper_negentropy_thresh <- quantile(ordering_df$negentropy_residual, probs=seq(0,1,0.01), na.rm=T)[upper_negentropy_bound]
  
  ordering_genes <- row.names(subset(ordering_df, 
                                     negentropy_residual >= lower_negentropy_thresh & 
                                       negentropy_residual <= upper_negentropy_thresh))
  ordering_genes
}


# TODO: we need to rename this function and its arguments.  What it actually
# does is very confusing.
####
#' Filter genes outside of a given range of expression
#'
#' @export
selectGenesInExpressionRange <- function(cds, 
                                         min_expression_threshold = -Inf, 
                                         max_expression_threshold = Inf, 
                                         detectionLimit=-Inf, 
                                         stat_fun=median, 
                                         relative_expr=TRUE)
{
  gene_nz_median = apply(exprs(cds), 1, function(x) { x <- x[x > detectionLimit]; stat_fun(x)})
  gene_nz_median<-esApply(cds,1,
                          function(x) { 
                            if (relative_expr && cds@expressionFamily@vfamily == "negbinomial"){
                              x <- x / Size_Factor
                            }
                            x <- x[x > detectionLimit]
                            stat_fun(x)
                          })
  
  #gene_nz_median
  names(gene_nz_median[is.na(gene_nz_median) == FALSE & gene_nz_median > min_expression_threshold & gene_nz_median < max_expression_threshold ])
}


# TODO: we need to rename this function and its arguments.  What it actually
# does is very confusing.
####

#' Compute a table of 
#'
#' @export
dispersionTable <- function(cds){
  disp_df<-data.frame(row.names=row.names(cds@dispFitInfo[["blind"]]$disp_table),
                      mean_expression=cds@dispFitInfo[["blind"]]$disp_table$mu, 
                      dispersion_fit=cds@dispFitInfo[["blind"]]$disp_func(cds@dispFitInfo[["blind"]]$disp_table$mu),
                      dispersion_empirical=cds@dispFitInfo[["blind"]]$disp_table$disp)
  return(disp_df)
}

#####
#' Sets the global expression detection threshold to be used with this CellDataSet.
#' Counts how many cells each feature in a CellDataSet object that are detectably expressed 
#' above a minimum threshold. Also counts the number of genes above this threshold are 
#' detectable in each cell.
#'
#' @param cds the CellDataSet upon which to perform this operation
#' @param min_expr the expression threshold 
#' @return an updated CellDataSet object
#' @export
#' @examples
#' \dontrun{
#' HSMM <- detectGenes(HSMM, min_expr=0.1)
#' }
detectGenes <- function(cds, min_expr=NULL){
  if (is.null(min_expr))
  {
    min_expr <- cds@lowerDetectionLimit
  }
#   FM_genes <- do.call(rbind, apply(FM, 1, 
#                                    function(x) {
#                                      return(data.frame(
#                                        num_cells_expressed=sum(unlist(as.list(x)) >= min_expr)
#                                      )
#                                      )
#                                    })
#   )
#   
#   FM_cells <- do.call(rbind, apply(FM, 2, 
#                                    function(x) {
#                                      return(data.frame(
#                                        num_genes_expressed=sum(unlist(as.list(x)) >= min_expr)
#                                      )
#                                      )
#                                    })
#   )
#   
#   
#   
#   fData(cds)$num_cells_expressed <-  FM_genes[row.names(fData(cds)),]
#   
#   pData(cds)$num_genes_expressed <-  FM_cells[row.names(pData(cds)),]
#   
  fData(cds)$num_cells_expressed <- Matrix::rowSums(exprs(cds) > min_expr)
  pData(cds)$num_genes_expressed <- Matrix::colSums(exprs(cds))

  cds
}

# Convert a slam matrix to a sparseMatrix
#' @import slam
#' @import Matrix
asSparseMatrix = function (simpleTripletMatrix) {
  retVal = sparseMatrix(i=simpleTripletMatrix[["i"]],
                        j=simpleTripletMatrix[["j"]],
                        x=simpleTripletMatrix[["v"]],
                        dims=c(simpleTripletMatrix[["nrow"]],
                               simpleTripletMatrix[["ncol"]]))
  if (!is.null(simpleTripletMatrix[["dimnames"]]))
    dimnames(retVal) = simpleTripletMatrix[["dimnames"]]
  return(retVal)
}

# Convert a sparseMatrix from Matrix package to a slam matrix
#' @import slam
asSlamMatrix = function (sp_mat) {
  sp <- Matrix::summary(sp_mat)
  simple_triplet_matrix(sp[,"i"], sp[,"j"], sp[,"x"], ncol=ncol(sp_mat), nrow=nrow(sp_mat), dimnames=dimnames(sp_mat))
}

# Convert a sparseMatrix from Matrix package to a slam matrix
#' @import Matrix
isSparseMatrix <- function(x){
  class(x) %in% c("dgCMatrix", "dgTMatrix")
}

# Estimate size factors for each column, given a sparseMatrix from the Matrix
# package
#' @import slam
estimateSizeFactorsForSparseMatrix <- function(counts, 
                                               locfunc = median, 
                                               round_exprs=TRUE, 
                                               pseudocount=0.0, 
                                               method="mean-geometric-mean-total"){
  CM <- counts
  if (round_exprs)
    CM <- round(CM)
  CM <- CM + pseudocount
  CM <- asSlamMatrix(CM)
  
  if (method == "weighted-median"){

    log_medians <- rowapply_simple_triplet_matrix(CM, function(cell_expr) { 
      log(locfunc(cell_expr))
    })
    
    weights <- rowapply_simple_triplet_matrix(CM, function(cell_expr) {
      num_pos <- sum(cell_expr > 0)
      num_pos / length(cell_expr)
    })
    
    sfs <- colapply_simple_triplet_matrix(CM, function(cnts) {
      norm_cnts <-  weights * (log(cnts) -  log_medians)
      norm_cnts <- norm_cnts[is.nan(norm_cnts) == FALSE]
      norm_cnts <- norm_cnts[is.finite(norm_cnts)]
      #print (head(norm_cnts))
      exp( mean(norm_cnts) )
    })
  }else if (method == "median-geometric-mean"){
    log_geo_means <- rowapply_simple_triplet_matrix(CM, function(x) { mean(log(CM)) })
    
    sfs <- colapply_simple_triplet_matrix(CM, function(cnts) {
      norm_cnts <- log(cnts) -  log_geo_means
      norm_cnts <- norm_cnts[is.nan(norm_cnts) == FALSE]
      norm_cnts <- norm_cnts[is.finite(norm_cnts)]
      #print (head(norm_cnts))
      exp( locfunc( norm_cnts ))
    })
  }else if(method == "median"){
    stop("Error: method 'median' not yet supported for sparse matrices")
  }else if(method == 'mode'){
    stop("Error: method 'mode' not yet supported for sparse matrices")
  }else if(method == 'geometric-mean-total') {
    cell_total <- col_sums(CM)
    sfs <- log(cell_total) / mean(log(cell_total))
  }else if(method == 'mean-geometric-mean-total') {
    cell_total <- col_sums(CM)
    sfs <- cell_total / exp(mean(log(cell_total)))
  } 
  
  sfs[is.na(sfs)] <- 1 
  sfs   
}

estimateSizeFactorsForDenseMatrix <- function(counts, locfunc = median, round_exprs=TRUE, pseudocount=0.0, method="mean-geometric-mean-total"){
  
  CM <- counts
  if (round_exprs)
    CM <- round(CM)
  CM <- CM + pseudocount
  if (method == "weighted-median"){
    log_medians <- apply(CM, 1, function(cell_expr) { 
      log(locfunc(cell_expr))
    })
    
    weights <- apply(CM, 1, function(cell_expr) {
      num_pos <- sum(cell_expr > 0)
      num_pos / length(cell_expr)
    })
    
    sfs <- apply( CM, 2, function(cnts) {
      norm_cnts <-  weights * (log(cnts) -  log_medians)
      norm_cnts <- norm_cnts[is.nan(norm_cnts) == FALSE]
      norm_cnts <- norm_cnts[is.finite(norm_cnts)]
      #print (head(norm_cnts))
      exp( mean(norm_cnts) )
    })
  }else if (method == "median-geometric-mean"){
    log_geo_means <- rowMeans(log(CM))
    
    sfs <- apply( CM, 2, function(cnts) {
      norm_cnts <- log(cnts) -  log_geo_means
      norm_cnts <- norm_cnts[is.nan(norm_cnts) == FALSE]
      norm_cnts <- norm_cnts[is.finite(norm_cnts)]
      #print (head(norm_cnts))
      exp( locfunc( norm_cnts ))
    })
  }else if(method == "median"){
    row_median <- apply(CM, 1, median)
    sfs <- apply(Matrix::t(Matrix::t(CM) - row_median), 2, median)
  }else if(method == 'mode'){
    sfs <- estimate_t(CM)
  }else if(method == 'geometric-mean-total') {
    cell_total <- apply(CM, 2, sum)
    sfs <- log(cell_total) / mean(log(cell_total))
  }else if(method == 'mean-geometric-mean-total') {
    cell_total <- apply(CM, 2, sum)
    sfs <- cell_total / exp(mean(log(cell_total)))
  } 
  
  sfs[is.na(sfs)] <- 1 
  sfs  
}



#' Function to calculate the size factor for the single-cell RNA-seq data
#'  
#' @param counts The matrix for the gene expression data, either read counts or FPKM values or transcript counts
#' @param locfunc The location function used to find the representive value 
#' @param round_exprs A logic flag to determine whether or not the expression value should be rounded
#' @param pseudocount Pseudo count added to the expression data counts 
#' @param method A character to specify the size factor calculation appraoches. It can be either "mean-geometric-mean-total" (default), 
#' "weighted-median", "median-geometric-mean", "median", "mode", "geometric-mean-total". 
#' @export
#'
estimateSizeFactorsForMatrix <- function(counts, locfunc = median, round_exprs=TRUE, pseudocount=0.0, method="mean-geometric-mean-total")
{
  if (isSparseMatrix(counts)){
    estimateSizeFactorsForSparseMatrix(counts, locfunc = median, round_exprs=TRUE, pseudocount=0.0, method="mean-geometric-mean-total")
  }else{
    estimateSizeFactorsForDenseMatrix(counts, locfunc = median, round_exprs=TRUE, pseudocount=0.0, method="mean-geometric-mean-total")
  }
  
}

################

# Some convenience functions for loading the HSMM data

#' Return the names of classic muscle genes
#' @export
get_classic_muscle_markers <- function(){
  c("MEF2C", "MEF2D", "MYF5", "ANPEP", "PDGFRA",
    "MYOG", "TPM1", "TPM2", "MYH2", "MYH3", "NCAM1", "TNNT1", "TNNT2", "TNNC1",
    "CDK1", "CDK2", "CCNB1", "CCNB2", "CCND1", "CCNA1", "ID1")
}

#' Build a CellDataSet from the HSMMSingleCell package
#' 
#' @export
load_HSMM <- function(){
  
  data(HSMM_expr_matrix, envir = environment())
  data(HSMM_gene_annotation, envir = environment())
  data(HSMM_sample_sheet, envir = environment())
  pd <- new("AnnotatedDataFrame", data = HSMM_sample_sheet)
  fd <- new("AnnotatedDataFrame", data = HSMM_gene_annotation)
  HSMM <- newCellDataSet(as.matrix(HSMM_expr_matrix), phenoData = pd, featureData = fd)
  HSMM
}

#' Return a CellDataSet of classic muscle genes
#' @export
load_HSMM_markers <- function(){
  HSMM <- load_HSMM()
  marker_names <- get_classic_muscle_markers()
  HSMM[row.names(subset(fData(HSMM), gene_short_name %in% marker_names)),]
}

#' Build a CellDataSet from the data stored in inst/extdata directory
#' 
#' @export
load_lung <- function(){
  baseLoc <- system.file(package="monocle")
  #baseLoc <- './inst'
  extPath <- file.path(baseLoc, "extdata")
  load(file.path(extPath, "lung_phenotype_data.RData"))
  load(file.path(extPath, "lung_exprs_data.RData"))
  load(file.path(extPath, "lung_feature_data.RData"))
  lung_exprs_data <- lung_exprs_data[,row.names(lung_phenotype_data)]

  pd <- new("AnnotatedDataFrame", data = lung_phenotype_data)
  fd <- new("AnnotatedDataFrame", data = lung_feature_data)

  # Now, make a new CellDataSet using the RNA counts
  lung <- newCellDataSet(as.matrix(lung_exprs_data), 
                         phenoData = pd, 
                         featureData = fd,
                         lowerDetectionLimit=1,
                         expressionFamily=negbinomial())

  lung <- estimateSizeFactors(lung)
  lung <- estimateDispersions(lung)

  pData(lung)$Total_mRNAs <- colSums(exprs(lung))
  lung <- detectGenes(lung, min_expr = 1)
  expressed_genes <- row.names(subset(fData(lung), num_cells_expressed >= 5))
  ordering_genes <- expressed_genes
  lung <- setOrderingFilter(lung, ordering_genes)
  lung <- reduceDimension(lung, use_vst = F, pseudo_expr = 1)
  lung <- orderCells(lung)
  lung <- orderCells(lung, root_state=3)

  lung
}
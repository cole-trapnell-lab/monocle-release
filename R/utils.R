
#' Creates a new CellDateSet object.
#'
#' @param cellData expression data matrix for an experiment
#' @param phenoData data frame containing attributes of individual cells
#' @param featureData data frame containing attributes of features (e.g. genes)
#' @param lowerDetectionLimit the minimum expression level that consistitutes true expression
#' @param expressionFamily the VGAM family function to be used for expression response variables
#' @return a new CellDataSet object
#' @import VGAM
<<<<<<< 7afe2e5a61fe398816056758e84b77f3d44ceb50
=======
#' @importFrom Biobase annotatedDataFrameFrom assayDataNew
>>>>>>> First commit
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
  
  if (class(cellData) != "matrix" && isSparseMatrix(cellData) == FALSE){
    stop("Error: argument cellData must be a matrix (either sparse from the Matrix package or dense)")
  }
  
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

sparseApply <- function(Sp_X, MARGIN, FUN, convert_to_dense, ...){
  if (convert_to_dense){
    if (MARGIN == 1){
      Sp_X <- Matrix::t(Sp_X)
      res <- lapply(colnames(Sp_X), function(i, FUN, ...) {
        FUN(as.matrix(Sp_X[,i]), ...) 
      }, FUN, ...)
    }else{
      res <- lapply(colnames(Sp_X), function(i, FUN, ...) {
        FUN(as.matrix(Sp_X[,i]), ...) 
      }, FUN, ...)
    }
  }else{
    if (MARGIN == 1){
      Sp_X <- Matrix::t(Sp_X)
      res <- lapply(colnames(Sp_X), function(i, FUN, ...) {
        FUN(Sp_X[,i], ...) 
      }, FUN, ...)
    }else{
      res <- lapply(colnames(Sp_X), function(i, FUN, ...) {
        FUN(Sp_X[,i], ...) 
      }, FUN, ...)
    }
  }

  return(res)
  
}

<<<<<<< 7afe2e5a61fe398816056758e84b77f3d44ceb50
=======
#' @importFrom parallel splitIndices
>>>>>>> First commit
splitRows <- function (x, ncl) {
  lapply(splitIndices(nrow(x), ncl), function(i) x[i, , drop = FALSE])
}

<<<<<<< 7afe2e5a61fe398816056758e84b77f3d44ceb50
=======
#' @importFrom parallel splitIndices
>>>>>>> First commit
splitCols <- function (x, ncl) {
  lapply(splitIndices(ncol(x), ncl), function(i) x[, i, drop = FALSE])
}

#' @importFrom BiocGenerics clusterApply
sparseParRApply <- function (cl, x, FUN, convert_to_dense, ...) 
{
  par_res <- do.call(c, clusterApply(cl = cl, x = splitRows(x, length(cl)), 
                          fun = sparseApply, MARGIN = 1L, FUN = FUN, convert_to_dense=convert_to_dense, ...), quote = TRUE)
  names(par_res) <- row.names(x)
  par_res
}

#' @importFrom BiocGenerics clusterApply
sparseParCApply <- function (cl = NULL, x, FUN, convert_to_dense, ...) 
{
  par_res <- do.call(c, clusterApply(cl = cl, x = splitCols(x, length(cl)), 
                          fun = sparseApply, MARGIN = 2L, FUN = FUN, convert_to_dense=convert_to_dense, ...), quote = TRUE)
  names(par_res) <- colnames(x)
  par_res
}


#' Multicore apply-like function for CellDataSet
#' 
#' mcesApply computes the row-wise or column-wise results of FUN, just like esApply.
#' Variables in pData from X are available in FUN. 
#'
#' @param X a CellDataSet object
#' @param MARGIN The margin to apply to, either 1 for rows (samples) or 2 for columns (features)
#' @param FUN Any function
#' @param required_packages A list of packages FUN will need. Failing to provide packages needed by FUN will generate errors in worker threads.
#' @param convert_to_dense Whether to force conversion a sparse matrix to a dense one before calling FUN
#' @param ... Additional parameters for FUN
#' @param cores The number of cores to use for evaluation
#' 
#' @return The result of with(pData(X) apply(exprs(X)), MARGIN, FUN, ...))
#' @importFrom parallel makeCluster stopCluster
#' @importFrom BiocGenerics clusterCall parRapply parCapply
<<<<<<< 7afe2e5a61fe398816056758e84b77f3d44ceb50
=======
#' @importFrom Biobase pData exprs multiassign
>>>>>>> First commit
#' @export
mcesApply <- function(X, MARGIN, FUN, required_packages, cores=1, convert_to_dense=TRUE, ...) {
  parent <- environment(FUN)
  if (is.null(parent))
    parent <- emptyenv()
  e1 <- new.env(parent=parent)
  multiassign(names(pData(X)), pData(X), envir=e1)
  environment(FUN) <- e1
  
  # Note: use outfile argument to makeCluster for debugging
  platform <- Sys.info()[['sysname']]
  if (platform == "Windows")
    cl <- makeCluster(cores)
  if (platform %in% c("Linux", "Darwin")) 
    cl <- makeCluster(cores)
  
  cleanup <- function(){
    stopCluster(cl)
  }
  on.exit(cleanup)
  
  if (is.null(required_packages) == FALSE){
    clusterCall(cl, function(pkgs) {
      for (req in pkgs) {
        library(req, character.only=TRUE)
      }
    }, required_packages)
  }
  #clusterExport(cl, ls(e1), e1)
  #force(exprs(X))
  if (MARGIN == 1){
    suppressWarnings(res <- sparseParRApply(cl, exprs(X), FUN, convert_to_dense, ...))
  }else{
    suppressWarnings(res <- sparseParCApply(cl, exprs(X), FUN, convert_to_dense, ...))
  }
  
  res
}

<<<<<<< 7afe2e5a61fe398816056758e84b77f3d44ceb50
=======
#' @importFrom Biobase multiassign
>>>>>>> First commit
smartEsApply <- function(X, MARGIN, FUN, convert_to_dense, ...) {
  parent <- environment(FUN)
  if (is.null(parent))
    parent <- emptyenv()
  e1 <- new.env(parent=parent)
  multiassign(names(pData(X)), pData(X), envir=e1)
  environment(FUN) <- e1
  
  if (isSparseMatrix(exprs(X))){
    res <- sparseApply(exprs(X), MARGIN, FUN, convert_to_dense, ...)
  }else{
    res <- apply(exprs(X), MARGIN, FUN, ...)
  }
  
  if (MARGIN == 1)
  {
    names(res) <- row.names(X)
  }else{
    names(res) <- colnames(X)
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
<<<<<<< 7afe2e5a61fe398816056758e84b77f3d44ceb50
=======
#' @importFrom stats quantile
>>>>>>> First commit
#' @export
#' @examples
#' \dontrun{
#' reasonableNegentropy <- selectNegentropyGenes(HSMM, "07%", "95%", 1, 100)
#' }
selectNegentropyGenes <- function(cds, lower_negentropy_bound="0%",
                                  upper_negentropy_bound="99%", 
                                  expression_lower_thresh=0.1,
                                  expression_upper_thresh=Inf){
  .Deprecated("dispersionTable")
  log_expression <- NULL
  
  FM <- exprs(cds)
  if (cds@expressionFamily@vfamily %in% c("negbinomial", "negbinomial.size"))
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



#' Retrieve a table of values specifying the mean-variance relationship
#' 
#' Calling estimateDispersions computes a smooth function describing how variance
#' in each gene's expression across cells varies according to the mean. This 
#' function only works for CellDataSet objects containing count-based expression
#' data, either transcripts or reads.
#' 
#' @param cds The CellDataSet from which to extract a dispersion table.
#' @return A data frame containing the empirical mean expression, 
#' empirical dispersion, and the value estimated by the dispersion model. 
#'
#' @export
dispersionTable <- function(cds){
  
  if (is.null(cds@dispFitInfo[["blind"]]))
    stop("Error: no dispersion model found. Please call estimateDispersions() before calling this function.")
  
  disp_df<-data.frame(gene_id=cds@dispFitInfo[["blind"]]$disp_table$gene_id,
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
<<<<<<< 7afe2e5a61fe398816056758e84b77f3d44ceb50
=======
#' @importFrom Biobase fData fData<- exprs pData pData<- 
>>>>>>> First commit
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
  pData(cds)$num_genes_expressed <- Matrix::colSums(exprs(cds) > min_expr)

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
<<<<<<< 7afe2e5a61fe398816056758e84b77f3d44ceb50
=======
#' @importFrom stats median
>>>>>>> First commit
estimateSizeFactorsForSparseMatrix <- function(counts, 
                                               locfunc = median, 
                                               round_exprs=TRUE, 
                                               method="mean-geometric-mean-total"){
  CM <- counts
  if (round_exprs)
    CM <- round(CM)
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

<<<<<<< 7afe2e5a61fe398816056758e84b77f3d44ceb50
=======
#' @importFrom stats median
>>>>>>> First commit
estimateSizeFactorsForDenseMatrix <- function(counts, locfunc = median, round_exprs=TRUE, method="mean-geometric-mean-total"){
  
  CM <- counts
  if (round_exprs)
    CM <- round(CM)
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
<<<<<<< 7afe2e5a61fe398816056758e84b77f3d44ceb50
=======
#'  @importFrom stats median
>>>>>>> First commit
#' @param counts The matrix for the gene expression data, either read counts or FPKM values or transcript counts
#' @param locfunc The location function used to find the representive value 
#' @param round_exprs A logic flag to determine whether or not the expression value should be rounded
#' @param method A character to specify the size factor calculation appraoches. It can be either "mean-geometric-mean-total" (default), 
#' "weighted-median", "median-geometric-mean", "median", "mode", "geometric-mean-total". 
#'
estimateSizeFactorsForMatrix <- function(counts, locfunc = median, round_exprs=TRUE,  method="mean-geometric-mean-total")
{
  if (isSparseMatrix(counts)){
    estimateSizeFactorsForSparseMatrix(counts, locfunc = locfunc, round_exprs=round_exprs, method=method)
  }else{
    estimateSizeFactorsForDenseMatrix(counts, locfunc = locfunc, round_exprs=round_exprs,  method=method)
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
<<<<<<< 7afe2e5a61fe398816056758e84b77f3d44ceb50
#' @export
load_HSMM <- function(){
  
  data(HSMM_expr_matrix, envir = environment())
  data(HSMM_gene_annotation, envir = environment())
  data(HSMM_sample_sheet, envir = environment())
=======
#' @import HSMMSingleCell
#' @importFrom utils data
#' @export
load_HSMM <- function(){
  
  HSMM_sample_sheet <- NA
  HSMM_gene_annotation <- NA
  HSMM_expr_matrix <- NA
  gene_short_name <- NA
  
>>>>>>> First commit
  pd <- new("AnnotatedDataFrame", data = HSMM_sample_sheet)
  fd <- new("AnnotatedDataFrame", data = HSMM_gene_annotation)
  HSMM <- newCellDataSet(as.matrix(HSMM_expr_matrix), phenoData = pd, featureData = fd)
  HSMM
}

#' Return a CellDataSet of classic muscle genes
<<<<<<< 7afe2e5a61fe398816056758e84b77f3d44ceb50
#' @return A CellDataSet object
#' @export
load_HSMM_markers <- function(){
=======
#' @importFrom Biobase fData
#' @export
load_HSMM_markers <- function(){
  gene_short_name <- NA
>>>>>>> First commit
  HSMM <- load_HSMM()
  marker_names <- get_classic_muscle_markers()
  HSMM[row.names(subset(fData(HSMM), gene_short_name %in% marker_names)),]
}

#' Build a CellDataSet from the data stored in inst/extdata directory
<<<<<<< 7afe2e5a61fe398816056758e84b77f3d44ceb50
#' @export
load_lung <- function(){
=======
#' @importFrom Biobase pData pData<- exprs fData
#' @importFrom DESeq2 estimateSizeFactors estimateDispersions
#' @export
load_lung <- function(){
  lung_phenotype_data <- NA
  lung_feature_data <- NA
  num_cells_expressed <- NA
  
>>>>>>> First commit
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
  lung <- newCellDataSet(lung_exprs_data, 
                         phenoData = pd, 
                         featureData = fd,
                         lowerDetectionLimit=1,
                         expressionFamily=negbinomial.size())

  lung <- estimateSizeFactors(lung)
  pData(lung)$Size_Factor <- lung_phenotype_data$Size_Factor

  lung <- estimateDispersions(lung)

  pData(lung)$Total_mRNAs <- colSums(exprs(lung))
  lung <- detectGenes(lung, min_expr = 1)
  expressed_genes <- row.names(subset(fData(lung), num_cells_expressed >= 5))
  ordering_genes <- expressed_genes
  lung <- setOrderingFilter(lung, ordering_genes)
  
  # DDRTree based ordering:
  lung <- reduceDimension(lung, norm_method="log", pseudo_expr = 1) #
  lung <- orderCells(lung)
  E14_state = as.numeric(pData(lung)['SRR1033936_0', 'State'])
  if(E14_state != 1)
    lung <- orderCells(lung, root_state=E14_state)

  lung
}

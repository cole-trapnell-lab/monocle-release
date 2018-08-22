
#' Creates a new CellDateSet object.
#'
#' @param cellData expression data matrix for an experiment
#' @param phenoData data frame containing attributes of individual cells
#' @param featureData data frame containing attributes of features (e.g. genes)
#' @param lowerDetectionLimit the minimum expression level that consistitutes true expression
#' @param expressionFamily the VGAM family function to be used for expression response variables
#' @return a new CellDataSet object
#' @importFrom VGAM negbinomial.size
#' @importFrom Rcpp compileAttributes
#' @importFrom Biobase annotatedDataFrameFrom assayDataNew
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
                            expressionFamily=VGAM::negbinomial.size())
{
  #cellData <- as.matrix( cellData )
  
  if(!('gene_short_name' %in% colnames(featureData))) {
    warning("Warning: featureData must contain a column verbatim named 'gene_short_name' for certain functions")
  }
  
  if (class(cellData) != "matrix" && isSparseMatrix(cellData) == FALSE){
    stop("Error: argument cellData must be a matrix (either sparse from the Matrix package or dense)")
  }
  
  if(!('gene_short_name' %in% colnames(featureData))) {
   warning("Warning: featureData must contain a column verbatim named 'gene_short_name' for certain functions") 
  }
  
  sizeFactors <- rep( NA_real_, ncol(cellData) )
  
  
  if( is.null( phenoData ) )
    phenoData <- annotatedDataFrameFrom( cellData, byrow=FALSE )
  if( is.null( featureData ) ) 
    featureData <- annotatedDataFrameFrom(cellData, byrow=TRUE)
  
  if(!('gene_short_name' %in% colnames(featureData))) {
    warning("Warning: featureData must contain a column verbatim named 'gene_short_name' for certain functions")
  }
  
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

#' @importFrom parallel splitIndices
splitRows <- function (x, ncl) {
  lapply(splitIndices(nrow(x), ncl), function(i) x[i, , drop = FALSE])
}

#' @importFrom parallel splitIndices
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
#' @importFrom Biobase pData exprs multiassign
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

#' @importFrom Biobase multiassign
#' @importFrom pbapply pbapply
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
    res <- pbapply(exprs(X), MARGIN, FUN, ...)
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
#' @description Examines all the genes in the CellDataSet passed in and removes
#' all the genes that contain extremely high or low negentropies. You can specify
#' which genes to filter out based on the boundaries you can set for expression levels
#' and the boundaries you set for which centile to include. the function "dispersionTable"
#' is a better form of this function.
#'
#' @param cds a CellDataSet object upon which to perform this operation
#' @param lower_negentropy_bound the centile below which to exclude to genes 
#' @param upper_negentropy_bound the centile above which to exclude to genes
#' @param expression_lower_thresh the expression level below which to exclude genes used to determine negentropy
#' @param expression_upper_thresh the expression level above which to exclude genes used to determine negentropy
#' @return a vector of gene names
#' @importFrom stats quantile
#' @importFrom VGAM gaussianff
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
  
  if (is.null(cds@dispFitInfo[["blind"]])){
    warning("Warning: estimateDispersions only works, and is only needed, when you're using a CellDataSet with a negbinomial or negbinomial.size expression family")
    stop("Error: no dispersion model found. Please call estimateDispersions() before calling this function")
  }
  
  #if(!(('negbinomial()' == cds@expressionFamily) || ('negbinomial.size()' == cds@expressionFamily))){
    
  #}
  disp_df<-data.frame(gene_id=cds@dispFitInfo[["blind"]]$disp_table$gene_id,
                      mean_expression=cds@dispFitInfo[["blind"]]$disp_table$mu, 
                      dispersion_fit=cds@dispFitInfo[["blind"]]$disp_func(cds@dispFitInfo[["blind"]]$disp_table$mu),
                      dispersion_empirical=cds@dispFitInfo[["blind"]]$disp_table$disp)
  return(disp_df)
}

#####
#' Detects genes above minimum threshold.
#'
#' @description Sets the global expression detection threshold to be used with this CellDataSet.
#' Counts how many cells each feature in a CellDataSet object that are detectably expressed 
#' above a minimum threshold. Also counts the number of genes above this threshold are 
#' detectable in each cell.
#'
#' @param cds the CellDataSet upon which to perform this operation
#' @param min_expr the expression threshold 
#' @return an updated CellDataSet object
#' @importFrom Biobase fData fData<- exprs pData pData<- 
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
#' @importFrom stats median
estimateSizeFactorsForSparseMatrix <- function(counts, 
                                               locfunc = median, 
                                               round_exprs=TRUE, 
                                               method="mean-geometric-mean-total"){
  #counts <- DelayedArray(counts)
  if (round_exprs)
    counts <- round(counts)
  
  if(method == 'mean-geometric-mean-total') {
    cell_total <- Matrix::colSums(counts)
    sfs <- cell_total / exp(mean(log(cell_total)))
  }else if(method == 'mean-geometric-mean-log-total') {
    cell_total <- Matrix::colSums(counts)
    sfs <- log(cell_total) / exp(mean(log(log(cell_total))))
  }
  
  sfs[is.na(sfs)] <- 1 
  sfs   
}

#' @importFrom stats median
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
#'  @importFrom stats median
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
#' 
#' @description Returns a list of classic muscle genes. Used to
#' add conveinence for loading HSMM data.
#' 
#' @export
get_classic_muscle_markers <- function(){
  c("MEF2C", "MEF2D", "MYF5", "ANPEP", "PDGFRA",
    "MYOG", "TPM1", "TPM2", "MYH2", "MYH3", "NCAM1", "TNNT1", "TNNT2", "TNNC1",
    "CDK1", "CDK2", "CCNB1", "CCNB2", "CCND1", "CCNA1", "ID1")
}

#' Build a CellDataSet from the HSMMSingleCell package
#' 
#' @description Creates a cellDataSet using the data from the
#' HSMMSingleCell package.
#' 
#' @import HSMMSingleCell
#' @importFrom utils data
#' @export
load_HSMM <- function(){
  HSMM_sample_sheet <- NA
  HSMM_gene_annotation <- NA
  HSMM_expr_matrix <- NA
  gene_short_name <- NA
  data(HSMM_expr_matrix, envir = environment())
  data(HSMM_gene_annotation, envir = environment())
  data(HSMM_sample_sheet, envir = environment())
  pd <- new("AnnotatedDataFrame", data = HSMM_sample_sheet)
  fd <- new("AnnotatedDataFrame", data = HSMM_gene_annotation)
  HSMM <- newCellDataSet(as.matrix(HSMM_expr_matrix), phenoData = pd, featureData = fd)
  HSMM <- estimateSizeFactors(HSMM)
  HSMM <- estimateSizeFactors(HSMM)
  
  HSMM
}

#' Return a CellDataSet of classic muscle genes.
#' @importFrom Biobase fData
#' @return A CellDataSet object
#' @export
load_HSMM_markers <- function(){
  gene_short_name <- NA
  HSMM <- load_HSMM()
  marker_names <- get_classic_muscle_markers()
  HSMM[row.names(subset(fData(HSMM), gene_short_name %in% marker_names)),]
}

#' A helper function to identify the root principal points
#' @param cds CellDataSet
#' @param cell_phenotype A column in the pData. This describes the characteristics of the cells in the cds
#' @param root_type A value in from the cell_phenotype column that corresponds to the predicted starting cell
#' @importFrom IRanges which.max
#' @export
get_correct_root_state <- function(cds, cell_phenotype, root_type){
  cell_ids <- which(pData(cds)[, cell_phenotype] == root_type)
  
  closest_vertex <- cds@auxOrderingData[[cds@rge_method]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])  
  root_pr_nodes <- V(cds@minSpanningTree)$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}

#' Build a CellDataSet from the data stored in inst/extdata directory.
#' @importFrom Biobase pData pData<- exprs fData
#' @export
load_lung <- function(){
  lung_phenotype_data <- NA
  lung_feature_data <- NA
  num_cells_expressed <- NA
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
  lung <- estimateDispersions(lung)

  pData(lung)$Total_mRNAs <- colSums(exprs(lung))
  lung <- detectGenes(lung, min_expr = 1)
  expressed_genes <- row.names(subset(fData(lung), num_cells_expressed >= 5))
  ordering_genes <- expressed_genes
  lung <- setOrderingFilter(lung, ordering_genes)
  
  # DDRTree based ordering:
  lung <- preprocessCDS(lung, num_dim = 5)
  lung <- reduceDimension(lung, norm_method="log", method = 'DDRTree', pseudo_expr = 1) #
  lung <- orderCells(lung, root_pr_nodes = get_correct_root_state(lung, 
                                                           cell_phenotype = 'Time', 
                                                           "E14.5"))
  lung
}

#' Principal Components Analysis
#'
#' Efficient computation of a truncated principal components analysis of a given data matrix
#' using an implicitly restarted Lanczos method from the \code{\link{irlba}} package.
#'
#' @param x a numeric or complex matrix (or data frame) which provides
#'          the data for the principal components analysis.
#' @param retx a logical value indicating whether the rotated variables should be returned.
#' @param center a logical value indicating whether the variables should be
#'          shifted to be zero centered. Alternately, a centering vector of length
#'          equal the number of columns of \code{x} can be supplied.
#' @param scale. a logical value indicating whether the variables should be
#'          scaled to have unit variance before the analysis takes place.
#'          The default is \code{FALSE} for consistency with S, but scaling is often advisable.
#'          Alternatively, a vector of length equal the number of columns of \code{x} can be supplied.
#'
#'          The value of \code{scale} determines how column scaling is performed
#'          (after centering).  If \code{scale} is a numeric vector with length
#'          equal to the number of columns of \code{x}, then each column of \code{x} is
#'          divided by the corresponding value from \code{scale}.  If \code{scale} is
#'          \code{TRUE} then scaling is done by dividing the (centered) columns of
#'          \code{x} by their standard deviations if \code{center=TRUE}, and the
#'          root mean square otherwise.  If \code{scale} is \code{FALSE}, no scaling is done.
#'          See \code{\link{scale}} for more details.
#' @param n integer number of principal component vectors to return, must be less than
#' \code{min(dim(x))}.
#' @param ... additional arguments passed to \code{\link{irlba}}.
#'
#' @return
#' A list with class "prcomp" containing the following components:
#' \itemize{
#'    \item{sdev} {the standard deviations of the principal components (i.e.,
#'          the square roots of the eigenvalues of the
#'          covariance/correlation matrix, though the calculation is
#'          actually done with the singular values of the data matrix).}
#'   \item{rotation} {the matrix of variable loadings (i.e., a matrix whose columns
#'          contain the eigenvectors).}
#'   \item {x} {if \code{retx} is \code{TRUE} the value of the rotated data (the centred
#'          (and scaled if requested) data multiplied by the \code{rotation}
#'         matrix) is returned.  Hence, \code{cov(x)} is the diagonal matrix
#'          \code{diag(sdev^2)}.}
#'   \item{center, scale} {the centering and scaling used, or \code{FALSE}.}
#' }
#'
#' @note
#' The signs of the columns of the rotation matrix are arbitrary, and
#' so may differ between different programs for PCA, and even between
#' different builds of R.
#'
#' NOTE DIFFERENCES WITH THE DEFAULT \code{\link{prcomp}} FUNCTION!
#' The \code{tol} truncation argument found in \code{prcomp} is not supported.
#' In place of the truncation tolerance in the original function, the
#' \code{prcomp_irlba}  function has the argument \code{n} explicitly giving the
#' number of principal components to return. A warning is generated if the
#' argument \code{tol} is used, which is interpreted differently between
#' the two functions.
#'
#' @examples
#' set.seed(1)
#' x  <- matrix(rnorm(200), nrow=20)
#' p1 <- prcomp_irlba(x, n=3)
#' summary(p1)
#'
#' # Compare with
#' p2 <- prcomp(x, tol=0.7)
#' summary(p2)
#'
#'
#' @seealso \code{\link{prcomp}}
#' @import Matrix
#' @importFrom stats rnorm prcomp sd var
#' @importFrom methods slotNames slot
#' @importFrom DelayedArray DelayedArray
#' @export
sparse_prcomp_irlba <- function(x, n = 3, retx = TRUE, center = TRUE, scale. = FALSE, ...)
{
  a <- names(as.list(match.call()))
  ans <- list(scale=scale.)
  if ("tol" %in% a)
    warning("The `tol` truncation argument from `prcomp` is not supported by
            `prcomp_irlba`. If specified, `tol` is passed to the `irlba` function to
            control that algorithm's convergence tolerance. See `?prcomp_irlba` for help.")
  # Try to convert to a matrix...
  #if (!is.matrix(x)) x <- as.matrix(x)
  orig_x = x
  if (class(x) != "DelayedMatrix")
    x = DelayedArray(x)
  
  args <- list(A=orig_x, nv=n)
  if (is.logical(center))
  {
    if (center) args$center <- DelayedMatrixStats::colMeans2(x)
  } else args$center <- center
  if (is.logical(scale.))
  {
    if (is.numeric(args$center))
    {
      #f <- function(i) sqrt(sum((x[, i] - args$center[i]) ^ 2) / (nrow(x) - 1L))
      #scale. <- vapply(seq(ncol(x)), f, pi, USE.NAMES=FALSE)
      scale. = sqrt(DelayedMatrixStats::colVars(x))
      if (ans$scale) ans$totalvar <- ncol(x)
      else ans$totalvar <- sum(scale. ^ 2)
    } else
    {
      if (ans$scale)
      {
        #scale. <- apply(x, 2L, function(v) sqrt(sum(v ^ 2) / max(1, length(v) - 1L)))
        #f <- function(i) sqrt(sum((x[, i] / scale.[i]) ^ 2) / (nrow(x) - 1L))
        #ans$totalvar <- sum(vapply(seq(ncol(x)), f, pi, USE.NAMES=FALSE) ^ 2)
        
        scale. = sqrt(DelayedMatrixStats::colSums2(x ^ 2) / (max(1, nrow(x) - 1L)))
        ans$totalvar = sum(sqrt(DelayedMatrixStats::colSums2(t(t(x)/scale.) ^ 2) / (nrow(x) - 1L)))
      } else
      {
        #f <- function(i) sum(x[, i] ^ 2) / (nrow(x) - 1L)
        #ans$totalvar <- sum(vapply(seq(ncol(x)), f, pi, USE.NAMES=FALSE))
        ans$totalvar = sum(DelayedMatrixStats::colSums2(x ^ 2) / (nrow(x) - 1L))
      }
    }
    if (ans$scale) args$scale <- scale.
  } else
  {
    args$scale <- scale.
    #f <- function(i) sqrt(sum((x[, i] / scale.[i]) ^ 2) / (nrow(x) - 1L))
    #ans$totalvar <- sum(vapply(seq(ncol(x)), f, pi, USE.NAMES=FALSE))
    ans$totalvar = sum(sqrt(DelayedMatrixStats::colSums2(t(t(x)/scale.) ^ 2) / (nrow(x) - 1L)))
  }
  if (!missing(...)) args <- c(args, list(...))
  
  #args$A = as(args$A, "sparseMatrix") 
  s <- do.call(irlba, args=args)
  ans$sdev <- s$d / sqrt(max(1, nrow(x) - 1))
  ans$rotation <- s$v
  colnames(ans$rotation) <- paste("PC", seq(1, ncol(ans$rotation)), sep="")
  ans$center <- args$center
  if (retx)
  {
    ans <- c(ans, list(x = sweep(s$u, 2, s$d, FUN=`*`)))
    colnames(ans$x) <- paste("PC", seq(1, ncol(ans$rotation)), sep="")
  }
  class(ans) <- c("irlba_prcomp", "prcomp")
  ans
}

#' Select cells in an RGL scene
#' 
#' This function provides an interactive session to allow the users to select multiple cells or region for downsampling analysis  
#' @param cds CellDataSet that you'd like to select cells from
#' @param show_pp_only whether to show only the principal points. You should always turn this on for large datasets to avoid computational overloading.
#' @param return_index whether to return the index or return the cds. Default to be FALSE
#' @export
selectCells <- function(cds, show_pp_only = FALSE, return_index = FALSE) {
  if(show_pp_only) {
    data_matrix <- reducedDimK(cds)
  }
  else {
    data_matrix <- reducedDimS(cds)
  }
  
  if(nrow(data_matrix) >= 3) {
    data_df <- data.frame(t(data_matrix[1:3,]))
  } else if(nrow(data_matrix) == 2) {
    data_df <- data.frame(cbind(t(data_matrix[1:2,]), 0))
  }
  
  radius <- min(diff(apply(data_matrix, 1, range))) * 0.05
  ids <- plot3d(data_df)
  id <- selectpoints3d(ids["data"], value = FALSE,
                       multiple = function(ids) {
                         spheres3d(data_df[ids[, "index"], , drop = FALSE], color = "red", 
                                   alpha = 0.3, radius = radius)
                         TRUE
                       })
  if(return_index) {
    return(id[, 'index'])
  } else { 
    pData(cds)$Marked <- FALSE
    pData(cds)$Marked[id[, 'index']] <- TRUE      
  }
  return(cds)
}

#' This function reads in a list of 10X pipeline output directories and output a Monocle cell dataset for downstream analysis
#'
#' @description Takes a list of 10X pipeline output directories and generates a cellDataSet containing all cells in these experiments. 
#' This function is originally from Andrew Hill. 
#'
#' @param pipeline_dirs Directory name or list of directory names of the top level 10X output directory for an experiment(s)
#' @param genome String with genome name specified for 10X run (such as hg19)
#' @param lowerDetectionLimit A number that signifies the minimum expression level required in order for a gene to be truly expressed in a cell
#' @param include_analysis A boolean that signifies whether or not an analysis path can be found in the 10X output and if you wish to include it
#' @return A cellDataSet object containing data from all experiments.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' cds <- tenx_to_cds(c("/net/trapnell/vol1/home/xqiu/aggregated_samples_1", 
#'                         "/net/trapnell/vol1/home/xqiu/aggregated_samples_1", 'hg19'))
#' }
#' @importFrom utils read.table read.delim
tenx_to_cds = function(pipeline_dirs, 
                       genome="hg19", 
                       include_analysis = F,
                       lowerDetectionLimit = 1) {
  
  # Outer scope variables for collecting data
  expression_matrices <- vector(mode = 'list', length = length(pipeline_dir))
  metadata_dfs <- vector(mode = 'list', length = length(pipeline_dir))
  gene_table <- NULL # will keep first gene table so can match ordering for each dataset
  current_i <- 1
  
  lapply(pipeline_dirs, function(pipeline_dir) {
    # Check initial user input
    if(!file.exists(pipeline_dir)) {
      stop(paste("Specified 10X output directory does not exist:", pipeline_dir)) 
    }
    
    # Construct paths for pipeline output files and check that they exist
    base_path = file.path(pipeline_dir, "outs", "filtered_gene_bc_matrices_mex", genome)
    
    if(!file.exists(base_path)) { 
      stop(paste("Specified genome does not appear in 10X output:", base_path)) 
    }
    
    matrix_path <- file.path(base_path, "matrix.mtx")
    genes_path <- file.path(base_path, "genes.tsv")
    barcodes_path <- file.path(base_path, "barcodes.tsv")
    
    if(include_analysis) {
      analysis_path <- file.path(pipeline_dir, "outs", "analysis")
    }
    
    if(!file.exists(matrix_path)) { 
      stop(paste("Expression matrix not found in 10X output:", matrix_path)) 
    }
    if(!file.exists(genes_path) ) { 
      stop(paste("Genes file not found in 10X output:", genes_path)) 
    }
    if(!file.exists(barcodes_path) ) { s
      top(paste("Barcodes file not found in 10X output:", barcodes_path)) 
    }
    if(include_analysis) {
      if(!file.exists(analysis_path) ) { 
        stop(paste("Analysis  path not found in 10X output:", analysis_path)) 
      }
    }
    
    # All files exist, read them in
    matrix <- Matrix::readMM(matrix_path)
    
    barcodes <- read.table(barcodes_path, header=F, as.is=T)[,1]
    
    current_gene_table <- read.delim(genes_path, stringsAsFactors = FALSE, sep = "\t", header = FALSE)
    
    if(include_analysis) {
      tsne <- read.delim(file.path(analysis_path, "tsne", "/2_components/projection.csv"), sep=',')
    }
    
    if(is.null(gene_table)) { 
      gene_table <<- current_gene_table 
    } ## store the first gene table so can match ordering between all experiments
    
    genes <- current_gene_table[, 1]
    
    # Add gene and sample names to expression matrix (adding dataset post-fix in case barcodes appear in multiple samples)
    sample <- basename(pipeline_dir)
    row.names(matrix) <- genes
    colnames(matrix) <- paste(barcodes, sample, sep="_") ## adds dataset post-fix
    matrix <- matrix[gene_table[, 1], ] ## ensures order of genes matches between experiments
    
    # Construct metadata table that includes directory samples come from and other stats
    total_umis <- Matrix::colSums(matrix)
    
    if(include_analysis) {
      metadata_df <- data.frame(
        cell = colnames(matrix),
        total_umis = total_umis,
        sample = sample,
        TSNE.1 = tsne$TSNE.1,
        TSNE.2 = tsne$TSNE.2
      )
    } else {
      metadata_df <- data.frame(
        cell = colnames(matrix),
        total_umis = total_umis,
        sample = sample)      
    }
    
    # Add both matrices to the running list
    expression_matrices[[current_i]] <<- matrix
    metadata_dfs[[current_i]] <<- metadata_df
    
    current_i <- current_i + 1
  })
  
  # Now combine all the dataframes into one and make CDS
  combined_expression_matrix <- do.call(Matrix::cBind, expression_matrices)
  row.names(combined_expression_matrix) <- gene_table[, 1]
  
  combined_metadata_df <- do.call(rbind, metadata_dfs)
  row.names(combined_metadata_df) <- combined_metadata_df$cell
  
  colnames(gene_table) <- c("id", "gene_short_name")
  row.names(gene_table) <- gene_table$id
  
  pd <- new("AnnotatedDataFrame", data = combined_metadata_df)
  fd <- new("AnnotatedDataFrame", data = gene_table)
  cds <- newCellDataSet(combined_expression_matrix,
                        phenoData=pd,
                        featureData=fd,
                        expressionFamily=negbinomial.size(),
                        lowerDetectionLimit=lowerDetectionLimit)
  return(cds)
}

#' Subset a cds which only includes cells provided with the argument cells
#'
#' @param cds a cell dataset after trajectory reconstruction
#' @param cells a vector contains all the cells you want to subset
#' @return a new cds containing only the cells from the cells argument
#' @importFrom igraph graph.adjacency
#' @importFrom igraph induced_subgraph
#' @export
#' @examples 
#' \dontrun{
#' lung <- load_lung()
#' tmp <- subset_cds(lung, cells = row.names(subset(pData(lung), State == 1)))
#' plot_cell_trajectory(tmp)
#' }

subsetCDS <- function(cds, cells, principal_nodes = NULL){
  # Note that this function doesn't really `subset` every element in the cds as some of them are not subsettable (for example, the kNN graph, etc.)
  cells <- unique(intersect(cells, colnames(cds)))
  if(length(cells) == 0) {
    stop("Cannot find any cell from the cds matches with the cell name from the cells argument! Please make sure the cell name you input is correct.")
  }
  
  if(is.null(principal_nodes)) {
    corresponding_cells <- which(paste0('Y_', cds@auxOrderingData[[cds@rge_method]]$pr_graph_cell_proj_closest_vertex) %in% principal_nodes)
    cells_tmp <- row.names(cds@auxOrderingData[[cds@rge_method]]$pr_graph_cell_proj_closest_vertex)[corresponding_cells]

    cells <- c(cells, cells_tmp)
  }

  exprs_mat <- exprs(cds[, cells])
  cds_subset <- newCellDataSet(exprs_mat,
                                 phenoData = new("AnnotatedDataFrame", data = pData(cds)[cells, ]),
                                 featureData = new("AnnotatedDataFrame", data = fData(cds)),
                                 lowerDetectionLimit=cds@lowerDetectionLimit, 
                                 expressionFamily=cds@expressionFamily)
  pData(cds_subset)$Size_Factor <- pData(cds)[cells, 'Size_Factor']
  cds_subset@dispFitInfo <- cds@dispFitInfo
  
  if(ncol(cds@reducedDimS) == ncol(cds)) {
    cds_subset@reducedDimS <- cds@reducedDimS[, cells, drop = F]
  } else {
    cds_subset@reducedDimS <- cds@reducedDimS
  }
  if(ncol(cds@reducedDimW) == ncol(cds)) {
    cds_subset@reducedDimW <- cds@reducedDimW[, cells, drop = F]
  } else {
    cds_subset@reducedDimW <- cds@reducedDimW
  }
  if(ncol(cds@reducedDimA) == ncol(cds)) {
    cds_subset@reducedDimA <- cds@reducedDimA[, cells, drop = F]
  } else {
    cds_subset@reducedDimA <- cds@reducedDimA
  }

  if(nrow(cds@normalized_data_projection) == ncol(cds)) {
    cds_subset@normalized_data_projection <- cds@normalized_data_projection[cells, , drop = F]
  } else {
    cds_subset@normalized_data_projection <- cds@normalized_data_projection
  }
  
  cds_subset@auxOrderingData <- new.env( hash=TRUE ) # cds@auxOrderingData
  # we may also subset results from any RGE methods 
  if('DDRTree' %in% names(cds@auxOrderingData)) {
    principal_node_ids <- c(principal_nodes, paste0('Y_', cds@auxOrderingData$DDRTree$pr_graph_cell_proj_closest_vertex[cells, 1]))
    #cds_subset@auxOrderingData$DDRTree$stree <- cds_subset@auxOrderingData$DDRTree$stree[principal_node_ids, principal_node_ids]
    #cds_subset@auxOrderingData$DDRTree$R <- cds_subset@auxOrderingData$DDRTree$R[cells, principal_node_ids, drop = F]
    cds_subset@auxOrderingData$DDRTree$pr_graph_cell_proj_closest_vertex <- cds@auxOrderingData$DDRTree$pr_graph_cell_proj_closest_vertex[cells, , drop = F]
  }
  if('SimplePPT' %in% names(cds@auxOrderingData)) {
    principal_node_ids <- c(principal_nodes, paste0('Y_', cds@auxOrderingData$SimplePPT$pr_graph_cell_proj_closest_vertex[cells, 1]))
    #cds_subset@auxOrderingData$SimplePPT$stree <- cds_subset@auxOrderingData$SimplePPT$stree[principal_node_ids, principal_node_ids]
    #cds_subset@auxOrderingData$SimplePPT$R <- cds_subset@auxOrderingData$SimplePPT$R[principal_node_ids, principal_node_ids, drop = F]
    cds_subset@auxOrderingData$SimplePPT$pr_graph_cell_proj_closest_vertex <- cds@auxOrderingData$SimplePPT$pr_graph_cell_proj_closest_vertex[cells, , drop = F]
  }
  if('L1graph' %in% names(cds@auxOrderingData)) {
    principal_node_ids <- c(principal_nodes, paste0('Y_', cds@auxOrderingData$L1graph$pr_graph_cell_proj_closest_vertex[cells, 1]))
    #cds_subset@auxOrderingData$L1graph$stree <- cds_subset@auxOrderingData$L1graph$stree[principal_node_ids, principal_node_ids]
    #cds_subset@auxOrderingData$L1graph$R <- cds_subset@auxOrderingData$L1graph$R[cells, principal_node_ids, drop = F]
    cds_subset@auxOrderingData$L1graph$pr_graph_cell_proj_closest_vertex <- cds@auxOrderingData$L1graph$pr_graph_cell_proj_closest_vertex[cells, , drop = F]
  }
  
  if(ncol(cds@auxOrderingData$normalize_expr_data) == ncol(cds)) {
    cds_subset@auxOrderingData$normalize_expr_data <- cds@auxOrderingData$normalize_expr_data[, cells]
  }

  cds_subset@minSpanningTree <- cds@minSpanningTree
  # find the corresponding principal graph nodes for those selected cells, followed by subseting the trajectories 
  if(length(cds@rge_method) > 0) {
    principal_graph_points <- c(principal_nodes, paste0('Y_', cds_subset@auxOrderingData[[cds@rge_method]]$pr_graph_cell_proj_closest_vertex[cells, 1, drop = F]))
    tmp <- induced_subgraph(cds@minSpanningTree, principal_graph_points)
    cds_subset@minSpanningTree <- tmp
    V(cds_subset@minSpanningTree)$name <- paste0('Y_', 1:vcount(cds_subset@minSpanningTree))
  }
  
  # find the corresponding principal graph nodes for those selected cells, followed by subseting the trajectories 
  cds_subset@auxClusteringData <- new.env( hash=TRUE ) # cds@auxClusteringData
  # browser()
  if('louvain' %in% names(cds@auxClusteringData)) {
    cds_subset@auxClusteringData$louvain$louvain_res$g <- induced_subgraph(cds@auxClusteringData$louvain$louvain_res$g, cells)
    cds_subset@auxClusteringData$louvain$louvain_res$relations <- subset(cds@auxClusteringData$louvain$louvain_res$relations, from %in% cells & to %in% cells)
    cds_subset@auxClusteringData$louvain$louvain_res$distMatrix <- cds@auxClusteringData$louvain$louvain_res$distMatrix[match(cells, colnames(cds)), , drop = F]
    cds_subset@auxClusteringData$louvain$louvain_res$optim_res$membership <- cds@auxClusteringData$louvain$louvain_res$optim_res$membership[match(cells, colnames(cds))]
  } else if('partitionCells' %in% names(cds@auxClusteringData)) {
    cds_subset@auxClusteringData$partitionCells$g <- induced_subgraph(cds@auxClusteringData$partitionCells$g, cells)
    cds_subset@auxClusteringData$partitionCells$relations <- subset(cds@auxClusteringData$partitionCells$relations, from %in% cells & to %in% cells)
    cds_subset@auxClusteringData$partitionCells$distMatrix <- cds@auxClusteringData$partitionCells$distMatrix[match(cells, colnames(cds)), , drop = F]
    cds_subset@auxClusteringData$partitionCells$optim_res$membership <- cds@auxClusteringData$partitionCells$optim_res$membership[match(cells, colnames(cds))]
  }
  
  cds_subset@reducedDimK <- cds@reducedDimK[, V(tmp)$name, drop = F]
  colnames(cds_subset@reducedDimK) <- paste0('Y_', 1:ncol(cds_subset@reducedDimK))
  cds_subset@dim_reduce_type <- cds@dim_reduce_type
  cds_subset@rge_method <- cds@rge_method
  cds_subset 
}

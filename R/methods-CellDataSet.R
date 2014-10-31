#use the deconvoluated linear regression parameters to normalize the log_fpkm
norm_kb <- function(kb, exprs_cds) {
    k <- kb[1]
    b <- kb[2]
    tmp <- k * log10(exprs_cds) + b
    norm_exprs <- 10^tmp
    
    norm_exprs
}

dmode <- function(x, breaks="Sturges") {
    den <- density(x, kernel=c("gaussian"))
    ( den$x[den$y==max(den$y)] )
}

#' Find the most abundant FPKM estimator of single cell FPKM matrix
#' @param fpkm_matrix a matrix of fpkm for single cell RNA-seq values with each row and column representing genes/isoforms and cells. Row and column names should be included
#' @param return_all parameter for the intended return results. If setting TRUE, matrix of dmode as well as max mu and min mu of two gaussian distribution mixture will be returned
#' @param fpkm_thresh the threshold for FPKM value used in detecting the most abundant FPKM value
#' @return an vector of most abundant fpkm value corresponding to the transcrpt copy 1. If setting return_all = TRUE, the mode based on gaussian density function and the max or min
#' mode from the mixture gaussian model
#' @detail This function estimates the most abundant FPKM (t^*) using gaussian kernel density function. It can also optionally output the t^* based on a two gaussian mixture model
#' based on the smsn.mixture from mixsmsn package
#' @export
#' @example
#' \dontrun{
#' HSMM_fpkm_matrix <- exprs(HSMM)
#' t_estimate = estimate_t(HSMM_fpkm_matrix)
#'}
estimate_t <- function(fpkm_matrix, return_all = F, fpkm_thresh = 0.1, ...) {
    #peak finder (using mixture Gauissan model for the FPKM distribution fitting, choose the minial peaker as default)
    
    smsn_mode_test <- function(fpkm, g = 2, fpkm_thresh = 0.1) {
        log_fpkms <- log10(fpkm[fpkm > fpkm_thresh]) #only look the transcripts above a certain threshold
        sm_2 <- smsn.mix(log_fpkms, nu = 3, g = 2, get.init = TRUE, criteria = TRUE, iter.max=1000,calc.im = FALSE, family="Normal")
        #   print (sm_2)
        
        sm_1 <- smsn.mix(log_fpkms, nu = 3, g = 1, get.init = TRUE, criteria = TRUE, iter.max=1000, calc.im = FALSE, family="Normal")
        #   print (sm_1)
        
        if (sm_1$aic >= sm_2$aic){
            trick_location_max <- 10^max(sm_2$mu)
            trick_location_min <- 10^min(sm_2$mu) #use the min / max
        }else{
            trick_location_max <- 10^max(sm_1$mu)
            trick_location_min <- 10^min(sm_1$mu)
        }
        
        #best_cov <- 10^sm_analysis$mu[best_location]
        best_cov_dmode <- 10^(dmode(log_fpkms))
        
        best_cov_max <- trick_location_max
        best_cov_min <- trick_location_min
        sd_fpkm <- 0
        #   print (c(best_cov_max, best_cov_min, best_cov_dmode))
        data.frame(best_cov_dmode = best_cov_dmode,
        best_cov_max = best_cov_max,
        best_cov_min = best_cov_min
        )
    }
    
    #apply each column
    if(return_all){
        do.call(rbind, apply(fpkm_matrix, 2, function(fpkm) smsn_mode_test(fpkm, ...)))
    }
    else{
        apply(fpkm_matrix, 2, function(fpkm) 10^dmode(log10(fpkm[fpkm > fpkm_thresh]))) #best coverage estimate}
    }
}

#linear transform by t* and m, c
opt_norm_t <- function(t, fpkm, m, c, return_norm = FALSE) {
    a_matrix <- matrix(c(log10(t), 1, m,
    -1), ncol = 2, nrow = 2, byrow = T)
    colnames(a_matrix) <- c("k", "b")
    b_matrix <- matrix(c(0, -c), nrow = 2, byrow = T)
    kb <- t(solve(a_matrix, b_matrix))
    
    k <- kb[1]
    b <- kb[2]
    tmp <- k * log10(fpkm) + b
    abs_cnt <- 10^tmp
    
    if(return_norm) return(abs_cnt)
    10^dmode(log10(abs_cnt[abs_cnt > 0]))
}

#linear transform by kb
opt_norm_kb <- function(fpkm, kb) {
    
    k <- kb[1]
    b <- kb[2]
    tmp <- k * log10(fpkm) + b
    abs_cnt <- 10^tmp
    
    10^dmode(log10(abs_cnt[abs_cnt > 0]))
}

#rmse between the dmode from t estimate based linear transformation and the spike-dmode
t_rmse_abs_cnt <- function (par, t_estimate, fpkm_mat, split_fpkm, alpha = 1, cores = 1, ...) {
  cell_num <- ncol(fpkm_mat)
  #t_estimate <- par[1:cell_num] #t*: the estimates for the best coverage
  names(t_estimate) <- colnames(fpkm_mat)
  split_t <- split(t(t_estimate), col(as.matrix(t(t_estimate)), as.factor = T))
  print(paste("t_estimate is: ", paste(as.character(t_estimate), sep = '', collapse = ' '), sep = '', collapse = ''))
  
  #mc_guess <- par[(cell_num + 1):(cell_num + 2)] #m, c parameters: b = m k + c
  mc_guess <- par
  print(paste("mc_guess is", mc_guess[1], mc_guess[2], sep = ' '))
  
  m_val <- mc_guess[1]
  c_val <- mc_guess[2]
  cell_dmode <- tryCatch({
    if(cores > 1){
      cell_dmode <- mcmapply(opt_norm_t, split_t, split_fpkm, m = m_val, c = c_val, mc.cores = cores)
      adj_est_std_cds <- mcmapply(opt_norm_t, split_t, split_fpkm, m = m_val, c = c_val, return_norm = T, mc.cores = cores)
    }
    else {
      cell_dmode <- mapply(opt_norm_t, split_t, split_fpkm, m = m_val, c = c_val)
      adj_est_std_cds <- mapply(opt_norm_t, split_t, split_fpkm, m = m_val, c = c_val, return_norm = T)
    }
    cell_dmode},
    error = function(e) {print(e); t_estimate} #return what is better?
  )
  print(paste("cell_dmode is: ", paste(as.character(cell_dmode), sep = '', collapse = ' '), sep = '', collapse = ''))
  
  rmse <- sqrt(mean((cell_dmode - alpha)^2)) #rmse between the estimated cell_dmode and the 1 copy of transcripts
  print(paste('rmse is:', rmse, sep = ' '))
  
  sum_total_cells_rna <- colSums(adj_est_std_cds)
  #sum_total_cells_rna[is.infinite(sum_total_cells_rna)] <- 7861584 * 2
  print(paste('sum of all total RNA is', sum_total_cells_rna))
  sqrt(mean(((cell_dmode - alpha) / sum_total_cells_rna)^2)) #new 
}

#' Transform FPKM matrix to absolute transcript matrix based on the decomposed linear regression parameters from most abundant isoform fpkm value.
#'
#' @param fpkm_matrix an matrix of fpkm for single cell RNA-seq values with each row and column representing genes/isoforms and cells. Row and column names should be included
#' @param t_estimate an vector for the estimated most abundant FPKM value of isoform for a single cell. Estimators based on gene FPKM can also give good approximation but estimators
#' based on isoform FPKM will give better results in general
#' @param global_scaling parameter for globaling scaling. Not perform global scaling by default (global_scaling = 1)
#' @param return_all parameter for the intended return results. If setting TRUE, matrix of k^*, b^* and vector of global_scaling as well the transformed absolute cds will be returned
#' in a list format
#' @param total_fragment: vector of total fragment sequenced for each cell, the element should matched with fpkm_matrix columns
#' @param alpha_v: vector of most abundant transcripts number estimates in each cell, which are dependent on the total_fragment
#' @return an matrix of absolute count for isoforms or genes after the transformation. For more details on other output, please refer to detail
#' @detail This function takes a FPKM matrix and a vector of estimated most abundant FPKM value from the isoform FPKM matrix and transform it into absolute transcript number.
#' It is based on the fact that most isoforms of gene in a single cell only express one copy so that the most abundant isoform FPKM (t^*) will corresponding to 1 copy transcript. The
#' function takes the the vector t^* and then decomposes it into two parameters vectors k^* and b^* (corresponding to each cell) which correspond to the slope and intercept when
#' we perform the robust linear regression for the spikein data. This decomposition is based on an observed relationship between k and b in terms of b = -3.652201 k + 2.263576. The
#' function will then apply a linear transformation to convert the FPKM to estimated absolute transcript counts based on the the k^* and b^*. The function can also apply a global
#' adjustment if setting global_scaling = TRUE. The k*, b* parameters vectors and the global scaling factor can be output in a list format (together with norm_cds) if setting return_all
#' == TRUE
#' @export
#' @example
#' \dontrun{
#' HSMM_fpkm_matrix <- exprs(HSMM)
#' HSMM_abs_matrix <- fpkm2abs(HSMM_fpkm_matrix, t_estimate = estimate_t(fpkm_matrix))
#'}
#use the deconvoluated linear regression parameters to normalize the log_fpkm

#accept an log fpkm matrix and transform it into the absolute transcript matrix
#row and column names of the matrix are genes and cells respectively
fpkm2abs <- function(fpkm_matrix, t_estimate = estimate_t(fpkm_matrix), total_fragment = 1.5e6, alpha_v = 1, m = -3.652201, c = 2.263576, global_scaling = FALSE, return_all = FALSE, num_cores = 1) {
    #alpha_v <- exp(log10(1.5 * 10e6  / sample_sheet$Mapped.Fragments)) #we can estimate alpha_v by certain function
    print('optimizing t_estimates and m and c...') #silence the optimization output

    split_fpkms <- split(as.matrix(fpkm_matrix), col(fpkm_matrix, as.factor = T)) #ensure the split dataset is matrix

    optim_para <- optim(par = c(m, c), 
                    t_rmse_abs_cnt, gr = NULL, t_estimate = t_estimate, alpha = alpha_v, cores = num_cores, fpkm_mat = fpkm_matrix, split_fpkm = split_fpkms,
                    method = c("L-BFGS-B"),
                    lower = c(rep(as.vector(t_estimate) - 0, 0), -10, 0.1), #search half low or up of the t_estimate
                    upper = c(rep(as.vector(t_estimate) + 0, 0), -0.1, 10), #m, c is between (-0.1 to -10 and 0.1 to 10)
                    control = list(factr = 1e12, pgtol = 1e-3, trace = 1, ndeps = c(1e-3, 1e-3) ), #as.vector(t_estimate) / 1000,
                    hessian = FALSE)
    print('optimization is done!')
    
    #t_estimate <- optim_para$par[1:length(t_estimate)]
    #regression line between b^* = m * k^* + c
    #m <- optim_para$par[length(t_estimate) + 1]
    #c <- optim_para$par[length(t_estimate) + 2]
    m <- optim_para$par[1]
    c <- optim_para$par[2]
    
    #estimate the t^* by smsn two gaussian model, choose minial peak  1
    names(t_estimate) <- colnames(fpkm_matrix)
    
    total_rna_df <- data.frame(Cell = colnames(fpkm_matrix), t_estimate = t_estimate)
    
    #solve k and b for t by matrix formulation (B = Ax)
    k_b_solution <- ddply(total_rna_df, .(Cell), function(x){
        a_matrix <- matrix(c(log10(x[, "t_estimate"]), 1, m, -1), ncol = 2, nrow = 2, byrow = T)
        
        colnames(a_matrix) <- c("k", "b")
        b_matrix <- matrix(c(0, -c), nrow = 2, byrow = T)
        k_b_solution <- t(solve(a_matrix, b_matrix))
    })
    
    print(k_b_solution[1:5, ])
    #predict the absolute isoform based on the adjustment
    
    rownames(k_b_solution) <- k_b_solution$Cell
    k_b_solution <- t(k_b_solution[, c(2, 3)]) #ddply give Cell, k, b columns, take the last two
    split_kb <- split(k_b_solution, col(k_b_solution, as.factor =  T))
    
    adj_split_fpkm <- mcmapply(norm_kb, split_kb, split_fpkms, mc.cores = num_cores)
    
    #confirm our adjustment
    #genes only have one isoform
    #qplot(x = adj_est_iso_cds[1,], y = exprs(isoform_trick_exprs)[1, ])
    
    total_rna_df$estimate_k <- k_b_solution[1, ]
    total_rna_df$estimate_b <- k_b_solution[2, ]
    
    #global_scaling seems unneccessary, remove this part in future
    #global normalization
    if(global_scaling == TRUE) {
        median_intercept <- median(total_rna_df$estimate_b)
        scaling_factor <- 10^(median_intercept - total_rna_df$estimate_b) #ensure the cells with intercept close to median_intercept don't scale
        split_scaling_factor <- split(scaling_factor, factor(1:length(scaling_factor)))
        split_gene_cds <- split(adj_split_fpkm, col(adj_split_fpkm, as.factor = T))
        norm_cds <- mcmapply(function(x, y) x *y, split_gene_cds, split_scaling_factor, mc.cores = num_cores)
    }
    else {
        scaling_factor = 1
        norm_cds <- adj_split_fpkm
    }
    
    row.names(norm_cds) <- row.names(fpkm_matrix)
    colnames(norm_cds) <- colnames(fpkm_matrix)
    
    if(return_all == T) { #also return the trick cds for genes and isoform if return_trick_cds is true otherwise only return total_rna_df
        return (list(norm_cds = norm_cds, k_b_solution = k_b_solution, scaling_factor = scaling_factor, optim_para = optim_para))
    }
    norm_cds
}

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
  cellData <- as.matrix( cellData )
  
  if( is.null( phenoData ) )
    phenoData <- annotatedDataFrameFrom( cellData, byrow=FALSE )
  if( is.null( featureData ) ) 
    featureData <- annotatedDataFrameFrom( cellData, byrow=TRUE )
  
  cds <- new( "CellDataSet",
              assayData = assayDataNew( "environment", exprs=cellData ),
              phenoData=phenoData, 
              featureData=featureData, 
              lowerDetectionLimit=lowerDetectionLimit,
              expressionFamily=expressionFamily )
  
  validObject( cds )
  cds
}

#######

setValidity( "CellDataSet", function( object ) {
#   if( any( counts(object) < 0 ) )
#     return( "the count data contains negative values" )
  TRUE
} )

#' Retrieves the coordinates of each cell in the reduced-dimensionality space generated by calls to 
#' reduceDimension.
#'
#' @param cds A CellDataSet object.
#' @return A matrix, where rows are cell coordinates and columns correspond to dimensions of the 
#' reduced space.
#' @export
#' @examples
#' data(HSMM)
#' S <- reducedDimS(HSMM)
reducedDimS <- function( cds ) {
  stopifnot( is( cds, "CellDataSet" ) )
  cds@reducedDimS
}   

#' Sets the coordinates of each cell in the reduced-dimensionality space.  Not intended to be called directly.
#'
#' @param cds A CellDataSet object.
#' @param value A matrix of coordinates specifying each cell's position in the reduced-dimensionality space.
#' @return An update CellDataSet object
#' @export
`reducedDimS<-` <- function( cds, value ) {
  stopifnot( is( cds, "CellDataSet" ) )
  cds@reducedDimS <- value
  validObject( cds )
  cds
}   

#' Retrieves the expression values for each cell (as a matrix) after whitening during independent component analysis.
#'
#' @param cds A CellDataSet object.
#' @return A matrix, where each row is a set of whitened expression values for a feature and columns are cells.
#' @export
#' @examples
#' data(HSMM)
#' W <- reducedDimW(HSMM)
reducedDimW <- function( cds ) {
  stopifnot( is( cds, "CellDataSet" ) )
  cds@reducedDimW
}   

#' Sets the whitened expression values for each cell prior to independent component analysis. Not intended to be called directly.
#'
#' @param cds A CellDataSet object.
#' @param value A whitened expression data matrix
#' @return An updated CellDataSet object
#' @export
`reducedDimW<-` <- function( cds, value ) {
  stopifnot( is( cds, "CellDataSet" ) )
  cds@reducedDimW <- value
  validObject( cds )
  cds
}   

#' Retrieves the weights that transform the cells' coordinates in the reduced dimension space back to the full (whitened) space.
#'
#' @param cds A CellDataSet object.
#' @return A matrix that when multiplied by a reduced-dimension set of coordinates for the CellDataSet, 
#' recovers a matrix in the full (whitened) space
#' @export
#' @examples
#' data(HSMM)
#' A <- reducedDimA(HSMM)
reducedDimA <- function( cds ) {
  stopifnot( is( cds, "CellDataSet" ) )
  cds@reducedDimA
}   

#' Sets the weights transform the cells' coordinates in the reduced dimension space back to the full (whitened) space.
#'
#' @param cds A CellDataSet object.
#' @param value A whitened expression data matrix
#' @return An updated CellDataSet object
#' @export
`reducedDimA<-` <- function( cds, value ) {
  stopifnot( is( cds, "CellDataSet" ) )
  cds@reducedDimA <- value
  validObject( cds )
  cds
}   

#' Retrieves the variance-stabilized expression matrix for the CellDataSet.
#'
#' @param cds A CellDataSet object.
#' @return A matrix of variance stabilized expression values for the CellDataSet
#' @export
#' @examples
#' data(HSMM)
#' A <- reducedDimA(HSMM)
vstExprs <- function( cds ) {
  stopifnot( is( cds, "CellDataSet" ) )
  cds@vstExprs[row.names(fData(cds)), row.names(pData(cds))]
}   

#' Sets the variance-stabilized expression matrix for the CellDataSet.
#'
#' @param cds A CellDataSet object.
#' @param value A variance-stabilized expression data matrix
#' @return An updated CellDataSet object
#' @export
`vstExprs<-` <- function( cds, value ) {
  stopifnot( is( cds, "CellDataSet" ) )
  
  #TODO: check row and column consistency with exprs(cds)
  
  cds@vstExprs <- value
  validObject( cds )
  cds
}   

#' Retrieves the size factors for the CellDataSet.
#'
#' @param cds A CellDataSet object.
#' @return A vector of length equal to the number of cells in the CellDataSet
#' @export
#' @examples
#' data(HSMM)
#' A <- reducedDimA(HSMM)
sizeFactors <- function( cds ) {
  stopifnot( is( cds, "CellDataSet" ) )
  pData(cds)$Size_Factor
}   

#' Sets the size factors for the CellDataSet.
#'
#' @param cds A CellDataSet object.
#' @param value A vector of length equal to the number of cells in the CellDataSet
#' @return An updated CellDataSet object
#' @export
`sizeFactors<-` <- function( cds, value ) {
  stopifnot( is( cds, "CellDataSet" ) )
  
  #TODO: check column consistency with exprs(cds)
  
  pData(cds)$Size_Factor <- value
  validObject( cds )
  cds
}   

#' Retrieves the minimum spanning tree generated by Monocle during cell ordering.
#'
#' @param cds expression data matrix for an experiment
#' @return An igraph object representing the CellDataSet's minimum spanning tree.
#' @export
#' @examples
#' data(HSMM)
#' T <- minSpanningTree(HSMM)
minSpanningTree <- function( cds ) {
  stopifnot( is( cds, "CellDataSet" ) )
  cds@minSpanningTree
}   

#' Sets the minimum spanning tree used by Monocle during cell ordering. Not intended to be called directly.
#'
#' @param cds A CellDataSet object.
#' @param value an igraph object describing the minimum spanning tree.
#' @return An updated CellDataSet object
#' @export
`minSpanningTree<-` <- function( cds, value ) {
  stopifnot( is( cds, "CellDataSet" ) )
  cds@minSpanningTree <- value
  validObject( cds )
  cds
}   

#' Retrieves a matrix capturing distances between each cell in the reduced-dimensionality space
#'
#' @param cds expression data matrix for an experiment
#' @return A square, symmetric matrix containing the distances between each cell in the reduced-dimensionality space.
#' @export
#' @examples
#' data(HSMM)
#' D <- cellPairwiseDistances(HSMM)
cellPairwiseDistances <- function( cds ) {
  stopifnot( is( cds, "CellDataSet" ) )
  cds@cellPairwiseDistances
}   

#' Sets the matrix containing distances between each pair of cells used by Monocle during cell ordering. Not intended to be called directly.
#'
#' @param cds A CellDataSet object.
#' @param value a square, symmetric matrix containing pairwise distances between cells.
#' @return An updated CellDataSet object
#' @export
`cellPairwiseDistances<-` <- function( cds, value ) {
  stopifnot( is( cds, "CellDataSet" ) )
  cds@cellPairwiseDistances <- value
  validObject( cds )
  cds
}   

#####

# Methods for PQ-tree based ordering

get_next_node_id <- function()
{
  next_node <<- next_node + 1
  return (next_node) 
}

# Recursively builds and returns a PQ tree for the MST
pq_helper<-function(mst, use_weights=TRUE, root_node=NULL)
{
  new_subtree <- graph.empty()
  
  root_node_id <- paste("Q_", get_next_node_id(), sep="")
  
  new_subtree <- new_subtree + vertex(root_node_id, type="Q", color="black")
  
  if (is.null(root_node) == FALSE){
    sp <- get.all.shortest.paths(mst, from=V(mst)[root_node])
    #print(sp)
    sp_lengths <- sapply(sp$res, length)
    target_node_idx <- which(sp_lengths == max(sp_lengths))[1]
    #print(unlist(sp$res[target_node_idx]))
    diam <- V(mst)[unlist(sp$res[target_node_idx])]
    #print(diam)
  }else{
    if (use_weights){
      diam <- V(mst)[get.diameter(mst)]
    }else{
      diam <- V(mst)[get.diameter(mst, weights=NA)]
    }
  }
  
  
  #print (diam)
  
  V(new_subtree)[root_node_id]$diam_path_len = length(diam)
  
  diam_decisiveness <- igraph::degree(mst, v=diam) > 2
  ind_nodes <- diam_decisiveness[diam_decisiveness == TRUE]
  
  first_diam_path_node_idx <- head(as.vector(diam), n=1)
  last_diam_path_node_idx <- tail(as.vector(diam), n=1)
    if (sum(ind_nodes) == 0 || 
          (igraph::degree(mst, first_diam_path_node_idx) == 1 && 
             igraph::degree(mst, last_diam_path_node_idx) == 1))
    {
      ind_backbone <- diam
    }
    else 
    {
      last_bb_point <- names(tail(ind_nodes, n=1))[[1]]
      first_bb_point <- names(head(ind_nodes, n=1))[[1]]	
      #diam_path_vertex_names <- as.vector()
      #print (last_bb_point)
      #print (first_bb_point)
      diam_path_names <- V(mst)[as.vector(diam)]$name
      last_bb_point_idx <- which(diam_path_names == last_bb_point)[1]
      first_bb_point_idx <- which(diam_path_names == first_bb_point)[1]
      ind_backbone_idxs <- as.vector(diam)[first_bb_point_idx:last_bb_point_idx]
      #print (ind_backbone_idxs)
      ind_backbone <- V(mst)[ind_backbone_idxs]
      
      #ind_backbone <- diam[first_bb_point:last_bb_point]
    }
  
  
  
  mst_no_backbone <- mst - ind_backbone
  #print (V(mst_no_backbone)$name)
  
  for (backbone_n in ind_backbone)
  {
    #print (n)
    #backbone_n <- ind_backbone[[i]]
    
    if (igraph::degree(mst, v=backbone_n) > 2)
    {
      new_p_id <- paste("P_", get_next_node_id(), sep="")
      #print(new_p_id)
      new_subtree <- new_subtree + vertex(new_p_id, type="P", color="grey")
      new_subtree <- new_subtree + vertex(V(mst)[backbone_n]$name, type="leaf", color="white")
      new_subtree <- new_subtree + edge(new_p_id, V(mst)[backbone_n]$name)
      new_subtree <- new_subtree + edge(root_node_id, new_p_id)
      
      nb <- graph.neighborhood(mst, 1, nodes=backbone_n)[[1]]
      
      #print (E(nb))
      #print (V(nb))
      for (n_i in V(nb))
      {
        n <- V(nb)[n_i]$name			
        if (n %in% V(mst_no_backbone)$name)
        {	
          #print (n)
          
          sc <- subcomponent(mst_no_backbone, n)
          
          sg <- induced.subgraph(mst_no_backbone, sc, impl="copy_and_delete")
          
          
          if (ecount(sg) > 0)
          {
            #print (E(sg))	
            sub_pq <- pq_helper(sg, use_weights)
            
            
            # Works, but slow:
            for (v in V(sub_pq$subtree))
            {
              new_subtree <- new_subtree + vertex(V(sub_pq$subtree)[v]$name, type=V(sub_pq$subtree)[v]$type, color=V(sub_pq$subtree)[v]$color, diam_path_len=V(sub_pq$subtree)[v]$diam_path_len)
            }
            
            edge_list <- get.edgelist(sub_pq$subtree)
            for (i in 1:nrow(edge_list))
            {
              new_subtree <- new_subtree + edge(V(sub_pq$subtree)[edge_list[i, 1]]$name, V(sub_pq$subtree)[edge_list[i, 2]]$name)
            }   					
            #plot (new_subtree)
            
            new_subtree <- new_subtree + edge(new_p_id, V(sub_pq$subtree)[sub_pq$root]$name)  
          }
          else
          {
            new_subtree <- new_subtree + vertex(n, type="leaf", color="white")
            new_subtree <- new_subtree + edge(new_p_id, n)
          }
        }
        
      }
      #print ("##########################")
    }
    else
    {
      new_subtree <- new_subtree + vertex(V(mst)[backbone_n]$name, type="leaf", color="white")
      new_subtree <- new_subtree + edge(root_node_id, V(mst)[backbone_n]$name)
    }
  }
  # else
  # {
  #     for (backbone_n in diam)
  #     {
  #           new_subtree <- new_subtree + vertex(backbone_n, type="leaf")
  #           new_subtree <- new_subtree + edge(root_node_id, backbone_n)
  #     }
  # }
  
  return (list(root=root_node_id, subtree=new_subtree))
}

make_canonical <-function(pq_tree)
{
  canonical_pq <- pq_tree
  
  V(canonical_pq)[type == "P" & igraph::degree(canonical_pq, mode="out") == 2]$color="black"
  V(canonical_pq)[type == "P" &  igraph::degree(canonical_pq, mode="out")== 2]$type="Q"
  
  single_child_p <- V(canonical_pq)[type == "P" & igraph::degree(canonical_pq, mode="out")== 1]
  V(canonical_pq)[type == "P" & igraph::degree(canonical_pq, mode="out")==1]$color="blue"
  
  for (p_node in single_child_p)
  {
    child_of_p_node <- V(canonical_pq) [ nei(p_node, mode="out") ]
    parent_of_p_node <- V(canonical_pq) [ nei(p_node, mode="in") ]
    
    for (child_of_p in child_of_p_node)
    {
      canonical_pq[parent_of_p_node, child_of_p] <- TRUE
      # print (p_node)
      # print (child_of_p)
      # print (parent_of_p_node)
      # print ("*********")
    }
  }
  
  canonical_pq <- delete.vertices(canonical_pq, V(canonical_pq)[type == "P" & igraph::degree(canonical_pq, mode="out")==1])
  #print (V(canonical_pq)[type == "Q" & igraph::degree(canonical_pq, mode="in")==0])
  return (canonical_pq)
}

extract_ordering <- function(pq_tree, curr_node)
{
  if (V(pq_tree)[curr_node]$type == "leaf")
  {
    return (V(pq_tree)[curr_node]$name)
  }
  else if (V(pq_tree)[curr_node]$type == "P")
  {
    p_level <- list()
    for (child in V(pq_tree) [ nei(curr_node, mode="out") ])
    {
      p_level[[length(p_level)+1]] <- extract_ordering(pq_tree, child)
    }
    p_level <- p_level[sample(length(p_level))]
    p_level <- unlist(p_level)
    #print (p_level)
    return (p_level)
  }
  else if(V(pq_tree)[curr_node]$type == "Q")
  {
    q_level <- list()
    for (child in V(pq_tree) [ nei(curr_node, mode="out") ])
    {
      q_level[[length(q_level)+1]] <- extract_ordering(pq_tree, child)
    }
    if (runif(1) >= 0.5)
    {
      q_level <- rev(q_level)
    }
    q_level <- unlist(q_level)
    #print (q_level)
    return (q_level)
  }
}

extract_fixed_ordering <- function(pq_tree, curr_node)
{
  if (V(pq_tree)[curr_node]$type == "leaf")
  {
    return (V(pq_tree)[curr_node]$name)
  }
  else if (V(pq_tree)[curr_node]$type == "P")
  {
    p_level <- list()
    for (child in V(pq_tree) [ nei(curr_node, mode="out") ])
    {
      p_level[[length(p_level)+1]] <- extract_ordering(pq_tree, child)
    }
    #p_level <- p_level[sample(length(p_level))]
    p_level <- unlist(p_level)
    #print (p_level)
    return (p_level)
  }
  else if(V(pq_tree)[curr_node]$type == "Q")
  {
    q_level <- list()
    for (child in V(pq_tree) [ nei(curr_node, mode="out") ])
    {
      q_level[[length(q_level)+1]] <- extract_ordering(pq_tree, child)
    }
    # if (runif(1) >= 0.5)
    # {
    #     q_level <- rev(q_level)
    # }
    q_level <- unlist(q_level)
    #print (q_level)
    return (q_level)
  }
}

weight_of_p_node_order<-function(p_level_list, dist_matrix)
{
  cost <- 0
  if (length(p_level_list) <= 1)
  {
    return (0)
  }
  #print (p_level_list)
  #print (length(p_level_list))
  for (i in (1:(length(p_level_list)-1)))
  {
    #print (i)
    #print ("***")
    #print (p_level_list[[i]])
    #print (paste("...", p_level_list[[i]][length(p_level_list[[i]])]))
    #print (p_level_list[[i+1]])
    #print (paste("...", p_level_list[[i+1]][1]))
    #print (dist_matrix[p_level_list[[i]][length(p_level_list[[i]])], p_level_list[[i+1]][1]])
    cost <- cost + dist_matrix[p_level_list[[i]][length(p_level_list[[i]])], p_level_list[[i+1]][1]]
  }
  return(cost)
}

order_p_node <- function(q_level_list, dist_matrix)
{ 
  q_order_res <- permn(q_level_list, fun=order_q_node, dist_matrix)
  #print (q_order_res)
  all_perms <- lapply(q_order_res, function(x) { x$ql } )
  #print ("perm ql:")
  #print(all_perms)
  all_perms_weights <- unlist(lapply(q_order_res, function(x) { x$wt }))
  #print ("perm weights:")
  #print (all_perms_weights)
  
  opt_perm_idx <- head((which(all_perms_weights == min(all_perms_weights))), 1)
  opt_perm <- all_perms[[opt_perm_idx]]
  
  #print ("opt_path:")
  #print (opt_perm)
  #print ("opt_all_weight:")
  #print (min(all_perms_weights))
  #print ("weights:")
  #print (all_perms_weights)
  # print ("q_level_list:")
  # print (q_level_list)
  stopifnot (length(opt_perm) == length(q_level_list))
  
  return(opt_perm)
}

order_q_node <- function(q_level_list, dist_matrix)
{
  new_subtree <- graph.empty()
  
  if (length(q_level_list) == 1)
  {
    return (list(ql=q_level_list, wt=0))
  }
  for (i in 1:length(q_level_list))
  {
    new_subtree <- new_subtree + vertex(paste(i,"F"), type="forward")
    new_subtree <- new_subtree + vertex(paste(i,"R"), type="reverse")
  }
  
  for (i in (1:(length(q_level_list)-1)))
  {
    cost <- dist_matrix[q_level_list[[i]][length(q_level_list[[i]])], q_level_list[[i+1]][1]]
    new_subtree <- new_subtree + edge(paste(i,"F"), paste(i+1,"F"), weight=cost)
    
    cost <- dist_matrix[q_level_list[[i]][length(q_level_list[[i]])], q_level_list[[i+1]][length(q_level_list[[i+1]])]]
    new_subtree <- new_subtree + edge(paste(i,"F"), paste(i+1,"R"), weight=cost)
    
    cost <- dist_matrix[q_level_list[[i]][1], q_level_list[[i+1]][1]]
    new_subtree <- new_subtree + edge(paste(i,"R"), paste(i+1,"F"), weight=cost)
    
    cost <- dist_matrix[q_level_list[[i]][1], q_level_list[[i+1]][length(q_level_list[[i+1]])]]
    new_subtree <- new_subtree + edge(paste(i,"R"), paste(i+1,"R"), weight=cost)
  }
  
  first_fwd = V(new_subtree)[paste(1,"F")]
  first_rev = V(new_subtree)[paste(1,"R")]
  last_fwd = V(new_subtree)[paste(length(q_level_list),"F")]
  last_rev = V(new_subtree)[paste(length(q_level_list),"R")]
  
  FF_path <- unlist(get.shortest.paths(new_subtree, from=as.vector(first_fwd), to=as.vector(last_fwd), mode="out", output="vpath")$vpath)
  FR_path <- unlist(get.shortest.paths(new_subtree, from=as.vector(first_fwd), to=as.vector(last_rev), mode="out", output="vpath")$vpath)
  RF_path <- unlist(get.shortest.paths(new_subtree, from=as.vector(first_rev), to=as.vector(last_fwd), mode="out", output="vpath")$vpath)
  RR_path <- unlist(get.shortest.paths(new_subtree, from=as.vector(first_rev), to=as.vector(last_rev), mode="out", output="vpath")$vpath)
  
  # print (FF_path)
  # print (FR_path)
  # print (RF_path)
  # print (RR_path)
  
  FF_weight <- sum(E(new_subtree, path=FF_path)$weight)
  FR_weight <- sum(E(new_subtree, path=FR_path)$weight)
  RF_weight <- sum(E(new_subtree, path=RF_path)$weight)
  RR_weight <- sum(E(new_subtree, path=RR_path)$weight)
  
  # print (FF_weight)
  # print (FR_weight)
  # print (RF_weight)
  # print (RR_weight)
  
  paths <- list(FF_path, FR_path, RF_path, RR_path)
  path_weights <- c(FF_weight, FR_weight, RF_weight, RR_weight)
  opt_path_idx <- head((which(path_weights == min(path_weights))), 1)
  opt_path <- paths[[opt_path_idx]]
  
  # print ("opt_path:")
  # print (opt_path)
  # print ("q_level_list:")
  # print (q_level_list)
  stopifnot (length(opt_path) == length(q_level_list))
  
  directions <- V(new_subtree)[opt_path]$type
  #print (directions)
  q_levels <- list()
  for (i in 1:length(directions))
  {
    if (directions[[i]] == "forward"){
      q_levels[[length(q_levels)+1]] <- q_level_list[[i]]
    }else{
      q_levels[[length(q_levels)+1]] <- rev(q_level_list[[i]])
    }
  }
  
  return(list(ql=q_levels, wt=min(path_weights)))
}

#order_q_node(q_level, dp)

extract_good_ordering <- function(pq_tree, curr_node, dist_matrix)
{
  if (V(pq_tree)[curr_node]$type == "leaf")
  {
    #print ("ordering leaf node")
    return (V(pq_tree)[curr_node]$name)
  }else if (V(pq_tree)[curr_node]$type == "P"){
    #print ("ordering P node")
    p_level <- list()
    for (child in V(pq_tree) [ nei(curr_node, mode="out") ])
    {
      p_level[[length(p_level)+1]] <- extract_good_ordering(pq_tree, child, dist_matrix)
    }
    p_level <- order_p_node(p_level, dist_matrix)
    p_level <- unlist(p_level)
    #print (p_level)
    return (p_level)
  }else if(V(pq_tree)[curr_node]$type == "Q"){
    #print ("ordering Q node")
    q_level <- list()
    for (child in V(pq_tree) [ nei(curr_node, mode="out") ])
    {
      q_level[[length(q_level)+1]] <- extract_good_ordering(pq_tree, child, dist_matrix)
    }
    q_level <- order_q_node(q_level, dist_matrix)
    q_level <- q_level$ql
    q_level <- unlist(q_level)
    #print (q_level)
    return (q_level)
  }
}

count_leaf_descendents <- function(pq_tree, curr_node, children_counts)
{
  if (V(pq_tree)[curr_node]$type == "leaf")
  {
    #V(pq_tree)[curr_node]$num_desc = 0
    children_counts[curr_node] = 0
    return(children_counts)
  } else {
    children_count = 0
    for (child in V(pq_tree) [ nei(curr_node, mode="out") ])
    {
      children_counts <- count_leaf_descendents(pq_tree, child, children_counts)
      if (V(pq_tree)[child]$type == "leaf")
      {
        children_count <- children_count + 1
      }
      else
      {
        children_count <- children_count + children_counts[child]
      }
    }
    #print (curr_node)
    children_counts[curr_node] = children_count
    return(children_counts)
  }
}

measure_diameter_path <- function(pq_tree, curr_node, path_lengths)
{
  if (V(pq_tree)[curr_node]$type != "Q")
  {
    #V(pq_tree)[curr_node]$num_desc = 0
    path_lengths[curr_node] = 0
    return(path_lengths)
  } else {
    
    children_count = 0
    for (child in V(pq_tree) [ nei(curr_node, mode="out") ])
    {
      children_counts <- count_leaf_descendents(pq_tree, child, children_counts)
      if (V(pq_tree)[child]$type == "leaf")
      {
        children_count <- children_count + 1
      }
      else
      {
        children_count <- children_count + children_counts[child]
      }
    }
    
    
    path_lengths[curr_node] = children_count
    return(children_counts)
  }
}

# Assign leaf nodes reachable in pq_tree from curr_node to assigned_state
assign_cell_lineage <- function(pq_tree, curr_node, assigned_state, node_states)
{
  if (V(pq_tree)[curr_node]$type == "leaf")
  {
    #V(pq_tree)[curr_node]$num_desc = 0
    #print (curr_node)
    node_states[V(pq_tree)[curr_node]$name] = assigned_state
    return(node_states)
  } else {
    for (child in V(pq_tree) [ nei(curr_node, mode="out") ])
    {
      node_states <- assign_cell_lineage(pq_tree, child, assigned_state, node_states)
    }
    return(node_states)
  }
}

extract_good_branched_ordering <- function(orig_pq_tree, curr_node, dist_matrix, num_branches, reverse_main_path=FALSE)
{
  pq_tree <- orig_pq_tree
  
  # children_counts <- rep(0, length(as.vector(V(pq_tree))))
  #     names(children_counts) <- V(pq_tree)$name
  # children_counts <- count_leaf_descendents(pq_tree, curr_node, children_counts)
  # 
  # branch_node_counts <- children_counts[V(res$subtree)[type == "P"]]
  # branch_node_counts <- sort(branch_node_counts, decreasing=TRUE)
  # print (branch_node_counts)
  
  
  branch_node_counts <- V(pq_tree)[type == "Q"]$diam_path_len
  names(branch_node_counts) <- V(pq_tree)[type == "Q"]$name
  branch_node_counts <- sort(branch_node_counts, decreasing=TRUE)
  #print (branch_node_counts)
  
  
  cell_states <- rep(NA, length(as.vector(V(pq_tree)[type=="leaf"])))
  names(cell_states) <- V(pq_tree)[type=="leaf"]$name
  
  cell_states <- assign_cell_lineage(pq_tree, curr_node, 1, cell_states)
  
  branch_point_roots <- list()
  
  # Start building the ordering tree. Each pseudo-time segment will be a node.
  branch_tree <- graph.empty()
  #root_branch_id <- "Q_1"
  #branch_tree <- branch_tree + vertex(root_branch_id)
  
  for (i in 1:num_branches)
  {
    #cell_states <- assign_cell_lineage(pq_tree, names(branch_node_counts)[i], i+1, cell_states)
    #print (head(cell_states))
    #print(names(branch_node_counts)[i])
    
    branch_point_roots[[length(branch_point_roots) + 1]] <- names(branch_node_counts)[i]
    branch_id <- names(branch_node_counts)[i]
    #print (branch_id)
    branch_tree <- branch_tree + vertex(branch_id)
    parents <- V(pq_tree)[nei(names(branch_node_counts)[i], mode="in")]
    if (length(parents) > 0 && parents$type == "P")
    {
      p_node_parent <- V(pq_tree)[nei(names(branch_node_counts)[i], mode="in")]
      parent_branch_id <- V(pq_tree)[nei(p_node_parent, mode="in")]$name
      #print (parent_branch_id)
      #print (branch_id)
      branch_tree <- branch_tree + edge(parent_branch_id, branch_id)
    }
    pq_tree[V(pq_tree) [ nei(names(branch_node_counts)[i], mode="in") ], names(branch_node_counts)[i] ] <- FALSE
  }
  
  #branch_point_roots[[length(branch_point_roots) + 1]] <- curr_node
  #branch_point_roots <- rev(branch_point_roots)
  branch_pseudotimes <- list()
  
  for (i in 1:length(branch_point_roots))
  {
    branch_ordering <- extract_good_ordering(pq_tree, branch_point_roots[[i]], dist_matrix)
    branch_ordering_time <- weight_of_ordering(branch_ordering, dist_matrix)
    names(branch_ordering_time) <- branch_ordering
    branch_pseudotimes[[length(branch_pseudotimes) + 1]] = branch_ordering_time
    names(branch_pseudotimes)[length(branch_pseudotimes)] = branch_point_roots[[i]]
  }
  
  cell_ordering_tree <- graph.empty()
  curr_branch <- "Q_1"
  
  extract_branched_ordering_helper <- function(branch_tree, curr_branch, cell_ordering_tree, branch_pseudotimes, dist_matrix, reverse_ordering=FALSE)
  {
    curr_branch_pseudotimes <- branch_pseudotimes[[curr_branch]]
    #print (curr_branch_pseudotimes)
    curr_branch_root_cell <- NA
    for (i in 1:length(curr_branch_pseudotimes))
    {
      cell_ordering_tree <- cell_ordering_tree + vertex(names(curr_branch_pseudotimes)[i])
      if (i > 1)
      {
        if (reverse_ordering == FALSE){
          cell_ordering_tree <- cell_ordering_tree + edge(names(curr_branch_pseudotimes)[i-1], names(curr_branch_pseudotimes)[i])
        }else{
          cell_ordering_tree <- cell_ordering_tree + edge(names(curr_branch_pseudotimes)[i], names(curr_branch_pseudotimes)[i-1])
        }
      }
    }
    
    if (reverse_ordering == FALSE)
    {
      curr_branch_root_cell <- names(curr_branch_pseudotimes)[1]
    }else{
      curr_branch_root_cell <- names(curr_branch_pseudotimes)[length(curr_branch_pseudotimes)]
    }
    
    for (child in V(branch_tree) [ nei(curr_branch, mode="out") ])
    {
      child_cell_ordering_subtree <- graph.empty()
      
      child_head <- names(branch_pseudotimes[[child]])[1]
      child_tail <- names(branch_pseudotimes[[child]])[length(branch_pseudotimes[[child]])]
      
      # find the closest cell in the parent branch for each of the head and the tail
      
      curr_branch_cell_names <- names(branch_pseudotimes[[curr_branch]])
      head_dist_to_curr <- dist_matrix[child_head, curr_branch_cell_names]
      closest_to_head <- names(head_dist_to_curr)[which(head_dist_to_curr == min(head_dist_to_curr))]
      
      head_dist_to_anchored_branch = NA
      branch_index_for_head <- NA
      
      head_dist_to_anchored_branch <- dist_matrix[closest_to_head, child_head]
      
      tail_dist_to_curr <- dist_matrix[child_tail, curr_branch_cell_names]
      closest_to_tail <- names(tail_dist_to_curr)[which(tail_dist_to_curr == min(tail_dist_to_curr))]
      
      tail_dist_to_anchored_branch = NA
      branch_index_for_tail <- NA
      
      tail_dist_to_anchored_branch <- dist_matrix[closest_to_tail, child_tail]
      
      if (tail_dist_to_anchored_branch < head_dist_to_anchored_branch)
      {
        reverse_child <- TRUE
      }else{
        reverse_child <- FALSE
      }
      
      res <- extract_branched_ordering_helper(branch_tree, child, child_cell_ordering_subtree, branch_pseudotimes, dist_matrix, reverse_child)
      child_cell_ordering_subtree <- res$subtree
      child_subtree_root <- res$root
      
      # Works, but slow:
      for (v in V(child_cell_ordering_subtree))
      {
        cell_ordering_tree <- cell_ordering_tree + vertex(V(child_cell_ordering_subtree)[v]$name)
      }
      
      edge_list <- get.edgelist(child_cell_ordering_subtree)
      for (i in 1:nrow(edge_list))
      {
        cell_ordering_tree <- cell_ordering_tree + edge(V(cell_ordering_tree)[edge_list[i, 1]]$name, V(cell_ordering_tree)[edge_list[i, 2]]$name)
      }   					
      
      if (tail_dist_to_anchored_branch < head_dist_to_anchored_branch)
      {
        cell_ordering_tree <- cell_ordering_tree + edge(closest_to_tail, child_subtree_root)
      }else{
        cell_ordering_tree <- cell_ordering_tree + edge(closest_to_head, child_subtree_root)
      }
      
    }
    
    return (list(subtree=cell_ordering_tree, root=curr_branch_root_cell, last_cell_state=1, last_cell_pseudotime=0.0))
  }
  
  res <- extract_branched_ordering_helper(branch_tree, curr_branch, cell_ordering_tree, branch_pseudotimes, dist_matrix, reverse_main_path)
  cell_ordering_tree <- res$subtree
  
  curr_state <- 1
  
  assign_cell_state_helper <- function(ordering_tree_res, curr_cell)
  {
    cell_tree <- ordering_tree_res$subtree
    V(cell_tree)[curr_cell]$cell_state = curr_state
    
    children <- V(cell_tree) [ nei(curr_cell, mode="out") ]
    ordering_tree_res$subtree <- cell_tree
    
    if (length(children) == 1){
      ordering_tree_res <- assign_cell_state_helper(ordering_tree_res, V(cell_tree)[children]$name)
    }else{
      for (child in children)	{
        curr_state <<- curr_state + 1
        ordering_tree_res <- assign_cell_state_helper(ordering_tree_res, V(cell_tree)[child]$name)
      }
    }
    return (ordering_tree_res)
  }
  
  res <- assign_cell_state_helper(res, res$root)
  
  assign_pseudotime_helper <- function(ordering_tree_res, dist_matrix, last_pseudotime, curr_cell)
  {
    cell_tree <- ordering_tree_res$subtree
    curr_cell_pseudotime <- last_pseudotime
    V(cell_tree)[curr_cell]$pseudotime = curr_cell_pseudotime
    #print (curr_cell_pseudotime)
    
    ordering_tree_res$subtree <- cell_tree
    children <- V(cell_tree) [ nei(curr_cell, mode="out") ]
    
    for (child in children)	{
      next_node <- V(cell_tree)[child]$name
      delta_pseudotime <- dist_matrix[curr_cell, next_node]
      ordering_tree_res <- assign_pseudotime_helper(ordering_tree_res, dist_matrix, last_pseudotime + delta_pseudotime, next_node)
    }
    
    return (ordering_tree_res)
  }
  
  res <- assign_pseudotime_helper(res, dist_matrix, 0.0, res$root)
  
  cell_names <- V(res$subtree)$name
  cell_states <- V(res$subtree)$cell_state
  cell_pseudotime <- V(res$subtree)$pseudotime
  # print (cell_names)
  # print (cell_states)
  # print (cell_pseudotime)
  ordering_df <- data.frame(sample_name = cell_names,
                            cell_state = factor(cell_states),
                            pseudo_time = cell_pseudotime)
  
  ordering_df <- arrange(ordering_df, pseudo_time)
  return(ordering_df)
}

reverse_ordering <- function(pseudo_time_ordering)
{
  pt <- pseudo_time_ordering$pseudo_time
  names(pt) <- pseudo_time_ordering$sample_name
  
  rev_pt <- -((pt - max(pt)))
  rev_df <- pseudo_time_ordering
  rev_df$pseudo_time <- rev_pt
  return(rev_df)
}


weight_of_ordering <- function(ordering, dist_matrix)
{
  time_delta <- c(0)
  curr_weight <- 0
  ep <- 0.01
  for (i in 2:length(ordering))
  {
    d <- dist_matrix[ordering[[i]], ordering[[i-1]]]
    curr_weight <- curr_weight + d + ep
    time_delta <- c(time_delta, curr_weight)
  }
  
  return(time_delta)
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
#' data(HSMM)
#' HSMM <- detectGenes(HSMM, min_expr=0.1)
detectGenes <- function(cds, min_expr=NULL){
  FM <- exprs(cds)
  if (is.null(min_expr))
  {
    min_expr <- cds@lowerDetectionLimit
  }
  FM_genes <- do.call(rbind, apply(FM, 1, 
                                       function(x) {
                                         return(data.frame(
                                           num_cells_expressed=sum(unlist(as.list(x)) >= min_expr)
                                         )
                                         )
                                       })
  )
  
  FM_cells <- do.call(rbind, apply(FM, 2, 
                                       function(x) {
                                         return(data.frame(
                                           num_genes_expressed=sum(unlist(as.list(x)) >= min_expr)
                                         )
                                         )
                                       })
  )
 
 
  fData(cds)$num_cells_expressed <-  FM_genes[row.names(fData(cds)),]
  
  pData(cds)$num_genes_expressed <-  FM_cells[row.names(pData(cds)),]
  
  cds
}

#' Sets the features (e.g. genes) to be used for ordering cells in pseudotime.
#' @param cds the CellDataSet upon which to perform this operation
#' @param ordering_genes a vector of feature ids (from the CellDataSet's featureData) used for ordering cells
#' @return an updated CellDataSet object
#' @export
setOrderingFilter <- function(cds, ordering_genes){
  fData(cds)$use_for_ordering <- row.names(fData(cds)) %in% ordering_genes
  cds
}

ica_helper <- function(X, n.comp, alg.typ = c("parallel", "deflation"), fun = c("logcosh", "exp"), alpha = 1, 
                       row.norm = TRUE, maxit = 200, tol = 1e-4, verbose = FALSE, w.init = NULL, use_irlba=TRUE){
  dd <- dim(X)
  d <- dd[dd != 1L]
  if (length(d) != 2L) 
    stop("data must be matrix-conformal")
  X <- if (length(d) != length(dd)) 
    matrix(X, d[1L], d[2L])
  else as.matrix(X)
  if (alpha < 1 || alpha > 2) 
    stop("alpha must be in range [1,2]")
  alg.typ <- match.arg(alg.typ)
  fun <- match.arg(fun)
  n <- nrow(X)
  p <- ncol(X)
  if (n.comp > min(n, p)) {
    message("'n.comp' is too large: reset to ", min(n, p))
    n.comp <- min(n, p)
  }
  if (is.null(w.init)) 
    w.init <- matrix(rnorm(n.comp^2), n.comp, n.comp)
  else {
    if (!is.matrix(w.init) || length(w.init) != (n.comp^2)) 
      stop("w.init is not a matrix or is the wrong size")
  }
  
  if (verbose) 
    message("Centering")
  X <- scale(X, scale = FALSE)
  X <- if (row.norm) 
    t(scale(X, scale = row.norm))
  else t(X)
  if (verbose) 
    message("Whitening")
  V <- X %*% t(X)/n
  
  if (verbose) 
    message("Finding SVD")
  if (use_irlba)
  {
    s <- irlba(V, min(n,p), min(n,p))  
    svs <- s$d  
  }
  else
  {
    s <- La.svd(V)
    svs <- s$d  
  }
  
  D <- diag(c(1/sqrt(s$d)))
  K <- D %*% t(s$u)
  K <- matrix(K[1:n.comp, ], n.comp, p)
  X1 <- K %*% X
  
  if (verbose) 
    message("Running ICA")
  if (alg.typ == "deflation") {
    a <- ica.R.def(X1, n.comp, tol = tol, fun = fun, 
                   alpha = alpha, maxit = maxit, verbose = verbose, 
                   w.init = w.init)
  }
  else if (alg.typ == "parallel") {
    a <- ica.R.par(X1, n.comp, tol = tol, fun = fun, 
                   alpha = alpha, maxit = maxit, verbose = verbose, 
                   w.init = w.init)
  }
  w <- a %*% K
  S <- w %*% X
  A <- t(w) %*% solve(w %*% t(w))
  return(list(X = t(X), K = t(K), W = t(a), A = t(A), S = t(S), svs=svs))
}

## This function was swiped from DESeq (Anders and Huber) and modified for our purposes
parametricDispersionFit <- function( means, disps )
{
  coefs <- c( .1, 1 )
  iter <- 0
  while(TRUE) {
    residuals <- disps / ( coefs[1] + coefs[2] / means )
    good <- which( (residuals > 1e-4) & (residuals < 15) )
    fit <- glm( disps[good] ~ I(1/means[good]),
                family=Gamma(link="identity"), start=coefs )
    oldcoefs <- coefs
    coefs <- coefficients(fit)
    if( !all( coefs > 0 ) )
      stop( "Parametric dispersion fit failed. Try a local fit and/or a pooled estimation. (See '?estimateDispersions')" )
    if( sum( log( coefs / oldcoefs )^2 ) < 1e-6 )
      break
    iter <- iter + 1
    if( iter > 10 ) {
      warning( "Dispersion fit did not converge." )
      break }
  }
  
  names( coefs ) <- c( "asymptDisp", "extraPois" )
  ans <- function( q )
    coefs[1] + coefs[2] / q
  #ans
  coefs
}

## This function was swiped from DESeq (Anders and Huber) and modified for our purposes
getVarianceStabilizedData <- function( cds, coefs ) {
  ncounts <- exprs(cds)
  ncounts <- t(t(ncounts) / sizeFactors(cds))
  ncounts <- round(ncounts)
  vst <- function( q )
    log( (1 + coefs["extraPois"] + 2 * coefs["asymptDisp"] * q +
            2 * sqrt( coefs["asymptDisp"] * q * ( 1 + coefs["extraPois"] + coefs["asymptDisp"] * q ) ) )
         / ( 4 * coefs["asymptDisp"] ) ) / log(2)
  vst( ncounts )
}

# FIXME: Maybe we're supposed to be extending BiocGenerics here?
estimateSizeFactors <- function( cds, locfunc = median )
{
  sizeFactors(cds) <- estimateSizeFactorsForMatrix(exprs(cds), locfunc=locfunc)
  cds
}

estimateSizeFactorsForMatrix <- function( counts, locfunc = median )
{
  CM <- round(counts)
  loggeomeans <- rowMeans( log(CM) )
  apply( CM, 2, function(cnts)
    exp( locfunc( ( log(cnts) - loggeomeans )[ is.finite(loggeomeans) ] ) ) )
}

vst_helper <- function(x, modelFormulaStr, expressionFamily){
  if (expressionFamily@vfamily == "negbinomial"){
    x <- x / Size_Factor
    expression <- round(x)
  }else if (expressionFamily@vfamily %in% c("gaussianff")){
    expression <- x
  }else{
    expression <- log10(x)
  }
  
  tryCatch({
    fitted_model <-  suppressWarnings(vgam(as.formula(modelFormulaStr), family=expressionFamily))
    disp_vals <- as.data.frame(predict(fitted_model))
    colnames(disp_vals) <- c("mu", "disp")
    disp_vals$disp <- signif(disp_vals$disp)
    disp_vals <- distinct(disp_vals, disp)
    
    disp_vals$mu <- exp(disp_vals$mu)
    disp_vals$disp <- 1.0/exp(disp_vals$disp)
    disp_vals
    #ddply(disp_vals, .(log_size), summarize, log_mu=mean(log_mu))
  }, 
  #warning = function(w) { FM_fit },
  error = function(e) { print (e); NULL }
  )
  
 }

# TODO: make this function return a function with coefs as an attr, 
# set the function as an attribute of cds, export the function
estimateDispersions <- function(cds, modelFormulaStr, relative_expr, cores=1)
{
  if (cores > 1){
    disp_table<-mcesApply(cds, 1, vst_helper, cores=cores, 
                          modelFormulaStr=modelFormulaStr, 
                          expressionFamily=cds@expressionFamily)
  }else{
    disp_table<-esApply(cds,1,vst_helper, 
                        modelFormulaStr=modelFormulaStr, 
                        expressionFamily=cds@expressionFamily)
  }
  disp_table <- do.call(rbind.data.frame, disp_table)
  
  coefs <- parametricDispersionFit(disp_table$mu, disp_table$disp)
  #print (coefs)
  
  coefs
}

checkSizeFactors <- function(cds)
{
  if (cds@expressionFamily@vfamily == "negbinomial")
  {
    if (is.null(sizeFactors(cds))){
      stop("Error: you must call estimateSizeFactors() before calling this function.")
    }  
  }
}

computeVarianceStabilizedValues <- function(cds, 
                                            model_row_names=NULL, 
                                            modelFormulaStr="expression~1", 
                                            relative_expr=TRUE, 
                                            cores=1) {
  if (is.null(model_row_names)){
    cds_subset <- cds
  }else{
    cds_subset <- cds[model_row_names,]
  }
  
  if (relative_expr){
    checkSizeFactors(cds_subset)
  }
  
  coefs <- estimateDispersions(cds_subset, 
                               modelFormulaStr=modelFormulaStr, 
                               relative_expr=relative_expr, 
                               cores=cores)
  
  vstExprs(cds) <-  getVarianceStabilizedData(cds, coefs)
  cds
}

#' Computes a projection of a CellDataSet object into a lower dimensional space
#' @param cds the CellDataSet upon which to perform this operation
#' @param max_components the dimensionality of the reduced space
#' @param use_irlba Whether to use the IRLBA package for ICA reduction.
#' @param pseudo_expr amount to increase expression values before dimensionality reduction
#' @param batch a vector of labels specifying batch for each cell, the effects of which will be removed prior to dimensionality reduction.
#' @param covariates a numeric vector or matrix specifying continuous effects to be removed prior to dimensionality reduction
#' @param ... additional arguments to pass to the dimensionality reduction function
#' @return an updated CellDataSet object
#' @details Currently, Monocle supports dimensionality reduction with Independent Component Analysis (ICA).
#' @export
reduceDimension <- function(cds, max_components=2, use_irlba=TRUE, pseudo_expr=1, batch=NULL, covariates=NULL, use_vst=F, ...){
  FM <- exprs(cds)
  
  # If we aren't using VST, then normalize the expression values by size factor
  if (use_vst == FALSE && cds@expressionFamily@vfamily == "negbinomial")
  {
    checkSizeFactors(cds)
    size_factors <- sizeFactors(cds)
    #print (size_factors)
    FM <- t(t(FM) / size_factors)
    #FM <- log2(FM)
  }
  
  if (is.null(fData(cds)$use_for_ordering) == FALSE)
    FM <- FM[fData(cds)$use_for_ordering,]
  
  FM <- FM + pseudo_expr
  FM <- FM[rowSds(FM) > 0,]
    
  if (use_vst){
    if (is.null(vstExprs(cds)) == FALSE){
      FM <- vstExprs(cds)
      FM <- FM[fData(cds)$use_for_ordering,]
    }else{
      stop("Error: set the variance-stabilized value matrix with vstExprs(cds) <- computeVarianceStabilizedValues() before calling this function with use_vst=TRUE")
    }
    #
  }else{
    FM <- log2(FM)
  }

  # TODO: get rid of this if possible.  Would be good just to subtract batch effects through the VST
  if (is.null(batch) == FALSE || is.null(covariates) == FALSE)
  {
    message("Removing batch effects")
    #FM <- log2(FM)
    FM <- removeBatchEffect(FM, batch=batch, covariates=covariates)
    FM <- 2^FM
  }
  
  #FM <- log2(FM)
  
  message("Reducing to independent components")
  
  #FM <- t(scale(t(FM)))
  #FM <- FM[rowSds(FM) > 0,]
  init_ICA <- ica_helper(t(FM), max_components, use_irlba=use_irlba, ...)
  
  x_pca <- t(t(FM) %*% init_ICA$K)
  W <- t(init_ICA$W)
  
  weights <- W
  
  # print(dim (init_ICA$K))
  # print(dim (solve(weights)))
  
  A <- t(solve(weights) %*% t(init_ICA$K))
  
  colnames(A) <- colnames(weights)
  rownames(A) <- rownames(FM)
  
  S <- weights %*% x_pca
  
  rownames(S) <- colnames(weights)
  colnames(S) <- colnames(FM) 
  
  reducedDimW(cds) <- W
  reducedDimA(cds) <- A
  reducedDimS(cds) <- S
  
  cds
}

#' Orders cells according to progress through a learned biological process.
#' @param cds the CellDataSet upon which to perform this operation
#' @param num_paths the number of end-point cell states to allow in the biological process.
#' @param reverse whether to reverse the beginning and end points of the learned biological process.
#' @param root_cell the name of a cell to use as the root of the ordering tree.
#' @return an updated CellDataSet object, in which phenoData contains values for State and Pseudotime for each cell
#' @export
orderCells <- function(cds, num_paths=1, reverse=FALSE, root_cell=NULL){
  
  adjusted_S <- t(cds@reducedDimS)
  
  dp <- as.matrix(dist(adjusted_S))
  
  #print (sum(rowSums(dp)))
  #dp <- as.matrix(dist(dp))
  #dp <- as.matrix(dist(adjusted_S))
  cellPairwiseDistances(cds) <- as.matrix(dist(adjusted_S))
  # Build an MST of the cells in ICA space.
  gp <- graph.adjacency(dp, mode="undirected", weighted=TRUE)
  dp_mst <- minimum.spanning.tree(gp)
  minSpanningTree(cds) <- dp_mst
  # Build the PQ tree
  next_node <<- 0
  res <- pq_helper(dp_mst, use_weights=FALSE, root_node=root_cell)
  #stopifnot(length(V(res$subtree)[type == "leaf"]) == nrow(pData(cds)))
  
  cc_ordering <- extract_good_branched_ordering(res$subtree, res$root, cellPairwiseDistances(cds), num_paths, reverse)
  row.names(cc_ordering) <- cc_ordering$sample_name
  
  pData(cds)$Pseudotime <-  cc_ordering[row.names(pData(cds)),]$pseudo_time
  pData(cds)$State <-  cc_ordering[row.names(pData(cds)),]$cell_state
  cds
}

fit_model_helper <- function(x, modelFormulaStr, expressionFamily, relative_expr){
  if (expressionFamily@vfamily == "negbinomial"){
    if (relative_expr)
    {
      x <- x / Size_Factor
    }
    expression <- round(x)
  }else if (expressionFamily@vfamily %in% c("gaussianff")){
    expression <- x
  }else{
    expression <- log10(x)
  }
  
  tryCatch({
    FM_fit <-  suppressWarnings(vgam(as.formula(modelFormulaStr), family=expressionFamily))
    FM_fit
  }, 
  #warning = function(w) { FM_fit },
  error = function(e) { print (e); NULL }
  )
}

diff_test_helper <- function(x, fullModelFormulaStr, reducedModelFormulaStr, expressionFamily, relative_expr){
  if (expressionFamily@vfamily == "negbinomial"){
    if (relative_expr)
    {
      x <- x / Size_Factor
    }
    expression <- round(x)
  }else if (expressionFamily@vfamily %in% c("gaussianff")){
    expression <- x
  }else{
    expression <- log10(x)
  }
  
  test_res <- tryCatch({
    full_model_fit <- suppressWarnings(vgam(as.formula(fullModelFormulaStr), family=expressionFamily))
    reduced_model_fit <- suppressWarnings(vgam(as.formula(reducedModelFormulaStr), family=expressionFamily))
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


mcesApply <- function(X, MARGIN, FUN, cores=1, ...) {
  parent <- environment(FUN)
  if (is.null(parent))
    parent <- emptyenv()
  e1 <- new.env(parent=parent)
  multiassign(names(pData(X)), pData(X), envir=e1)
  environment(FUN) <- e1
  cl <- makeCluster(cores)
  clusterEvalQ(cl, {require(VGAM);})
  if (MARGIN == 1){
    res <- parRapply(cl, exprs(X), FUN, ...)
  }else{
    res <- parCapply(cl, exprs(X), FUN, ...)
  }

  stopCluster(cl)
  res
}

#' Fits a model for each gene in a CellDataSet object.
#' @param cds the CellDataSet upon which to perform this operation
#' @param modelFormulaStr a formula string specifying the model to fit for the genes.
#' @param cores the number of processor cores to be used during fitting.
#' @return a list of VGAM model objects
#' @export
#' @details
#' 
#' This function fits a Tobit-family vector generalized additive model (VGAM) from the VGAM package for each gene in a CellDataSet. 
#' The default formula string speficies that the (log transformed) expression values follow a Tobit distribution with upper and lower bounds
#' specificed by max_expr and min_expr, respectively. By default, expression levels are modeled as smooth functions of the Pseudotime value of each 
#' cell. That is, expression is a function of progress through the biological process.  More complicated formulae can be provided to account for
#' additional covariates (e.g. day collected, genotype of cells, media conditions, etc).
fitModel <- function(cds,
                     modelFormulaStr="expression~sm.ns(Pseudotime, df=3)",
                     relative_expr=TRUE,
                     cores=1){
  if (cores > 1){
    f<-mcesApply(cds, 1, fit_model_helper, cores=cores, 
                 modelFormulaStr=modelFormulaStr, 
                 expressionFamily=cds@expressionFamily,
                 relative_expr=relative_expr)
    f
  }else{
    f<-esApply(cds,1,fit_model_helper, 
               modelFormulaStr=modelFormulaStr, 
               expressionFamily=cds@expressionFamily,
               relative_expr=relative_expr)
    f
  }
}

#' Response values
#' 
#' Generates a matrix of response values for a set of fitted models
#' @param models a list of models, e.g. as returned by fitModels()
#' @return a matrix where each row is a vector of response values for a particular feature's model, and columns are cells.
#' @export
responseMatrix <- function(models){
  res_list <- lapply(models, function(x) { 
    if (is.null(x)) { NA } else { 
      if (x@family@vfamily  == "negbinomial"){
        predict(x, type="response") 
      }else if (x@family@vfamily %in% c("gaussianff")){
        predict(x, type="response")
      }else{
        10^predict(x, type="response") 
      }
    } 
  } )
  res_list_lengths <- lapply(res_list[is.na(res_list) == FALSE], length)
  
  stopifnot(length(unique(res_list_lengths)) == 1)
  num_na_fits <- length(res_list[is.na(res_list)])
  if (num_na_fits > 0){
    na_matrix<- matrix(rep(rep(NA, res_list_lengths[[1]]), num_na_fits),nrow=num_na_fits) 
    row.names(na_matrix) <- names(res_list[is.na(res_list)])
    
    non_na_matrix <- t(do.call(cbind, lapply(res_list[is.na(res_list) == FALSE], unlist)))
    row.names(non_na_matrix) <- names(res_list[is.na(res_list) == FALSE])
    res_matrix <- rbind(non_na_matrix, na_matrix)
    res_matrix <- res_matrix[names(res_list),]
  }else{
    res_matrix <- t(do.call(cbind, lapply(res_list, unlist)))
    row.names(res_matrix) <- names(res_list[is.na(res_list) == FALSE])
  }
  
  res_matrix
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
#' @return a data frame containing the p values and q-values from the likelihood ratio tests on the parallel arrays of models.
#' @export
differentialGeneTest <- function(cds, 
                                 fullModelFormulaStr="expression~sm.ns(Pseudotime, df=3)",
                                 reducedModelFormulaStr="expression~1", cores=1, relative_expr=TRUE){
  if (cds@expressionFamily@vfamily == "negbinomial"){
    if (is.null(sizeFactors(cds))){
      stop("Error: to call this function with relative_expr==TRUE, you must first call estimateSizeFactors() on the CellDataSet.")
    }
  }
  
  if (cores > 1){
    diff_test_res<-mcesApply(cds, 1, diff_test_helper, cores=cores, 
                             fullModelFormulaStr=fullModelFormulaStr,
                             reducedModelFormulaStr=reducedModelFormulaStr,
                             expressionFamily=cds@expressionFamily,
                             relative_expr=relative_expr)
    diff_test_res
  }else{
    diff_test_res<-esApply(cds,1,diff_test_helper, 
               fullModelFormulaStr=fullModelFormulaStr,
               reducedModelFormulaStr=reducedModelFormulaStr, 
               expressionFamily=cds@expressionFamily, 
               relative_expr=relative_expr)
    diff_test_res
  }
  
  diff_test_res <- do.call(rbind.data.frame, diff_test_res)
  diff_test_res$qval <- p.adjust(diff_test_res$pval, method="BH")
  diff_test_res
}

# differentialGeneTest <- function(cds, 
#                                  fullModelFormulaStr="expression~sm.ns(Pseudotime, df=3)",
#                                  reducedModelFormulaStr="expression~1", cores=1){
#   full_model_fits <- fitModel(cds,  modelFormulaStr=fullModelFormulaStr, cores=cores)
#   reduced_model_fits <- fitModel(cds, modelFormulaStr=reducedModelFormulaStr, cores=cores)
#   test_res <- compareModels(full_model_fits, reduced_model_fits)
#   test_res
# }

#' Plots the minimum spanning tree on cells.
#'
#' @param expr_matrix a matrix of expression values to cluster together
#' @param k how many clusters to create
#' @param method the distance function to use during clustering
#' @param ... extra parameters to pass to pam() during clustering
#' @return a pam cluster object
#' @export
#' @examples
#' \dontrun{
#' full_model_fits <- fitModel(HSMM[sample(nrow(fData(HSMM_filtered)), 100),],  modelFormulaStr="expression~sm.ns(Pseudotime)")
#' expression_curve_matrix <- responseMatrix(full_model_fits)
#' clusters <- clusterGenes(expression_curve_matrix, k=4)
#' plot_clusters(HSMM_filtered[ordering_genes,], clusters)
#' }
clusterGenes<-function(expr_matrix, k, method=function(x){as.dist((1 - cor(t(x)))/2)}, ...){
  expr_matrix <- expr_matrix[rowSums(is.na(expr_matrix)) == 0,] 
  #expr_matrix <- t(scale(t(log10(expr_matrix))))
  expr_matrix <- expr_matrix[is.nan(rowSums(expr_matrix)) == FALSE,] 
  expr_matrix[is.na(expr_matrix)] <- 0
  n<-method(expr_matrix)
  clusters<-pam(n,k, ...)
  class(clusters)<-"list"
  clusters$exprs<-expr_matrix
  clusters
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
  
  FM <- exprs(cds)
  if (cds@expressionFamily@vfamily == "negbinomial")
  {
    expression_lower_thresh <- expression_lower_thresh / colSums(FM)
    expression_upper_thresh <- expression_upper_thresh / colSums(FM)
    FM <- t(t(FM)/colSums(FM))
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


selectGenesInExpressionRange <- function(cds, min_expression_threshold = -Inf, max_expression_threshold = Inf, detectionLimit=-Inf, stat_fun=median){
  gene_nz_median = apply(exprs(cds), 1, function(x) { x <- x[x > detectionLimit]; stat_fun(x)})
  #gene_nz_median
  names(gene_nz_median[is.na(gene_nz_median) == FALSE & gene_nz_median > min_expression_threshold & gene_nz_median < max_expression_threshold ])
}


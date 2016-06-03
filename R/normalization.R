# Calculate the probability vector 
makeprobsvec<-function(p){
  phat<-p/sum(p)
  phat[is.na(phat)] = 0
  phat
}

# Calculate the probability matrix for a relative abundance matrix
makeprobs<-function(a){
  colSums<-apply(a,2,sum)
  b<-Matrix::t(Matrix::t(a)/colSums)
  b[is.na(b)] = 0
  b
}
# makeprobs <- function(a){
#     colSums<-apply(a,2, function(x) sum(x[x > 0 & is.finite(x)]))
#     b<-Matrix::t(Matrix::t(a)/colSums)
#     b[is.na(b)] = 0; b[b <= 0] = 0
#     b
# }

# Calculate the Shannon entropy based on the probability vector
shannon.entropy <- function(p) {
  if (min(p) < 0 || sum(p) <=0)
    return(Inf)
  p.norm<-p[p>0]/sum(p)
  -sum( log10(p.norm)*p.norm)
}

# Calculate the Jessen-Shannon distance for two probability distribution 
JSdistVec <- function (p, q) 
{
  JSdiv <- shannon.entropy((p + q)/2) - (shannon.entropy(p) + 
                                           shannon.entropy(q)) * 0.5
  JSdiv[is.infinite(JSdiv)] <- 1
  JSdiv[JSdiv < 0] <- 0
  JSdist <- sqrt(JSdiv)
  JSdist
}

# Recover the absolute transcript counts based on m,c and t estimate for a single cell (Used in the optimization function)
opt_norm_t <- function(t, fpkm, mRNAs_for_mode, kb_slope, kb_intercept, expr_thresh = 0.1, pseudocnt = NULL, return_norm = FALSE) {
  a_matrix <- matrix(c(log10(t), 1, kb_slope,
                       -1), ncol = 2, nrow = 2, byrow = T)
  colnames(a_matrix) <- c("k", "b")
  b_matrix <- matrix(c(log10(mRNAs_for_mode), -kb_intercept), nrow = 2, byrow = T)
  kb <- Matrix::t(solve(a_matrix, b_matrix))
  
  k <- kb[1]
  b <- kb[2]
  #print(kb)
  tmp <- k * log10(fpkm) + b
  abs_cnt <- 10^tmp
  
  if(return_norm) return(abs_cnt)
  
  selected <- abs_cnt[abs_cnt > expr_thresh & is.na(abs_cnt) == FALSE]
  if (length(selected) == 0) 
    return (0);
  
  if(!is.null(pseudocnt)){
    10^dmode(log10(selected + pseudocnt)) #keep the original scale
  }
  else
    10^dmode(log10(selected))
}

#linear transform by kb
opt_norm_kb <- function(relative_expr, kb) {
  
  k <- kb[1]
  b <- kb[2]
  tmp <- k * log10(relative_expr) + b
  abs_cnt <- 10^tmp
  
  selected <- abs_cnt[abs_cnt > 0 & is.na(abs_cnt) == FALSE]
  if (length(selected) == 0) 
    return (0);
  
  10^dmode(log10(selected))
}

#use the deconvoluated linear regression parameters to normalize the log_relative_expr
norm_kb <- function(kb, exprs_cds) {
  k <- kb[1]
  b <- kb[2]
  tmp <- k * log10(exprs_cds) + b
  norm_exprs <- 10^tmp
  
  norm_exprs
}

#use gaussian kernel to calculate the mode of transcript counts
dmode <- function(x, breaks="Sturges") {
  if (length(x) < 2) return (0);
  den <- density(x, kernel=c("gaussian"))
  ( den$x[den$y==max(den$y)] )
}

# Calculate the optimization function based on mode of transcript counts, Jessen-Shannon distance as well as the hypothetical total RNA counts
optim_mc_func_fix_c <- function (kb_slope_intercept, kb_intercept = NULL, t_estimate = estimate_t(TPM_isoform_count_cds),
          relative_expr_matrix = relative_expr_mat, split_relative_expr_matrix = split_relative_exprs,
          alpha = rep(1, ncol(relative_expr_matrix)), total_RNAs = rep(150000, ncol(relative_expr_matrix)),
          cores = 1, weight_mode=0.17, weight_relative_expr=0.50, weight_total_rna=0.33, verbose = F,  ...) {
  data('spike_df') #add the spikein dataset

  if(is.null(spike_df$log_numMolecule))
    spike_df$log_numMolecules <- log10(spike_df$numMolecules)
  
  if(is.null(kb_intercept)) {
    kb_slope_val <- kb_slope_intercept[1]
    kb_intercept_val <- kb_slope_intercept[2]
  }
  else {
    kb_slope_val <- kb_slope_intercept[1]
    kb_intercept_val <- kb_intercept 
  }
  
  cell_num <- ncol(relative_expr_matrix)
  names(t_estimate) <- colnames(relative_expr_matrix)
  split_t <- split(Matrix::t(t_estimate), col(as.matrix(Matrix::t(t_estimate)), as.factor = T))
  
  total_rna_df <- data.frame(Cell = colnames(relative_expr_matrix), t_estimate = t_estimate, alpha_v = alpha)
  
  t_k_b_solution <- tryCatch({
    k_b_solution <- plyr::ddply(total_rna_df, .(Cell), function(x) {
      a_matrix <- matrix(c(log10(x[, "t_estimate"]), 1,
                           kb_slope_val, -1), ncol = 2, nrow = 2, byrow = T)
      colnames(a_matrix) <- c("k", "b")
      b_matrix <- matrix(c(log10(x[, "alpha_v"]), -kb_intercept_val), nrow = 2, byrow = T)
      k_b_solution <- Matrix::t(solve(a_matrix, b_matrix))
    })
    k_b_solution},
    error = function(e) {print(e); c(NA, NA)}
  )
  
  if(any(is.na(t_k_b_solution)))
    return(NA)
  
  cell_dmode <- tryCatch({
    if(cores > 1){
      cell_dmode <- unlist(mcmapply(opt_norm_t, split_t, split_relative_expr_matrix, alpha, kb_slope = kb_slope_val, kb_intercept = kb_intercept_val, pseudocnt = 0.01, mc.cores = cores))
      adj_est_std_cds <- unlist(mcmapply(opt_norm_t, split_t, split_relative_expr_matrix, alpha, kb_slope = kb_slope_val, kb_intercept = kb_intercept_val, pseudocnt = 0.01, return_norm = T, mc.cores = cores))
    }
    else {
      cell_dmode <- unlist(mapply(opt_norm_t, split_t, split_relative_expr_matrix, alpha, kb_slope = kb_slope_val, kb_intercept = kb_intercept_val, pseudocnt = 0.01))
      adj_est_std_cds <- unlist(mapply(opt_norm_t, split_t, split_relative_expr_matrix,  alpha, kb_slope = kb_slope_val, kb_intercept = kb_intercept_val, pseudocnt = 0.01, return_norm = T))
    }
    cell_dmode},
    error = function(e) {print(e); NA}
  )
  
  if(any(is.na(cell_dmode)))
    return(NA)
  
  sum_total_cells_rna <- colSums(adj_est_std_cds)
  
  #Jenn-Shannon divergence:
  p_df <- makeprobs(relative_expr_matrix) #relative expression
  p_list <- split(Matrix::t(p_df), 1:ncol(p_df))
  q_df <- makeprobs(adj_est_std_cds) #no rounding
  q_list <- split(Matrix::t(q_df), 1:ncol(q_df))
  
  dist_divergence <- mcmapply(function(x, y) {
    JSdistVec(x, y)
  }, p_list, q_list, mc.cores = cores)
  
  
  gm_dist_divergence <- exp(mean(log(dist_divergence)))
  
  #total RNA MSE
  sum_total_cells_rna_finite <- sum_total_cells_rna[is.finite(sum_total_cells_rna)]
  total_RNAs_finite <- total_RNAs[is.finite(sum_total_cells_rna)]
  total_rna_obj <- exp(mean(log(((sum_total_cells_rna_finite -  total_RNAs_finite)/total_RNAs_finite)^2))) #use geometric mean to avoid outlier cells

  #mode MSE
  mode_obj <- exp(mean(log(((cell_dmode[is.finite(cell_dmode)] - alpha[is.finite(cell_dmode)])/alpha[is.finite(cell_dmode)])^2)))

  relative_expr_obj <- gm_dist_divergence * 10 #give more weight on this
  
  #objective
  res <- weight_mode * mode_obj + weight_relative_expr * relative_expr_obj + weight_total_rna * total_rna_obj
   
  if(verbose){
    message('current m, c values are ', paste(kb_slope_val, kb_intercept_val, sep = ', '))
    message('objective is:')
    print(res)

    message('\n\ntotal_rna_obj is ', total_rna_obj)
    message('mode_obj is ', mode_obj)
    message('relative_expr_obj is ', relative_expr_obj)

    message('\n\nmean modes are:')
    print (mean(cell_dmode))
    message('mean target modes are:')
    print (mean(alpha))
    message('mean mode delta is:')
    print (mean(cell_dmode - alpha))
    message('mean total RNA delta is:')
    print (mean(sum_total_cells_rna_finite -  total_RNAs_finite))

  }
  #   return(list(m = m_val, c = c_val, dmode_rmse_weight_total = dmode_rmse_weight_total, gm_dist_divergence = gm_dist_divergence, dist_divergence_round = dist_divergence_round,
  #               cell_dmode = cell_dmode, t_k_b_solution = t_k_b_solution, sum_total_cells_rna = sum_total_cells_rna, optim_res = res))
  #
  if(is.finite(res))
    return(res)
  else
    return(1e10 * runif(1)) #Census should not run this part since non-finite values are removed 
}

#' Find the most commonly occuring relative expression value in each cell
#' 
#' Converting relative expression values to mRNA copies per cell requires
#' knowing the most commonly occuring relative expression value in each cell
#' This value typically corresponds to an RPC value of 1. This function 
#' finds the most commonly occuring (log-transformed) relative expression value
#' for each column in the provided expression matrix. 
#' 
#' @param relative_expr_matrix a matrix of relative expression values for 
#' values with each row and column representing genes/isoforms and cells, 
#' respectively. Row and column names should be included. 
#' Expression values should not be log-transformed.
#' @param relative_expr_thresh Relative expression values below this threshold 
#' are considered zero.
#' @return an vector of most abundant relative_expr value corresponding to the 
#' RPC 1. 
#' @details This function estimates the most abundant relative expression value 
#' (t^*) using a gaussian kernel density function. It can also optionally 
#' output the t^* based on a two gaussian mixture model
#' based on the smsn.mixture from mixsmsn package
#' @export
#' @examples
#' \dontrun{
#' HSMM_fpkm_matrix <- exprs(HSMM)
#' t_estimate = estimate_t(HSMM_fpkm_matrix)
#'}

estimate_t <- function(relative_expr_matrix, relative_expr_thresh = 0.1) {
  #apply each column
  unlist(apply(relative_expr_matrix, 2, function(relative_expr) 10^mean(dmode(log10(relative_expr[relative_expr > relative_expr_thresh]))))) #avoid multiple output
}


#' Make an educated guess on the spike-based slope regression parameters
#' @importFrom plyr ldply
#' @importFrom MASS rlm
calibrate_mc <- function(total_mRNA, capture_rate, ladder, total_ladder_transcripts, reads, trials=100){
  if(length(capture_rate) > 1) capture_rate <- dmode(capture_rate) 
  if(length(reads) > 1) reads <- 10^dmode(log10(reads)) 
  if(length(total_mRNA) > 1) total_mRNA <- 10^dmode(log10(total_mRNA))

  kb_df <- ldply (seq(0,trials, by=1), function(i){
    
    hypothetical_ladder <- rmultinom(1, (total_ladder_transcripts * capture_rate), (ladder / sum(ladder)))
    asympt_proportions <- hypothetical_ladder / (sum(hypothetical_ladder) + (total_mRNA * capture_rate))
    asympt_proportions <- c(asympt_proportions, 1 - sum(asympt_proportions))
    
    trial_reads <- rmultinom(1, reads, asympt_proportions)
    trial_tpm <- 1e6 * trial_reads / sum(trial_reads)
    
    hypothetical_ladder_tpm <- trial_tpm[1:length(hypothetical_ladder)]
    
    ladder_df <- data.frame(hypothetical_ladder_tpm=hypothetical_ladder_tpm, 
                            hypothetical_ladder=hypothetical_ladder,
                            ladder = ladder)
    ladder_df <- subset(ladder_df, hypothetical_ladder_tpm > 0 & ladder > 0)
    
    ladder_reg <- MASS::rlm (log10(ladder) ~ log10(hypothetical_ladder_tpm), data=ladder_df)
    b <- coefficients(ladder_reg)[1]
    k <- coefficients(ladder_reg)[2]
    data.frame(k=k,b=b)
  })
  
  kb_reg <- MASS::rlm (b ~ k, data=kb_df)
  return (list(m=coefficients(kb_reg)[2], c=coefficients(kb_reg)[1], kb_df = kb_df))
}

calibrate_mode <- function(ind, tpm_distribution, ladder, total_ladder_transcripts, total_mRNA, capture_rate, reads, trials=100){
  tpm_distribution <- tpm_distribution[[ind]] / sum(tpm_distribution[[ind]]) * 1e6
  total_mRNA <- total_mRNA[ind] 
  capture_rate <- capture_rate[ind]
  reads <- reads[ind]

  mode_df <- ldply (seq(0,trials, by=1), function(i){
    hypothetical_ladder <- rmultinom(1, (total_ladder_transcripts * capture_rate), (ladder / sum(ladder)))
    
    asympt_proportions <- hypothetical_ladder / (sum(hypothetical_ladder) + (total_mRNA * capture_rate))
    asympt_proportions <- c(asympt_proportions, 1 - sum(asympt_proportions))
    
    trial_reads <- rmultinom(1, reads, asympt_proportions)
    trial_tpm <- 1e6 * trial_reads / sum(trial_reads)
    
    hypothetical_ladder_tpm <- trial_tpm[1:length(hypothetical_ladder)]
    tpm_distribution <-  tpm_distribution * trial_reads[length(hypothetical_ladder) + 1] / sum(trial_reads) #put the tpm for the spike-in and the endogenous RNA at the same space
    
    ladder_df <- data.frame(asympt_proportions_tpm=asympt_proportions[1:length(hypothetical_ladder)] * 10e6,
                            hypothetical_ladder_tpm=hypothetical_ladder_tpm, 
                            hypothetical_ladder=hypothetical_ladder,
                            ladder = ladder)
    ladder_df <- subset(ladder_df, hypothetical_ladder_tpm > 0 & ladder > 0)
    
    ladder_reg <- tryCatch({
       ladder_reg <-  MASS::rlm (log10(ladder) ~ log10(hypothetical_ladder_tpm), data=ladder_df)

        ladder_reg
      }, error = function(e) {
        print(e)
        NULL
      })

    if(is.null(ladder_reg))
      return(data.frame(hypothetical_mode=NULL, k = NULL, b = NULL, 
      ladder = NULL, hypothetical_mode = NULL))

    b <- coefficients(ladder_reg)[1]
    k <- coefficients(ladder_reg)[2]

    fpkm_hypothetical_mode <- dmode(log10(tpm_distribution[tpm_distribution > 0]))

    df <- data.frame(hypothetical_mode=10^(k * fpkm_hypothetical_mode + b), 
      hypothetical_ladder_tpm = hypothetical_ladder_tpm, ladder = ladder 
      )
    df
  })

  ladder_df <- subset(mode_df, hypothetical_ladder_tpm > 0 & ladder > 0)
  ladder_reg <-  MASS::rlm (log10(ladder) ~ log10(hypothetical_ladder_tpm), data=ladder_df)
  b <- coefficients(ladder_reg)[1]
  k <- coefficients(ladder_reg)[2]

  return(data.frame(mean_hypotetical_mode = dmode(mode_df$hypothetical_mode), 
    k = k, b = b))
}


#' Transform relative expression values into absolute transcript counts.
#' 
#' Transform a relative expression matrix to absolute transcript matrix based on the inferred linear regression parameters from most abundant isoform relative expression value.
#' This function takes a relative expression matrix and a vector of estimated most abundant expression value from the isoform-level matrix and transform it into absolute transcript number.
#' It is based on the observation that the recovery efficient of the single-cell RNA-seq is relative low and that most expressed isoforms of gene in a single cell therefore only sequenced one copy so that the 
#' most abundant isoform log10-FPKM (t^*) will corresponding to 1 copy transcript. It is also based on the fact that the spikein regression parameters k/b for each cell will fall on a line because of the 
#' intrinsic properties of spikein experiments. We also assume that if we perform the same spikein experiments as Treutlein et al. did, the regression parameters should also fall on a line in the same way. The
#' function takes the the vector t^* and the detection limit as input, then it uses the t^* and the m/c value corresponding to the detection limit to calculate two parameters vectors k^* and b^* (corresponding to each cell)
#' which correspond to the slope and intercept for the linear conversion function between log10 FPKM and log10 transcript counts. The function will then apply a linear transformation 
#' to convert the FPKM to estimated absolute transcript counts based on the the k^* and b^*. The default m/c values used in the algoritm are 3.652201, 2.263576, respectively.
#' 
#' @param relative_cds the cds object of relative expression values for single cell RNA-seq with each row and column representing genes/isoforms and cells. Row and column names should be included
#' @param t_estimate an vector for the estimated most abundant FPKM value of isoform for a single cell. Estimators based on gene-level relative expression can also give good approximation but estimators
#' based on isoform FPKM will give better results in general
#' @param modelFormulaStr modelformula used to grouping cells for transcript counts recovery. Default is "~ 1", which means to recover the transcript counts from all cells.
#' @param m the initial guess of the slope for the regression line of b_i (intercept of spikein regression in i-th cell) and k_i (slope of spikein regression in i-th cell)
#' @param c the initial guess of the intercept for the regression line of b_i (intercept of spikein regression in i-th cell) and k_i (slope of spikein regression in i-th cell). Note that this value can be approximated by calculation based on the spikein data (See method section in the paper).  
#' @param m_rng the range of m values used by the optimization function to optimize. By default, it is between -10 and -0.1
#' @param c_rng the range of c values. Since we can approximate this value based on spikein data. By default, it is fixed. Under certain cases, we can provide a small range for the optimization function to optimize. 
#' @param ERCC_controls the FPKM matrix for each ERCC spike-in transcript in the cells if user wants to perform the transformation based on their spike-in data. Note that the row and column names should match up with the ERCC_annotation and relative_exprs_matrix respectively. 
#' @param ERCC_annotation the ERCC_annotation matrix from illumina USE GUIDE which will be ued for calculating the ERCC transcript copy number for performing the transformation. 
#' @param volume the approximate volume of the lysis chamber (nanoliters). Default is 10
#' @param dilution the dilution of the spikein transcript in the lysis reaction mix. Default is 40, 000. The number of spike-in transcripts per single-cell lysis reaction was calculated from
#' @param mixture_type the type of spikein transcripts from the spikein mixture added in the experiments. By default, it is mixture 1. Note that m/c we inferred are also based on mixture 1. 
#' @param detection_threshold the lowest concentration of spikein transcript considered for the regression. Default is 800 which will ensure (almost) all included spike transcripts expressed in all the cells. Also note that the value of c is based on this concentration. 
#' @param alpha_v the hypothesized mode of transcript counts for each cell. Default is 1. 
#' @param total_RNAs the guess of total transcript counts in a single cell. Default is 150000. 
#' @param weight the weight for the first term associate with the mode of transcript in the optimization function (See the method section in the paper for more details)
#' @param return_all parameter for the intended return results. If setting TRUE, matrix of m, c, k^*, b^* as well as the transformed absolute cds will be returned
#' in a list format
#' @param cores number of cores to perform the recovery. The recovery algorithm is very efficient so multiple cores only needed when we have very huge number of cells or genes.
#' @param verbose a logical flag to determine whether or not we should print all the optimization details 
#' @param optim_num The number of rounds of optimization to perform.
#' @return an matrix of absolute count for isoforms or genes after the transformation. It can also be a list including the m, c values, the dataframe for the k_i/b_i in each cell as well as the recovered absolute transcript counts if return_all is set to be TRUE. 
#' @export
#' @importFrom plyr ddply
#' @examples
#' \dontrun{
#' HSMM_relative_expr_matrix <- exprs(HSMM)
#' HSMM_abs_matrix <- relative2abs(HSMM_relative_expr_matrix, t_estimate = estimate_t(HSMM_relative_expr_matrix))
#'}

relative2abs <- function(relative_cds, 
  t_estimate = estimate_t(exprs(relative_cds)),
  modelFormulaStr = "~1", 
  #kb_slope = -3.652201, 
  #kb_intercept = 2.263576, 
  #kb_slope_rng = c(-10, -0.1), 
  #kb_intercept_rng = c(kb_intercept, kb_intercept), 
  kb_slope = NULL,
  kb_intercept = NULL,
  kb_slope_rng = NULL,
  kb_intercept_rng = NULL,
  use_fixed_intercept=TRUE,
  ERCC_controls = NULL, 
  ERCC_annotation = NULL, 
  volume = 10, 
  dilution = 40000, 
  mixture_type = 1,
  detection_threshold = 800, 
  expected_mRNA_mode = NULL, 
  expected_total_mRNAs = 100000, #based on lung endogenous RNA
  calibrate_total_mRNA = T,
  calculation_trials = 100, 
  reads_per_cell = 1e6,
  expected_capture_rate = 0.25,
  weight_mode=0.17, 
  weight_relative_expr=0.5, 
  weight_total_rna=0.33,
  verbose = FALSE, 
  return_all = FALSE, 
  cores = 1, 
  optim_num = 1, ...) {
  relative_expr_matrix <- exprs(relative_cds)
  # relative_expr_matrix <- apply(relative_expr_matrix, 2, function(x) x / sum(x) * 10^6) #convert to TPM

  parameters <- c(t_estimate, volume, dilution, detection_threshold,  weight_mode, weight_relative_expr, weight_total_rna,  optim_num)
  if(any(c(!is.finite(parameters), is.null(parameters))))
    stop('Your input parameters should not contain either null or infinite numbers')

  if (detection_threshold < 0.01430512 | detection_threshold >
        7500) 
    stop("concentration detection limit should be between 0.01430512 and 7500")

  if (mixture_type == 1) 
    mixture_name = "conc_attomoles_ul_Mix1"
  else mixture_name = "conc_attomoles_ul_Mix2"
  Cell <- NULL
  if (!is.null(ERCC_controls) | !is.null(ERCC_annotation)) {
    if (is.null(ERCC_controls) | is.null(ERCC_annotation)) 
      stop("If you want to transform the data to copy number with your spikein data, please provide both of ERCC_controls and ERCC_annotation data frame...")
    valid_ids <- which(ERCC_annotation[, mixture_name] >= 
                         detection_threshold)
    if (verbose) 
      message("Performing robust linear regression for each cell based on the spikein data...")
    molModels <- apply(ERCC_controls, 2, function(cell_exprs, 
                                                  input.ERCC.annotation, valid_ids) {
      spike_df <- input.ERCC.annotation
      spike_df <- cbind(spike_df, cell_exprs[row.names(spike_df)])
      colnames(spike_df)[length(colnames(spike_df))] <- "FPKM"
      spike_df$numMolecules <- spike_df[, mixture_name] * 
        (volume * 10^(-3) * 1/dilution * 10^(-18) * 6.02214129 * 
           10^(23))
      spike_df$rounded_numMolecules <- round(spike_df$numMolecules)
      if (is.null(valid_ids)) 
        spike_df <- subset(spike_df, FPKM >= 1e-10)
      else {
        spike_df <- spike_df[valid_ids, ]
        spike_df <- subset(spike_df, FPKM >= 1e-10)
      }
      spike_df$log_fpkm <- log10(spike_df$FPKM)
      spike_df$log_numMolecules <- log10(spike_df$numMolecules)
      molModel <- tryCatch({
        molModel <- MASS::rlm(log_numMolecules ~ log_fpkm, 
                        data = spike_df)
        molModel
      }, error = function(e) {
        print(e)
        NULL
      })
      molModel
    }, ERCC_annotation, valid_ids)
    if (verbose) 
      message("Apply the fitted robust linear regression model to recovery the absolute copy number for all transcripts each cell...")
    norm_fpkms <- mcmapply(function(cell_exprs, molModel) {
      tryCatch({
        norm_df <- data.frame(log_fpkm = log10(cell_exprs))
        res <- 10^predict(molModel, type = "response", 
                          newdata = norm_df)
      }, error = function(e) {
        rep(NA, length(cell_exprs))
      })
    }, split(as.matrix(relative_expr_matrix), rep(1:ncol(relative_expr_matrix), 
                                                  each = nrow(relative_expr_matrix))), molModels, mc.cores = cores)
    k_b_solution <- data.frame(b = unlist(lapply(molModels, 
                                                 FUN = function(x) {
                                                   intercept = x$coefficients[1]
                                                 })), k = unlist(lapply(molModels, FUN = function(x) {
                                                   slope = x$coefficients[2]
                                                 })))
    kb_model <- MASS::rlm(b ~ k, data = k_b_solution)
    kb_slope <- kb_model$coefficients[2]
    kb_intercept <- kb_model$coefficients[1]
    if (return_all == T) {
      return(list(norm_cds = norm_fpkms, kb_slope = kb_slope, kb_intercept = kb_intercept, 
                  k_b_solution = k_b_solution))
    }
    norm_fpkms
  }
  else {
      formula_all_variables <- all.vars(as.formula(modelFormulaStr))
      
      names(t_estimate) <- colnames(relative_expr_matrix)
          
      pd <- pData(relative_cds)
      pd$Cell = row.names(pd) #use pData instead of the merged data.frame
      norm_cds_list <- plyr::dlply(pd, c(formula_all_variables), function(x){
        relative_expr_matrix_subsets <- relative_expr_matrix[, x$Cell]
        
        split_relative_exprs <- split(as.matrix(relative_expr_matrix_subsets), col(relative_expr_matrix_subsets, as.factor = T))
 
        #calibrate the mode/mc by groups: 
        if(length(expected_total_mRNAs) == 1)
          expected_total_mRNAs <- rep(expected_total_mRNAs, length(split_relative_exprs))
        if(length(expected_capture_rate) == 1)
          expected_capture_rate <- rep(expected_capture_rate, length(split_relative_exprs))
        if(length(reads_per_cell) == 1)
          reads_per_cell <- rep(reads_per_cell, length(split_relative_exprs))

        calibrated_mc <- NULL 
        calibrated_modes_df <- NULL

        ladder_df <- subset(spike_df, mixture_name > detection_threshold) 
        ladder_df$numMolecules <- ladder_df[, mixture_name] * 
          (volume * 10^(-3) * 1/dilution * 10^(-18) * 6.02214129 * 
             10^(23))
        ladder_df$rounded_numMolecules <- round(ladder_df$numMolecules)
        
        calibrated_modes <- mclapply(1:length(split_relative_exprs), 
                       calibrate_mode, 
                       tpm_distribution = split_relative_exprs, 
                       ladder = ladder_df$numMolecules, 
                       total_ladder_transcripts = sum(ladder_df$numMolecules),
                       total_mRNA = expected_total_mRNAs, 
                       capture_rate = expected_capture_rate,
                       reads = reads_per_cell,
                       trials = calculation_trials, 
                       mc.cores = cores)
        calibrated_modes_df <- do.call(rbind.data.frame, calibrated_modes)
        if(verbose)
          message('Calibrating mean total_mRNAs is ...')

        if(calibrate_total_mRNA) {
          num_gene_expressed <- apply(relative_expr_matrix, 2, function(x) sum(x > 1))
          mean_relative_expression <- apply(relative_expr_matrix, 2, function(x) mean((x[x > 1])))
          expected_total_mRNAs <- mean(num_gene_expressed * (mean_relative_expression / t_estimate) * calibrated_modes_df$mean_hypotetical_mode)
          expected_total_mRNAs <- rep(expected_total_mRNAs, length(split_relative_exprs))
        }

        if(verbose)
          message(paste('The calibrated mean total_mRNAs is', expected_total_mRNAs[1]))

        if(is.null(kb_slope) || is.null(kb_intercept)){
          # expected_mRNA_mode <- ceiling(calibrated_modes)
          calibrated_mc <- coef(MASS::rlm(b ~ k, data = calibrated_modes_df))

          kb_slope <- calibrated_mc[2]

          if(is.null(kb_intercept))
            kb_intercept <- calibrated_mc[1]
        }
      
        if (is.null(kb_slope_rng)){
          if(kb_slope > 0)
            kb_slope_rng = c(0.2 * kb_slope, 2 * kb_slope) #note that m is a negative value
          else 
            kb_slope_rng = c(2 * kb_slope, 0.2 * kb_slope) 
        }
        
        if (is.null(kb_intercept_rng)){
          if(kb_intercept > 0)
            kb_intercept_rng = c(0.75 * kb_intercept, 1.25 * kb_intercept)
          else
            kb_intercept_rng = c(1.25 * kb_intercept, 0.75 * kb_intercept)
        }

        if(verbose){
          message('the m/c values are: ')
          print (paste(kb_slope, kb_intercept, sep = ', '))
          message('the range for m is: ')
          print (kb_slope_rng)          
          message('the range for c is: ')
          print (kb_intercept_rng)
        }
        
        if (is.null(expected_mRNA_mode)){
            expected_mRNA_mode <- calibrated_modes_df$mean_hypotetical_mode

            if(verbose){
              message('the calibrated modes are: ')
              print (expected_mRNA_mode)
            }
        }
          
        t_estimate_subset <- t_estimate[colnames(relative_expr_matrix_subsets)]
        if (verbose)
          message("optimizating mc values...")
        for (optim_iter in 1:optim_num) {
            if (verbose)
                message(paste("optimization cycle", optim_iter,
                "..."))
           # only optimize both m and c if the user provided us with a value for c
           # otherwise just use the fixed c, which is easy to calibrate correctly
            if (use_fixed_intercept == FALSE && kb_intercept_rng[1] != kb_intercept_rng[2]) {
                if (verbose)
                  message("optimization m and c values (NOTE that c range should not be huge)")
               if(length(expected_total_mRNAs) == 1)
                expected_total_mRNAs <- rep(expected_total_mRNAs, nrow(relative_expr_matrix))

                optim_para <- optim(par = c(kb_slope=kb_slope, kb_intercept=kb_intercept), optim_mc_func_fix_c,
                  gr = NULL, t_estimate = t_estimate_subset,
                  verbose = verbose, 
                  alpha = expected_mRNA_mode, 
                  total_RNAs = expected_total_mRNAs, # * expected_capture_rate
                  cores = cores, pseudocnt = 0.01,
                  relative_expr_matrix = relative_expr_matrix_subsets,
                  split_relative_expr_matrix = split_relative_exprs,
                  weight_mode=weight_mode, 
                  weight_relative_expr=weight_relative_expr, 
                  weight_total_rna=weight_total_rna,
                  method = c("L-BFGS-B"), 
                  # maxit = max_iterations,
                  lower = c(kb_slope_rng[1], kb_intercept_rng[1]), 
                  upper = c(kb_slope_rng[2], kb_intercept_rng[2]), 
                  control = list(factr = 1, pgtol = 0.1,trace = 1, ndeps = c(0.001, 0.001)), 
                  hessian = FALSE)
            }
            else {
                if (verbose){
                  message("optimization m and fix c as discussed in the method")
                  message("current m value before optimization is, ", kb_slope)
                  message("fixed c value before optimization is, ", kb_intercept)
                }
                optim_para <- optim(par = c(kb_slope = kb_slope), optim_mc_func_fix_c,
                  gr = NULL, kb_intercept = kb_intercept, t_estimate = t_estimate_subset,
                  alpha = expected_mRNA_mode, 
                  total_RNAs = expected_total_mRNAs, #* expected_capture_rate
                  cores = cores,
                  weight_mode=weight_mode, 
                  weight_relative_expr=weight_relative_expr, 
                  weight_total_rna=weight_total_rna,
                  verbose=verbose,
                  pseudocnt = 0.01, relative_expr_matrix = relative_expr_matrix_subsets,
                  split_relative_expr_matrix = split_relative_exprs,
                  method = c("Brent"), 
                  # maxit = max_iterations,
                  lower = c(kb_slope_rng[1]), 
                  upper = c(kb_slope_rng[2]), 
                  control = list(reltol=1e-1, trace = 1), ndeps = c(0.001), 
                  hessian = FALSE)

            }
            if (verbose)
                message("optimization is done!")
            kb_slope <- optim_para$par[1]
            if (use_fixed_intercept == FALSE && kb_intercept_rng[1] != kb_intercept_rng[2])
              kb_intercept <- optim_para$par[2]
            
            if (verbose){
              message("current m value after optimization is, ", kb_slope)
              message("current c value after optimization is, ", kb_intercept)
            }
            
            total_rna_df <- data.frame(Cell = colnames(relative_expr_matrix_subsets),
                t_estimate = t_estimate_subset, alpha_v = expected_mRNA_mode)
            if (verbose)
                message("Estimating the slope and intercept for the linear regression between relative expression value and copy number...")
            # save(file = 'debug_relative2abs', total_rna_df, split_relative_exprs)
            k_b_solution <- plyr::ddply(total_rna_df, .(Cell),
              function(x) {
                  a_matrix <- matrix(c(log10(x[, "t_estimate"]),
                  1, kb_slope, -1), ncol = 2, nrow = 2, byrow = T)
                  colnames(a_matrix) <- c("k", "b")
                  b_matrix <- matrix(c(log10(x[, "alpha_v"]), -kb_intercept), nrow = 2, byrow = T)
                  k_b_solution <- t(solve(a_matrix, b_matrix))
              })
            # print (k_b_solution)
            rownames(k_b_solution) <- k_b_solution$Cell
            k_b_solution <- t(k_b_solution[, c(2, 3)])
            split_kb <- split(k_b_solution, col(k_b_solution,
                as.factor = T))
            if (verbose)
              message("Apply the estimated linear regression model to recovery the absolute copy number for all transcripts each cell...")
            adj_split_relative_expr <- mcmapply(norm_kb, split_kb,
                split_relative_exprs, mc.cores = cores)
            total_rna_df$estimate_k <- k_b_solution[1, ]
            total_rna_df$estimate_b <- k_b_solution[2, ]
            norm_cds <- adj_split_relative_expr
            row.names(norm_cds) <- row.names(relative_expr_matrix_subsets)
            colnames(norm_cds) <- colnames(relative_expr_matrix_subsets)
            t_estimate_subset <- 10^(-(kb_slope + kb_intercept/total_rna_df$estimate_k))
            alpha_v <- estimate_t(norm_cds)
            total_RNAs <- apply(norm_cds, 2, sum)
        }

        return(list(norm_cds = norm_cds, kb_slope = kb_slope, kb_intercept = kb_intercept, k_b_solution = k_b_solution, 
                    expected_mRNA_mode = expected_mRNA_mode, calibrated_mc = calibrated_mc, calibrated_modes_df = calibrated_modes_df)) #
      }
      )
      norm_cds <- do.call(cbind.data.frame, lapply(norm_cds_list, function(x) x$norm_cds))
      colnames(norm_cds) <- as.character(unlist(lapply(norm_cds_list, function(x) colnames(x$norm_cds)))) #set colnames
      norm_cds <- norm_cds[, colnames(relative_cds)] #switch back to the original order

      if (verbose)
        message("Return results...")
      if (return_all == T) {
        kb_slope_vec <- do.call(cbind.data.frame, lapply(norm_cds_list, function(x) x$kb_slope))
        kb_intercept_vec <- do.call(cbind.data.frame, lapply(norm_cds_list, function(x) x$kb_intercept))
        
        k_b_solution <- do.call(cbind.data.frame, lapply(norm_cds_list, function(x) x$k_b_solution))
        colnames(k_b_solution) <- as.character(unlist(lapply(norm_cds_list, function(x) colnames(x$k_b_solution)))) #colnames
        norm_cds <- norm_cds[, colnames(relative_cds)]
        
        # return all estimated values
        expected_mRNA_mode <- NULL
        if(!is.null(norm_cds_list[[1]]$expected_mRNA_mode)){
          expected_mRNA_mode <- unlist(lapply(norm_cds_list, function(x) x$expected_mRNA_mode))
          names(expected_mRNA_mode) <- colnames(k_b_solution)
        }

        calibrated_mc <- NULL
        if(!is.null(norm_cds_list[[1]]$calibrated_mc))
          calibrated_mc <- do.call(cbind.data.frame, lapply(norm_cds_list, function(x) x$calibrated_mc))

        calibrated_modes_df <- NULL
        if(!is.null(norm_cds_list[[1]]$calibrated_modes_df)){
          calibrated_modes_df <- do.call(rbind.data.frame, lapply(norm_cds_list, function(x) x$calibrated_modes_df))
          rownames(calibrated_modes_df) <- colnames(k_b_solution) #colnames
        }
        
        calibrated_total_mRNAs <- NULL
        if(calibrate_total_mRNA)
          calibrated_total_mRNAs <- expected_total_mRNAs

        return(list(norm_cds = norm_cds, kb_slope = t(kb_slope_vec), kb_intercept = kb_intercept_vec, k_b_solution = k_b_solution, 
          expected_mRNA_mode = expected_mRNA_mode, calibrated_mc = calibrated_mc, calibrated_total_mRNAs = calibrated_total_mRNAs,
          calibrated_modes_df = calibrated_modes_df))
    }
    norm_cds
  }
}
#' Spike-in transcripts data.
#'
#' A dataset containing the information for the 92 ERCC spikein transcripts (This dataset is based on the data from the
#'  Nature paper from Stephen Quake group)
#'
#' @format A data frame with 92 rows and 9 variables:
#' \describe{
#'   \item{ERCC_ID}{ID for ERCC transcripts}
#'   \item{subgroup}{Subgroup for ERCC transcript}
#'   \item{conc_attomoles_ul_Mix1}{Contration of Mix 1 (attomoles / ul)}
#'   \item{conc_attomoles_ul_Mix2}{Contration of Mix 2 (attomoles / ul)}
#'   \item{exp_fch_ratio}{expected fold change between mix 1 over mix 2}
#'   \item{numMolecules}{number of molecules calculated from concentration and volume}
#'   \item{rounded_numMolecules}{number in rounded digit of molecules calculated from concentration and volume}
#' }
"spike_df"

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
  
  if(!is.null(pseudocnt)){
    10^dmode(log10(abs_cnt[fpkm > expr_thresh] + pseudocnt)) #keep the original scale
    #k * dmode(log10(fpkm[fpkm > expr_thresh])) + b
  }
  else
    10^dmode(log10(abs_cnt[abs_cnt > expr_thresh]))
}

#linear transform by kb
opt_norm_kb <- function(relative_expr, kb) {
  
  k <- kb[1]
  b <- kb[2]
  tmp <- k * log10(relative_expr) + b
  abs_cnt <- 10^tmp
  
  10^dmode(log10(abs_cnt[abs_cnt > 0]))
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
          relative_expr_matrix = relative_expr_matrix, split_relative_expr_matrix = split_relative_exprs,
          alpha = rep(1, ncol(relative_expr_matrix)), total_RNAs = rep(50000, ncol(relative_expr_matrix)),
          cores = 1, weight_mode=0.33, weight_relative_expr=0.33, weight_total_rna=0.33, verbose = F,  ...) {
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
  
  total_rna_df <- data.frame(Cell = colnames(relative_expr_matrix), t_estimate = t_estimate)
  
  t_k_b_solution <- tryCatch({
    k_b_solution <- plyr::ddply(total_rna_df, .(Cell), function(x) {
      a_matrix <- matrix(c(log10(x[, "t_estimate"]), 1,
                           kb_slope_val, -1), ncol = 2, nrow = 2, byrow = T)
      colnames(a_matrix) <- c("k", "b")
      b_matrix <- matrix(c(0, -kb_intercept_val), nrow = 2, byrow = T)
      k_b_solution <- Matrix::t(solve(a_matrix, b_matrix))
    })
    k_b_solution},
    error = function(e) {print(e); c(NA, NA)}
  )
  
  if(any(is.na(t_k_b_solution)))
    return(NA)
  
  cell_dmode <- tryCatch({
    if(cores > 1){
      cell_dmode <- mcmapply(opt_norm_t, split_t, split_relative_expr_matrix, alpha, kb_slope = kb_slope_val, kb_intercept = kb_intercept_val, pseudocnt = 0.01, mc.cores = cores)
      adj_est_std_cds <- mcmapply(opt_norm_t, split_t, split_relative_expr_matrix, alpha, kb_slope = kb_slope_val, kb_intercept = kb_intercept_val, pseudocnt = 0.01, return_norm = T, mc.cores = cores)
    }
    else {
      cell_dmode <- mapply(opt_norm_t, split_t, split_relative_expr_matrix, alpha, kb_slope = kb_slope_val, kb_intercept = kb_intercept_val, pseudocnt = 0.01)
      adj_est_std_cds <- mapply(opt_norm_t, split_t, split_relative_expr_matrix,  alpha, kb_slope = kb_slope_val, kb_intercept = kb_intercept_val, pseudocnt = 0.01, return_norm = T)
    }
    cell_dmode},
    error = function(e) {print(e); NA}
  )
  
  if(any(is.na(cell_dmode)))
    return(NA)
  
  #adj_est_std_cds <- mapply(opt_norm_t, split_t, split_fpkm, m = m_val, c = c_val, return_norm = T)
  sum_total_cells_rna <- colSums(adj_est_std_cds)
  
  #minimization function:
  #8.
  #dmode_rmse_weight_total <- mean(weight*((cell_dmode - alpha)/alpha)^2 + (1 - weight)*((sum_total_cells_rna -  total_RNAs)/total_RNAs)^2)
  #add the JS distance measure:
  split_relative_expr_matrix <- split(Matrix::t(adj_est_std_cds), 1:ncol(adj_est_std_cds))
  round_split_relative_expr_matrix <- split(Matrix::t(round(adj_est_std_cds)), 1:ncol(adj_est_std_cds))
  
  p_df <- makeprobs(relative_expr_matrix) #relative expression
  p_list <- split(Matrix::t(p_df), 1:ncol(p_df))
  q_df_round <- makeprobs(round(adj_est_std_cds)) #round
  q_df <- makeprobs(adj_est_std_cds) #no rounding
  q_list <- split(Matrix::t(q_df), 1:ncol(q_df))
  q_list_round <- split(Matrix::t(q_df_round), 1:ncol(q_df_round))
  
  dist_divergence <- mcmapply(function(x, y) {
    JSdistVec(x, y)
  }, p_list, q_list, mc.cores = cores)
  
  dist_divergence_round <- mcmapply(function(x, y) {
    JSdistVec(x, y)
  }, p_list, q_list_round, mc.cores = cores)
  
  gm_dist_divergence <- exp(mean(log(dist_divergence)))
  
  total_rna_obj <- mean(((sum_total_cells_rna -  total_RNAs)/total_RNAs)^2)
  mode_obj <- mean(((cell_dmode - alpha)/alpha)^2)
  relative_expr_obj <- gm_dist_divergence
  
  res <- weight_mode * mode_obj + weight_relative_expr * relative_expr_obj + weight_total_rna * total_rna_obj
  
  # if(add_kl_divergence)
  #   res <- 0.25 * log10(dmode_rmse_weight_total + 1) + 0.75 * gm_dist_divergence
  # else
  #   res <- log10(dmode_rmse_weight_total + 1)
  # 
  # #use the algorithm:
  # if(add_kl_divergence)
  #   res <- (weight * (mean(((cell_dmode - alpha)/alpha)^2) - 0)) + (1 - weight) * (gm_dist_divergence - 0.0) + dmode_rmse_weight_total
  # else
  #   res <- log10(dmode_rmse_weight_total + 1)
  
  if(verbose){
    message('current m, c values are ', paste(kb_slope_val, kb_intercept_val, sep = ', '))
    message('total_rna_obj is ', total_rna_obj)
    message('mode_obj is ', mode_obj)
    message('relative_expr_obj is ', relative_expr_obj)
    message('modes are:')
    print (cell_dmode)
    message('target modes are:')
    print (alpha)
    message('mode delta is:')
    print (cell_dmode - alpha)
    message('total RNA delta is:')
    print (sum_total_cells_rna -  total_RNAs)
  }
  #   return(list(m = m_val, c = c_val, dmode_rmse_weight_total = dmode_rmse_weight_total, gm_dist_divergence = gm_dist_divergence, dist_divergence_round = dist_divergence_round,
  #               cell_dmode = cell_dmode, t_k_b_solution = t_k_b_solution, sum_total_cells_rna = sum_total_cells_rna, optim_res = res))
  #
  if(is.finite(res))
    return(res)
  else
    return(10)
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
  #peak finder (using mixture Gauissan model for the FPKM distribution fitting, choose the minial peaker as default)
  
  #apply each column
  apply(relative_expr_matrix, 2, function(relative_expr) 10^dmode(log10(relative_expr[relative_expr > relative_expr_thresh]))) #best coverage estimate}
}


#' Make an educated guess on the spike-based slope regression parameters
#' @importFrom plyr ldply
#' @importFrom MASS rlm
calibrate_mc <- function(total_mRNA, capture_rate, ladder, total_ladder_transcripts, reads, trials=100){
  kb_df <- ldply (seq(0,trials, by=1), function(i){
    #print (i)
    
    hypothetical_ladder <- rmultinom(1, (total_ladder_transcripts * capture_rate), (ladder / sum(ladder)))
    #print (hypothetical_ladder)
    asympt_proportions <- hypothetical_ladder / (sum(hypothetical_ladder) + (total_mRNA * capture_rate))
    asympt_proportions <- c(asympt_proportions, 1 - sum(asympt_proportions))
    
    trial_reads <- rmultinom(1, reads, asympt_proportions)
    trial_tpm <- 1e6 * trial_reads / sum(trial_reads)
    
    #hypothetical_ladder <- hypothetical_ladder[hypothetical_ladder > 0]
    #hypothetical_ladder_tpm <- 1e6 * hypothetical_ladder / ((total_ladder_transcripts + total_mRNA) * capture_rate)
    hypothetical_ladder_tpm <- trial_tpm[1:length(hypothetical_ladder)]
    #print (hypothetical_ladder_tpm)
    
    ladder_df <- data.frame(hypothetical_ladder_tpm=hypothetical_ladder_tpm, 
                            hypothetical_ladder=hypothetical_ladder)
    ladder_df <- subset(ladder_df, hypothetical_ladder_tpm > 0 & hypothetical_ladder > 0)
    
    ladder_reg <- rlm (log10(hypothetical_ladder) ~ log10(hypothetical_ladder_tpm), data=ladder_df)
    #print (summary(ladder_reg))
    b <- coefficients(ladder_reg)[1]
    k <- coefficients(ladder_reg)[2]
    #print (qplot(hypothetical_ladder_tpm, hypothetical_ladder, log="xy") + geom_abline(slope=k, intercept=b))
    data.frame(k=k,b=b)
  })
  #print (kb_df)
  
  kb_reg <- rlm (b ~ k, data=kb_df)
  #print (summary(kb_reg))
  return (list(m=coefficients(kb_reg)[2], c=coefficients(kb_reg)[1]))
}

calibrate_mode <- function(tpm_distribution, total_mRNA, capture_rate, reads, trials=100){
  mode_df <- ldply (seq(0,trials, by=1), function(i){
    proportions <- tpm_distribution / sum(tpm_distribution)
    
    hypothetical_endo <- rmultinom(1, (total_mRNA * capture_rate), proportions)
    hypothetical_proportions <- hypothetical_endo / sum(hypothetical_endo)
    hypothetical_library <- rmultinom(1, reads, hypothetical_proportions)
    
    hypothetical_mode <- 10^dmode(log10(hypothetical_endo[hypothetical_library > 0]))
    data.frame(hypothetical_mode=hypothetical_mode)
  })
  return(mean(mode_df$hypothetical_mode))
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
#' @param total_RNAs the guess of total transcript counts in a single cell. Default is 50000. 
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
  ERCC_controls = NULL, 
  ERCC_annotation = NULL, 
  volume = 10, 
  dilution = 40000, 
  mixture_type = 1,
  detection_threshold = 800, 
  expected_mRNA_mode = NULL, 
  expected_total_mRNAs = 500000, 
  reads_per_cell = 1e6,
  expected_capture_rate = 0.1,
  weight_mode=0.33, 
  weight_relative_expr=0.33, 
  weight_total_rna=0.33,
  verbose = FALSE, 
  return_all = FALSE, 
  cores = 1, 
  optim_num = 1) {
  relative_expr_matrix <- exprs(relative_cds)

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
      spike_df$rounded_numMolecules <- round(spike_df[, 
                                                      mixture_name] * (volume * 10^(-3) * 1/dilution * 
                                                                         10^(-18) * 6.02214129 * 10^(23)))
      if (is.null(valid_ids)) 
        spike_df <- subset(spike_df, FPKM >= 1e-10)
      else {
        spike_df <- spike_df[valid_ids, ]
        spike_df <- subset(spike_df, FPKM >= 1e-10)
      }
      spike_df$log_fpkm <- log10(spike_df$FPKM)
      spike_df$log_numMolecules <- log10(spike_df$numMolecules)
      molModel <- tryCatch({
        molModel <- rlm(log_numMolecules ~ log_fpkm, 
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
    kb_model <- rlm(b ~ k, data = k_b_solution)
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
      
      user_provided_kb_intercept = kb_intercept
      
      if(is.null(kb_slope) || is.null(kb_intercept)){
        ladder_df <- subset(spike_df, mixture_name > detection_threshold) 
        ladder_df$numMolecules <- ladder_df[, mixture_name] * 
          (volume * 10^(-3) * 1/dilution * 10^(-18) * 6.02214129 * 
             10^(23))
        ladder_df$rounded_numMolecules <- round(ladder_df[, 
                                                        mixture_name] * (volume * 10^(-3) * 1/dilution * 
                                                                           10^(-18) * 6.02214129 * 10^(23)))
        calibrated_mc <- calibrate_mc(expected_total_mRNAs, 
                                      expected_capture_rate, 
                                      ladder_df$numMolecules, 
                                      sum(ladder_df$numMolecules), 
                                      reads_per_cell,
                                      trials=100) 
        kb_slope <- calibrated_mc[["m"]]

        if(is.null(kb_intercept))
          kb_intercept <- calibrated_mc[["c"]]
      }
      
      if (is.null(kb_slope_rng)){
        kb_slope_rng = c(1.2 * kb_slope, 0.8 * kb_slope) #note that m is a negative value
      }
      
      if (is.null(kb_intercept_rng)){
        kb_intercept_rng = c(0.8 * kb_intercept, 1.2 * kb_intercept)
      }
      
      pd <- pData(relative_cds)
      pd$Cell = row.names(pd) #use pData instead of the merged data.frame
      norm_cds_list <- plyr::dlply(pd, c(formula_all_variables), function(x){
          relative_expr_matrix_subsets <- relative_expr_matrix[, x$Cell]
          
          split_relative_exprs <- split(as.matrix(relative_expr_matrix_subsets), col(relative_expr_matrix_subsets, as.factor = T))
          
          if (is.null(expected_mRNA_mode)){
            calibrated_modes <- lapply(split_relative_exprs, 
                                      calibrate_mode, 
                                      expected_total_mRNAs, 
                                      expected_capture_rate,
                                      reads_per_cell)
            calibrated_modes <- as.vector(unlist(calibrated_modes))
            expected_mRNA_mode <- calibrated_modes
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
              if (is.null(user_provided_kb_intercept) == FALSE && kb_intercept_rng[1] != kb_intercept_rng[2]) {
                  if (verbose)
                    message("optimization m and c values (NOTE that c range should not be huge)")
                  optim_para <- optim(par = c(kb_slope=kb_slope, kb_intercept=kb_intercept), optim_mc_func_fix_c,
                    gr = NULL, t_estimate = t_estimate_subset,
                    verbose = verbose, 
                    alpha = expected_mRNA_mode, 
                    total_RNAs = expected_total_mRNAs * expected_capture_rate,
                    cores = cores, pseudocnt = 0.01,
                    relative_expr_matrix = relative_expr_matrix_subsets,
                    split_relative_expr_matrix = split_relative_exprs,
                    weight_mode=weight_mode, 
                    weight_relative_expr=weight_relative_expr, 
                    weight_total_rna=weight_total_rna,
                    method = c("L-BFGS-B"), 
                    lower = c(kb_slope_rng[1], kb_intercept_rng[1]), 
                    upper = c(kb_slope_rng[2], kb_intercept_rng[2]), 
                    control = list(factr = 1e+12,
                      pgtol = 0.1, trace = 1, ndeps = c(0.001,
                        0.001),
                      factr=1e14), hessian = FALSE)
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
                    total_RNAs = expected_total_mRNAs * expected_capture_rate,
                    cores = cores,
                    weight_mode=weight_mode, 
                    weight_relative_expr=weight_relative_expr, 
                    weight_total_rna=weight_total_rna,
                    verbose=verbose,
                    pseudocnt = 0.01, relative_expr_matrix = relative_expr_matrix_subsets,
                    split_relative_expr_matrix = split_relative_exprs,
                    method = c("Brent"), 
                    lower = c(kb_slope_rng[1]), 
                    upper = c(kb_slope_rng[2]), 
                    control = list(factr = 1e+12, pgtol = 0.1, reltol=1e-1,
                       trace = 1),
                  hessian = FALSE)
              }
              if (verbose)
                  message("optimization is done!")
              kb_slope <- optim_para$par[1]
              if (kb_intercept_rng[1] != kb_intercept_rng[2])
                kb_intercept <- optim_para$par[2]
              
              if (verbose){
                message("current m value after optimization is, ", kb_slope)
                message("current c value after optimization is, ", kb_intercept)
              }
              
              total_rna_df <- data.frame(Cell = colnames(relative_expr_matrix_subsets),
                  t_estimate = t_estimate_subset, alpha_v = expected_mRNA_mode)
              if (verbose)
                  message("Estimating the slope and intercept for the linear regression between relative expression value and copy number...")
              k_b_solution <- plyr::ddply(total_rna_df, .(Cell),
                function(x) {
                    a_matrix <- matrix(c(log10(x[, "t_estimate"]),
                    1, kb_slope, -1), ncol = 2, nrow = 2, byrow = T)
                    colnames(a_matrix) <- c("k", "b")
                    b_matrix <- matrix(c(log10(x[, "alpha_v"]), -kb_intercept), nrow = 2, byrow = T)
                    k_b_solution <- t(solve(a_matrix, b_matrix))
                })
              print (k_b_solution)
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
          return(list(norm_cds = norm_cds, kb_slope = kb_slope, kb_intercept = kb_intercept, k_b_solution = k_b_solution)) #
      }
      )
      
      norm_cds <- do.call(cbind.data.frame, lapply(norm_cds_list, function(x) x$norm_cds))
      colnames(norm_cds) <- as.character(unlist(lapply(norm_cds_list, function(x) colnames(x$norm_cds)))) #set colnames
      norm_cds <- norm_cds[, colnames(relative_cds)] #switch back to the original order

      kb_slope_vec <- do.call(cbind.data.frame, lapply(norm_cds_list, function(x) x$kb_slope))
      kb_intercept_vec <- do.call(cbind.data.frame, lapply(norm_cds_list, function(x) x$kb_intercept))
      
      k_b_solution <- do.call(cbind.data.frame, lapply(norm_cds_list, function(x) x$k_b_solution))
      colnames(k_b_solution) <- as.character(unlist(lapply(norm_cds_list, function(x) colnames(x$k_b_solution)))) #colnames
      norm_cds <- norm_cds[, colnames(relative_cds)]
      
      if (verbose)
        message("Return results...")
      if (return_all == T) {
        return(list(norm_cds = norm_cds, kb_slope = kb_slope_vec, kb_intercept = kb_intercept_vec, k_b_solution = k_b_solution))
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
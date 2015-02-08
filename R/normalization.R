

#use the deconvoluated linear regression parameters to normalize the log_relative_expr
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

#' Find the most abundant relative expression value for each cell
#' @param relative_expr_matrix a matrix of relative expression values for values with each row and column representing genes/isoforms and cells, respectively. Row and column names should be included
#' @param return_all parameter for the intended return results. If setting TRUE, matrix of dmode as well as max mu and min mu of two gaussian distribution mixture will be returned
#' @param relative_expr_thresh Relative expression values below this threshold are considered zero.
#' @return an vector of most abundant relative_expr value corresponding to the transcrpt copy 1. If setting return_all = TRUE, the mode based on gaussian density function and the max or min
#' mode from the mixture gaussian model
#' @details This function estimates the most abundant relative expression value (t^*) using a gaussian kernel density function. It can also optionally output the t^* based on a two gaussian mixture model
#' based on the smsn.mixture from mixsmsn package
#' @export
#' @examples
#' \dontrun{
#' HSMM_fpkm_matrix <- exprs(HSMM)
#' t_estimate = estimate_t(HSMM_fpkm_matrix)
#'}

estimate_t <- function(relative_expr_matrix, return_all = F, relative_expr_thresh = 0.1, ...) {
  #peak finder (using mixture Gauissan model for the FPKM distribution fitting, choose the minial peaker as default)
  
  smsn_mode_test <- function(relative_expr, g = 2, relative_expr_thresh = 0.1) {
    log_relative_exprs <- log10(relative_expr[relative_expr > relative_expr_thresh]) #only look the transcripts above a certain threshold
    sm_2 <- smsn.mix(log_relative_exprs, nu = 3, g = 2, get.init = TRUE, criteria = TRUE, iter.max=1000,calc.im = FALSE, family="Normal")
    #   print (sm_2)
    
    sm_1 <- smsn.mix(log_relative_exprs, nu = 3, g = 1, get.init = TRUE, criteria = TRUE, iter.max=1000, calc.im = FALSE, family="Normal")
    #   print (sm_1)
    
    if (sm_1$aic >= sm_2$aic){
      trick_location_max <- 10^max(sm_2$mu)
      trick_location_min <- 10^min(sm_2$mu) #use the min / max
    }else{
      trick_location_max <- 10^max(sm_1$mu)
      trick_location_min <- 10^min(sm_1$mu)
    }
    
    #best_cov <- 10^sm_analysis$mu[best_location]
    best_cov_dmode <- 10^(dmode(log_relative_exprs))
    
    best_cov_max <- trick_location_max
    best_cov_min <- trick_location_min
    sd_relative_expr <- 0
    #   print (c(best_cov_max, best_cov_min, best_cov_dmode))
    data.frame(best_cov_dmode = best_cov_dmode,
               best_cov_max = best_cov_max,
               best_cov_min = best_cov_min
    )
  }
  
  #apply each column
  if(return_all){
    do.call(rbind, apply(relative_expr_matrix, 2, function(relative_expr) smsn_mode_test(relative_expr, ...)))
  }
  else{
    apply(relative_expr_matrix, 2, function(relative_expr) 10^dmode(log10(relative_expr[relative_expr > relative_expr_thresh]))) #best coverage estimate}
  }
}

#linear transform by t* and m, c
opt_norm_t <- function(t, relative_expr, m, c, return_norm = FALSE) {
  a_matrix <- matrix(c(log10(t), 1, m,
                       -1), ncol = 2, nrow = 2, byrow = T)
  colnames(a_matrix) <- c("k", "b")
  b_matrix <- matrix(c(0, -c), nrow = 2, byrow = T)
  kb <- t(solve(a_matrix, b_matrix))
  
  k <- kb[1]
  b <- kb[2]
  tmp <- k * log10(relative_expr) + b
  abs_cnt <- 10^tmp
  
  if(return_norm) return(abs_cnt)
  10^dmode(log10(abs_cnt[abs_cnt > 0]))
}

#linear transform by kb
opt_norm_kb <- function(relative_expr, kb) {
  
  k <- kb[1]
  b <- kb[2]
  tmp <- k * log10(relative_expr) + b
  abs_cnt <- 10^tmp
  
  10^dmode(log10(abs_cnt[abs_cnt > 0]))
}

#rmse between the dmode from t estimate based linear transformation and the spike-dmode
t_rmse_abs_cnt <- function (par, t_estimate, relative_expr_mat, split_relative_expr, alpha = 1, cores = 1, verbose = T, ...) {
  cell_num <- ncol(relative_expr_mat)
  #t_estimate <- par[1:cell_num] #t*: the estimates for the best coverage
  names(t_estimate) <- colnames(relative_expr_mat)
  split_t <- split(t(t_estimate), col(as.matrix(t(t_estimate)), as.factor = T))
  
  if(verbose)
    print(paste("t_estimate is: ", paste(as.character(t_estimate), sep = '', collapse = ' '), sep = '', collapse = ''))
  
  #mc_guess <- par[(cell_num + 1):(cell_num + 2)] #m, c parameters: b = m k + c
  mc_guess <- par
  
  if(verbose)
    print(paste("mc_guess is", mc_guess[1], mc_guess[2], sep = ' '))
  
  m_val <- mc_guess[1]
  c_val <- mc_guess[2]
  cell_dmode <- tryCatch({
    if(cores > 1){
      cell_dmode <- mcmapply(opt_norm_t, split_t, split_relative_expr, m = m_val, c = c_val, mc.cores = cores)
      adj_est_std_cds <- mcmapply(opt_norm_t, split_t, split_relative_expr, m = m_val, c = c_val, return_norm = T, mc.cores = cores)
    }
    else {
      cell_dmode <- mapply(opt_norm_t, split_t, split_relative_expr, m = m_val, c = c_val)
      adj_est_std_cds <- mapply(opt_norm_t, split_t, split_relative_expr, m = m_val, c = c_val, return_norm = T)
    }
    cell_dmode},
    error = function(e) {print(e); t_estimate} #return what is better?
  )
  
  if(verbose)
    print(paste("cell_dmode is: ", paste(as.character(cell_dmode), sep = '', collapse = ' '), sep = '', collapse = ''))
  
  rmse <- sqrt(mean((cell_dmode - alpha)^2)) #rmse between the estimated cell_dmode and the 1 copy of transcripts
  
  if(verbose)
    print(paste('rmse is:', rmse, sep = ' '))
  
  sum_total_cells_rna <- colSums(adj_est_std_cds)
  #sum_total_cells_rna[is.infinite(sum_total_cells_rna)] <- 7861584 * 2
 
  if(verbose)
    print(paste('sum of all total RNA is', sum_total_cells_rna))
 
  sqrt(mean(((cell_dmode - alpha) / sum_total_cells_rna)^2)) #new 
}

#' Transform relative expression values into absolute transcript counts
#' 
#' Transform a relative expression matrix to absolute transcript matrix based on the decomposed linear regression parameters from most abundant isoform relative expression value.
#' This function takes a relative expression matrix and a vector of estimated most abundant expression value from the isoform-level matrix and transform it into absolute transcript number.
#' It is based on the fact that most isoforms of gene in a single cell only express one copy so that the most abundant isoform FPKM (t^*) will corresponding to 1 copy transcript. The
#' function takes the the vector t^* and then decomposes it into two parameters vectors k^* and b^* (corresponding to each cell) which correspond to the slope and intercept when
#' we perform the robust linear regression for the spikein data. This decomposition is based on an observed relationship between k and b in terms of b = -3.652201 k + 2.263576. The
#' function will then apply a linear transformation to convert the FPKM to estimated absolute transcript counts based on the the k^* and b^*. The function can also apply a global
#' adjustment if setting global_scaling = TRUE. The k*, b* parameters vectors and the global scaling factor can be output in a list format (together with norm_cds) if setting return_all
#' == TRUE
#' @param relative_expr_matrix an matrix of relative expression values for single cell RNA-seq values with each row and column representing genes/isoforms and cells. Row and column names should be included
#' @param t_estimate an vector for the estimated most abundant FPKM value of isoform for a single cell. Estimators based on gene-level relative expression can also give good approximation but estimators
#' based on isoform FPKM will give better results in general
#' @param detecthion_threshold the lowest ERCC concentration used as a sequencing detection limit, by default is 0.01430512 attomole / Ul. Note that by default we use Mix 1 in ERCC spike-in kit. For all other concentrations, please refer to the illumina ERCC spikein USE GUIDE.   
#' @param ERCC_controls the FPKM matrix for each ERCC spike-in transcript in the cells if user wants to perform the transformation based on their spike-in data. Note that the row and column names should match up with the ERCC_annotation and relative_exprs_matrix respectively. 
#' @param ERCC_annotation the ERCC_annotation matrix from illumina USE GUIDE which will be ued for calculating the ERCC transcript copy number for performing the transformation. 
#' @param volume the approximate volume of the lysis chamber (nl). Default is 10
#' @param dilution the dilution of the spikein transcript in the lysis reaction mix. Default is 40, 000The number of spike-in transcripts per single-cell lysis reaction was calculated from
#' @param return_all parameter for the intended return results. If setting TRUE, matrix of m, c, k^*, b^* as well as the transformed absolute cds will be returned
#' in a list format
#' @param verbose: a logic flag to determine whether or not we should print all the optimization details 
#' @return an matrix of absolute count for isoforms or genes after the transformation. For more details on other output, please refer to detail
#' @export
#' @importFrom plyr ddply
#' @examples
#' \dontrun{
#' HSMM_relative_expr_matrix <- exprs(HSMM)
#' HSMM_abs_matrix <- relative2abs(HSMM_relative_expr_matrix, t_estimate = estimate_t(HSMM_relative_expr_matrix))
#'}
relative2abs <- function(relative_expr_matrix, 
                         t_estimate = estimate_t(relative_expr_matrix), 
                         detecthion_threshold = 0.01430512, #c(1.430512e-02, 2.861023e-02, 5.722046e-02, 1.144409e-01, 2.288818e-01, 4.577637e-01, 9.155273e-01, 1.831055e+00, 3.662109e+00, 7.324219e+00, 1.464844e+01, 2.929688e+01, 5.859375e+01, 1.171875e+02, 2.343750e+02, 4.687500e+02, 9.375000e+02, 1.875000e+03, 3.750000e+03, 7.500000e+03) 
                         ERCC_controls = NULL, ERCC_annotation = NULL, volume = 10, dilution = 40000, 
                         return_all = FALSE, 
                         cores = 1, 
                         verbose = TRUE){
  if(detecthion_threshold < 0.01430512 | detecthion_threshold > 7500)
    stop('concentration detection limit should be between 0.01430512 and 7500')
  else 
   mc_id <- round(detecthion_threshold / 0.01430512)

  Cell <- NULL

  if(!is.null(ERCC_controls) | !is.null(ERCC_annotation)){
    if(is.null(ERCC_controls) | is.null(ERCC_annotation))
      stop('If you want to transform the data to copy number with your spikein data, please provide both of ERCC_controls and ERCC_annotation data frame...')
    
    valid_ids <- which(ERCC_annotation[, 'conc_attomoles_ul_Mix1'] > detecthion_threshold)

    if(verbose)
      message('Performing robust linear regression for each cell based on the spikein data...')

    molModels <- apply(ERCC_controls, 2, function(cell_exprs, input.ERCC.annotation, valid_ids) {
      
      spike_df <- input.ERCC.annotation 
      spike_df <- cbind(spike_df, cell_exprs[row.names(spike_df)])
      colnames(spike_df)[length(colnames(spike_df))] <- "FPKM"
      spike_df$numMolecules <- spike_df$conc_attomoles_ul_Mix1*(volume*10^(-3)*1/dilution*10^(-18)*6.02214129*10^(23))
      spike_df$rounded_numMolecules <- round(spike_df$conc_attomoles_ul_Mix1*(volume*10^(-3)*1/dilution*10^(-18)*6.02214129*10^(23)))
      
      if(is.null(valid_ids))
        spike_df <- subset(spike_df, FPKM >= 1e-10)
      else{
        spike_df <- spike_df[valid_ids, ]
        spike_df <- subset(spike_df, FPKM >= 1e-10)
      }
      
      spike_df$log_fpkm <- log10(spike_df$FPKM)
      spike_df$log_numMolecules <- log10(spike_df$numMolecules)
      
      molModel <- tryCatch({
        molModel <- rlm(log_numMolecules ~ log_fpkm, data=spike_df)
        
          #
          #         qp <- qplot(FPKM, numMolecules, data=spike_df, log="xy") +
          #       geom_abline(color="green") +
          # geom_smooth(aes(y=10^log_numMolecules), method="lm") +
          # geom_smooth(aes(y=10^log_numMolecules), method="rlm", color="red") +
          # geom_line(aes(y=10^predict(molModel,type="response")))
          # print(qp)
          molModel
        }, 
        error = function(e) { print(e); NULL })
      molModel
    }, ERCC_annotation, valid_ids)

    if(verbose)
      message('Apply the fitted robust linear regression model to recovery the absolute copy number for all transcripts each cell...')

    norm_fpkms <- mapply(function(cell_exprs, molModel) {
      tryCatch({
        norm_df <- data.frame(log_fpkm=log10(cell_exprs))
        res <- 10^predict(molModel, type="response", newdata=norm_df)
      }, 
      error = function(e) {
        rep(NA, length(cell_exprs))
        })
    }, 
    split(relative_expr_matrix, rep(1:ncol(relative_expr_matrix), each = nrow(standard_cds))), 
    molModels)

    k_b_solution <- data.frame(b = unlist(lapply(molModels, FUN=function(x) { intercept=x$coefficients[1] })), k = unlist(lapply(molModels, FUN=function(x) { slope=x$coefficients[2] })))
    kb_model <- rlm(b ~ k, data = k_b_solution)
    m <- kb_model$coefficients[2]
    c <- kb_model$coefficients[1]
  
   if(return_all == T) { #also return the trick cds for genes and isoform if return_trick_cds is true otherwise only return total_rna_df
      return (list(norm_cds = norm_cds, m = m, c = c, k_b_solution = k_b_solution))
    }
    norm_cds
  }
  else {
    if(verbose)
      message('Estimating the slope and intercept for the linear regression between relative expression value and copy number based on the lowest detection limits...')
    mc_list <- get_mc_list()
    
    m <- mc_list[mc_id, 1]
    c <- mc_list[mc_id, 2]
 
    split_relative_exprs <- split(as.matrix(relative_expr_matrix), col(relative_expr_matrix, as.factor = T)) #ensure the split dataset is matrix
     
    #estimate the t^* by smsn two gaussian model, choose minial peak  1
    names(t_estimate) <- colnames(relative_expr_matrix)
    
    total_rna_df <- data.frame(Cell = colnames(relative_expr_matrix), t_estimate = t_estimate)
    
    #solve k and b for t by matrix formulation (B = Ax)
    k_b_solution <- plyr::ddply(total_rna_df, .(Cell), function(x){
      a_matrix <- matrix(c(log10(x[, "t_estimate"]), 1, m, -1), ncol = 2, nrow = 2, byrow = T)
      
      colnames(a_matrix) <- c("k", "b")
      b_matrix <- matrix(c(0, -c), nrow = 2, byrow = T)
      k_b_solution <- t(solve(a_matrix, b_matrix))
    })
       
    rownames(k_b_solution) <- k_b_solution$Cell
    k_b_solution <- t(k_b_solution[, c(2, 3)]) #ddply give Cell, k, b columns, take the last two
    split_kb <- split(k_b_solution, col(k_b_solution, as.factor =  T))

    if(verbose)
      message('Apply the estimated linear regression model to recovery the absolute copy number for all transcripts each cell...')

    adj_split_relative_expr <- mcmapply(norm_kb, split_kb, split_relative_exprs, mc.cores = cores)
    
    total_rna_df$estimate_k <- k_b_solution[1, ]
    total_rna_df$estimate_b <- k_b_solution[2, ]
    
    norm_cds <- adj_split_relative_expr
    
    row.names(norm_cds) <- row.names(relative_expr_matrix)
    colnames(norm_cds) <- colnames(relative_expr_matrix)
    
    if(verbose)
      message('Return results...')

    if(return_all == T) { #also return the trick cds for genes and isoform if return_trick_cds is true otherwise only return total_rna_df
      return (list(norm_cds = norm_cds, m = m, c = c, k_b_solution = k_b_solution))
    }
    norm_cds
  }
}

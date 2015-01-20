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
t_rmse_abs_cnt <- function (par, t_estimate, relative_expr_mat, split_relative_expr, alpha = 1, cores = 1, ...) {
  cell_num <- ncol(relative_expr_mat)
  #t_estimate <- par[1:cell_num] #t*: the estimates for the best coverage
  names(t_estimate) <- colnames(relative_expr_mat)
  split_t <- split(t(t_estimate), col(as.matrix(t(t_estimate)), as.factor = T))
  print(paste("t_estimate is: ", paste(as.character(t_estimate), sep = '', collapse = ' '), sep = '', collapse = ''))
  
  #mc_guess <- par[(cell_num + 1):(cell_num + 2)] #m, c parameters: b = m k + c
  mc_guess <- par
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
  print(paste("cell_dmode is: ", paste(as.character(cell_dmode), sep = '', collapse = ' '), sep = '', collapse = ''))
  
  rmse <- sqrt(mean((cell_dmode - alpha)^2)) #rmse between the estimated cell_dmode and the 1 copy of transcripts
  print(paste('rmse is:', rmse, sep = ' '))
  
  sum_total_cells_rna <- colSums(adj_est_std_cds)
  #sum_total_cells_rna[is.infinite(sum_total_cells_rna)] <- 7861584 * 2
  print(paste('sum of all total RNA is', sum_total_cells_rna))
  sqrt(mean(((cell_dmode - alpha) / sum_total_cells_rna)^2)) #new 
}

#' Transform relative expression values into absolute transcript counts
#' 
#' Transform a relative expression matrix to absolute transcript matrix based on the decomposed linear regression parameters from most abundant isoform relative expression value.
#'
#' This function takes a relative expression matrix and a vector of estimated most abundant expression value from the isoform-level matrix and transform it into absolute transcript number.
#' It is based on the fact that most isoforms of gene in a single cell only express one copy so that the most abundant isoform FPKM (t^*) will corresponding to 1 copy transcript. The
#' function takes the the vector t^* and then decomposes it into two parameters vectors k^* and b^* (corresponding to each cell) which correspond to the slope and intercept when
#' we perform the robust linear regression for the spikein data. This decomposition is based on an observed relationship between k and b in terms of b = -3.652201 k + 2.263576. The
#' function will then apply a linear transformation to convert the FPKM to estimated absolute transcript counts based on the the k^* and b^*. The function can also apply a global
#' adjustment if setting global_scaling = TRUE. The k*, b* parameters vectors and the global scaling factor can be output in a list format (together with norm_cds) if setting return_all
#' == TRUE
#'
#' @param relative_expr_matrix an matrix of relative expression values for single cell RNA-seq values with each row and column representing genes/isoforms and cells. Row and column names should be included
#' @param t_estimate an vector for the estimated most abundant FPKM value of isoform for a single cell. Estimators based on gene-level relative expression can also give good approximation but estimators
#' based on isoform FPKM will give better results in general
#' @param global_scaling parameter for globaling scaling. Not perform global scaling by default (global_scaling = 1)
#' @param return_all parameter for the intended return results. If setting TRUE, matrix of k^*, b^* and vector of global_scaling as well the transformed absolute cds will be returned
#' in a list format
#' @param total_fragment: vector of total fragment sequenced for each cell, the element should matched with relative_expr_matrix columns
#' @param alpha_v: vector of most abundant transcripts number estimates in each cell, which are dependent on the total_fragment
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
                         total_fragment = 1.5e6, 
                         alpha_v = 1, 
                         m = -3.652201, c = 2.263576, 
                         global_scaling = FALSE, 
                         return_all = FALSE, 
                         cores = 1) 
  {
  Cell <- NULL
  
  #alpha_v <- exp(log10(1.5 * 10e6  / sample_sheet$Mapped.Fragments)) 
  print('optimizing t_estimates and m and c...') #silence the optimization output
  
  split_relative_exprs <- split(as.matrix(relative_expr_matrix), col(relative_expr_matrix, as.factor = T)) #ensure the split dataset is matrix
  
  optim_para <- optim(par = c(m, c), 
                      t_rmse_abs_cnt, 
                      gr = NULL, 
                      t_estimate = t_estimate, 
                      alpha = alpha_v, 
                      cores = cores, 
                      relative_expr_mat = relative_expr_matrix, 
                      split_relative_expr = split_relative_exprs,
                      method = c("L-BFGS-B"),
                      lower = c(rep(as.vector(t_estimate) - 0, 0), -10, 0.1), #search half low or up of the t_estimate
                      upper = c(rep(as.vector(t_estimate) + 0, 0), -0.1, 10), #m, c is between (-0.1 to -10 and 0.1 to 10)
                      control = list(factr = 1e12, pgtol = 1e-3, trace = 1, ndeps = c(1e-3, 1e-3) ), #as.vector(t_estimate) / 1000,
                      hessian = FALSE)
  message('optimization is done!')
  
  #t_estimate <- optim_para$par[1:length(t_estimate)]
  #regression line between b^* = m * k^* + c
  #m <- optim_para$par[length(t_estimate) + 1]
  #c <- optim_para$par[length(t_estimate) + 2]
  m <- optim_para$par[1]
  c <- optim_para$par[2]
  
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
  
  #predict the absolute isoform based on the adjustment
  
  rownames(k_b_solution) <- k_b_solution$Cell
  k_b_solution <- t(k_b_solution[, c(2, 3)]) #ddply give Cell, k, b columns, take the last two
  split_kb <- split(k_b_solution, col(k_b_solution, as.factor =  T))
  
  adj_split_relative_expr <- mcmapply(norm_kb, split_kb, split_relative_exprs, mc.cores = cores)
  
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
    split_gene_cds <- split(adj_split_relative_expr, col(adj_split_relative_expr, as.factor = T))
    norm_cds <- mcmapply(function(x, y) x *y, split_gene_cds, split_scaling_factor, mc.cores = num_cores)
  }
  else {
    scaling_factor = 1
    norm_cds <- adj_split_relative_expr
  }
  
  row.names(norm_cds) <- row.names(relative_expr_matrix)
  colnames(norm_cds) <- colnames(relative_expr_matrix)
  
  if(return_all == T) { #also return the trick cds for genes and isoform if return_trick_cds is true otherwise only return total_rna_df
    return (list(norm_cds = norm_cds, k_b_solution = k_b_solution, scaling_factor = scaling_factor, optim_para = optim_para))
  }
  norm_cds
}

library(monocle)
context("relative2abs")

test_that("relative2abs() recovers the absolute transcript counts from the relative abundance", {

	HSMM <- load_HSMM()
	rpc_matrix <- relative2abs(HSMM, cores = detectCores())
	expect_equal(round(as.numeric(rpc_matrix[1, 1:5])), 
             c(7, 0, 2, 0, 1))

	#test relative2abs by using model formula and other non-default parameters
    rpc_matrix_2 <- relative2abs(HSMM,
	  modelFormulaStr = "~Hours", 
	  t_estimate = estimate_t(exprs(HSMM)),
	  m = -3.652201, 
	  c = 2.263576, 
	  m_rng = c(-10, -0.1), 
	  c_rng = c(2.263576, 2.263576), 
	  ERCC_controls = NULL, 
	  ERCC_annotation = NULL, 
	  volume = 10, 
	  dilution = 40000, 
	  mixture_type = 1,
	  detection_threshold = 800, 
	  alpha_v = 1, 
	  total_RNAs = 50000, 
	  weight = 0.01, 
	  verbose = T, 
	  return_all = T, 
	  cores = detectCores(), 
	  optim_num = 2)

	 expect_equal(names(rpc_matrix_2), 
		             c('norm_cds', 'm', 'c', 'k_b_solution'))

})
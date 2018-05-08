def umap(i, j, val, dim, n, n_c, metric, negative_sample_rate, alpha, init, mdist, spread, set_op_mix_ratio, local_connectivity, gamma, bandwidth, angular_rp_forest, verbose):
  import umap
  import numpy
  from scipy.sparse import csc_matrix
  data = csc_matrix((val, (i, j)), shape = dim)
  res = umap.UMAP(n_neighbors = n, 
				  	n_components = n_c, 
				  	metric = metric, 
				  	negative_sample_rate = negative_sample_rate,
				  	alpha = alpha,
				  	init = init,
				  	min_dist = mdist,
				  	spread = spread,
				  	set_op_mix_ratio = set_op_mix_ratio,
				  	local_connectivity = local_connectivity,
				  	gamma = gamma,
				  	bandwidth = bandwidth,
				  	angular_rp_forest = angular_rp_forest,
				  	random_state=0,
				  	verbose = verbose).fit(data)
  return res

def umap(i, j, val, dim, n, n_c, metric, n_epochs, negative_sample_rate, alpha, init, mdist, spread, set_op_mix_ratio, local_connectivity, bandwidth, gamma, a, b, random_state, metric_kwds, angular_rp_forest, verbose):
  import umap
  import numpy
  from scipy.sparse import csc_matrix
  data = csc_matrix((val, (i, j)), shape = dim)
  res = umap.UMAP(n_neighbors = n, 
				  	n_components = n_c, 
				  	metric = metric, 
				  	n_epochs = n_epochs, 
				  	negative_sample_rate = negative_sample_rate,
				  	alpha = alpha,
				  	init = init,
				  	min_dist = mdist,
				  	spread = spread,
				  	set_op_mix_ratio = set_op_mix_ratio,
				  	local_connectivity = local_connectivity,
				  	bandwidth = bandwidth, 
				  	gamma = gamma,
				  	a = a, 
				  	b = b, 
				  	random_state = random_state,
				  	metric_kwds = metric_kwds, 
				  	angular_rp_forest = angular_rp_forest,
				  	verbose = verbose).fit(data) 
  return res

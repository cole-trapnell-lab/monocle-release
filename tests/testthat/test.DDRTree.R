library(monocle)
library('R.matlab')
library(igraph)
library(ggplot2)
library(DDRTree)
library(irlba)
library(Rcpp)
library(RcppEigen)
context("DDRTRee")

#note that dim (the inherent dimension the data reduced to and where the tree is constructed) can be changed (currently it is 2)

test_that("DDRTRee() perform the DDRTree construction", {

#ko data:
ko_exprs <- readMat('/Users/xqiu/Dropbox (Personal)/bifurcation_path/simplePPT/data/ko_exprs.mat')
ko_exprs <- as.matrix(ko_exprs$ko.exprs)
storage.mode(ko_exprs) <- 'numeric'
X <- ko_exprs
X <- log2(X + 1)
X <- apply(X, 1, function(x) x - mean(x))
X <- t(X)

pca_res <- pca_projection(X[, ] %*% t(X [, ]), 2)
sqdist_res <- sqdist(X, X)
params <- list(maxIter = 20, eps = 1e-3, dim = 2, lambda = 2435, sigma = 1e-3, gamma = 10)
DDRTree_res <- DDRTree(X = X, params = params, T)
qplot(x = DDRTree_res$Y[1, ], y = DDRTree_res$Y[2, ])

#cpp implementation
system.time(DDRTree_res_cpp <- DDRTree_cpp(X = X, params = params, verbose = T))
qplot(x = DDRTree_res_cpp$Y[1, ], y = DDRTree_res_cpp$Y[2, ])

#ko tree:
ko_DDRTree_res <- DDRTree_res

#lung data:
lung_exprs <- readMat('/Users/xqiu/Dropbox (Personal)/bifurcation_path/simplePPT/data/lung_exprs_mat.mat')
lung_exprs <- as.matrix(lung_exprs$x)
storage.mode(lung_exprs) <- 'numeric'
X <- lung_exprs
X <- log2(X + 1)
X <- apply(X, 1, function(x) x - mean(x))
X <- t(X)

pca_res <- pca_projection(X[, ] %*% t(X [, ]), 2)
sqdist_res <- sqdist(X, X)
params <- list(maxIter = 20, eps = 1e-3, dim = 2, lambda = 2435, sigma = 1e-3, gamma = 10)
DDRTree_res <- DDRTree(X = X, params = params, T)
qplot(x = DDRTree_res$Y[1, ], y = DDRTree_res$Y[2, ])

#lung tree:
lung_DDRTree_res <- DDRTree_res

#golgi data:
golgi_exprs <- readMat('/Users/xqiu/Dropbox (Personal)/bifurcation_path/simplePPT/data/golgi_exprs.mat')
golgi_exprs <- as.matrix(golgi_exprs$golgi.exprs)
storage.mode(golgi_exprs) <- 'numeric'
X <- golgi_exprs
X <- log2(X + 1)
X <- apply(X, 1, function(x) x - mean(x))
X <- t(X)

pca_res <- pca_projection(X[, ] %*% t(X [, ]), 2)
sqdist_res <- sqdist(X, X)
params <- list(maxIter = 20, eps = 1e-3, dim = 2, lambda = 2435, sigma = 1e-3, gamma = 10)
DDRTree_res <- DDRTree_cpp(X = X, lambda = 2435, T)
qplot(x = DDRTree_res$Y[1, ], y = DDRTree_res$Y[2, ])

#golgi tree:
golgi_DDRTree_res <- DDRTree_res

#Cell data:
valid_subset_GSE72857_exprs <- readMat('/Users/xqiu/Dropbox (Personal)/bifurcation_path/simplePPT/data/valid_subset_GSE72857.mat')
valid_subset_GSE72857_exprs <- as.matrix(valid_subset_GSE72857_exprs$valid.subset.GSE72857)
storage.mode(valid_subset_GSE72857_exprs) <- 'numeric'
X <- valid_subset_GSE72857_exprs
X <- log2(X + 1)
X <- apply(X, 1, function(x) x - mean(x))
X <- t(X)

pca_res <- pca_projection(X[, ] %*% t(X [, ]), 2)
sqdist_res <- sqdist(X, X)
params <- list(maxIter = 20, eps = 1e-3, dim = 2, lambda = 5 * ncol(X), sigma = 1e-3, gamma = 10)
DDRTree_res <- DDRTree(X = X, params = params, T)
qplot(x = DDRTree_res$Y[1, ], y = DDRTree_res$Y[2, ])

#cell tree:
valid_subset_GSE72857_DDRTree_res <- DDRTree_res
qplot(x = DDRTree_res$Y[1, ], y = DDRTree_res$Y[2, ], color = MAP_cells_exprs)

#cell tree:
cell_DDRTree_res <- DDRTree_res

#HSMM data:
HSMM_exprs <- readMat('/Users/xqiu/Dropbox (Personal)/bifurcation_path/simplePPT/data/HSMM_exprs.mat')
HSMM_exprs <- as.matrix(HSMM_exprs$HSMM.exprs)
storage.mode(HSMM_exprs) <- 'numeric'
X <- HSMM_exprs
X <- log2(X + 1)
X <- apply(X, 1, function(x) x - mean(x))
X <- t(X)

pca_res <- pca_projection(X[, ] %*% t(X [, ]), 2)
sqdist_res <- sqdist(X, X)
params <- list(maxIter = 20, eps = 1e-3, dim = 2, lambda = 2435, sigma = 1e-3, gamma = 10)
DDRTree_res <- DDRTree_cpp(X, verbose = T)
qplot(x = DDRTree_res$Y[1, ], y = DDRTree_res$Y[2, ])

#HSMM tree:
HSMM_DDRTree_res <- DDRTree_res

load('/Users/xqiu/Dropbox (Personal)/bifurcation_path/simplePPT/data/analysis_shalek_data.RData')

shalek_custom_color_scale_plus_states= c(shalek_custom_color_scale, c('1'='#40A43A', '2'='#CB1B1E', '3'='#3660A5', 'Unstimulated_Replicate.' = 'gray'))

Shalek_golgi_update <- reduceDimension(Shalek_golgi_update, max_components = 2, use_vst = T, use_irlba=F, pseudo_expr = 0, covariates = as.vector(pData(Shalek_golgi_update)$num_genes_expressed) )
Shalek_golgi_update <- orderCells(Shalek_golgi_update, num_paths=2, root_state = NULL)
plot_spanning_tree(Shalek_golgi_update, color_by="interaction(experiment_name, time)", cell_size=2) + 
  scale_color_manual(values=shalek_custom_color_scale_plus_states)


Shalek_abs_subset_ko_LPS <- reduceDimension(Shalek_abs_subset_ko_LPS, max_components = 2, use_vst = T, use_irlba=F, pseudo_expr = 0, covariates = as.vector(pData(Shalek_abs_subset_ko_LPS)$num_genes_expressed) )
Shalek_abs_subset_ko_LPS <- orderCells(Shalek_abs_subset_ko_LPS, num_paths=2, root_state = NULL)
plot_spanning_tree(Shalek_abs_subset_ko_LPS, color_by="interaction(experiment_name, time)", cell_size=2) + 
  scale_color_manual(values=shalek_custom_color_scale_plus_states)

load('/Users/xqiu/Dropbox (Personal)/Infer_GRN/analysis_HSMM_data.RData')

HSMM_myo <- reduceDimension(HSMM_myo, max_components = 2, use_irlba=T, use_vst=T, pseudo_expr=0)
HSMM_myo <- orderCells(HSMM_myo, num_paths=1, root_state = NULL)
plot_spanning_tree(HSMM_myo, color_by="Time", cell_size=2) 

#all genes: 
use_for_ordering_ori <- fData(HSMM_myo)$use_for_ordering
fData(HSMM_myo)$use_for_ordering <- T
HSMM_myo <- reduceDimension(HSMM_myo, max_components = 2, use_irlba=T, use_vst=T, pseudo_expr=0)
HSMM_myo <- orderCells(HSMM_myo, num_paths=1, root_state = NULL)
plot_spanning_tree(HSMM_myo, color_by="Time", cell_size=2) 

#cell paper datasets: 
MAP_cells_clusters <- readMat('/Users/xqiu/Dropbox (Personal)/bifurcation_path/simplePPT/data/MAP_cells_clusters.mat')
MAP_cells_clusters <- MAP_cells_clusters$MAP.cells.clusters
storage.mode(MAP_cells_clusters) <- 'numeric'
MAP_cells_clusters <- read.csv('/Users/xqiu/Downloads/MAP.csv', header = F)

# valid_subset_GSE72857_exprs <- read.table('/Users/xqiu/Downloads/GSE72857_umitab.txt', header = T, row.names = 1)
rownames(valid_subset_GSE72857_exprs) <- paste('g', 1:nrow(valid_subset_GSE72857_exprs), sep = '')
colnames(valid_subset_GSE72857_exprs) <- MAP_cells_clusters$V1

fd <- new("AnnotatedDataFrame", data = data.frame(gene_short_name = row.names(valid_subset_GSE72857_exprs), row.names = row.names(valid_subset_GSE72857_exprs)))
pd <- new("AnnotatedDataFrame", data = data.frame(clusters = MAP_cells_clusters$V2, row.names = MAP_cells_clusters$V1))

# Now, make a new CellDataSet using the RNA counts
valid_subset_GSE72857_exprs <- newCellDataSet(as.matrix(valid_subset_GSE72857_exprs), 
                       phenoData = pd, 
                       featureData = fd,
                       lowerDetectionLimit=1,
                       expressionFamily=negbinomial())

valid_subset_GSE72857_exprs <- estimateSizeFactors(valid_subset_GSE72857_exprs)
valid_subset_GSE72857_exprs <- estimateDispersions(valid_subset_GSE72857_exprs)

options(expressions=500000)
fData(valid_subset_GSE72857_exprs)$use_for_ordering <- T
valid_subset_GSE72857_exprs <- reduceDimension(valid_subset_GSE72857_exprs, max_components = 2, use_irlba=T, use_vst=T)
valid_subset_GSE72857_exprs <- orderCells(valid_subset_GSE72857_exprs, num_paths=1, root_state = NULL)
plot_spanning_tree(valid_subset_GSE72857_exprs, color_by="clusters", cell_size=2) 

})



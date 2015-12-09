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
DDRTree_res <- DDRTree_cpp(X, verbose = T))
qplot(x = DDRTree_res$Y[1, ], y = DDRTree_res$Y[2, ])

#HSMM tree:
HSMM_DDRTree_res <- DDRTree_res

#create a direct graph from the stree:

curr_state <- 1

#stree <- DDRTree_res$stree +  t(DDRTree_res$stree) != 0
stree <- as.matrix(DDRTree_res$stree)
# stree <- stree +  t(stree)
tmp <- stree > 0
stree[tmp == T] <- 1
stree[tmp == F] <- 0

#stree[upper.tri(stree)] <- 0
dimnames(stree) <- list(as.character(1:nrow(stree)), as.character(1:nrow(stree)))
stree_g <- graph.adjacency(stree, mode = "undirected", diag = F, weighted = NULL)

res <- list(subtree = stree_g, root = "6")

load('/Users/xqiu/Dropbox (Personal)/bifurcation_path/simplePPT/data/analysis_shalek_data.RData')

Shalek_golgi_update@reducedDimS[1:2, ] <- DDRTree_res$Y
res <- assignPseudotimePT(Shalek_golgi_update, 'LPS_4h_GolgiPlug_2h_S78_0', plotting = F, scale_pseudotime =
                            F)

qplot(reducedDimS(res)[1, ], reducedDimS(res)[2, ], color = as.character(pData(res)$State), size = pData(res)$Pseudotime)

Shalek_golgi_update <- reduceDimension(Shalek_golgi_update, max_components = 2, use_vst = T, use_irlba=F, pseudo_expr = 0, covariates = as.vector(pData(Shalek_golgi_update)$num_genes_expressed) )
#HSMM_myo <- reduceDimension(HSMM_myo, max_components = 2, use_irlba=T, use_vst=F, pseudo_expr=0.5, fun="exp")

Shalek_golgi_update <- orderCells(Shalek_golgi_update, reverse=F, num_paths=2, root_cell = NULL)
plot_spanning_tree()
load('/Users/xqiu/Dropbox (Personal)/Infer_GRN/analysis_HSMM_data.RData')

HSMM_myo <- reduceDimension(HSMM_myo, max_components = 2, use_irlba=T, use_vst=T, pseudo_expr=0, fun="exp")
#HSMM_myo <- reduceDimension(HSMM_myo, max_components = 2, use_irlba=T, use_vst=F, pseudo_expr=0.5, fun="exp")

HSMM_myo <- orderCells(HSMM_myo, reverse=F, num_paths=1, root_cell = NULL)

HSMM_myo@reducedDimS[1:2, ] <- DDRTree_res$Y
res <- assignPseudotimePT(HSMM_myo, 'T0_CT_D03_0', plotting = F, scale_pseudotime =
                              F)

qplot(reducedDimS(res)[1, ], reducedDimS(res)[2, ], color = as.character(pData(res)$State), size = pData(res)$Pseudotime)


})

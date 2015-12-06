library(monocle)
library('R.matlab')
library(igraph)
library(ggplot2)
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
DDRTree_res <- DDRTree(X = X, params = params, T)
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
DDRTree_res <- DDRTree(X = X, params = params, T)
qplot(x = DDRTree_res$Y[1, ], y = DDRTree_res$Y[2, ])

#HSMM tree: 
HSMM_DDRTree_res <- DDRTree_res

})

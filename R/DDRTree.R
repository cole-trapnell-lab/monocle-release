#' perform PCA projection
#' solve the problem size(C) = NxN, size(W) = NxL
#' max_W trace( W' C W ) : W' W = I
#' @param C a matrix with N x N dimension
#' @param L a dataframe used to generate new data for interpolation of time points
#' @return a matrix (W) with N x L dimension
#' @export
#' 
pca_projection <- function(C, L) {
	eigen_res <- eigen(C)

	U <- eigen_res$vector
	V <- eigen_res$value
	eig_sort <- sort(V, decreasing = T, index.return = T)
	eig_idx <- eig_sort$ix

	W <- U[, eig_idx[1:L]]
}

#' perform PCA projection
#' solve the problem size(C) = NxN, size(W) = NxL
#' max_W trace( W' C W ) : W' W = I
#' @param a a matrix with D x N dimension
#' @param b a matrix with D x N dimension
#' @return a numeric value for the different between a and b
#' @export
#' 
sqdist <- function(a, b) {
	aa <- colSums(a^2)
	bb <- colSums(b^2)
	ab <- t(a) %*% b

	aa_repmat <- matrix(rep(aa, times = ncol(b)), ncol = ncol(b), byrow = F)
	bb_repmat <- matrix(rep(bb, times = ncol(a)), nrow = ncol(a), byrow = T)
	dist <- abs(aa_repmat + bb_repmat - 2 * ab)
}

# X : DxN data matrix
# params.
#       maxIter : maximum iterations
#       eps     : relative objective difference
#       dim     : reduced dimension
#       lambda  : regularization parameter for inverse graph embedding
#       sigma   : bandwidth parameter
#       gamma   : regularization parameter for k-means
#' Perform DDRTree construction
#' @param X a matrix with D x N dimension which is needed to perform DDRTree construction
#' @param params a list with the following parameters: 
#' maxIter : maximum iterations
#' eps     : relative objective difference
#' dim     : reduced dimension
#' lambda  : regularization parameter for inverse graph embedding
#' sigma   : bandwidth parameter
#' gamma   : regularization parameter for k-means
#' @return a list with W, Z, stree, Y, history
#' @export
#' gamma   : regularization parameter for k-means
#' 
DDRTree <- function(X, params, verbose = F) {

	D <- nrow(X) 
	N <- ncol(X)

	#initialization
	W <- pca_projection(X %*% t(X), params$dim)
	Z <- t(W) %*% X

	if(!('ncenter' %in% names(params))) {
		K <- N 
		Y <- Z[, 1:K]
	}
	else {
		K <- params$ncenter
		kmean_res <- kmeans(t(Z), K)
		Y <- kmean_res$centers
		Y <- t(Y)
	}

	#main loop: 
	objs <- c()
	history <- list()
	for(iter in 1:params$maxIter) {

# 		#Kruskal method to find optimal B (use RBGL algorithm: http://stackoverflow.com/questions/16605825/minimum-spaning-tree-with-kruskal-algorithm)
		distsqMU <- sqdist(Y, Y)
# 		#convert with graph packagege to BAM class of graph an calculate mst
# 		mstKruskalBAM <- mstree.kruskal(graphBAM(as.data.frame(distsqMU)))
# 		#build new data frame with resut
# 		stree <- data.frame(cbind(t(mstKruskalBAM$edgeList),
# 		                                 t(mstKruskalBAM$weight)))

	  ##########################use mst from igraph: ##########################
	  g <- graph.adjacency(distsqMU, mode = 'lower', diag = T, weighted = T)
	  g_mst <- mst(g)		
	  stree <- get.adjacency(g_mst, attr = 'weight', type = 'lower')
	  
	  #convert to matrix: 
	  stree <- as.matrix(stree)
	  stree <- stree + t(stree)
		B_tmp <- stree != 0
		B <- B_tmp
		B[B_tmp == FALSE] <- 0
		B[B_tmp == TRUE] <- 1
		L <- diag(colSums(B)) - B

		# #convert back to igraph package
		# stree <- graph.data.frame(mstKruskalDF, directed=FALSE)
		
		#compute R usingmean-shift update rule 
		distZY <- sqdist(Z, Y)
		min_dist <- matrix(rep(apply(distZY, 1, min), times = K), ncol = K, byrow = F)
		tmp_distZY <- distZY - min_dist
		tmp_R <- exp(-tmp_distZY / params$sigma)
		R <- tmp_R / matrix(rep(rowSums(tmp_R), times = K), byrow = F, ncol = K)
		Gamma <- matrix(rep(0, ncol(R) ^ 2), nrow = ncol(R))
		diag(Gamma) <- colSums(R)

		#termination condition 
		obj1 <- - params$sigma * sum(log(rowSums(exp(-tmp_distZY / params$sigma))) 
				- min_dist[, 1] /params$sigma)
		objs[iter] <- (norm(X - W %*% Z, '2'))^2 + params$lambda * sum(diag(Y %*% L %*% t(Y))) + params$gamma * obj1 #sum(diag(A))

		if(verbose)
			message('iter = ', iter, ' ', objs[iter])

		history$W[iter] <- W
		history$Z[iter] <- Z
		history$Y[iter] <- Y
		history$stree[iter] <- stree
		history$R[iter] <- R

		if(iter > 1) {
			if(abs(objs[iter] - objs[iter - 1]) / abs(objs[iter - 1]) < params$eps) {
				break
			}

		}

		#compute low dimension projection matrix
		tmp <- t(solve((((params$gamma + 1) / params$gamma) * ((params$lambda / params$gamma) * L + Gamma) - t(R) %*% R), t(R))) 
		Q <- 1 / (params$gamma + 1) * (diag(1, N) + tmp %*% t(R))
		C <- X %*% Q
		tmp1 <- C %*% t(X)
		W <- pca_projection((tmp1 + t(tmp1)) / 2, params$dim)
		Z <- t(W) %*% C
		Y <- t(solve((params$lambda / params$gamma * L + Gamma), t(Z %*% R)))
	}

	history$objs <- objs
	  
	return(list(W = W, Z = Z, stree = stree, Y = Y, history = history))
}
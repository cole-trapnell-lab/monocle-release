#pca_projection function: 
# X : DxN data matrix
# params.
#       maxIter : maximum iterations
#       eps     : relative objective difference
#       dim     : reduced dimension
#       lambda  : regularization parameter for inverse graph embedding
pca_projection <- function(C, L) {
	eigen_res <- eigen(C)

	U <- eigen_res$vector
	V <- eigen_res$value
	eig_sort <- sort(V, decreasing = T, index.return = T)
	eig_idx <- eig_sort$ix

	W <- U[, eig_idx[1:L]]
}

# a, b : [D, N] = size(a)
# calculate the square distance between a, b
sqdist <- function(a, b) {
	aa <- rowSums(a^2)
	bb <- rowSums(b^2)
	ab <- t(a) * b

	aa_repmat <- matrix(rep(aa, times = ncol(bb)), ncol = ncol(bb), byrow = T)
	bb_repmat <- matrix(rep(bb, times = ncol(aa)), ncol = ncol(aa), byrow = F)
	dist <- abs(aa_repmat + bb_repmat - 2 * a * b)
}

# X : DxN data matrix
# params.
#       maxIter : maximum iterations
#       eps     : relative objective difference
#       dim     : reduced dimension
#       lambda  : regularization parameter for inverse graph embedding
#       sigma   : bandwidth parameter
#       gamma   : regularization parameter for k-means
DDRTree <- function(X, params) {

	D <- nrow(X) 
	N <- ncol(X)

	#initialization
	W <- pca_projection(X * t(X), params$dim)
	Z <- t(W) * W

	if('ncenter' %in% names(params)) {
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
	for(iter in 1:params$maxIter) {

		#Kruskal method to find optimal B (use RBGL algorithm: http://stackoverflow.com/questions/16605825/minimum-spaning-tree-with-kruskal-algorithm)
		distsqMU <- sqdist(Y, Y)
		#convert with graph packagege to BAM class of graph an calculate mst
		mstKruskalBAM <- mstree.kruskal(graphBAM(Y))
		#build new data frame with resut
		stree <- data.frame(cbind(t(mstKruskalBAM$edgeList),
		                                 t(mstKruskalBAM$weight)))
		stree <- stree + t(stree)
		B_tmp <- stree != 0
		B_tmp[B_tmp == 0] <- 0
		B_tmp[B_tmp == 1] <- 1
		L <- diag(colSums(B)) - B

		# #convert back to igraph package
		# stree <- graph.data.frame(mstKruskalDF, directed=FALSE)

		#compute R usingmean-shift update rule 
		distZY <- sqdist(Z, Y)
		min_dist <- matrix(rep(apply(distZY, 2, min), times = K), ncol = K, byrow = T)
		tmp_distZY <- distZY - min_dist
		tmp_R <- exp(-tmp_distZY ./ params$sigma)
		R <- tmp_R ./ matrix(rep(rowSums(tmp_R), times = K), byrow = T)
		Gamma <- diag(sum(R))

		#termination condition 
		obj1 <- - params$sigma * sum(log(rowSums(exp(-tmp_distZY ./ params$sig))) 
				- min_dist[, 1] ./params$sigma)
		objs[iter] <- (norm(X - W * Z))^2 + params$lambda .* trace(Y * L * t(Y)) + params$gamma * obj1

		if(verbose)
			message('iter = ', obj, ' = ', iter, ' ', objs[iter])

		histor$W[iter] <- W
		histor$Z[iter] <- Z
		histor$Y[iter] <- Y
		histor$stree[iter] <- stree
		histor$R[iter] <- R

		if(iter > 1) {
			if(objs[iter] - objs[iter - 1] / abs(objs[iter - 1])) < params$eps {
				break
			}

		}

		#compute low dimension projection matrix
		tmp <- R / (((params$gamma + 1) / params$gamma) .* ((params$lambda / params$gamma) .* L + gamma) - t(R) * R)
		Q <- 1 / (params$gamma + 1) .* (eye(N, N) + tmp * t(R))
		C <- X * Q
		tmp1 <- C * t(X)
		W <- pca_projection((tmp1 + t(tmp1)) ./ 2, params$dim)
		Z <- t(W) * C
		Y <- Z * R / (params$lambda / params$gamma .* L + Gamma)
	}

	history$objs <- objs
}



























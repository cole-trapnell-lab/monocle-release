
#' a function to assign pseudotime for the MST
assign_cell_state_helper <- function(ordering_tree_res, curr_cell, visited_node = curr_cell)
{
    nei <- NULL
    
    cell_tree <- ordering_tree_res$subtree
    V(cell_tree)[curr_cell]$cell_state = curr_state
    
    children <- V(cell_tree) [ suppressWarnings(nei(curr_cell, mode="all")) ]
    children <- setdiff(children, visited_node)
    
    ordering_tree_res$subtree <- cell_tree
    message('curr_cell: ', curr_cell)
    message('children: ', children)
    
    if (length(children) == 1){
        visited_node <- union(children, visited_node)
        
        ordering_tree_res <- assign_cell_state_helper(ordering_tree_res, V(cell_tree)[children]$name, visited_node)
    }else{
        for (child in children)	{
            visited_node <- union(child, visited_node)
            
            curr_state <<- curr_state + 1
            ordering_tree_res <- assign_cell_state_helper(ordering_tree_res, V(cell_tree)[child]$name, visited_node)
        }
    }
    return (ordering_tree_res)
}

#' a function to assign pseudotime and states based on the projected coordinates from the DDRTree algorithm, the stree matrix maybe used later
#' also: fix the bug when the scale_pseudotime can be used
#' @param cds a matrix with N x N dimension
#' @param root_cell a dataframe used to generate new data for interpolation of time points
#' @param scale_pseudotime a matrix with N x N dimension
#' @param verbose a matrix with N x N dimension
#' @return a cds object with the states and the branch assigned correctly
#' @export
#' 
assignPseudotimeBranchPT <- function(cds, root_cell = NULL, scale_pseudotime = F, verbose = F) {
    adjusted_S <- t(cds@reducedDimS)
    dp <- as.matrix(dist(adjusted_S))
    cellPairwiseDistances(cds) <- as.matrix(dist(adjusted_S))
    gp <- graph.adjacency(dp, mode = "undirected", weighted = T)
    dp_mst <- minimum.spanning.tree(gp)
    minSpanningTree(cds) <- dp_mst

    terminal_cell_ids <- V(dp_mst)[which(degree(dp_mst, mode = 'total') == 1)]
    if(verbose) {
        print('the cells on the end of MST: ')
        print((degree(dp_mst, mode = 'total') == 1)[terminal_cell_ids])
    }

    Pseudotime <- rep(0, ncol(cds))
    names(Pseudotime) <- V(dp_mst)

    if(is.null(root_cell))
        root_cell = terminal_cell_ids[1]

    Pseudotime <- shortest.paths(dp_mst, v=root_cell, to=V(dp_mst))
    pData(cds)$Pseudotime <- as.vector(Pseudotime)

    curr_state <- 1

    res <- list(subtree = dp_mst, root = root_cell)
    res <- assign_cell_state_helper(res, res$root)
    pData(cds)$State <- V(res$subtree)[colnames(cds)]$cell_state

    if (scale_pseudotime) {
        cds <- scale_pseudotime(cds)
    }

    cds
}

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
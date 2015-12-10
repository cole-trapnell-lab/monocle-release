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
plot_spanning_tree(Shalek_golgi_update, color_by="interaction(experiment_name, time)", cell_size=1) + 
  scale_color_manual(values=shalek_custom_color_scale_plus_states)

load('/Users/xqiu/Dropbox (Personal)/Infer_GRN/analysis_HSMM_data.RData')

HSMM_myo <- reduceDimension(HSMM_myo, max_components = 2, use_irlba=T, use_vst=T, pseudo_expr=0, fun="exp")
#HSMM_myo <- reduceDimension(HSMM_myo, max_components = 2, use_irlba=T, use_vst=F, pseudo_expr=0.5, fun="exp")

HSMM_myo <- orderCells(HSMM_myo, reverse=F, num_paths=1, root_cell = NULL)

HSMM_myo@reducedDimS[1:2, ] <- DDRTree_res$Y
res <- assignPseudotimePT(HSMM_myo, 'T0_CT_D03_0', plotting = F, scale_pseudotime =
                              F)

qplot(reducedDimS(res)[1, ], reducedDimS(res)[2, ], color = as.character(pData(res)$State), size = pData(res)$Pseudotime)


})

function(cds, 
         x=1, 
         y=2, 
         color_by="State", 
         show_tree=TRUE, 
         show_backbone=TRUE, 
         backbone_color="black", 
         markers=NULL, 
         show_cell_names=FALSE, 
         cell_size=1.5,
         cell_link_size=0.75,
         cell_name_size=2,
         show_all_lineages = F){
  gene_short_name <- NULL
  sample_name <- NULL
  
  #TODO: need to validate cds as ready for this plot (need mst, pseudotime, etc)
  lib_info_with_pseudo <- pData(cds)
  
  #print (lib_info_with_pseudo)
  S_matrix <- reducedDimS(cds)
  
  if (is.null(S_matrix)){
    stop("You must first call reduceDimension() before using this function")
  }
  
  ica_space_df <- data.frame(t(S_matrix[c(x,y),]))
  colnames(ica_space_df) <- c("ICA_dim_1", "ICA_dim_2")
  
  ica_space_df$sample_name <- row.names(ica_space_df)
  ica_space_with_state_df <- merge(ica_space_df, lib_info_with_pseudo, by.x="sample_name", by.y="row.names")
  #print(ica_space_with_state_df)
  dp_mst <- minSpanningTree(cds)
  
  if (is.null(dp_mst)){
    stop("You must first call orderCells() before using this function")
  }
  
  edge_list <- as.data.frame(get.edgelist(dp_mst))
  colnames(edge_list) <- c("source", "target")
  
  edge_df <- merge(ica_space_with_state_df, edge_list, by.x="sample_name", by.y="source", all=TRUE)
  
  edge_df <- plyr::rename(edge_df, c("ICA_dim_1"="source_ICA_dim_1", "ICA_dim_2"="source_ICA_dim_2"))
  edge_df <- merge(edge_df, ica_space_with_state_df[,c("sample_name", "ICA_dim_1", "ICA_dim_2")], by.x="target", by.y="sample_name", all=TRUE)
  edge_df <- plyr::rename(edge_df, c("ICA_dim_1"="target_ICA_dim_1", "ICA_dim_2"="target_ICA_dim_2"))
  
  diam <- as.data.frame(as.vector(V(dp_mst)[get.diameter(dp_mst, weights=NA)]$name))
  colnames(diam) <- c("sample_name")
  diam <- plyr::arrange(merge(ica_space_with_state_df,diam, by.x="sample_name", by.y="sample_name"), Pseudotime)
  
  if(show_all_lineages) {
    pro_state_pseudotime <- diam[as.numeric(diam$State) == min(as.numeric(diam$State)),
                                 "Pseudotime"]
    bifurcation_sample <- diam[which(diam$Pseudotime == max(pro_state_pseudotime)),
                               ]
    bifurcation_sample$State <- diam$State[which(diam$Pseudotime ==
                                                   max(pro_state_pseudotime)) + 1]
    diam <- rbind(diam[1:which(diam$Pseudotime == max(pro_state_pseudotime)),
                       ], bifurcation_sample, diam[(which(diam$Pseudotime ==
                                                            max(pro_state_pseudotime)) + 1):nrow(diam), ])
    no_diam_states <- setdiff(lib_info_with_pseudo$State, lib_info_with_pseudo[diam[,
                                                                                    1], "State"])
    for (state in no_diam_states) {
      state_sample <- ica_space_with_state_df[ica_space_with_state_df$State ==
                                                state, "sample_name"]
      subset_dp_mst <- induced.subgraph(dp_mst, state_sample,
                                        impl = "auto")
      subset_diam <- as.data.frame(as.vector(V(subset_dp_mst)[get.diameter(subset_dp_mst,
                                                                           weights = NA)]$name))
      colnames(subset_diam) <- c("sample_name")
      subset_diam <- plyr::arrange(merge(ica_space_with_state_df,
                                         subset_diam, by.x = "sample_name", by.y = "sample_name"),
                                   Pseudotime)
      subset_diam$State <- state
      bifurcation_sample$State <- state
      diam <- rbind(diam, bifurcation_sample, subset_diam)
    }
  }
  
  markers_exprs <- NULL
  if (is.null(markers) == FALSE){
    markers_fData <- subset(fData(cds), gene_short_name %in% markers)
    if (nrow(markers_fData) >= 1){
      markers_exprs <- reshape2::melt(exprs(cds[row.names(markers_fData),]))
      markers_exprs <- merge(markers_exprs, markers_fData, by.x = "Var1", by.y="row.names")
      #print (head( markers_exprs[is.na(markers_exprs$gene_short_name) == FALSE,]))
      markers_exprs$feature_label <- as.character(markers_exprs$gene_short_name)
      markers_exprs$feature_label[is.na(markers_exprs$feature_label)] <- markers_exprs$Var1
    }
  }
  if (is.null(markers_exprs) == FALSE && nrow(markers_exprs) > 0){
    edge_df <- merge(edge_df, markers_exprs, by.x="sample_name", by.y="Var2")
    #print (head(edge_df))
    g <- ggplot(data=edge_df, aes(x=source_ICA_dim_1, y=source_ICA_dim_2, size=log10(value + 0.1))) + facet_wrap(~feature_label)
  }else{
    g <- ggplot(data=edge_df, aes(x=source_ICA_dim_1, y=source_ICA_dim_2)) 
  }
  if (show_tree){
    g <- g + geom_segment(aes_string(xend="target_ICA_dim_1", yend="target_ICA_dim_2", color=color_by), size=.3, linetype="solid", na.rm=TRUE)
  }
  
  if (is.null(markers_exprs) == FALSE && nrow(markers_exprs) > 0){
    for(i in unique(diam[, 'State'])) { #plot each lineage separately
      g <- g + geom_path(aes(x = ICA_dim_1, y = ICA_dim_2),
                         color = I(backbone_color), size = I(cell_link_size),
                         data = subset(diam, State == i), na.rm = TRUE)
    }
  }
  else {
    g <- g + geom_point(aes_string(color = color_by), size = I(cell_size),
                        na.rm = TRUE)
  }
  
  if (show_backbone){
    #print (diam)
    if(backbone_color %in% colnames(diam))
      g <- g +geom_path(aes(x=ICA_dim_1, y=ICA_dim_2), color=diam[, backbone_color], size=I(cell_link_size), data=diam, na.rm=TRUE)
    else
      g <- g +geom_path(aes(x=ICA_dim_1, y=ICA_dim_2), color=I(backbone_color), size=I(cell_link_size), data=diam, na.rm=TRUE)
  }
  
  if (show_cell_names){
    g <- g +geom_text(aes(label=sample_name), size=cell_name_size)
  }
  g <- g + 
    #scale_color_brewer(palette="Set1") +
    theme(panel.border = element_blank(), axis.line = element_line()) +
    theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) +
    theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) + 
    ylab("Component 1") + xlab("Component 2") +
    theme(legend.position="top", legend.key.height=grid::unit(0.35, "in")) +
    #guides(color = guide_legend(label.position = "top")) +
    theme(legend.key = element_blank()) +
    theme(panel.background = element_rect(fill='white'))
  g
}

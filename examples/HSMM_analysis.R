library(monocle)

# This example performs a greatly simplified version of the single-cell RNA-Seq
# analysis of skeletal myoblast differentiation
# described in Trapnell, Cacchiarelli et al (Nature Biotechnology, 2014).

# Count how many cells each gene is expressed in, and how many
# genes are expressed in each cell
HSMM <- detectGenes(HSMM, min_expr = 0.1)

# Get a list of genes expressed in at least 50 cells
expressed_genes <- row.names(subset(fData(HSMM), num_cells_expressed >= 50))

# Test the above genes for differential expression in response from switch from GM to DM
# Note: this step can take several hours on a single core, so you might want to parallelize it
# with the 'cores' argument 
diff_test_res <- differentialGeneTest(HSMM[expressed_genes,], fullModelFormulaStr="~Media", cores=24)

# Use the differentially expressed genes as the basis for ordering the cells
# by progress through differentiation

# First: collet the list of genes that are significantly differentially expressed (at FDR < 1%)
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))

ordering_genes <- intersect(ordering_genes, row.names(subset(fData(HSMM), biotype %in% c("protein_coding"))))

# Second: mark those genes as the ones used for ordering
HSMM <- setOrderingFilter(HSMM, ordering_genes)

# Third: perform dimensionality reduction using ICA
HSMM <- reduceDimension(HSMM, method="ICA")

# Fourth: compute the minimum spanning tree in the reduced space and use it to order the cells.
# Note that we're allowing a branch with two outcomes in the biological process
HSMM <- orderCells(HSMM, num_paths=2, reverse=T)

# Done! Let's see what the tree and its ordering backbone look like
pdf("mst_plot.pdf")
plot_spanning_tree(HSMM)
dev.off()

# Let's remove all the cells that Monocle tagged as being interstitial 
# mesenchymal cells (e.g. fibroblast contamination) so they don't interfere 
# with downstream analysis
HSMM_filtered <- HSMM[expressed_genes, pData(HSMM)$State != 3]
#HSMM_filtered <- HSMM

# Grab a few specific genes so we can plot their kinetics in pseudotime
cds_subset <- HSMM_filtered[row.names(subset(fData(HSMM_filtered), gene_short_name %in% c("CDK1", "MEF2C", "MYH3"))),]
pdf("muscle_marker_pseudotime_plot.pdf")
plot_genes_in_pseudotime(cds_subset, color_by="State")
dev.off()

# Test the above genes for differential expression across pseudotime
# Note: Again, these steps can take several hours on a single core, so you might want to parallelize it
# with the 'cores' argument 

# We could just run the code below to get differential expression results, but in this 
# example, we'll perform the analysis in pieces, because we'll need some of the model
# information for the gene clustering analysis we'll perform later.
#diff_test_res <- differentialGeneTest(HSMM_filtered, fullModelFormulaStr="~Pseudotime", cores=24)

# Fit the full model for each genefirst
full_model_fits <- fitModel(HSMM_filtered,  modelFormulaStr="~VGAM::bs(Pseudotime)", min_expr = 0.1, cores=24)

# Now fit the reduced models
reduced_model_fits <- fitModel(HSMM_filtered, modelFormulaStr="~1", min_expr = 0.1, cores=24)

# Compare them with an approximate likelihood ratio test
pseudotime_test_res <- compareModels(full_model_fits, reduced_model_fits)

# Merge the test results with the metadata for each gene, including HUGO symbol, etc.
pseudotime_test_res <- merge(fData(HSMM), pseudotime_test_res, by="row.names")

# Collect the response curves in a matrix for use with clustering.
expression_curve_matrix <- responseMatrix(full_model_fits)

# Use k-means clustering to group genes with common pseudotemporal kinetics.
# We're using 6 clusters here, because it's the smallest value that doesn't 
# produce redundant clusters.
clusters <- clusterGenes(expression_curve_matrix, k=6)

pdf("clustering_by_pseudotime.pdf")
plot_clusters(HSMM_filtered, clusters)
dev.off()



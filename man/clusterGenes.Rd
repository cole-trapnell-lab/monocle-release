\name{clusterGenes}
\alias{clusterGenes}
\title{Plots the minimum spanning tree on cells.}
\usage{
clusterGenes(expr_matrix, k, method = function(x) {     as.dist((1 -
  cor(t(x)))/2) }, ...)
}
\arguments{
  \item{expr_matrix}{a matrix of expression values to
  cluster together}

  \item{k}{how many clusters to create}

  \item{method}{the distance function to use during
  clustering}

  \item{...}{extra parameters to pass to pam() during
  clustering}
}
\value{
a pam cluster object
}
\description{
Plots the minimum spanning tree on cells.
}
\examples{
\dontrun{
full_model_fits <- fitModel(HSMM[sample(nrow(fData(HSMM_filtered)), 100),],  modelFormulaStr="expression~VGAM::bs(Pseudotime)")
expression_curve_matrix <- responseMatrix(full_model_fits)
clusters <- clusterGenes(expression_curve_matrix, k=4)
plot_clusters(HSMM_filtered[ordering_genes,], clusters)
}
}


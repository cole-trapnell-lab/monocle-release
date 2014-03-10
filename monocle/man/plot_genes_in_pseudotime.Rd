\name{plot_genes_in_pseudotime}
\alias{plot_genes_in_pseudotime}
\title{Plots expression for one or more genes as a function of pseudotime}
\usage{
plot_genes_in_pseudotime(cds_subset, min_expr = 0.1, cell_size = 0.75,
  nrow = NULL, ncol = 1, panel_order = NULL, color_by = "State")
}
\arguments{
  \item{cds_subset}{CellDataSet for the experiment}

  \item{min_expr}{the minimum (untransformed) expression
  level to use in plotted the genes.}

  \item{cell_size}{the size (in points) of each cell used
  in the plot}

  \item{nrow}{the number of rows used when laying out the
  panels for each gene's expression}

  \item{ncol}{the number of columns used when laying out
  the panels for each gene's expression}

  \item{panel_order}{the order in which genes should be
  layed out (left-to-right, top-to-bottom)}

  \item{color_by}{the cell attribute (e.g. the column of
  pData(cds)) to be used to color each cell}
}
\value{
a ggplot2 plot object
}
\description{
Plots expression for one or more genes as a function of
pseudotime
}
\examples{
my_genes <- row.names(subset(fData(HSMM), gene_short_name \%in\% c("CDK1", "MEF2C", "MYH3")))
cds_subset <- HSMM[my_genes,]
plot_genes_in_pseudotime(cds_subset, color_by="Time")
}


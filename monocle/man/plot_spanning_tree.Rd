\name{plot_spanning_tree}
\alias{plot_spanning_tree}
\title{Plots the minimum spanning tree on cells.}
\usage{
plot_spanning_tree(cds, x = 1, y = 2, color_by = "State", show_tree = T,
  show_backbone = T, backbone_color = "black", marker = NULL)
}
\arguments{
  \item{cds}{CellDataSet for the experiment}

  \item{x}{the column of reducedDimS(cds) to plot on the
  horizontal axis}

  \item{y}{the column of reducedDimS(cds) to plot on the
  vertical axis}

  \item{color_by}{the cell attribute (e.g. the column of
  pData(cds)) to map to each cell's color}

  \item{show_tree}{whether to show the links between cells
  connected in the minimum spanning tree}

  \item{show_backbone}{whether to show the diameter path of
  the MST used to order the cells}

  \item{backbone_color}{the color used to render the
  backbone.}

  \item{marker}{a gene name or gene id to use for setting
  the size of each cell in the plot}
}
\value{
a ggplot2 plot object
}
\description{
Plots the minimum spanning tree on cells.
}
\examples{
\dontrun{
data(HSMM)
plot_spanning_tree(HSMM)
plot_spanning_tree(HSMM, color_by="Pseudotime", show_backbone=F)
plot_spanning_tree(HSMM, marker="MYH3")
}
}


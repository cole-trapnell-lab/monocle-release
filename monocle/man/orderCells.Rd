\name{orderCells}
\alias{orderCells}
\title{Orders cells according to progress through a learned biological process.}
\usage{
orderCells(cds, num_paths = 1, reverse = FALSE)
}
\arguments{
  \item{cds}{the CellDataSet upon which to perform this
  operation}

  \item{num_paths}{the number of end-point cell states to
  allow in the biological process.}

  \item{reverse}{whether to reverse the beginning and end
  points of the learned biological process.}
}
\value{
an updated CellDataSet object, in which phenoData contains
values for State and Pseudotime for each cell
}
\description{
Orders cells according to progress through a learned
biological process.
}


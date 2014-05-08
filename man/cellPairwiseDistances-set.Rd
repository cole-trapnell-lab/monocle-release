\name{cellPairwiseDistances<-}
\alias{cellPairwiseDistances<-}
\title{Sets the matrix containing distances between each pair of cells used by Monocle during cell ordering. Not intended to be called directly.}
\usage{
cellPairwiseDistances(cds) <- value
}
\arguments{
  \item{cds}{A CellDataSet object.}

  \item{value}{a square, symmetric matrix containing
  pairwise distances between cells.}
}
\value{
An updated CellDataSet object
}
\description{
Sets the matrix containing distances between each pair of
cells used by Monocle during cell ordering. Not intended to
be called directly.
}


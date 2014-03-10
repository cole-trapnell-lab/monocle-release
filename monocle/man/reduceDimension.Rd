\name{reduceDimension}
\alias{reduceDimension}
\title{Computes a projection of a CellDataSet object into a lower dimensional space}
\usage{
reduceDimension(cds, max_components = 2, use_irlba = TRUE)
}
\arguments{
  \item{cds}{the CellDataSet upon which to perform this
  operation}

  \item{max_components}{the dimensionality of the reduced
  space}

  \item{use_irlba}{Whether to use the IRLBA package for ICA
  reduction.}
}
\value{
an updated CellDataSet object
}
\description{
Computes a projection of a CellDataSet object into a lower
dimensional space
}
\details{
Currently, Monocle supports dimensionality reduction with
Independent Component Analysis (ICA).
}


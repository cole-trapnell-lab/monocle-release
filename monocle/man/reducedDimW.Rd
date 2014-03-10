\name{reducedDimW}
\alias{reducedDimW}
\title{Retrieves the expression values for each cell (as a matrix) after whitening during independent component analysis.}
\usage{
reducedDimW(cds)
}
\arguments{
  \item{cds}{A CellDataSet object.}
}
\value{
A matrix, where each row is a set of whitened expression
values for a feature and columns are cells.
}
\description{
Retrieves the expression values for each cell (as a matrix)
after whitening during independent component analysis.
}
\examples{
data(HSMM)
W <- reducedDimW(HSMM)
}


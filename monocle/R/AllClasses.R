setOldClass(c("igraph"), prototype=structure(list(), class="igraph"))

#' The CellDataSet class
#'
#' The main class used by Monocle to hold single cell expression data. CellDataSet extends the basic Bioconductor ExpressionSet class.
#' 
#' This class is initialized from a matrix of expression values Methods that operate on CellDataSet objects constitute the basic Monocle workflow.
#'
#'
#'@section Slots: 
#'  \describe{
#'    \item{\code{reducedDimS}:}{Matrix of class \code{"numeric"}, containing the source values computed by Independent Components Analysis.}
#'    \item{\code{reducedDimW}:}{Matrix of class \code{"numeric"}, containing the whitened expression values computed during Independent Components Analysis.}
#'    \item{\code{reducedDimA}:}{Matrix of class \code{"numeric"}, containing the weight values computed by Independent Components Analysis.}
#'    \item{\code{minSpanningTree}:}{Object of class \code{"igraph"}, containing the minimum spanning tree used by Monocle to order cells according to progress through a biological process.}
#'    \item{\code{cellPairwiseDistances}:}{Matrix of class \code{"numeric"}, containing the pairwise distances between cells in the reduced dimension space.}
#'  }
#'
#' @name CellDataSet 
#' @rdname CellDataSet
#' @aliases CellDataSet-class
#' @exportClass CellDataSet
setClass( "CellDataSet", 
          contains = "ExpressionSet",
          slots = c(reducedDimS = "matrix",
                    reducedDimW = "matrix",
                    reducedDimA = "matrix",
                    minSpanningTree="igraph",
                    cellPairwiseDistances="matrix"),
          prototype = prototype( new( "VersionedBiobase",
                                      versions = c( classVersion("ExpressionSet"), CellDataSet = "1.0.0" ) ))
)


#' The HSMM single-cell RNA-Seq timeseries from Trapnell, Cacchiarelli et al. 
#'
#' @name HSMM
#' @docType HSMM
#' @author Cole Trapnell \email{cole@@broadinstitute.org}
#' @references \url{http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52529}
#' @keywords data
NULL


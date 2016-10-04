setOldClass(c("igraph"), prototype=structure(list(), class="igraph"))

#' The CellDataSet class
#'
#' The main class used by Monocle to hold single cell expression data. 
#' CellDataSet extends the basic Bioconductor ExpressionSet class.
#' 
#' This class is initialized from a matrix of expression values Methods that 
#' operate on CellDataSet objects constitute the basic Monocle workflow.
#'
#'
#'@section Slots: 
#'  \describe{
#'    \item{\code{reducedDimS}:}{Matrix of class \code{"numeric"}, containing the source values computed by Independent Components Analysis.}
#'    \item{\code{reducedDimW}:}{Matrix of class \code{"numeric"}, containing the whitened expression values computed during Independent Components Analysis.}
#'    \item{\code{reducedDimA}:}{Matrix of class \code{"numeric"}, containing the weight values computed by Independent Components Analysis.}
#'    \item{\code{reducedDimK}:}{Matrix of class \code{"numeric"}, containing the pre-whitening matrix computed by Independent Components Analysis.}
#'    \item{\code{minSpanningTree}:}{Object of class \code{"igraph"}, containing the minimum spanning tree used by Monocle to order cells according to progress through a biological process.}
#'    \item{\code{cellPairwiseDistances}:}{Matrix of class \code{"numeric"}, containing the pairwise distances between cells in the reduced dimension space.}
#'    \item{\code{expressionFamily}:}{Object of class \code{"vglmff"}, specifying the VGAM family function used for expression responses.}
#'    \item{\code{lowerDetectionLimit}:}{A \code{"numeric"} value specifying the minimum expression level considered to be true expression.}
#'    \item{\code{dispFitInfo}:}{An \code{environment} containing lists, one for each set of estimated dispersion values. See \code{\link{estimateDispersions}}}
#'    \item{\code{dim_reduce_type}:}{A \code{"character"} A string encoding how this CellDataSet has been reduced in dimensionality}
#'    \item{\code{auxOrderingData}:}{A \code{environment} Auxilliary data structures used by various steps in Monocle. Not to be accessed by users directly.}
#'  }
#'
#' @name CellDataSet 
#' @rdname CellDataSet
#' @aliases CellDataSet-class
#' @exportClass CellDataSet
#' @importFrom Biobase ExpressionSet
<<<<<<< 7afe2e5a61fe398816056758e84b77f3d44ceb50
#' @import BiocGenerics
=======
>>>>>>> First commit
#' 
setClass( "CellDataSet", 
          contains = "ExpressionSet",
          slots = c(reducedDimS = "matrix",
                    reducedDimW = "matrix",
                    reducedDimA = "matrix",
                    reducedDimK = "matrix",
                    minSpanningTree="igraph",
                    cellPairwiseDistances="matrix",
                    expressionFamily="vglmff",
                    lowerDetectionLimit="numeric",
                    dispFitInfo = "environment",
                    dim_reduce_type="character",
                    auxOrderingData = "environment"),
          prototype = prototype( new( "VersionedBiobase",
                                      versions = c( classVersion("ExpressionSet"), CellDataSet = "1.1.0" ) ))
)


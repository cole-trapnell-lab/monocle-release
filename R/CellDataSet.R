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
#' @field reducedDimS Matrix of class numeric, containing the source values computed by Independent Components Analysis.
#' @field reducedDimW Matrix of class numeric, containing the whitened expression values computed during Independent Components Analysis.
#' @field reducedDimA Matrix of class numeric, containing the weight values computed by Independent Components Analysis.
#' @field reducedDimK A Matrix of class numeric, containing the pre-whitening matrix computed by Independent Components Analysis.
#' @field minSpanningTree An Object of class igraph, containing the minimum spanning tree used by Monocle to order cells according to progress through a biological process.
#' @field cellPairwiseDistances A Matrix of class numeric, containing the pairwise distances between cells in the reduced dimension space.
#' @field expressionFamily An Object of class vglmff, specifying the VGAM family function used for expression responses.
#' @field lowerDetectionLimit A numeric value specifying the minimum expression level considered to be true expression.
#' @field dispFitInfo An environment containing lists, one for each set of estimated dispersion values. See estimateDispersions.
#' @field dim_reduce_type A string encoding how this CellDataSet has been reduced in dimensionality
#' @field auxOrderingData An environment of auxilliary data structures used by various steps in Monocle. Not to be accessed by users directly.
#' @name CellDataSet 
#' @rdname CellDataSet
#' @aliases CellDataSet-class
#' @exportClass CellDataSet
#' @importFrom Biobase ExpressionSet
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
                    auxOrderingData = "environment", 
                    auxClusteringData = "environment"
                    ),
          prototype = prototype( new( "VersionedBiobase",
                                      versions = c( classVersion("ExpressionSet"), CellDataSet = "1.2.0" ) ))
)


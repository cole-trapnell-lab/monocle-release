#' The CellType class
#'
#' Classifies cells using a criterion function.
#' 
#' Classifies cells via a user-defined
#' gating function. The gating function accepts as input the entire matrix of
#' expression data from a \code{CellDataSet}, and return TRUE or FALSE for each
#' cell in it, depending on whether each meets the criteria in the gating 
#' function
#'
#'@section Slots: 
#'  \describe{
#'    \item{\code{classify_func}:}{A function that accepts a matrix of expression values as input, and returns a logical vector (of length equal to the number of columns in the matrix) as output}
#'  }
#'
#' @name CellType
#' @rdname CellType
#' @aliases CellType-class
#' @exportClass CellType
setClass( "CellType", 
          slots = c(classify_func="function"))


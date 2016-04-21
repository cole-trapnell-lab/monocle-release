setOldClass(c("igraph"), prototype=structure(list(), class="igraph"))

#' The CellTypeHierarchy class
#'
#' Classifies cells according to a hierarchy of types.
#' 
#' Classifies cells according to a hierarchy of types via user-defined
#' gating functions.
#'
#'
#'@section Slots: 
#'  \describe{
#'    \item{\code{classificationTree}:}{Object of class \code{"igraph"}}
#'  }
#'
#' @name CellTypeHierarchy 
#' @rdname CellTypeHierarchy
#' @aliases CellTypeHierarchy-class
#' @exportClass CellTypeHierarchy
#' @import igraph
setClass( "CellTypeHierarchy", 
          slots = c(classificationTree="igraph"))

#' The CellTypeHierarchy class
#'
#' Classifies cells according to a hierarchy of types.
#' 
#' Classifies cells according to a hierarchy of types via user-defined
#' gating functions.
#'
#'
#'@section Slots: 
#'  \describe{
#'    \item{\code{classificationTree}:}{Object of class \code{"igraph"}}
#'  }
#'
#' @name CellType
#' @rdname CellType
#' @aliases CellType-class
#' @exportClass CellType
#' @import igraph
setClass( "CellType", 
          slots = c(classify_func="function"))



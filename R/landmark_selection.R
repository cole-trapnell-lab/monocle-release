#' Downsample datasets with landmark algorithm to ensure full sampling of the entire data space
#'
#' @param cds a cell dataset for landmark selection
#' @param landmark_num number of landmark cells to retrieve
#' @return a list of three named vectors with the same number of cells from the cds: assign (which landmark the current cell closest to), 
#' dist (distance of current cell to its closest landmark) and flag (whether or not the current cell is a landmark, 1: yes, 0: no)
#' @export
landmark_selection <- function(cds, landmark_num) {
  if(landmark_num >= ncol(cds) | landmark_num < 1) {
    stop('number of landmarks should be smaller than the number of cells and larger than 1!')
  }
  
  data_class <- class(cds@assayData$exprs)
  if(data_class == "matrix" | data_class == 'data.frame') {
    cds@assayData$exprs <- as(cds@assayData$exprs, 'sparseMatrix')  
  }
  sp_data <- cds@assayData$exprs
  
  res <- monocle:::select_landmarks(sp_data@x, sp_data@i, sp_data@p, sp_data@Dim[2], sp_data@Dim[1], landmark_num)
  
  return(res)
}

#' Downsample datasets with landmark algorithm to ensure full sampling of the entire data space
#'
#' @param cds a cell dataset after trajectory reconstruction
#' @param landmark_num 
#' @return a new cds containing only the cells selected by the landmark selection algorithm 
#' @export
downsampleCDS <- function(cds, landmark_num) {
  # 
  res <- landmark_selection(cds, landmark_num)
  cds[, which(res$flag == 1)]
}
# this file contains functions to convert between Monocle cds and scran or Seurat object back and forth. 

#' Export a monocle CellDataSet object to other popular single cell analysis toolkit.
#' 
#' This function takes a monocle CellDataSet and converts it to another type of object used in another popular single cell analysis toolkit. It currently
#' supports Scran and Seurat packages.  
#' 
#' @param monocle_cds the Monocle CellDataSet you would like to export into a type used in another package 
#' @param export_to the object type you would like to export to, either Seurat or Scater 
#' @param export_all Whether or not to export all the slots in Monocle and keep in another object type. Default is FALSE (or only keep
#' minimal dataset). If export_all is setted to be true, the original monocle cds will be keeped in the other cds object too. 
#' This argument is also only applicable when export_to is Seurat.  
#' @return a new object in the format of another package, as described in the export_to argument. 
#' @export
#' @examples
#' \dontrun{
#' lung <- load_lung()
#' seurat_lung <- exportCDS(lung)
#' seurat_lung_all <- exportCDS(lung, export_all = T)
#' scater_lung <- exportCDS(lung, export_to = 'Scater')
#' scater_lung_all <- exportCDS(lung, export_to = 'Scater', export_all = T)
#' }
exportCDS <- function(monocle_cds, export_to = c('Seurat', 'Scater'), export_all = FALSE) {
  if(export_to == 'Seurat') {
    requireNamespace("Seurat")
    data <- exprs(monocle_cds)
    ident <- colnames(monocle_cds)
    
    if(export_all) {
      # mist_list <- list(Monocle = list(reducedDimS = monocle_cds@reducedDimS, 
      #                   reducedDimW = monocle_cds@reducedDimW,  
      #                   reducedDimA = monocle_cds@reducedDimA, 
      #                   reducedDimK = monocle_cds@reducedDimK, 
      #                   minSpanningTree = monocle_cds@minSpanningTree, 
      #                   cellPairwiseDistances = monocle_cds@cellPairwiseDistances, 
      #                   expressionFamily = monocle_cds@expressionFamily, 
      #                   dispFitInfo = monocle_cds@dispFitInfo, 
      #                   dim_reduce_type = monocle_cds@dim_reduce_type, 
      #                   auxOrderingData = monocle_cds@auxOrderingData, 
      #                   auxClusteringData = monocle_cds@auxClusteringData,
      #                   experimentData = monocle_cds@experimentData,
      #                   classVersion = monocle_cds@.__classVersion__, 
      #                   annotation = monocle_cds@annotation,
      #                   protocolData = monocle_cds@protocolData,
      #                   featureData = monocle_cds@featureData
      #                   ))
      # clean all conversion related slots 
      monocle_cds@auxClusteringData$seurat <- NULL
      monocle_cds@auxClusteringData$scran <- NULL
      mist_list <- monocle_cds
    } else {
      mist_list <- list()
    }
    if("use_for_ordering" %in% colnames(fData(monocle_cds))) {
      var.gene <- row.names(subset(fData(monocle_cds), use_for_ordering == TRUE)); 
    }
    
    export_cds <- Seurat::CreateSeuratObject(raw.data = data, 
    				  normalization.method = "LogNormalize",
    				  do.scale = TRUE, 
    				  do.center = TRUE,
                      is.expr = monocle_cds@lowerDetectionLimit,
                      project = "exportCDS",
                      # display.progress = FALSE, 
                      meta.data = pData(monocle_cds))
    
    export_cds@misc <- mist_list
    export_cds@meta.data <- pData(monocle_cds)
    
  } else if (export_to == 'Scater') {
    requireNamespace("scater")
    data <- log2(exprs(monocle_cds) + 1)
    pd <- new("AnnotatedDataFrame", data = pData(monocle_cds))
    fd <- new("AnnotatedDataFrame", data = fData(monocle_cds))
    experimentData = monocle_cds@experimentData
    
    export_cds <- SingleCellExperiment::SingleCellExperiment(exprsData = data, countData = NULL, tpmData = NULL,
              fpkmData = NULL, cpmData = NULL, phenoData = pd, featureData = fd,
              experimentData = experimentData, is_exprsData = NULL,
              cellPairwiseDistances = dist(vector()),
              featurePairwiseDistances = dist(vector()), 
              lowerDetectionLimit = monocle_cds@lowerDetectionLimit,
              logExprsOffset = 1)    
  } else {
    stop('the object type you want to export to is not supported yet')
  }
  
  return(export_cds)
}

#' Import a seurat or scatter/scran CellDataSet object and convert it to a monocle cds.
#' 
#' This function takes a monocle CellDataSet and converts it to another type of object used in another popular single cell analysis toolkit. It currently
#' supports Scran and Seurat packages.  
#' 
#' @param otherCDS the object you would like to convert into a monocle cds 
#' @param import_all Whether or not to import all the slots in seurat or scatter. Default is FALSE (or only keep
#' minimal dataset). 
#' @return a new monocle cell dataset object converted from other objects (Scatter or Seurat).  
#' @export
#' @examples
#' \dontrun{
#' lung <- load_lung()
#' seurat_lung <- exportCDS(lung)
#' seurat_lung_all <- exportCDS(lung, export_all = T)
#' scater_lung <- exportCDS(lung, export_to = 'Scater')
#' scater_lung_all <- exportCDS(lung, export_to = 'Scater', export_all = T)
#' 
#' importCDS(seurat_lung)
#' importCDS(seurat_lung, import_all = T)
#' importCDS(seurat_lung_all)
#' importCDS(seurat_lung_all, import_all = T)
#' importCDS(scater_lung)
#' importCDS(scater_lung, import_all = T)
#' importCDS(scater_lung_all)
#' importCDS(scater_lung_all, import_all = T)
#' }
#' @importFrom Matrix as.matrix
#' @importFrom VGAM uninormal
importCDS <- function(otherCDS, import_all = FALSE) {
  if(class(otherCDS)[1] == 'seurat') {
    requireNamespace("Seurat")
    data <- otherCDS@raw.data

    if(class(data) == "data.frame") {
      data <- as(as.matrix(data), "sparseMatrix")
    }
    
    pd <- tryCatch( {
      pd <- new("AnnotatedDataFrame", data = otherCDS@meta.data)
      pd
    }, 
    #warning = function(w) { },
    error = function(e) { 
      pData <- data.frame(cell_id = colnames(data), row.names = colnames(data))
      pd <- new("AnnotatedDataFrame", data = pData)
      
      message("This Seurat object doesn't provide any meta data");
      pd
    })
    
    # remove filtered cells from Seurat
    if(length(setdiff(colnames(data), rownames(pd))) > 0) {
      data <- data[, rownames(pd)]  
    }
    
    fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
    fd <- new("AnnotatedDataFrame", data = fData)
    lowerDetectionLimit <- otherCDS@is.expr
    
    if(all(data == floor(data))) {
      expressionFamily <- negbinomial.size()
    } else if(any(data < 0)){
      expressionFamily <- uninormal()
    } else {
      expressionFamily <- tobit()
    }
    
    valid_data <- data[, row.names(pd)]
    
    monocle_cds <- newCellDataSet(data,
                           phenoData = pd, 
                           featureData = fd,
                           lowerDetectionLimit=lowerDetectionLimit,
                           expressionFamily=expressionFamily)
    
    if(import_all) {
      if("Monocle" %in% names(otherCDS@misc)) {
        # if(slotNames(lung) == )
        # monocle_cds@reducedDimS = otherCDS@misc$Monocle@reducedDimS 
        # monocle_cds@reducedDimW = otherCDS@misc$Monocle@reducedDimW  
        # monocle_cds@reducedDimA = otherCDS@misc$Monocle@reducedDimA 
        # monocle_cds@reducedDimK = otherCDS@misc$Monocle@reducedDimK 
        # monocle_cds@minSpanningTree = otherCDS@misc$Monocle@minSpanningTree 
        # monocle_cds@cellPairwiseDistances = otherCDS@misc$Monocle@cellPairwiseDistances 
        # monocle_cds@expressionFamily = otherCDS@misc$Monocle@expressionFamily 
        # monocle_cds@dispFitInfo = otherCDS@misc$Monocle@dispFitInfo 
        # monocle_cds@dim_reduce_type = otherCDS@misc$Monocle@dim_reduce_type 
        # monocle_cds@auxOrderingData = otherCDS@misc$Monocle@auxOrderingData 
        # monocle_cds@auxClusteringData = otherCDS@misc$Monocle@auxClusteringData
        # monocle_cds@experimentData = otherCDS@misc$Monocle@experimentData
        # monocle_cds@classVersion = otherCDS@misc$Monocle@.__classVersion__ 
        # monocle_cds@annotation = otherCDS@misc$Monocle@annotation
        # monocle_cds@protocolData = otherCDS@misc$Monocle@protocolData
        # monocle_cds@featureData = otherCDS@misc$Monocle@featureData
        
        # clean all conversion related slots 
        otherCDS@misc$Monocle@auxClusteringData$seurat <- NULL
        otherCDS@misc$Monocle@auxClusteringData$scran <- NULL
        
        monocle_cds <- otherCDS@misc$Monocle
        mist_list <- otherCDS
        
      } else {
        # mist_list <- list(ident = ident, 
        #                   project.name = project.name,
        #                   dr = otherCDS@dr,
        #                   assay = otherCDS@assay,
        #                   hvg.info = otherCDS@hvg.info,
        #                   imputed = otherCDS@imputed,
        #                   cell.names = otherCDS@cell.names,
        #                   cluster.tree = otherCDS@cluster.tree,
        #                   snn = otherCDS@snn,
        #                   kmeans = otherCDS@kmeans,
        #                   spatial = otherCDS@spatial,
        #                   misc = otherCDS@misc
        # ) 
        mist_list <- otherCDS
      }
    } else {
      mist_list <- list()
    }
    
    if("var.genes" %in% slotNames(otherCDS)) {
      monocle_cds <- setOrderingFilter(monocle_cds, otherCDS@var.genes)
      
    }
    monocle_cds@auxClusteringData$seurat <- mist_list
    
  } else if (class(otherCDS)[1] == 'SCESet') {
    requireNamespace("scater")
    
    message('Converting the exprs data in log scale back to original scale ...')    
    data <- 2^otherCDS@assayData$exprs - otherCDS@logExprsOffset

    fd <- otherCDS@featureData
    pd <- otherCDS@phenoData
    experimentData = otherCDS@experimentData
    if("is.expr" %in% slotNames(otherCDS))
      lowerDetectionLimit <- otherCDS@is.expr
    else 
      lowerDetectionLimit <- 1
    
    if(all(data == floor(data))) {
      expressionFamily <- negbinomial.size()
    } else if(any(data < 0)){
      expressionFamily <- uninormal()
    } else {
      expressionFamily <- tobit()
    }
    
    if(import_all) {
      # mist_list <- list(iotherCDS@sc3,
      #                   otherCDS@reducedDimension)
      mist_list <- otherCDS 
                        
    } else {
      mist_list <- list()
    }

	monocle_cds <- newCellDataSet(data,
	                   phenoData = pd, 
	                   featureData = fd,
	                   lowerDetectionLimit=lowerDetectionLimit,
	                   expressionFamily=expressionFamily)
    # monocle_cds@auxClusteringData$sc3 <- otherCDS@sc3
    # monocle_cds@auxOrderingData$scran <- mist_list
    
    monocle_cds@auxOrderingData$scran <- mist_list
    
  } else {
    stop('the object type you want to export to is not supported yet')
  }
  
  return(monocle_cds)
}




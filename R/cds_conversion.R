# this file contains functions to convert between Monocle cds to scran or Seurat CDS back and forth. Note that those functionss   

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
#' @importFrom scater newSCESet
#' @import Seurat
#' @export
exportCDS <- function(monocle_cds, export_to = c('Seurat', 'Scater'), export_all = FALSE) {
  if(export_to == 'Seurat') {
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
      misc = list()
    }
    export_cds <- new("seurat", raw.data = data, 
                      data = log(data + 1), 
                      scale.data = t(scale(t(data))),
                      var.genes = row.names(subset(fData(monocle_cds), use_for_ordering == TRUE)), 
                      is.expr = monocle_cds@lowerDetectionLimit,
                      meta.data = pData(monocle_cds),
                      project.name = 'exportCds', 
                      misc = mist_list
                      )
    
  } else if (export_to == 'Scater') {
    data <- log2(exprs(monocle_cds) + 1)
    pd <- new("AnnotatedDataFrame", data = pData(monocle_cds))
    fd <- new("AnnotatedDataFrame", data = fData(monocle_cds))
    experimentData = monocle_cds@experimentData
    
    export_cds <- scater::newSCESet(exprsData = data, countData = NULL, tpmData = NULL,
              fpkmData = NULL, cpmData = NULL, phenoData = pd, featureData = fd,
              experimentData = experimentData, is_exprsData = NULL,
              cellPairwiseDistances = dist(vector()),
              featurePairwiseDistances = dist(vector()), 
              lowerDetectionLimit = monocle_cds@lowerDetectionLimit,
              logExprsOffset = 1, logged = TRUE, useForExprs = "exprs")    
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
#' 
importCDS <- function(otherCDS, import_all = FALSE) {
  if(class(otherCDS)[1] == 'seurat') {
    data <- otherCDS@raw.data

    pd <- new("AnnotatedDataFrame", data = otherCDS@meta.data)
    lowerDetectionLimit <- otherCDS@is.expr
    
    if(is.integer(data[1, 1])) {
      expressionFamily <- negbinomial.size()
    } else if(any(data < 0)){
      expressionFamily <- gaussianff()
    } else {
      expressionFamily <- tobit()
    }
    
    monocle_cds <- new("CellDataSet",
                       assayData = assayDataNew( "environment", exprs=data ),
                       phenoData=pd, 
                       lowerDetectionLimit=lowerDetectionLimit,
                       expressionFamily=expressionFamily,
                       dispFitInfo = new.env( hash=TRUE ))
    
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
    
    monocle_cds <- setOrderingFilter(monocle_cds, var.genes)
    monocle_cds@auxClusteringData$seurat <- mist_list
    
  } else if (class(otherCDS)[1] == 'scater') {
    if(otherCDS@logged) {
      data <- 2^otherCDS@assayData$exprs - otherCDS@logExprsOffset
    }
    else { 
      data <- otherCDS@assayData$exprs
    }
    
    pd <- otherCDS@featureData
    fd <- otherCDS@phenoData
    experimentData = monocle_cds@experimentData
    
    if(is.integer(data[1, 1])) {
      expressionFamily <- negbinomial.size()
    } else if(any(data < 0)){
      expressionFamily <- gaussianff()
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
    monocle_cds <- new("CellDataSet",
                        assayData = assayDataNew( "environment", exprs=data ),
                        phenoData=pd, 
                        featureData=fd, 
                        lowerDetectionLimit=lowerDetectionLimit,
                        expressionFamily=expressionFamily,
                        dispFitInfo = new.env( hash=TRUE ))
    
    # monocle_cds@auxClusteringData$sc3 <- otherCDS@sc3
    # monocle_cds@auxOrderingData$scran <- mist_list
    
    monocle_cds@auxOrderingData$scran <- mist_list
    
  } else {
    stop('the object type you want to export to is not supported yet')
  }
  
  return(monocle_cds)
}




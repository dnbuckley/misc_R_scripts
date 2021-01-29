library(QDNAseq)

# standard QDNAseq pipeline
callCopyNumbers <- function(readCounts, normRef = NULL, noisePlot = T, returnCGH = T, 
                            future.globals.maxSize.GB = NULL){
  if (!is.null(future.globals.maxSize.GB)){
    options(future.globals.maxSize = future.globals.maxSize.GB*1000*1024^2)
  }
  readCountsFiltered <- applyFilters(readCounts, residual = T, blacklist = T)
  readCountsFiltered <- estimateCorrection(readCountsFiltered)
  if (noisePlot) noisePlot(readCountsFiltered)
  copyNumbers <- correctBins(readCountsFiltered)
  copyNumbersNormalized <- normalizeBins(copyNumbers)
  copyNumbersSmooth <- smoothOutlierBins(copyNumbersNormalized)
  if (!is.null(normRef)){
    tumorVsNormal <- compareToReference(copyNumbersSmooth, normRef)
    copyNumbersSegmented <- segmentBins(tumorVsNormal, transformFun="sqrt")
  } else{
    copyNumbersSegmented <- segmentBins(copyNumbersSmooth, transformFun="sqrt")
  }
  copyNumbersSegmented <- normalizeSegmentedBins(copyNumbersSegmented)
  copyNumbersCalled <- callBins(copyNumbersSegmented)
  if (returnCGH){
    copyNumbersCalled <- makeCgh(copyNumbersCalled, filter = F)
  }
  return(copyNumbersCalled)
}

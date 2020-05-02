library(QDNAseq)

# standard QDNAseq pipeline
callCopyNumbers <- function(readCounts, normRef = NULL, noisePlot = T){
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
  copyNumbersCalled <- makeCgh(copyNumbersCalled, filter = F)
  return(copyNumbersCalled)
}
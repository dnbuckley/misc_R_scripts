library(QDNAseq)
library(GenomicRanges)
library(Biobase)
source("~/misc_R_scripts/combineInRange.R")

# attempt to find focal CNA from QDNAseq data
# CN cutoff: log2 value for cutoff
# columns: which columns of CNC to use
# aproxGR: GRanges of aprox location
# windowSize: how far apart can the points be and still be considered in the same CNA
# tweak all these params until you get the region you want
getFocalCNA <- function(copyNumbersCalled, CNcutoff, chr, columns, aproxGR = NULL, windowSize = 5){
  df <- data.frame(copynum = assayDataElement(copyNumbersCalled[, columns], "copynumber"),
                   segs = assayDataElement(copyNumbersCalled[, columns], "segmented"),
                   call = assayDataElement(copyNumbersCalled[, columns], "calls"))
  rownames(df) <- paste0("chr", rownames(df))
  df <- df[complete.cases(df), ]
  df <- df[grepl(paste0(chr, ":"), rownames(df)), ]
  if (CNcutoff > 1){
    ROI <- df[rowMeans(df[, grepl("copynum\\.", names(df))]) > CNcutoff, ]
  } else {
    ROI <- df[rowMeans(df[, grepl("copynum\\.", names(df))]) < CNcutoff, ]
  }
  ROIgr <- GRanges(rownames(ROI))
  ROIgr <- combineInRange(ROIgr, range = windowSize)
  if (!is.null(aproxGR)){
    ovlp <- findOverlaps(aproxGR, ROIgr)
    ROIgr <- ROIgr[to(ovlp)]
  }
  # nBins: number of bins in CNA
  # couldnt get this to work well
  # if (!is.null(nBins)){
  #   binSize <- width(GRanges(rownames(df)[1]))
  #   widthCutoff <- binSize * nBins
  #   # 1 bin forgivness window
  #   ROIgr <- ROIgr[width(ROIgr) > widthCutoff - binSize]
  #   ROIgr <- ROIgr[width(ROIgr) < widthCutoff + binSize]
  # }
  return(ROIgr)
}

library(rtracklayer)
library(GenomicRanges)
library(parallel)

# memory efficient version, reads in only relevant segments
# return mat will return a granges list or matrix of beta values 
# with granges as row names
summarizeBigWigFilesOver <- function(files, gr, cores = 4, returnMat = F){
  bws <- mclapply(files, import, which = gr, mc.cores = cores)
  bws <- mclapply(bws, .summarizeScoreOver, segs = gr, mc.cores = cores)
  if (returnMat) {
    source("~/misc_R_scripts/getGrangesScoreMatrix.R")
    bws <- getGrangesScoreMatrix(bws)
    rownames(bws) <- paste0(granges(gr))
  }
  return(bws)
}

#' helper function: gets median score across regions  
#'
#' @param   gr        a GRanges with scores to summarize
#' @param   segs      the regions to summarize over 
#' @return  matrix    a matrix of regional "beta values"
#'
#' @import  GenomicRanges
#' 
#' @export
.summarizeScoreOver <- function(gr, segs) { 
  gr <- subsetByOverlaps(gr, segs)
  segs <- granges(segs) # to score
  names(segs) <- as.character(segs)
  ol <- findOverlaps(gr, segs)
  segs$score <- rep(NA_real_, length(segs))
  bySeg <- sapply(split(score(gr)[queryHits(ol)], names(segs)[subjectHits(ol)]), mean, na.rm = T)
  segs[names(bySeg)]$score <- bySeg
  return(segs)  
}

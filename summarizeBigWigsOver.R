library(GenomicRanges)
library(rtracklayer)

#' like it says on the tin...
#'
#' @param   grl       a GRangesList of bigWig scores
#' @param   segs      the regions to summarize over 
#'
#' @return  matrix    a matrix of regional "beta values"
#'
#' @import  GenomicRanges
#' 
#' @export
summarizeBigWigsOver <- function(grl, segs) { 
  stopifnot(is(grl, "GRangesList"))
  mat <- matrix(NA_real_, nrow=length(segs), ncol=length(grl))
  rownames(mat) <- as.character(segs)
  colnames(mat) <- names(grl) 
  toSummarize <- GRangesList(lapply(grl, subsetByOverlaps, segs))
  for (i in colnames(mat)) {
    mat[, i] <- summarizeScoreOver(grl[[i]], segs)$score
  }
  return(mat)
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
summarizeScoreOver <- function(gr, segs) { 
  gr <- subsetByOverlaps(gr, segs)
  segs <- granges(segs) # to score
  names(segs) <- as.character(segs)
  ol <- findOverlaps(gr, segs)
  segs$score <- rep(NA_real_, length(segs))
  bySeg <- sapply(split(score(gr)[queryHits(ol)], names(segs)[subjectHits(ol)]), median)
  segs[names(bySeg)]$score <- bySeg
  return(segs)  
}

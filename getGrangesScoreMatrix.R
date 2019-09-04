library(GenomicRanges)

# converts a list of GRanges w/ score to score matrix
getGrangesScoreMatrix <- function(grs){
  # make  sure all granges are identical
  for (i in 1:length(grs)){
    stopifnot(identical(granges(grs[[1]]), granges(grs[[i]])))
  }
  scores <- lapply(grs, score)
  scores <- do.call(cbind, scores)
  rownames(scores) <- paste0(granges(grs[[1]]))
  colnames(scores) <- names(grs)
  return(scores)
}
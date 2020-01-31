library(GenomicRanges)

# the trick is to set the names of the
# granges list to NULL, stupid bug
combineListOfGR <- function(grl){
  names(grl) <- NULL
  gr <- do.call(c, grl)
  return(gr)
}

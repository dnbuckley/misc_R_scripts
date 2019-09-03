library(GenomicRanges)
library(bsseq)
# convert a list of DMRs to single granges obj
list2gr <- function(grl){
  gr <- GRanges()
  for (grx in grl) gr <- c(gr, granges(grx))
  return(gr)
}

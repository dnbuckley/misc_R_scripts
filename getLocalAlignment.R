library(GenomicRanges)
library(Biostrings)
source("~/misc_R_scripts/getBSgenomeObj.R")

getLocalAlignment <- function(gr, query.seq, genome, buffer = NULL,
                              max.mismatch = 0) {
  bsg <- getBSgenomeObj(genome = genome)
  if (!is.null(buffer)) {
    start(gr) <- start(gr)-buffer
    end(gr) <- end(gr)+buffer
  }
  ref.seq <- getSeq(bsg, gr)
  match.loc  <- unlist(vmatchPattern(pattern = query.seq, subject = ref.seq, max.mismatch = max.mismatch))
  if (length(match.loc) != 1) {
    return(length(match.loc))
  }
  seq.loc <- GRanges(seqnames = seqnames(gr), 
                     ranges = IRanges(start = start(gr)+start(match.loc)-1,
                                      end = start(gr)+end(match.loc)-1),
                     strand = strand(gr))
  return(seq.loc)
}

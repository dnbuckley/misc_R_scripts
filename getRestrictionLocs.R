library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Hsapiens.UCSC.hg19)
# get restriction locations for specified rescriction sequence
# default MSP1

.combineInRange <- function(gr, range){
  end(gr) <- end(gr) + range
  gr <- reduce(gr)
  end(gr) <- end(gr) - range
  return(gr)
}

getRestrictionLocs <- function(sequence = "CCGG", genome = "hg19", cores = 4, proximity = 1000){
  if (genome == "hg19") {
    Hsapiens <- BSgenome.Hsapiens.UCSC.hg19
  } else if (genome == "hg38") {
    Hsapiens <- BSgenome.Hsapiens.UCSC.hg38
  } else message("invalid genome selection")
  chrs <- names(Hsapiens)[1:24]
  locs <- mclapply(chrs, function(x) start(matchPattern("CCGG", 
                                              Hsapiens[[x]])), mc.cores = cores)
  gr <- do.call(c, lapply(1:24, function(x) GRanges(names(Hsapiens)[x], 
                                            IRanges(locs[[x]], width = length(sequence)))))
  seqinfo(gr) <- seqinfo(Hsapiens)[seqlevels(gr)]
  gr <- .combineInRange(gr, range = proximity)
  gr <- gr[width(gr) > length(sequence)]
  return(gr)
}


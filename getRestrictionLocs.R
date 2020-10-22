library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Hsapiens.UCSC.hg19)
library(BSgenome.Hsapiens.NCBI.GRCh38)
# get restriction locations for specified rescriction sequence
# default MSP1

.segmentCHR <- function(chr){
  starts <- start(chr)
  ends <- end(chr)
  x <- GRanges(cbind.data.frame(chr = unique(seqnames(chr)),
                                 start = ends[1:(length(ends) - 1)],
                                 end = starts[2:length(starts)]))
  return(x)
}
# default for RRBS regions determined from OVCA RRBS 9-19
# reduce option combines adjacent regions
getRestrictionLocs <- function(sequence = "CCGG", genome = "hg38", cores = 4,
                               minWidth = 30, maxWidth = 230, reduce = T){
  if (genome == "hg19") {
    Hsapiens <- BSgenome.Hsapiens.UCSC.hg19
  } else if (genome == "hg38") {
    Hsapiens <- BSgenome.Hsapiens.UCSC.hg38
  } else if (genome == "GRCh38") {
    Hsapiens <- BSgenome.Hsapiens.NCBI.GRCh38
  } else message("invalid genome selection")
  chrs <- names(Hsapiens)[1:24]
  locs <- mclapply(chrs, function(x) start(matchPattern(sequence, 
                                              Hsapiens[[x]])), mc.cores = cores)
  gr <- mclapply(1:24, function(x) GRanges(names(Hsapiens)[x], 
                                IRanges(locs[[x]], width = nchar(sequence))), mc.cores = cores)
  gr <- do.call("c", lapply(gr, .segmentCHR))
  seqinfo(gr) <- seqinfo(Hsapiens)[seqlevels(gr)]
  start(gr) <- start(gr) - 2
  end(gr) <- end(gr) + 2
  gr <- gr[width(gr) >= minWidth]
  gr <- gr[width(gr) <= maxWidth]
  if (reduce) gr <- reduce(gr)
  return(gr)
}


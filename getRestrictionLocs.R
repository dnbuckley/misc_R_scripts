library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Hsapiens.UCSC.hg19)
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

getRestrictionLocs <- function(sequence = "CCGG", genome = "hg19", cores = 4, minWidth = 20, maxWidth = 200){
  if (genome == "hg19") {
    Hsapiens <- BSgenome.Hsapiens.UCSC.hg19
  } else if (genome == "hg38") {
    Hsapiens <- BSgenome.Hsapiens.UCSC.hg38
  } else message("invalid genome selection")
  chrs <- names(Hsapiens)[1:24]
  locs <- mclapply(chrs, function(x) start(matchPattern(sequence, 
                                              Hsapiens[[x]])), mc.cores = cores)
  gr <- lapply(1:24, function(x) GRanges(names(Hsapiens)[x], 
                                IRanges(locs[[x]], width = nchar(sequence))))
  gr <- do.call(c, lapply(gr, .segmentCHR))
  seqinfo(gr) <- seqinfo(Hsapiens)[seqlevels(gr)]
  gr <- gr[width(gr) >= minWidth]
  gr <- gr[width(gr) <= maxWidth]
  export.bed(gr, "test.bed")
  return(gr)
}


library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Hsapiens.UCSC.hg19)
library(BSgenome.Hsapiens.NCBI.GRCh38)

findSequence <- function(sequence, genome = "hg38", cores = 4){
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
  gr <- do.call(c, gr)
  seqinfo(gr) <- seqinfo(Hsapiens)[seqlevels(gr)]
  return(gr)

}


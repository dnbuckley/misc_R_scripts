# get number of cpgs in a granges obj
library(GenomicRanges)

.getCpGs <- function (genome = "hg19") {
  require(GenomicRanges)
  require(paste0("BSgenome.Hsapiens.UCSC.", genome), character.only = TRUE)
  message("Building CpG positions for ", genome)
  chrs <- names(Hsapiens)[1:24]
  cgs <- lapply(chrs, function(x) start(matchPattern("CG", 
                                                     Hsapiens[[x]])))
  cpgr <- do.call(c, lapply(1:24, function(x) GRanges(names(Hsapiens)[x], 
                                                      IRanges(cgs[[x]], width = 2))))
  seqinfo(cpgr) <- seqinfo(Hsapiens)[seqlevels(cpgr)]
  return(cpgr)
}

getNumCpGs <- function(gr, cpGR = NULL, genome = "hg19"){
  if (is.null(cpGR)){
    cpGR <- .getCpGs(genome)
  }
  ovlp <- findOverlaps(cpGR, gr)
  return(as.vector(table(to(ovlp))))
}

library(GenomicRanges)

.calcLocs <- function(genome = "hg38", cores = 1) {
  if (genome == "hg19") {
    Hsapiens <- BSgenome.Hsapiens.UCSC.hg19
  } else if (genome == "hg38") {
    Hsapiens <- BSgenome.Hsapiens.UCSC.hg38
  } else if (genome == "GRCh38") {
    Hsapiens <- BSgenome.Hsapiens.NCBI.GRCh38
  } else message("invalid genome selection")
  chrs <- names(Hsapiens)[1:24]
  locs <- mclapply(chrs, function(x) start(matchPattern("CG", 
                                                        Hsapiens[[x]])), mc.cores = cores)
  gr <- mclapply(1:24, function(x) { GRanges(names(Hsapiens)[x], 
                                           IRanges(locs[[x]], width = 2))}, mc.cores = cores)
  gr <- do.call("c", gr)
  seqinfo(gr) <- seqinfo(Hsapiens)[seqlevels(gr)]
  seqlengths(gr) <- seqlengths(Hsapiens)[seqlevels(gr)]
  return(gr)
}

# does what it says on the label
getCpGLocs <- function(genome = "hg38", strand = NULL, cores = 1, recalc = F) {
  if (Sys.info()["user"] == "dbuckley" & !recalc){
    file <- paste0("~/Desktop/salhia_lab/useful_granges/cpg_locations_", genome, ".rds")
    cpgr <- readRDS(file)
  } else {
    source("~/misc_R_scripts/getRestrictionLocs.R")
    cpgr <- .calcLocs(cores = cores)
  }
  if (!is.null(strand)) {
    if (strand == "+") {
      end(cpgr) <- start(cpgr)
      strand(cpgr) <- "+"
    } else if (strand == "-") {
      start(cpgr) <- end(cpgr)
      strand(cpgr) <- "-"
    } else {
      message("Invalid strand, please select '+' || '-'")
    }
  }
  return(cpgr)
}

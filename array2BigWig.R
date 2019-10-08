library(GenomicRanges)
library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg19)
library(BSgenome.Hsapiens.UCSC.hg38)
source("~/misc_R_scripts/combineInRange.R")

array2BigWig <- function(input, outfile, outDir = "", genome = "hg19", rmDupCpG = F){
  if(class(input) == "character"){
    input <- readRDS(input)
  } 
  if (class(input) != "GRanges"){
    message("invalid input")
    return(NULL)
  }
  gr <- granges(input)
  score(gr) <- input$Beta_value
  gr <- gr[!is.na(score(gr))]
  if (genome == "hg19") {
    seqlevels(gr) <- seqlevels(BSgenome.Hsapiens.UCSC.hg19)
    seqlengths(gr) <- seqlengths(BSgenome.Hsapiens.UCSC.hg19)
  } else if (genome == "hg38") {
    seqlevels(gr) <- seqlevels(BSgenome.Hsapiens.UCSC.hg38)
    seqlengths(gr) <- seqlengths(BSgenome.Hsapiens.UCSC.hg38)
  } else {
    message("invalid genome")
    return(NULL)
  }
  if(rmDupCpG){
    # remove duplicates
    allCpGs <- table(paste0(seqnames(gr), start(gr), end(gr)))
    gr <- gr[!paste0(seqnames(gr), start(gr), end(gr)) %in% names(allCpGs[allCpGs > 1])]
    # remove adjacent
    badCG <- combineInRange(gr, range = 1)
    badCG <- badCG[width(badCG) > 2]
    ovlp <- findOverlaps(gr, badCG)
    gr <- gr[-from(ovlp)]
  }
  export.bw(gr, paste0(outDir, outfile))
}

library(GenomicRanges)
library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg19)
library(BSgenome.Hsapiens.UCSC.hg38)

array2BigWig <- function(input, outfile, outDir = "", genome = "hg19"){
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
  }
  else if (genome == "hg38") {
    seqlevels(gr) <- seqlevels(BSgenome.Hsapiens.UCSC.hg38)
    seqlengths(gr) <- seqlengths(BSgenome.Hsapiens.UCSC.hg38)
  }
  else {
    message("invalid genome")
    return(NULL)
  }
  export.bw(gr, paste0(outDir, outfile))
}

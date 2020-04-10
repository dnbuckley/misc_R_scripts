# returns bsgenome object for genome
getBSgenomeObj <- function(genome) {
  if (genome == "hg19") {
    library(BSgenome.Hsapiens.UCSC.hg19)
    return(BSgenome.Hsapiens.UCSC.hg19)
  }
  if (genome == "hg38") {
    library(BSgenome.Hsapiens.UCSC.hg38)
    return(BSgenome.Hsapiens.UCSC.hg38)
  }
  if (genome == "mm10"){
    library(BSgenome.Mmusculus.UCSC.mm10)
    return(BSgenome.Mmusculus.UCSC.mm10)
  }
  message("Invalid Genome")
  return()
}
# get chromosome lengths
getChrSeqLengths <- function(genome = "hg19", sex = F){
  if (genome == "hg19"){
    library(BSgenome.Hsapiens.UCSC.hg19)
    sl <- seqlengths(BSgenome.Hsapiens.UCSC.hg19)
  }
  if (genome == "hg38"){
    library(BSgenome.Hsapiens.UCSC.hg38)
    sl <- seqlengths(BSgenome.Hsapiens.UCSC.hg38)
  }
  if (sex) chrs <- paste0("chr", c(1:22, "X", "Y"))
  else chrs <- paste0("chr", 1:22)
  sl <- sl[names(sl) %in% chrs]
  return(sl)
}

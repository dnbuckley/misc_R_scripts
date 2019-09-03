# get txdb obj without having to download from UCSC
getTxDBOffline <- function(genome = "hg19"){
  if (genome == "hg19"){
    library(TxDb.Hsapiens.UCSC.hg19.knownGene)
    return(TxDb.Hsapiens.UCSC.hg19.knownGene)
  }
  else if (genome == "hg38"){
    library(TxDb.Hsapiens.UCSC.hg38.knownGene)
    return(TxDb.Hsapiens.UCSC.hg38.knownGene)
  }
  else if (genome == "mm10"){
    library(TxDb.Mmusculus.UCSC.mm10.knownGene)
    return(TxDb.Mmusculus.UCSC.mm10.knownGene)
  }
  else{
    message("Invalid genome option")
    return(NULL)
  }
}
library(GenomicFeatures)
library(org.Hs.eg.db)
library(org.Mm.eg.db)

# get genes in a certain range, use refGRObj for genomes not supported
getGenesInRange <- function(gr, refGRObj = NULL, genome = NULL,
                            withinOnly = F){
  if (is.null(refGRObj) & is.null(genome)){
    message("Provide either a reference GRanges ojbect with a \\'SYMBOL\\' mcol
             or specify a genome")
    return(NULL)
  }
  if (!is.null(refGRObj)){
    if(!any(grepl("^SYMBOL$", names(mcols(refGRObj))))){
      message("refGRObj must have \\'SYMBOL\\' mcol")
      return(NULL)
    }
  }
  if (!is.null(genome)){
    if (genome == "hg19"){
      library(TxDb.Hsapiens.UCSC.hg19.knownGene)
      genes <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
      geneList <- select(org.Hs.eg.db, 
                         keys = genes$gene_id,
                         columns = "SYMBOL")
    }
    else if (genome == "hg38"){
      library(TxDb.Hsapiens.UCSC.hg38.knownGene)
      genes <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
      geneList <- select(org.Hs.eg.db, 
                         keys = genes$gene_id,
                         columns = "SYMBOL")
    }
    else if (genome == "mm10"){
      library(TxDb.Mmusculus.UCSC.mm10.knownGene)
      genes <- genes(TxDb.Mmusculus.UCSC.mm10.knownGene)
      geneList <- select(org.Mm.eg.db, 
                         keys = genes$gene_id,
                         columns = "SYMBOL")
    }
    else{
      message("Invalid genome option")
      return(NULL)
    }
    refGRObj <- GRanges(granges(genes), mcols = geneList)
  }
  if (withinOnly) ovlp <- findOverlaps(refGRObj, gr, type = "within")
  else ovlp <- findOverlaps(refGRObj, gr)
  return(split(refGRObj$mcols.SYMBOL[from(ovlp)], paste0(granges(gr))[to(ovlp)]))
}










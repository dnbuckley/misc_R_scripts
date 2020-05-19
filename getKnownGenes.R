require(GenomicFeatures)
require(org.Hs.eg.db)
source("~/misc_R_scripts/getTxDBOffline.R")

getKnownGenes <- function(genome = "hg19") {
  refseq.db <- getTxDBOffline(genome)
  refseq.genes <- genes(refseq.db)
  df <- select(org.Hs.eg.db, 
               keys = refseq.genes$gene_id, 
               columns = c("SYMBOL"))
  mcols(refseq.genes) <- df
  return(refseq.genes)
}
library(GenomicFeatures)
library(org.Hs.eg.db)

### GetNearestGeneSymbolHg19
# Find nearest gene
# WARNING: ChrM does not map properly - please filter them out first!
#' @param granges       granges object containing regions of interest
#' @param refseq.genes  TxDb object generated from AnnotationDbi
#'                      This is NOT required but will speed up process
#' @return dataframe with gene_symbol and distance
# 
GetNearestGeneSymbol <- function (granges, refseq.genes=NULL, genome = "hg19") {
  if (is.null(refseq.genes)) {
    source("~/misc_R_scripts/getTxDBOffline.R")
    refseq.genes <- genes(getTxDBOffline(genome))
  }
  distance <- distanceToNearest(granges, refseq.genes)
  genelist <- select(org.Hs.eg.db, 
                     keys = values(refseq.genes)[to(distance),], 
                     columns = c("SYMBOL"))
  output <- data.frame(nearest_gene=genelist$SYMBOL, 
                       distance=values(distance))
  return(output)
}
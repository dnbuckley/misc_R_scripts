# converts ensembl ID to HGNC gene symbol
library(biomaRt)
ens2hgnc <- function(ensemblIDs) {
  library('biomaRt')
  mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  genes <- gsub("_.*", "", rownames(fpkm))
  G_list <- getBM(filters = "ensembl_gene_id", attributes = c("ensembl_gene_id","hgnc_symbol"),
                  values = ensemblIDs, mart = mart, uniqueRows = T)
  return(G_list)
}

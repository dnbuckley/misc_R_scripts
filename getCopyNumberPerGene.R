library(QDNAseq)

# copyNumbersCalled object from QDNA seq
# reference set of genes || a genome
getCopyNumberPerGene <- function(copyNumbersCalled, 
                                 genes = NULL, 
                                 genome = NULL,
                                 returnDF = F) {
  if (is.null(genes) & is.null(genome)) {
    stop("Must provide a gene granges object or a genome")
  } 
  if (is.null(genes)) {
    if (Sys.info()["user"] == "dbuckley") {
      file <- paste0("~/Desktop/salhia_lab/useful_granges/", genome, "_HUGO-genes_granges.rds")
      message("Loading ", file, "...")
      genes <- readRDS(file)
    } else {
      source("~/misc_R_scripts/getKnownGenes.R")
      genes <- getKnownGenes(genome)
    }
  }
  if (!any(grepl("SYMBOL", names(mcols(genes))))) {
    stop("Gene granges must contain SYMBOL column")
  }
  genes <- genes[order(genes)]
  genes <- genes[!is.na(genes$SYMBOL)]
  cn <- as.data.frame(assayDataElement(copyNumbersCalled, "copynumber"))
  if (width(GRanges(rownames(cn)[1])) > 5000) {
    warning("QDNAseq object has a binsize of > 5000, for gene level
copy number you may want to use a smaller binsize.")
  }
  rownames(cn) <- paste0("chr", rownames(cn))
  message("Finding overlaps between QDNAseq object and reference genes...")
  ovlp <- findOverlaps(GRanges(rownames(cn)), genes)
  cn <- split(cn[from(ovlp), ], genes$SYMBOL[to(ovlp)])
  cn <- mclapply(cn, colMeans, na.rm = T, mc.cores = 6)
  genes <- genes[genes$SYMBOL %in% names(cn)]
  genes <- genes[order(genes$SYMBOL)]
  df <- as.data.frame(do.call(rbind, cn))
  stopifnot(identical(rownames(df), genes$SYMBOL))
  if (returnDF) {
    return(df)
  }
  mcols(genes) <- cbind.data.frame(mcols(genes), df)
  genes <- genes[order(genes)]
  return(genes)
}



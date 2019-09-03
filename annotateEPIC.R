setwd("~/Desktop/salhia_lab/lenz_EPIC/")
library(minfi)
library(DMRcate)


# annotates CpG info to a gset object, groups are column groups in the gset object
annotateEPIC <- function(gset, groups, refAnnotations = NULL, 
                         dmpType = "categorical", reAnnotateGenes = F){
  if (is.null(refAnnotations)){
    refAnnotations <- readRDS("~/Desktop/salhia_lab/lenz_EPIC/EPIC.hg19.manifest.rds")
  }
  # in order to annotate custom
  refAnnotations <- refAnnotations[seqnames(refAnnotations) != "chrM"]
  refAnnotations <- refAnnotations[order(names(refAnnotations))]
  beta <- getBeta(gset)
  beta <- beta[order(rownames(beta)), ]
  dmp <- dmpFinder(beta,
                   pheno = groups,
                   type = dmpType,
                   shrinkVar = T)
  dmp <- dmp[order(rownames(dmp)), ]
  names(dmp)[1] <- "dmValue"
  CpGInt <- intersect(rownames(beta), names(refAnnotations))
  refAnnotations <- refAnnotations[names(refAnnotations) %in% CpGInt]
  beta <- beta[rownames(beta) %in% CpGInt, ]
  dmp <- dmp[rownames(dmp) %in% CpGInt, ]
  stopifnot(identical(rownames(beta), names(refAnnotations)))
  stopifnot(identical(rownames(beta), rownames(dmp)))
  colnames(beta) <- paste0(colnames(beta), ".meth")
  if (reAnnotateGenes){
    source("~/misc_R_scripts/GetNearestGeneSymbol.R")
    genes <- GetNearestGeneSymbol(granges(refAnnotations))
    df <- cbind(CpG_ID = rownames(beta),
                dmp,
                genes,
                beta,
                mcols(refAnnotations))
  } else {
    df <- cbind(CpG_ID = rownames(beta),
                dmp,
                beta,
                mcols(refAnnotations))
  }
  gr <- granges(refAnnotations)
  stopifnot(identical(df$CpG_ID, names(gr)))
  mcols(gr) <- df
  gr <- gr[order(gr)]
  return(gr)
}


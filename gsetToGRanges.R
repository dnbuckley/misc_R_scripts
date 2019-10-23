library(GenomicRanges)

gsetToGRanges <- function(gset, refGR = NULL){
  if (is.null(refGR)) {
    warning("no refGR object provided, using EPIC hg19")
    refGR <- readRDS("~/misc_R_scripts/script_files/EPIC.hg19.manifest.rds")
  }
  refGR <- granges(refGR)
  meth <- getBeta(gset)
  common <- intersect(names(refGR), rownames(meth))
  refGR <- refGR[names(refGR) %in% common]
  meth <- meth[rownames(meth) %in% common, ]
  refGR <- refGR[order(names(refGR))]
  meth  <- meth[order(rownames(meth)), ]
  if (any(names(pData(gset)) %in% "Sample.ID")) {
    colnames(meth) <- gset$Sample.ID
  }
  stopifnot(identical(names(refGR), rownames(meth)))
  gr <- GRanges(refGR, mcols = meth)
  return(gr[order(gr)])
}

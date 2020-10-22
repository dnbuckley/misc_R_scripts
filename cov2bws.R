library(rtracklayer)
library(bsseq)
source("~/misc_R_scripts/getBSgenomeObj.R")

cov2bws <- function(covFile, genome = "hg38", outFolder = "") {
  message(covFile)
  b <- read.bismark(covFile, verbose = T)
  bsg <- getBSgenomeObj(genome)
  g <- granges(b)
  seqlevels(g) <- seqlevels(bsg)
  seqinfo(g) <- seqinfo(bsg)
  m <- getMeth(b, type = "raw")
  c <- getCoverage(b)
  gm <- g
  gc <- g
  score(gm) <- m
  score(gc) <- c
  outFile <- gsub("\\.gz$", "", covFile)
  outFile <- gsub("\\.cov$", "", outFile)
  outFile <- gsub(".*\\/", "", outFile)
  export.bw(gm, paste0(outFolder, "/", outFile, ".meth.bw"))
  export.bw(gc, paste0(outFolder, "/", outFile, ".covg.bw"))
}

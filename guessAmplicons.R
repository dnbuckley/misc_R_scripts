# Guess the amplicon
source("~/misc_R_scripts/AnnotateNearestGene.R")

#' Guess amplicons from a bsseq object
#' Amplicon is defined as a cluster of 
#' @param bs         bsseq object 
#' @param minCov     minimum coverage across region (average across all CpGs)
#' @param minDist    minimum distance between two amplicons
#'                   also the maximum distance between CpGs within an amplicon
#' @param genome     genome build for finding nearest gene
#' @param manFormat  (TRUE - default) outputs in manifest format
#'                   (FALSE) outputs grange object
guessAmplicons <- function(bs, minCov=10, minDist=150, genome="hg19", manFormat=FALSE) {
  # Get all CpGs with coverage > minCov
  rowavg <- getCoverage(bs, type="Cov")
  rowavg[rowavg == 0] <- NA
  rowavg <- rowMeans(rowavg, na.rm = TRUE)
  bs2 <- bs[rowavg >= minCov]
  # Calculate distance between CpGs
  gr <- granges(bs2)
  dist <- end(gr)[-1] - start(gr)[-length(gr)]
  largedist <- which(abs(dist) > abs(minDist))
  starts <- c(1, (largedist + 1))
  ends <- c(largedist, length(gr))
  gr2 <- gr[starts]
  end(gr2) <- end(gr[ends])
  # Filter out width too short
  if (any(width(gr2) <= 1)) message("Warning: filtered out ", sum(width(gr2) <= 1), " amplicons with only 1 CpG!")
  gr2 <- gr2[width(gr2) > 1]
  # Width too long warning message
  if (any(width(gr2) > 1000)) message("Warning: found ", sum(width(gr2) > 1000), " amplicons over 1000 bp!")
  # Find gene of interest
  gr3 <- gr2; strand(gr3) <- "*" # unstrand amplicons
  gr2$genes <- GetNearestGeneSymbol(gr3, genome = genome)$genes
  # Find average amplicon coverage
  cov <- getCoverage(bs, regions = gr2, type = "Cov", what = "perRegionAverage")
  cov[cov == 0] <- NA
  cov <- rowMeans(cov, na.rm = TRUE)
  gr2$covg <- cov
  if (manFormat == TRUE) {
    res <- data.frame(Name = gr2$genes,
                      Strand = strand(gr2),
                      Chromosome = seqnames(gr2),
                      "Amplicon Start" = start(gr2),
                      "Amplicon End" = end(gr2))
    return(res)
  } else {
    return(gr2)
  }
}


setwd("~/Desktop/salhia_lab/komen_probes/")
library(openxlsx)
library(GenomicRanges)
library(rtracklayer)

.summarizeScoreOver <- function(gr, segs) { 
  gr <- subsetByOverlaps(gr, segs)
  segs <- granges(segs) # to score
  names(segs) <- as.character(segs)
  ol <- findOverlaps(gr, segs)
  segs$score <- rep(NA_real_, length(segs))
  bySeg <- sapply(split(score(gr)[queryHits(ol)], names(segs)[subjectHits(ol)]), mean, na.rm = T)
  segs[names(bySeg)]$score <- bySeg
  return(segs)  
}

beds <- list.files(".", ".bed$")
gc <- beds[grepl("gc", beds)]
ref <- beds[!grepl("gc", beds)]
names(gc) <- gsub("\\..*", "", gc)
names(ref) <- gsub("\\..*", "", ref)
gc <- lapply(gc, import.bed)
ref <- lapply(ref, import.bed)
stopifnot(identical(names(gc), names(ref)))

dfs <- list()
for (i in 1:length(ref)) {
  r <- ref[[i]]
  df <- data.frame(loc = paste0(granges(r)), 
                   ref = names(ref)[i])
  for (j in 1:length(gc)) {
    g <- gc[[j]]
    score(g) <- as.numeric(g$name)
    avgs <- score(.summarizeScoreOver(g, r))
    df[, j+2] <- avgs
  }
  dfs[[i]] <- df
  names(dfs)[i] <- names(ref)[i]
}







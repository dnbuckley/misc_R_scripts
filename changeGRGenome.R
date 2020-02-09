library(rtracklayer)

# convert genome positions from genome X to genome Y
# requires chain file in ~/misc_R_scripts/script_files/
# or provided file path
changeGRGenome <- function(gr, from = NULL, to = NULL, chainFile = NULL){
  if (all(is.null(c(from, to, chainFile)))) {
    message("ERROR: Please provide 'from' and 'to' genome names or a chainFile name")
    return(gr)
  }
  if (is.null(chainFile) & !is.null(from) & !is.null(to)){
    chainFile <- paste0("~/misc_R_scripts/script_files/", from, "To", to, ".over.chain")
  }
  chain <- import.chain(chainFile)
  x <- as.data.frame(liftOver(gr, chain))
  x <- x[, !names(x) %in% c("group", "group_name")]
  x <- GRanges(x)
  # remove identical ranges
  ovlp <- findOverlaps(x)
  ovlp <- ovlp[from(ovlp) != to(ovlp)]
  x <- x[-from(ovlp)]
  x <- x[-to(ovlp)]
  if(length(x) != length(gr)){
    message("\nWARNING: Conversion lost ", length(gr)-length(x), " ranges.\n")
  }
  return(x)
}


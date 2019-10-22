library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
library(BSgenome.Hsapiens.UCSC.hg38)

getRandomGranges <- function(n, genome = "hg19", meanW = 2e3, sdW = 750, minLen = 10){
  if (genome == "hg19"){
    lengths <- seqlengths(BSgenome.Hsapiens.UCSC.hg19)
    names(lengths) <- seqnames(BSgenome.Hsapiens.UCSC.hg19)
  } else if (genome == "hg38"){
    lengths <- seqlengths(BSgenome.Hsapiens.UCSC.hg38)
  } else{
    message("Invalid genome")
    return(NULL)
  }
  lengths <- lengths[names(lengths) %in% paste0("chr", c(1:22, "X", "Y"))]
  total <- sum(lengths)
  lenPct <- sapply(lengths, function(x, y){x/y}, total)
  n <- sapply(lenPct, function(x, y){floor(x*y)}, n)
  starts <- mapply(function(x, y){sample(1:y, x)}, n, lengths)
  widths <- lapply(n, function(x){floor(rnorm(x, mean = meanW, sd = sdW))})
  ends <- mapply(function(x, y){x+y}, starts, widths)
  grs <- list()
  for (chr in paste0("chr", c(1:22, "X", "Y"))){
    bad <- starts[[chr]] > ends[[chr]]
    starts[[chr]] <- starts[[chr]][!bad]
    ends[[chr]] <- ends[[chr]][!bad]
    gr <- GRanges(paste0(chr, ":", starts[[chr]], "-", end = ends[[chr]]))
    grs <- c(grs, gr)
  }
  gr <- Reduce(c, grs)
  gr <- reduce(gr)
  return(gr)
}

library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
library(BSgenome.Hsapiens.UCSC.hg38)

getRandomGranges <- function(n, genome = "hg19", meanW = 2e3,
                             sdW = 750, minLen = 10, chrs = paste0("chr", c(1:22, "X", "Y"))){
  if (genome == "hg19"){
    lengths <- seqlengths(BSgenome.Hsapiens.UCSC.hg19)
    names(lengths) <- seqnames(BSgenome.Hsapiens.UCSC.hg19)
  } else if (genome == "hg38"){
    lengths <- seqlengths(BSgenome.Hsapiens.UCSC.hg38)
  } else{
    message("Invalid genome")
    return(NULL)
  }
  lengths <- lengths[names(lengths) %in% chrs]
  total <- sum(lengths)
  lenPct <- sapply(lengths, function(x, y){x/y}, total)
  n <- sapply(lenPct, function(x, y){floor(x*y)}, n)
  starts <- mcmapply(function(x, y){sample(1:y, x)}, n, lengths, mc.cores = 10)
  widths <- mclapply(n, function(x){floor(rnorm(x, mean = meanW, sd = sdW))}, mc.cores = 10)
  ends <- mcmapply(function(x, y){x+y}, starts, widths, mc.cores = 10)
  grs <- list()
  for (chr in chrs){
    bad <- starts[[chr]] > ends[[chr]]
    starts[[chr]] <- starts[[chr]][!bad]
    ends[[chr]] <- ends[[chr]][!bad]
    gr <- GRanges(paste0(chr, ":", starts[[chr]], "-", end = ends[[chr]]))
    grs <- c(grs, gr)
  }
  gr <- Reduce(c, grs)
  gr <- reduce(gr)
  gr <- gr[width(gr) > minLen]
  return(gr)
}

library(bsseq)

# smooths in windows of size = window

# this is hopelessly inefficient but I'm too tired
# to figure out a better method at the moment
smoothGR <- function(gr, window = 10, cores = 10){
  gr <- gr[order(gr)]
  grs <- split(gr, as.character(seqnames(gr)))
  grs <- mclapply(grs, function(x, window){
    idx <- as.numeric(sapply(1:ceiling(length(x)/window), function(x){rep(x, window)}))
    idx <- idx[1:length(x)]
    x <- split(x, idx)
    return(x)
  }, window, mc.cores = cores)
  grs <- lapply(grs, function(x){
    x <- mclapply(x, function(y, seqname){
      z <- GRanges(data.frame(start = start(y)[1],
                              end = end(y)[length(y)],
                              seqnames = seqname,
                              strand = "*"))
      score(z) <- mean(score(y), na.rm = T)
      return(z)
      }, names(grs), mc.cores = cores)
  })
  grs <- lapply(grs, function(x){
    lapply(x, as.data.frame)
  })
  allRanges <- list()
  for (i in 1:length(grs)){
    allRanges <- c(allRanges, grs[[i]])
  }
  grs <- GRanges(do.call(rbind, allRanges))
  names(grs) <- paste0(granges(grs))
  return(grs)
}

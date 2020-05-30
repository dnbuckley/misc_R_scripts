source("~/misc_R_scripts/combineListOfGR.R")
source("~/misc_R_scripts/combineInRange.R")

# takes copy number data frame from QDNAseq and recalculates
# the bins, ie if binsize = 1kb and you want to increase
# the bin size to 5kb: reBinCopyNumbers(df, binN = 5)
reBinCopyNumbers <- function(df, binN, recalcRanges = T, cores = 10) {
  df <- as.data.frame(df)
  byChr <- split(df, gsub(":.*", "", rownames(df)))
  # split by Chr
  byChr <- lapply(byChr, function(x, binN){
    split(x, rep(1:ceiling(nrow(x)/binN), each = binN))
  }, binN)
  # get new ranges
  if (recalcRanges){
    chnks <- lapply(byChr, function(x){
      y <- mclapply(x, function(y){
        gr <- GRanges(rownames(y))
        gr <- combineInRange(gr, 1)
      }, mc.cores = cores)
      z <- combineListOfGR(y)
      return(z)
    })
    chnks <- combineListOfGR(chnks)
  }
  # average across
  byChr <- lapply(byChr, function(x){
    x <- mclapply(x, colMeans, mc.cores = cores, na.rm = T)
    x <- do.call(rbind, x)
    return(x)
  })
  df <- do.call(rbind, byChr)
  if (recalcRanges){
    rownames(df) <- paste0(chnks)
  }
  return(as.data.frame(df))
}
# helper function for each CHR
.chunkChr <- function(length, chr, width = 1000){
  starts <- seq(0, length, width)
  ends <- starts[2:length(starts)] - 1
  starts <- starts[1:length(starts) - 1]
  ends[length(ends)] <- length
  gr <- GRanges(cbind.data.frame(chr = chr,
                                 start = starts,
                                 end = ends))
  return(gr)
}

# returns whole genome non overlapping GRanges object of width n
getChrChunks <- function(chromLengths = NULL, genome = NULL, width = 1000, allosomes = F){
  if (is.null(chromLengths) && is.null(genome)){
    message("provide a seqLengths object or genome option... dumbass")
    return(NULL)
  }
  if (!is.null(genome)){
    if (genome == "hg38"){
      library(BSgenome.Hsapiens.UCSC.hg38)
      chromLengths <- seqlengths(BSgenome.Hsapiens.UCSC.hg38)
    }
    else if (genome == "hg19"){
      library(BSgenome.Hsapiens.UCSC.hg19)
      chromLengths <- seqlengths(BSgenome.Hsapiens.UCSC.hg19)
    }
    else {
      message("provide a valid genome: hg19 || hg38")
      return(NULL)
    }
  }
  if (allosomes){
    chromLengths <- chromLengths[names(chromLengths) %in% c(paste0("chr", 1:22), "chrX", "chrY")]
  } else chromLengths <- chromLengths[names(chromLengths) %in% paste0("chr", 1:22)]
  chrs <- names(chromLengths)
  chrs <- mapply(.chunkChr, chromLengths, chrs, width)
  gr <- GRanges()
  for (i in 1:length(chrs)) gr <- c(gr, chrs[[i]])
  return(gr)
}
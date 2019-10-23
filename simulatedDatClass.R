library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
library(BSgenome.Hsapiens.UCSC.hg38)

.getRandomGranges <- function(n, genome = "hg19", meanW = 2e3,
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
  starts <- mapply(function(x, y){sample(1:y, x)}, n, lengths)
  widths <- lapply(n, function(x){floor(rnorm(x, mean = meanW, sd = sdW))})
  ends <- mapply(function(x, y){x+y}, starts, widths)
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

.getNormalTracks <- function(covg, beta, cpgs){
  tracks <- lapply(1:covg, 
                   function(x, cpgs, beta){return(rbinom(length(cpgs), 1, beta))}, 
                   beta = beta, cpgs = cpgs)
  return(tracks)
}

.getDiseaseTracks <- function(num, betaRange = c(0.1, 0.5), cpgs){
  betas <- sample(seq(betaRange[1], betaRange[2], 0.05), num, replace = T)
  tracks <- lapply(betas, function(x, cpgs){rbinom(length(cpgs), 1, x)}, cpgs)
  names(tracks) <- betas
  return(tracks)
}

# def
setClass("fauxWGBS", representation(
  dTracks = "list",
  nTrakcs = "list",
  TPregions = "GRanges",
  randomRegions = "GRanges",
  loci = "GRanges"),
  contains = c("GRanges")
  )

covg <- 100
TPregions <- readRDS("hg19_HUGO-promoters_granges.rds")
cpgs <- cpgs[seqnames(cpgs) == "chr1"]
TPregions <- TPregions[seqnames(TPregions) == "chr1"]
# constructor
fauxWGBS <- function(covg, 
                     TPregions,
                     dBetaRange = c(0.1, 0.6), 
                     nBeta = 0.75, 
                     nRandomRegions = 1e3){
  nTracks <- .getNormalTracks(covg, beta = nBeta, cpgs = cpgs)
  dTracks <- .getDiseaseTracks(num = covg, betaRange = dBetaRange, cpgs = cpgs)
  randomRegions <- .getRandomGranges(n = 15e3)
  randomRegions <- randomRegions[seqnames(randomRegions) %in% unique(seqnames(cpgs))]
  randomRegions <- randomRegions[-from(findOverlaps(randomRegions, TPregions))]
  if (length(randomRegions) > nRandomRegions){
    randomRegions <- randomRegions[sample(1:length(randomRegions), nRandomRegions, replace = F)]
  }
  fWGBS <- new("fauxWGBS",
               dTracks = dTracks,
               nTrakcs = nTracks,
               TPregions = TPregions,
               randomRegions = randomRegions,
               loci = granges(cpgs))
}

.getNormal <- function(fWGBS, covg){
  tracks <- fWGBS@nTrakcs[sample(1:length(fWGBS@nTrakcs), covg, replace = F)]
  return(tracks)
}

.getDisease <- function(fWGBS, covg){
  tracks <- fWGBS@dTracks[sample(1:length(fWGBS@dTracks), covg, replace = F)]
  return(tracks)
}

getCOV <- function(fWGBS, covg, dFrac = NULL){
  gr <- fWGBS@loci
  if (is.null(dFrac)){
    n <- .getNormal(fWGBS, covg)
    n <- do.call(rbind, n)
    m <- colSums(n == 1)
    u <- colSums(n == 0)
    df <- data.frame("seqnames" = seqnames(gr),
                     "start" = start(gr),
                     "end" = end(gr),
                     "score" = m/(m+u),
                     "meth" = m,
                     "umeth" = u)
    return(df)
  } else{
    nD <- round(covg*dFrac)
    nN <- round(covg*(1-dFrac))
    d <- .getDisease(fWGBS, nD)
    n <- .getNormal(fWGBS, nN)
    nA <- .getNormal(fWGBS, covg)
    nA <- do.call(rbind, nA)
    c <- do.call(rbind, c(d, n))
    mc <- colSums(c == 1)
    uc <- colSums(c == 0)
    mn <- colSums(nA == 1)
    un <- colSums(nA == 0)
    grD <- GRanges(fWGBS@loci, mcols = cbind(mc, uc))
    grN <- GRanges(fWGBS@loci, mcols = cbind(mn, un))
    ovlp <- findOverlaps(grD, fWGBS@TPregions)
    grN[from(ovlp)] <- grD[from(ovlp)]
    df <- data.frame("seqnames" = seqnames(grN),
                     "start" = start(grN),
                     "end" = end(grN),
                     "score" = grN$mcols.mn/(grN$mcols.mn + grN$mcols.un),
                     "meth" = grN$mcols.mn,
                     "umeth" = grN$mcols.un)
    return(df)
  }
}

getTracks <- function(fWGBS, n, type){
  if (type == "d"){
    tracks <- .getDisease(fWGBS)
    tracks <- tracks[sample(1:length(tracks), size = n)]
  } 
  if (type == "n"){
    tracks <- .getNormal(fWGBS)
    tracks <- tracks[sample(1:length(tracks), size = n)]
  }
  return(tracks)
}







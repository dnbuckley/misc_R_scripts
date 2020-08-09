library(bsseq)
library(stringr)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Hsapiens.UCSC.hg19)
library(reshape2)
library(ggplot2)

# default for RRBS regions determined from OVCA RRBS 9-19
# reduce option combines adjacent regions
.getCpGLocs <- function(sequence = "CG", genome = "hg38", cores = 4){
  if (genome == "hg19") {
    Hsapiens <- BSgenome.Hsapiens.UCSC.hg19
  } else if (genome == "hg38") {
    Hsapiens <- BSgenome.Hsapiens.UCSC.hg38
  } else stop("invalid genome selection")
  chrs <- names(Hsapiens)[1:24]
  locs <- mclapply(chrs, function(x) start(matchPattern(sequence,
                                                        Hsapiens[[x]])), mc.cores = cores)
  gr <- lapply(1:24, function(x) GRanges(names(Hsapiens)[x],
                                         IRanges(locs[[x]], width = nchar(sequence))))
  gr <- do.call("c", gr)
  seqinfo(gr) <- seqinfo(Hsapiens)[seqlevels(gr)]
  return(gr)
}


.getPlexFromName <- function(name) {
  name <- gsub("MP", "plex", name)
  loc <- as.data.frame(str_locate(name, "(P|p)lex"))
  if (nrow(loc) != 1) {
    stop("Unable to determine location of 'plex' substring")
  }
  # if the name is plex1 or plex2
  plx <- as.numeric(substr(name, loc$end+1, loc$end+1))
  if (!is.na(plx)) {
    return(plx)
  }
  # if the name is plex-1 or plex_2
  plx <- as.numeric(substr(name, loc$end+2, loc$end+2))
  if (!is.na(plx)) {
    return(plx)
  }
  stop("Unable to determine plex number")
}

.getPlexFromName <- function(name) {
  plex <- gsub(".*plex(-|_|)", "", name, ignore.case = T)
  return((gsub("(_|-).*", "", plex)))
}

# .getSampleFromName <- function(name, is.std) {
#   message(name)
#   name <- gsub("MP", "plex", name)
#   loc <- as.data.frame(str_locate(name, "(P|p)lex"))
#   if (grepl("plex(_|-|\\.)", ignore.case = T, name)) {
#     sstr <- substr(name, loc$end+4, nchar(name))
#     message(sstr)
#   } else {
#     sstr <- substr(name, loc$end+3, nchar(name))
#     message(sstr)
#   }
#   return(gsub("_S[0-9]{1,3}.*", "", sstr))
# }

.getSampleFromName <- function(name, is.std) {
  name <- gsub("MP", "plex", name)
  if (is.std) {
    name <- gsub("(.*(plex|plex(_|-))|_S[0-9]{1,3}.*)", "", 
                 name, ignore.case = T)
  }
  else {
    name <- gsub("(.*(plex|plex(_|-))[0-9]{1,3}(_|-)|_S[0-9]{1,3}.*)", "", 
                 name, ignore.case = T)
  }
  return(name)
}

# attempt to get plex into a "plex" column
.formatGR <- function(gr) {
  reg <- granges(gr)
  mc <- mcols(gr)
  plexCol <- grep("plex", names(mc), ignore.case = T)
  if (length(plexCol) != 1) {
    stop("Unable to determine plex information from granges object")
  }
  reg$plex <- mc[, plexCol]
  n <- names(mc)[-plexCol]
  mcols(reg) <- cbind(mcols(reg), mc[, -plexCol])
  names(mcols(reg)) <- c("plex", n)
  return(reg)
}

.getCpGr <- function(genome) {
  if (Sys.info()[["user"]] == "dbuckley") {
    cpgr <- paste0("~/Desktop/salhia_lab/useful_granges/cpg_locations_", genome, ".rds")
    message("Reading ", cpgr)
    cpgr <- readRDS(cpgr)
  } else {
    message("finding CpG locations in ", genome, "...")
    cpgr <- .getCpGLocs(genome = genome, cores = detectCores()-2)
  }
  return(cpgr)
}

# .readPlex <- function(df, gr, cpgr) {
#   files <- as.character(df$file)
#   locs <- cpgr[from(findOverlaps(cpgr, gr))]
#   stopifnot(all(width(locs) == 2))
#   end(locs) <- start(locs)
#   bs <- read.bismark(files, loci = locs, nThread = detectCores()-2)
#   sampleNames(bs) <- df$sample
#   return(bs)
# }

.readPlex <- function(df, locs, is.std) {
  files <- as.character(df$file)
  bs <- read.bismark(files, loci = locs, nThread = detectCores()-2, 
                     strandCollapse = ifelse(all(strand(locs) == "*"), T, F))
  if (is.std) {
    sampleNames(bs) <- gsub(".*(_|-)", "", df$sample)
  } else {
    sampleNames(bs) <- df$sample
  }
  return(bs)
}


.getSampleManifest <- function(folder, is.std) {
  covs <- list.files(folder, "cov.gz$", full.names = T, recursive = T)
  tst <- gsub("MP", "plex", covs)
  if (!all(grepl("plex", ignore.case = T, tst))) {
    stop("Unable to determine plex from file name")
  }
  mfst <- data.frame(file = covs,
                     plex = sapply(covs, .getPlexFromName),
                     sample = sapply(covs, .getSampleFromName, is.std))
  t <- table(mfst$sample)
  if (length(unique(t)) != 1) {
    stop("Not all samples are represented an equal number of times")
  }
  return(mfst)
}

# function to stitch together a bsseq from multiple plexes.
# loci should be a non-overlapping granges object with
# at least one column specifying which plex each CpG is from 
# the folder should contain all .cov.gz files for every sample
# and each sample should have a file name with the following
# patterns:
# ....plex[0-9](_|-)sampleName_S....
# ....plex(-|_)[0-9](_|-)sampleName_S....
# ....MP[0-9](_|-)sampleName_S....
# ....MP(-|_)[0-9](_|-)sampleName_S....
# the function will attempt to guess the plex and sample, if 
# it's not doing that you can provide a manifest data.frame
# with columns 'file', 'plex' and 'sample' where file is the
# complete file path to the cov file, plex is the plex number
# and sample is the sampleName
stitchPlexes <- function(loci, folder = ".", mfst = NULL, 
                         qc = T, is.std = F) {
  # check for overlapping regions
  if (length(reduce(loci)) != length(loci)) {
    stop("Not all loci are unique")
  }
  if (!all(width(loci) == 1)){
    stop("Not all loci width = 1")
  }
  loci <- loci[order(loci)]
  sanity.check <- loci
  if (is.null(mfst)){
    mfst <- .getSampleManifest(folder, is.std)
  }
  mfst <- split(mfst, mfst$plex)
  loci <- .formatGR(loci)
  loci <- as.list(split(loci, loci$plex))
  loci <- loci[order(names(loci))]
  mfst <- mfst[order(names(mfst))]
  if (!identical(names(mfst), names(loci))) {
    stop("Plex numbers do not match")
  }
  message("Reading methylation information from ", length(loci), 
          " plexes and ", sum(sapply(mfst, nrow)), " files...")
  bsqs <- mapply(.readPlex, mfst, loci, is.std, SIMPLIFY = F)
  bs <- do.call(rbind, bsqs)
  bs <- realize(bs)
  bs <- bs[order(granges(bs))]
  # end(sanity.check) <- start(sanity.check)
  sanity.check <- sanity.check[order(sanity.check)]
  if (!identical(paste0(granges(bs)), paste0(granges(sanity.check)))) {
    stop("This should be an unreachable state...")
  }
  if (any(rowMeans(getCoverage(bs)) == 0)){
    warning("There were ", sum(rowMeans(getCoverage(bs)) == 0), " loci with 0 coverage in all samples.")
  }
  if (any(rowMedians(as.matrix(getCoverage(bs))) < 10)) {
    warning("There were ", sum(rowMedians(as.matrix(getCoverage(bs))) < 10), " loci with < 10x median coverage.")
  }
  if (qc) {
    df <- as.data.frame(getCoverage(bs))
    df$plex <- factor(granges(bs)$plex)
    df <- melt(df, id.vars = "plex")
    gp <- ggplot(df, aes(x = plex, y = value)) +
      geom_boxplot(fill = "grey50") +
      scale_y_continuous(trans = "log10") +
      ylab("log10 coverage") +
      theme_bw() +
      labs(title = "BSSEQ Coverage Boxplot")
    print(gp)
  }
  sampleNames(bs) <- gsub("\\.cov\\.gz", "", sampleNames(bs))
  return(bs)
}

















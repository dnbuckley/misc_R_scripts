library(openxlsx)
library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Hsapiens.UCSC.hg19)

# uses zymo bisulfite primer seeker to find primer loci
# for set of regions
getZymoPrimers <- function(gr, IDs, genome = "hg38",
                           bwsA = NULL, bwsB = NULL){
  source("~/misc_R_scripts/getNumCpGs.R")
  source("~/misc_R_scripts/getCpGLocs.R")
  if (genome == "hg19") BSgenome <- BSgenome.Hsapiens.UCSC.hg19
  if (genome == "hg38") BSgenome <- BSgenome.Hsapiens.UCSC.hg38
  if (any(width(gr) > 2000)) {
    message("Some regions are > 2000b, removing...")
    idx <- width(gr) < 2000
    gr <- gr[idx]
    IDs <- IDs[idx]
  }
  if (any(width(gr) < 115)) {
    message("Some regions are < 100b, removing...")
    idx <- width(gr) > 115
    gr <- gr[idx]
    IDs <- IDs[idx]
  }
  strand(gr) <- "*"
  pos <- gr
  neg <- gr
  strand(pos) <- "+"
  strand(neg) <- "-"
  posTab <- data.frame(seq = getSeq(BSgenome, pos, as.character = T),
                       IDs = IDs,
                       strand = strand(pos))
  negTab <- data.frame(seq = getSeq(BSgenome, neg, as.character = T),
                       IDs = IDs,
                       strand = strand(neg))
  table <- rbind.data.frame(posTab, negTab)
  write.table(table, "seqs.tsv", col.names = F, row.names = F, quote = F, sep = "\t")
  cmd <- "/opt/anaconda3/bin/python ~/Desktop/salhia_lab/automate_primer_design/zymo_design.py seqs.tsv"
  message("Launching Zymo Bisulfite Primer Seeker...")
  system(cmd)
  pri <- read.table("zymo_output.txt", sep = "\t", header = T)
  system("rm seqs.tsv")
  pri <- split(pri, pri$ID)
  grs <- as.list(split(gr, IDs))
  pri <- pri[order(names(pri))]
  grs <- grs[order(names(grs))]
  # remove rows with no primers found
  pri <- lapply(pri, function(x){x[rowSums(is.na(x)) == 0, ]})
  idx <- sapply(pri, function(x){nrow(x) == 0})
  if (any(idx)){
    message("No primers found for: \n", paste0(names(pri)[idx], collapse = "\n"))
    pri <- pri[!idx]
    grs <- grs[!idx]
  }
  stopifnot(identical(names(pri), names(grs)))
  message("Found ", sum(sapply(pri, nrow))/2, " potential primers.")
  message("Calculating genomic coordinates and CpG counts for ", genome, "...")
  cpgrpos <- getCpGLocs(strand = "+", genome = genome)
  cpgrneg <- getCpGLocs(strand = "-", genome = genome)
  all.primers <- mcmapply(.calculatePositions, grs, pri, 
                          MoreArgs = list(genome = genome, 
                                          cpgrpos = cpgrpos,
                                          cpgrneg = cpgrneg),
                          SIMPLIFY = F, mc.cores = detectCores()-2)
  all.primers <- do.call(rbind, all.primers)
  rownames(all.primers) <- NULL
  if (!is.null(bwsA) & !is.null(bwsB)) {
    all.primers$dmrvalue <- .getDMvalue(all.primers, bwsA, bwsB)
  }
  all.primers$mt.diff <- abs(all.primers$melting.temp.forward-all.primers$melting.temp.reverse)
  all.primers <- .classPrimers(all.primers)
  return(all.primers)
}

.getDMvalue <- function(primers, bwsA, bwsB) {
  source("~/misc_R_scripts/summarizeBigWigFilesOver.R")
  target.regions <- GRanges(primers$target.region)
  mA <- summarizeBigWigFilesOver(bwsA, target.regions)
  mA <- rowMeans(as.matrix(do.call(cbind, lapply(mA, mcols))))
  mB <- summarizeBigWigFilesOver(bwsB, target.regions)
  mB <- rowMeans(as.matrix(do.call(cbind, lapply(mB, mcols))))
  return(mA-mB)
}

getBestPrimer <- function(primers) {
  byID <- split(primers, primers$ID)
  byID <- lapply(byID, .getMostCpGs)
  byID <- lapply(byID, .getFewestPrimerCpGs)
  byID <- lapply(byID, .getHighestCpGDensity)
  byID <- lapply(byID, .getSmallest)
  df <- do.call(rbind, byID)
  return(df)
}

.classPrimers <- function(primers) {
  primers$tempid <- paste0(primers$ID, "_", primers$amp.number, ":", primers$origin.strand)
  classA <- .getClassA(primers)
  classB <- setdiff(.getClassB(primers), classA)
  classC <- setdiff(.getClassC(primers), c(classB, classA))
  classD <- setdiff(.getClassD(primers), c(classC, classB, classA))
  classF <- setdiff(primers$tempid, c(classD, classC, classB, classA))
  t <- sum(length(classA),
      length(classB),
      length(classC),
      length(classD),
      length(classF))
  stopifnot(t == nrow(primers))
  primers$class <- ifelse(primers$tempid %in% classA, "A",
                          ifelse(primers$tempid %in% classB, "B",
                          ifelse(primers$tempid %in% classC, "C",
                          ifelse(primers$tempid %in% classD, "D", 
                          ifelse(primers$tempid %in% classF, "F", NA)))))
  primers <- primers[, -grep("tempid", names(primers))]
  return(primers)
}


# no 3' Cpg, low mtdiff, high mt, high #CpG
.getClassA <- function(primers) {
  primers <- primers[!primers$cpg.3p.flag, ]
  primers <- primers[primers$mt.diff <= 2, ]
  primers <- primers[primers$melting.temp.forward > 60, ]
  primers <- primers[primers$melting.temp.reverse > 60, ]
  primers <- primers[primers$target.region.CpGs >= 5, ]
  return(primers$tempid)
}

# no 3' Cpg, low mtdiff, high #CpG
.getClassB <- function(primers){
  primers <- primers[!primers$cpg.3p.flag, ]
  primers <- primers[primers$mt.diff <= 2, ]
  primers <- primers[primers$target.region.CpGs >= 5, ]
  return(primers$tempid)
}

# no 3' Cpg, high #CpG
.getClassC <- function(primers){
  primers <- primers[!primers$cpg.3p.flag, ]
  primers <- primers[primers$target.region.CpGs >= 5, ]
  return(primers$tempid)
}

.getClassD <- function(primers) {
  primers <- primers[!primers$cpg.3p.flag, ]
  return(primers$tempid)
}

.getMostCpGs <- function(df) {
  df <- df[df$target.region.CpGs %in% max(df$target.region.CpGs), ]
  return(df)
}

.getFewestPrimerCpGs <- function(df) {
  total <- df$forward.primer.CpGs+df$reverse.primer.CpGs
  df <- df[total %in% min(total), ]
  return(df)
}

.getHighestCpGDensity <- function(df){
  density <- df$target.region.CpGs/width(GRanges(df$target.region))
  df <- df[density %in% max(density), ]
  return(df)
}

.getSmallest <- function(df) {
  df <- df[df$amplicon.width %in% min(df$amplicon.width), ]
  return(df)
}

.getAmpGRPos <- function(a, g, cpgr, ID) {
  stopifnot(nrow(a) == 2)
  f <- a[a$direction == "Forward", ]
  r <- a[a$direction == "Reverse", ]
  fP <- GRanges(seqnames = seqnames(g),
                ranges = IRanges(start = start(g) + f$start,
                                 width = f$size),
                strand = "+")
  rP <- GRanges(seqnames = seqnames(g),
                ranges = IRanges(start = start(g) + r$start - r$size + 1, 
                                 width = r$size),
                strand = "+")
  amp <- GRanges(seqnames = seqnames(g), 
                 ranges = IRanges(start = start(fP),
                                  end = end(rP)),
                 strand = "+")
  tgt <- GRanges(seqnames = seqnames(g),
                 ranges = IRanges(start = end(fP),
                                  end = start(rP)) - 1,
                 strand = "+")
  stopifnot(width(fP)+width(rP)+width(tgt) == width(amp))
  ret <- data.frame(ID = ID,
                    condition = unique(a$condition),
                    amp.number = unique(a$amp),
                    amplicon.width = width(amp),
                    origin.strand = "+",
                    amplicon = paste0(amp),
                    forward.primer = paste0(fP),
                    reverse.primer = paste0(rP),
                    target.region = paste0(tgt),
                    origin.region = paste0(g),
                    forward.sequence = f$seq,
                    reverse.sequence = r$seq,
                    forward.primer.CpGs = getNumCpGs(fP, cpgr),
                    reverse.primer.CpGs = getNumCpGs(rP, cpgr),
                    target.region.CpGs = getNumCpGs(tgt, cpgr),
                    melting.temp.forward = f$melting_temp,
                    melting.temp.reverse = r$melting_temp)
  return(ret)
}

.getAmpGRNeg <- function(a, g, cpgr, ID) {
  stopifnot(nrow(a) == 2)
  f <- a[a$direction == "Forward", ]
  r <- a[a$direction == "Reverse", ]
  fP <- GRanges(seqnames = seqnames(g),
                ranges = IRanges(end = end(g) - f$start,
                                 width = f$size),
                strand = "-")
  rP <- GRanges(seqnames = seqnames(g),
                ranges = IRanges(end = end(g) - r$start + r$size - 1, 
                                 width = r$size),
                strand = "-")
  amp <- GRanges(seqnames = seqnames(g), 
                 ranges = IRanges(end = end(fP),
                                  start = start(rP)),
                 strand = "-")
  tgt <- GRanges(seqnames = seqnames(g),
                 ranges = IRanges(end = start(fP),
                                  start = end(rP))-1,
                 strand = "-")
  stopifnot(width(fP)+width(rP)+width(tgt) == width(amp))
  ret <- data.frame(ID = ID,
                    condition = unique(a$condition),
                    amp.number = unique(a$amp),
                    amplicon.width = width(amp),
                    origin.strand = "-",
                    amplicon = paste0(amp),
                    forward.primer = paste0(fP),
                    reverse.primer = paste0(rP),
                    target.region = paste0(tgt),
                    origin.region = paste0(g),
                    forward.sequence = f$seq,
                    reverse.sequence = r$seq,
                    forward.primer.CpGs = getNumCpGs(fP, cpgr),
                    reverse.primer.CpGs = getNumCpGs(rP, cpgr),
                    target.region.CpGs = getNumCpGs(tgt, cpgr),
                    melting.temp.forward = f$melting_temp,
                    melting.temp.reverse = r$melting_temp)
  return(ret)
}

.calculatePositions <- function(g, p, genome, cpgrpos, cpgrneg){
  pos <- p[p$strand == "+", ]
  neg <- p[p$strand == "-", ]
  res <- NULL
  if (nrow(pos) > 0) {
    amps <- split(pos, as.character(pos$amp))
    pos <- do.call(rbind.data.frame, lapply(amps, .getAmpGRPos, g, cpgr = cpgrpos, ID = unique(pos$ID)))
    flag <- .flag3PrimeCpG(pos, cpgrpos, "+")
    pos$cpg.3p.flag <- flag$CpG.flag.either
    res[[length(res)+1]] <- pos
  } 
  if (nrow(neg) > 0) {
    amps <- split(neg, as.character(neg$amp))
    neg <- do.call(rbind.data.frame, lapply(amps, .getAmpGRNeg, g, cpgr = cpgrneg, ID = unique(neg$ID)))
    flag <- .flag3PrimeCpG(neg, cpgrneg, "-")
    neg$cpg.3p.flag <- flag$CpG.flag.either
    res[[length(res)+1]] <- neg
  }
  return(do.call(rbind.data.frame, res))
}

# get the 2/3 of the primer at the 3' end
# and flag if it has a cpg
.flag3PrimeCpG <- function(p, cpgr, strand) {
  f <- GRanges(p$forward.primer)
  r <- GRanges(p$reverse.primer)
  if (strand == "+") {
    start(f) <- start(f)+floor(width(f)*(1/3))
    start(r) <- start(r)+floor(width(r)*(1/3))
  }
  else if (strand == "-") {
    end(f) <- end(f)-floor(width(f)*(1/3))
    end(r) <- end(r)-floor(width(r)*(1/3))
  } else stop("Invalid strand")
  ovlpf <- findOverlaps(f, cpgr)
  ovlpr <- findOverlaps(r, cpgr)
  ret <- data.frame(forward.primer.CpG.flag = 1:length(f) %in% from(ovlpf),
                    reverse.primer.CpG.flag = 1:length(r) %in% from(ovlpr))
  ret$CpG.flag.either <- ret$forward.primer.CpG.flag | ret$reverse.primer.CpG.flag
  return(ret)
}















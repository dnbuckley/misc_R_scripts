library(openxlsx)
library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Hsapiens.UCSC.hg19)

# uses zymo bisulfite primer seeker to find primer loci
# for set of regions
# gr <- tnbc
# genome <- "hg19"
# nZymoDrivers <- 1
getZymoPrimers <- function(gr, IDs, genome = "hg38",
                           bwsA = NULL, bwsB = NULL,
                           nZymoDrivers = 1){
  source("~/misc_R_scripts/getNumCpGs.R")
  source("~/misc_R_scripts/getCpGLocs.R")
  # clear out any files from prior runs if they were canceled or failed
  system("rm zymo_input_* zymo_output_*")
  if (genome == "hg19") BSgenome <- BSgenome.Hsapiens.UCSC.hg19
  if (genome == "hg38") BSgenome <- BSgenome.Hsapiens.UCSC.hg38
  original.granges <- GRanges(gr)
  while (any(width(gr) > 2000)) {
    message("Some regions are > 2000b, shrinking by 20b...")
    idx <- width(gr) > 2000
    start(gr[idx]) <- start(gr[idx])+10
    end(gr[idx]) <- end(gr[idx])-10
  }
  while (any(width(gr) < 115)) {
    message("Some regions are < 115, expanding by 20b...")
    idx <- width(gr) < 115
    start(gr[idx]) <- start(gr[idx])-10
    end(gr[idx]) <- end(gr[idx])+10
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
  tabs <- split(table, .splitTable(table, nZymoDrivers))
  message("Launching Zymo Bisulfite Primer Seeker Using ", nZymoDrivers, " drivers...")
  mcmapply(.launchDriver, tabs, names(tabs), mc.cores = nZymoDrivers)
  zymoFiles <- list.files(".", "^zymo_output.*tsv")
  pri <- do.call(rbind, lapply(zymoFiles, read.table, sep = "\t", header = T))
  system("rm zymo_input_* zymo_output_*")
  pri$amp <- paste0(pri$amp, "-C", pri$condition)
  pri <- split(pri, pri$ID)
  pri <- lapply(pri, .removeIdenticalPrimers)
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
  message("Calculating genomic coordinates and CpG counts using ", genome, "...")
  cpgrpos <- getCpGLocs(strand = "+", genome = genome)
  cpgrneg <- getCpGLocs(strand = "-", genome = genome)
  all.primers <- mapply(.calculatePositions, grs, pri, 
                        MoreArgs = list(genome = genome, 
                                        cpgrpos = cpgrpos,
                                        cpgrneg = cpgrneg),
                        SIMPLIFY = F)
  all.primers <- do.call(rbind, all.primers)
  rownames(all.primers) <- NULL
  if (!is.null(bwsA) & !is.null(bwsB)) {
    all.primers$target.region.dmrvalue <- .getDMvalue(all.primers, bwsA, bwsB)
  }
  all.primers <- .classPrimers(all.primers)
  all.primers$method <- ifelse(width(GRanges(all.primers$origin.region)) >= 250, "tiling", "regular")
  all.primers$was.resized <- .checkAgainstOriginal(all.primers, original.granges)
  return(all.primers)
}

.launchDriver <- function(tab, n) {
  name <- paste0("zymo_input_", n, ".tsv")
  write.table(tab, name, col.names = F, row.names = F, quote = F, sep = "\t")
  cmd <- paste0("/Users/dbuckley/anaconda3/bin/python ~/misc_R_scripts/zymo_design.py ", name, " ", n)
  system(cmd)
}

.checkAgainstOriginal <- function(df, ref) {
  pr <- GRanges(df$origin.region)
  ovlp <- findOverlaps(pr, ref, type = "equal")
  was.resized <- rep(T, nrow(df))
  was.resized[from(ovlp)] <- F
  return(was.resized)
}

.splitTable <- function(table, nZymoDrivers) {
  df <- data.frame(id = unique(table$IDs))
  return(floor(as.numeric(rownames(df)) %% nZymoDrivers))
}

# remove identical primers found under different conditions
.removeIdenticalPrimers <- function(p) {
  amps <- split(p, paste0(p$amp, p$strand))
  stopifnot(all(sapply(amps, nrow) == 2))
  df <- data.frame(name = names(amps),
                   loc = sapply(amps, function(x){paste0(x$ID, x$seq, x$start, collapse = "")}),
                   condition = sapply(amps, function(x) unique(x$condition)))
  t <- as.data.frame(table(df$loc))
  rep <- df[df$loc %in% t$Var1[t$Freq > 1], ]
  rep <- split(rep, rep$loc)
  # get rid of the one with higher (worst) condition
  rep <- do.call(rbind, lapply(rep, function(x){x[x$condition %in% max(x$condition), ]}))
  amps <- amps[!names(amps) %in% rep$name]
  p <- do.call(rbind, amps)
  rownames(p) <- NULL
  return(p)
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

.classPrimers <- function(primers) {
  primers$tempid <- paste0(primers$ID, "_", primers$amp.number, ":", primers$origin.strand)
  classA <- .getClassA(primers)
  classB <- setdiff(.getClassB(primers), classA)
  classC <- setdiff(.getClassC(primers), c(classB, classA))
  classD <- setdiff(.getClassD(primers), c(classC, classB, classA))
  classE <- setdiff(.getClassE(primers), c(classD, classC, classB, classA))
  classF <- setdiff(primers$tempid, c(classE, classD, classC, classB, classA))
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
  primers <- primers[primers$melt.temp.diff <= 2, ]
  primers <- primers[primers$condition == 1, ]
  primers <- primers[primers$target.region.CpGs >= 5, ]
  return(primers$tempid)
}

# no 3' Cpg, low mtdiff, high #CpG
.getClassB <- function(primers){
  primers <- primers[!primers$cpg.3p.flag, ]
  primers <- primers[primers$melt.temp.diff <= 2, ]
  primers <- primers[primers$condition == 2, ]
  primers <- primers[primers$target.region.CpGs >= 5, ]
  return(primers$tempid)
}

# no 3' Cpg, high #CpG
.getClassC <- function(primers){
  primers <- primers[!primers$cpg.3p.flag, ]
  primers <- primers[primers$melt.temp.diff <= 2, ]
  primers <- primers[primers$condition == 3, ]
  primers <- primers[primers$target.region.CpGs >= 5, ]
  return(primers$tempid)
}

.getClassD <- function(primers) {
  primers <- primers[!primers$cpg.3p.flag, ]
  primers <- primers[primers$melt.temp.diff <= 2, ]
  primers <- primers[primers$condition == 4, ]
  primers <- primers[primers$target.region.CpGs >= 5, ]
  return(primers$tempid)
}

.getClassE <- function(primers) {
  primers <- primers[!primers$cpg.3p.flag, ]
  primers <- primers[primers$target.region.CpGs >= 3, ]
}

.getNumCpGforward <- function(fps) {
  fps <- strsplit(fps, "")
  if (any(grepl("R", fps))) {
    stop("Found 'R' in forward primer, should not happen.")
  }
  return(sum(grepl("Y", fps)))
}

.getNumCpGreverse <- function(rps) {
  rps <- strsplit(rps, "")
  if (any(grepl("Y", rps))) {
    stop("Found 'Y' in forward primer, should not happen.")
  }
  return(sum(grepl("R", rps)))
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
                    forward.primer.CpGs = .getNumCpGforward(f$seq),
                    reverse.primer.CpGs = .getNumCpGreverse(r$seq),
                    target.region.CpGs = getNumCpGs(tgt, cpgr),
                    melting.temp.forward = f$melting_temp,
                    melting.temp.reverse = r$melting_temp,
                    melt.temp.diff = abs(f$melting_temp-r$melting_temp))
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
                    forward.primer.CpGs = .getNumCpGforward(f$seq),
                    reverse.primer.CpGs = .getNumCpGreverse(r$seq),
                    target.region.CpGs = getNumCpGs(tgt, cpgr),
                    melting.temp.forward = f$melting_temp,
                    melting.temp.reverse = r$melting_temp,
                    melt.temp.diff = abs(f$melting_temp-r$melting_temp))
  return(ret)
}

.calculatePositions <- function(g, p, genome, cpgrpos, cpgrneg){
  message("Calculating genomic coordinates and CpG counts for ", unique(p$ID), "...")
  pos <- p[p$strand == "+", ]
  neg <- p[p$strand == "-", ]
  res <- NULL
  if (nrow(pos) > 0) {
    amps <- split(pos, as.character(pos$amp))
    pos <- do.call(rbind.data.frame, mclapply(amps, .getAmpGRPos, g, cpgr = cpgrpos, 
                                              ID = unique(pos$ID), mc.cores = detectCores()-2))
    # flag <- .flag3PrimeCpG(pos, cpgrpos, "+")
    # flag <-.flag3PrimeCpG(pos)
    # pos$cpg.3p.flag <- flag$CpG.flag.either
    res[[length(res)+1]] <- pos
  } 
  if (nrow(neg) > 0) {
    amps <- split(neg, as.character(neg$amp))
    neg <- do.call(rbind.data.frame, mclapply(amps, .getAmpGRNeg, g, cpgr = cpgrneg, 
                                              ID = unique(neg$ID), mc.cores = detectCores()-2))
    # flag <- .flag3PrimeCpG(neg, cpgrneg, "-")
    # flag <-.flag3PrimeCpG(neg)
    # neg$cpg.3p.flag <- flag$CpG.flag.either
    res[[length(res)+1]] <- neg
  }
  p <- do.call(rbind.data.frame, res)
  p$cpg.3p.flag <- .flag3PrimeCpG(p)
  return(p)
}

# get the 2/3 of the primer at the 3' end
# and flag if it has a cpg
.flag3PrimeCpG <- function(p) {
  f.2.3 <- sapply(p$forward.sequence, function(x){
    x <- substr(x, floor((nchar(x)*1)/3), nchar(x))
    grepl("Y", x)
  })
  r.2.3 <- sapply(p$reverse.sequence, function(x){
    x <- substr(x, floor((nchar(x)*1)/3), nchar(x))
    grepl("R", x)
  })
  return(f.2.3 | r.2.3)
}

# .flag3PrimeCpG <- function(p, cpgr, strand) {
#   f <- GRanges(p$forward.primer)
#   r <- GRanges(p$reverse.primer)
#   if (strand == "+") {
#     start(f) <- start(f)+floor(width(f)*(1/3))
#     start(r) <- start(r)+floor(width(r)*(1/3))
#   }
#   else if (strand == "-") {
#     end(f) <- end(f)-floor(width(f)*(1/3))
#     end(r) <- end(r)-floor(width(r)*(1/3))
#   } else stop("Invalid strand")
#   ovlpf <- findOverlaps(f, cpgr)
#   ovlpr <- findOverlaps(r, cpgr)
#   ret <- data.frame(forward.primer.CpG.flag = 1:length(f) %in% from(ovlpf),
#                     reverse.primer.CpG.flag = 1:length(r) %in% from(ovlpr))
#   ret$CpG.flag.either <- ret$forward.primer.CpG.flag | ret$reverse.primer.CpG.flag
#   return(ret)
# }















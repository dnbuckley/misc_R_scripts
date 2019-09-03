library(Gviz)
library(org.Hs.eg.db)
library(Homo.sapiens)
library(BSgenome.Hsapiens.UCSC.hg19)

# print Gviz of one granges
printGviz <- function(gr, bwFolder, bwRegex = NULL, gen = "hg19", return = F){
  bws <- list.files(bwFolder, "bw$")
  if (!is.null(bwRegex)) bws <- bws[grepl(bwRegex, bws)]
  gatrack <- GenomeAxisTrack()
  prevChr <- -1
  gr <- gr[order(seqnames(gr))]
  failedGRs <- list()
  ext <- end(gr) - start(gr)
  ext <- round(ext*(1/3))
  grE <- data.frame(chr = seqnames(gr),
                    start = start(gr) - ext,
                    end = end(gr) + ext,
                    strand = "*")
  grE <- GRanges(grE)
  dtracks <- list()
  chr <- as.character(seqnames(gr))
  message("Loading ", chr)
  # idk how this chunk of code works but it does so don't fuck with it
  ucscGenes <- UcscTrack(genome="hg19", table="ncbiRefSeq", track = 'NCBI RefSeq', trackType="GeneRegionTrack",
                         chromosome=chr, rstarts = "exonStarts", rends = "exonEnds",
                         gene = "name", symbol = 'name', transcript = "name", transcriptAnnotation = "symbol",
                         strand = "strand", showID = TRUE, geneSymbol = TRUE, name = "UCSC genes", fontsize = 12, 
                         from = start(grE), to = end(grE))
  z <- ranges(ucscGenes)
  # in case there are no genes within range
  if (ncol(mcols(z)) > 0){
    mcols(z)$symbol <- mapIds(org.Hs.eg.db, gsub("\\.[1-9]$", "", mcols(z)$symbol), "SYMBOL","REFSEQ")
    ucscGenes2 <- ucscGenes
    ranges(ucscGenes2) <- z
  }
  itrack <- IdeogramTrack(genome = gen, chromosome = chr)
  strack <- SequenceTrack(Hsapiens, chromosome = chr)
  if(exists("ucscGenes2")) dtracks <- c(itrack, gatrack, strack, ucscGenes2)
  else dtracks <- c(itrack, gatrack, strack)
  mtracks <- list()
  vals <- GRangesList(mclapply(paste0(bwFolder, bws), import, which=grE))
  names(vals) <- gsub(".bw", "", bws)
  for (j in 1:length(vals)){
    name <- names(vals)[j]
    mtrack <- DataTrack(vals[[j]], start = (start(gr)-ext), end = (end(gr)+ext),
                        genome = gen, chromosome = chr, name = name,
                        type = "histogram", fill = "black")
    mtracks <- c(mtracks, mtrack)
  }
  dtracks <- c(dtracks, mtracks)
  if (return) return(dtracks)
  else plotTracks(dtracks, cex = 1, from = start(grE), to = end(grE))
}

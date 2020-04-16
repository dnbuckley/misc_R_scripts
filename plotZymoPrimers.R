library(Gviz)
library(GenomicRanges)
library(gridExtra)
library(pdftools)
source("~/misc_R_scripts/getCpGLocs.R")
source("~/misc_R_scripts/getBSgenomeObj.R")

plotPrimers <- function(primers, pdfFileName, genome, bwsA, bwsB, simpleFilter = T){
  byID <- split(primers, as.character(primers$ID))
  byID <- lapply(byID, function(x){x[order(x$class), ]})
  if (simpleFilter) {
    byID <- lapply(byID, .simpleFilter)
  }
  cpgs <- getCpGLocs(genome = genome)
  plots <- lapply(byID, .gvizPrimerPlot, genome, cpgs, bwsA, bwsB)
  tables <- lapply(byID, .getPrimerTable)
  dir.create("all_pdfs/")
  for (i in 1:length(byID)) {
    df <- byID[[i]]
    kableExtra::kable(df, "latex")
    # grb <- arrangeGrob(plots[[i]], tables[[i]], nrow = 2,
    #                    top = textGrob(paste0(unique(df$ID), "\n", 
    #                                          unique(df$origin.region)), 
    #                                   gp=gpar(fontsize=20,font=8)))
    
  pdf(paste0("all_pdfs/", unique(df$ID), ".pdf"), height = 15, width = 26)
  # grid.arrange(grb)
  grid.arrange(plots[[i]])
  kableExtra::kable(tables[[i]], "latex")
  dev.off()
  }
  pdfs <- list.files("all_pdfs/", full.names = T)
  pdf_combine(pdfs, output = pdfFileName)
}

.getPrimerTable <- function(df) {
  df$ID <- seq(1:nrow(df))
  names(df)[names(df) == "ID"] <- "plot.amp.label"
  df <- df[, !names(df) %in% c("amp.number", "origin.regon")]
  # t <- tableGrob(df, theme = ttheme_default(base_size = 10), rows = NULL)
  # return(t)
  return(df)
}

.getDatTrackFromBWs <- function(b, g, title = "") {
  d <- lapply(b, import, which = g)
  m <- rowMeans(as.matrix(do.call(cbind, lapply(d, mcols))), na.rm = T)
  g <- paste0(do.call(c, lapply(d, granges)))
  g <- GRanges(unique(g))
  t <- DataTrack(g, 
                 data = m, 
                 name = title, 
                 type = "histogram", 
                 fill = "black", 
                 background.panel = "#FFFEDB")
  return(t)
}

.simpleFilter <- function(df) {
  if (any(df$class != "F")) {
    df <- df[df$class != "F", ]
  }
  return(df)
}

.gvizPrimerPlot <- function(df, genome, cpgs, bwsA, bwsB){
  HS <- getBSgenomeObj(genome)
  origin <- GRanges(unique(df$origin.region))
  message(paste0(origin))
  plotRange <- GRanges(seqnames = seqnames(origin), 
                       ranges = IRanges(start = start(origin)-floor(width(origin)*0.25),
                                        end = end(origin)+floor(width(origin)*0.25)))
  f <- GRanges(df$forward.primer)
  r <- GRanges(df$reverse.primer)
  g <- c(f, r)
  gatrack <- GenomeAxisTrack()
  pTrack <- AnnotationTrack(range = g, 
                            group = rep(1:nrow(df), 2),
                            showId = T, 
                            # stacking = "full",
                            collapse = F,
                            name = "Primers Found",
                            background.title = "black",
                            shape = "smallArrow")
  dmrTrack <- AnnotationTrack(origin, 
                              fill = "royalblue",
                              name = "DMR")
  itrack <- IdeogramTrack(chromosome = as.character(seqnames(origin)),
                          genome = genome)
  cpgs <- cpgs[from(findOverlaps(cpgs, plotRange))]
  cpgtrack <- AnnotationTrack(range = granges(cpgs),
                              fill = "red",
                              color = "red",
                              name = "CpGs", 
                              background.panel = "grey90")
  dTrackA <- .getDatTrackFromBWs(bwsA, plotRange, "Beta A")
  dTrackB <- .getDatTrackFromBWs(bwsB, plotRange, "Beta B")
  HS <- getBSgenomeObj(genome)
  strack <- SequenceTrack(HS, chromosome = seqnames(plotRange))
  plotTracks(c(itrack,
               gatrack,
               dmrTrack,
               dTrackA,
               dTrackB,
               strack,
               cpgtrack,
               pTrack),
             from = start(plotRange),
             to = end(plotRange))
  g <- grid.grab()
  return(g)
}



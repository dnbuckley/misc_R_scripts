library(ggplot2)
library(reshape2)
library(ggrepel)
library(ggpubr)
library(gridExtra)
library(Biobase)
source("~/misc_R_scripts/getKnownGenes.R")

.convertCN <- function(cn, plotPloidy) {
  if (median(as.numeric(unlist(cn)), na.rm = T) > 0.5) {
    warning("Median 'copynumber' field is high. Be sure you have called makeCgh(QDNAobj)")
  }
  if (plotPloidy) {
    cn <- (2^cn)*2
  }
  return(cn)
}

plotFocalCNA <- function(geneName, copyNumbersCalled,
                         expand = 5e5, plotPloidy = T,
                         geomSmooth = F, ylims = NULL,
                         genome = "hg19", returnGrob = T) {
  message("Plotting ", geneName, "...")
  cn <- assayDataElement(copyNumbersCalled, "copynumber")
  colnames(cn) <- gsub(" vs. .*", "", colnames(cn))
  genes <- getKnownGenes(genome = genome)
  gene <- genes[genes$SYMBOL %in% geneName]
  stopifnot(length(gene) == 1)
  plotRange <- GRanges(seqnames = seqnames(gene),
                       ranges = IRanges(start = start(gene) - expand,
                                        end = end(gene) + expand))
  cn <- cn[from(findOverlaps(GRanges(paste0("chr", rownames(cn))), plotRange)), ]
  cn <- .convertCN(cn, plotPloidy)
  yUnits <- ifelse(plotPloidy, "Copy Number", "CGH")
  cn <- melt(cn)
  df <- data.frame(loc = paste0("chr", cn$Var1), 
                   start = gsub("(.*:|-.*)", "", cn$Var1),
                   sample = cn$Var2,
                   cgh = cn$value)
  dfs <- split(df, df$sample)
  # plots
  cnPlot <- function(x, n, g, s, e, gn){
    o <- findOverlaps(GRanges(x$loc), g)
    x$color <- "black"
    x$color[from(o)] <- "red"
    x$label <- ""
    x$label[x$color == "red"][1] <- gn
    x$start <- as.numeric(x$start)
    xp <- x[x$start > s & x$start < e, ]
    w <- floor((max(xp$start)-min(xp$start))/1e3)
    gp <- ggplot(xp, aes(x = start, y = cgh, label = label)) +
      geom_point(color = xp$color, size = 0.7) +
      theme_pubr() +
      theme(axis.text.x = element_blank())
    labs(title = n) 
    return(gp)
  }
  plts <- mapply(cnPlot, dfs, names(dfs),
                 MoreArgs = list(g = gene, s = start(plotRange), e = end(plotRange), gn = geneName),
                 SIMPLIFY = F)
  plts <- lapply(plts, function(x){x + ylab(yUnits)})
  if (geomSmooth) {
    plts <- lapply(plts, function(x){x + geom_smooth(color = "blue")})
  }
  if (!is.null(ylims)) {
    plts <- lapply(plts, function(x, y){x + ylim(y)}, ylims)
  }
  if (returnGrob) {
    plts <- do.call(arrangeGrob, plts)
  }
  return(plts)
}

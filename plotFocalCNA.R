library(ggplot2)
library(reshape2)
library(ggrepel)
library(ggpubr)
library(gridExtra)
library(Biobase)
source("~/misc_R_scripts/getKnownGenes.R")
source("~/misc_R_scripts/combineListOfGR.R")
source("~/misc_R_scripts/combineInRange.R")

.reBin <- function(df, spl, recalcRanges = T) {
  df <- as.data.frame(df)
  byChr <- split(df, gsub(":.*", "", rownames(df)))
  # split by Chr
  byChr <- lapply(byChr, function(x, spl){
    split(x, rep(1:ceiling(nrow(x)/spl), each = spl))
  }, spl)
  # get new ranges
  if (recalcRanges){
    chnks <- lapply(byChr, function(x){
      y <- mclapply(x, function(y){
        gr <- GRanges(rownames(y))
        gr <- combineInRange(gr, 1)
      }, mc.cores = 10)
      z <- combineListOfGR(y)
      return(z)
    })
    chnks <- combineListOfGR(chnks)
  }
  # average across
  byChr <- lapply(byChr, function(x){
    x <- mclapply(x, colMeans, mc.cores = 10, na.rm = T)
    x <- do.call(rbind, x)
    return(x)
  })
  df <- do.call(rbind, byChr)
  if (recalcRanges){
    rownames(df) <- paste0(chnks)
  }
  return(as.data.frame(df))
}

.convertCN <- function(cn, plotPloidy) {
  if (median(as.numeric(unlist(cn)), na.rm = T) > 0.5) {
    warning("Median 'copynumber' field is high. Be sure you have called makeCgh(QDNAobj)")
  }
  if (plotPloidy) {
    cn <- (2^cn)*2
  }
  return(cn)
}

# geneName = HUGO gene symbol
# copyNumbersCalled = QDNAseq obj with makeCgh called
# expand = 1/2 width of plot window
# plotPloidy = plot actual copy number ie. normal = 2
# geomSmooth = should function add geom_smooth()
# genome = self explanatory
# returnGrob = should function combine plots using GridExtra
# combineNbins = combine (average) n # bins in an attempt to simplify/smooth output
# pointSize = size of geom_point()
# lineRange = plot a line at hight n showing gene location
plotFocalCNA <- function(geneName, copyNumbersCalled,
                         expand = 5e5, plotPloidy = T,
                         geomSmooth = F, ylims = NULL,
                         genome = "hg19", returnGrob = T,
                         combineNbins = NULL, pointSize = 1,
                         lineRange = NULL, lineRangeTbarHeight = 0.1,
                         lineRangeSize = 1) {
  message("Plotting ", geneName, "...")
  cn <- assayDataElement(copyNumbersCalled, "copynumber")
  colnames(cn) <- gsub(" vs. .*", "", colnames(cn))
  genes <- suppressMessages(getKnownGenes(genome = genome))
  gene <- genes[genes$SYMBOL %in% geneName]
  stopifnot(length(gene) == 1)
  plotRange <- GRanges(seqnames = seqnames(gene),
                       ranges = IRanges(start = start(gene) - expand,
                                        end = end(gene) + expand))
  cn <- cn[from(findOverlaps(GRanges(paste0("chr", rownames(cn))), plotRange)), ]
  if (!is.null(combineNbins)) {
    if (combineNbins <= 1) {
      warning("Invalid combineNbins value, skipping.")
    } else {
      w <- width(GRanges(rownames(cn)[1]))
      message("Converting bins of size ", w, "b to size ", w*combineNbins, "b.")
      cn <- .reBin(df = cn, spl = combineNbins)
    }
  }
  cn <- .convertCN(cn, plotPloidy)
  yUnits <- ifelse(plotPloidy, "Copy Number", "CGH")
  cn <- melt(as.matrix(cn))
  df <- data.frame(loc = paste0("chr", cn$Var1), 
                   start = gsub("(.*:|-.*)", "", cn$Var1),
                   sample = cn$Var2,
                   cgh = cn$value)
  dfs <- split(df, df$sample)
  # plots
  cnPlot <- function(x, n, g, s, e, gn, ps){
    o <- findOverlaps(GRanges(x$loc), g)
    x$color <- "black"
    x$color[from(o)] <- "red"
    x$label <- ""
    x$label[x$color == "red"][1] <- gn
    x$start <- as.numeric(x$start)
    xp <- x[x$start > s & x$start < e, ]
    w <- floor((max(xp$start)-min(xp$start))/1e3)
    gp <- ggplot(xp, aes(x = start, y = cgh)) +
      geom_point(color = xp$color, size = ps) +
      theme_pubr() +
      theme(axis.text.x = element_blank()) +
      labs(title = n)
    return(gp)
  }
  plts <- mapply(cnPlot, dfs, names(dfs),
                 MoreArgs = list(g = gene, 
                                 s = start(plotRange), 
                                 e = end(plotRange), 
                                 gn = geneName,
                                 ps = pointSize),
                 SIMPLIFY = F)
  plts <- lapply(plts, function(x){x + ylab(yUnits)})
  if (geomSmooth) {
    plts <- lapply(plts, function(x){x + geom_smooth(color = "blue")})
  }
  if (!is.null(ylims)) {
    plts <- lapply(plts, function(x, y){x + ylim(y)}, ylims)
  }
  if (!is.null(lineRange)) {
    plts <- lapply(plts, function(x, geneName) {
      x + 
        geom_segment(aes(x = start(gene), xend = end(gene),
                           y = lineRange, yend = lineRange),
                       color = "red", size = lineRangeSize) +
        geom_segment(aes(x = start(gene), xend = start(gene),
                         y = lineRange-lineRangeTbarHeight, yend = lineRange+lineRangeTbarHeight),
                     color = "red", size = 0.5) +
        geom_segment(aes(x = end(gene), xend = end(gene),
                         y = lineRange-lineRangeTbarHeight, yend = lineRange+lineRangeTbarHeight),
                     color = "red", size = 0.5) +
        geom_label(label = geneName, x = start(gene) + width(gene)/2, y = lineRange + lineRangeTbarHeight*3, 
                   inherit.aes = F)
    }, geneName)
  }
  if (returnGrob) {
    plts <- do.call(arrangeGrob, plts)
  }
  return(plts)
}

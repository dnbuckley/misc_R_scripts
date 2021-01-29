library(ggplot2)
library(Biobase)

# generic gplot of QDNAseq sample from copyNumersCalled
plotQDNAseq <- function(copyNumbersCalled, column, title = "CNV plot", 
                        ylim = c(0,2), plotLine = F, return = T, print = F){
  df <- data.frame(copynum = assayDataElement(copyNumbersCalled, "copynumber")[, column],
                   segs = assayDataElement(copyNumbersCalled, "segmented")[, column],
                   call = assayDataElement(copyNumbersCalled, "calls")[, column])
  df$loc <- as.numeric(gsub(".*-", "", rownames(df)))
  df$chr <- as.numeric(gsub(":.*", "", rownames(df)))
  df <- df[!is.na(df$chr), ]
  df$call[is.na(df$call)] <- 0
  pointCols <- c("-2" = "blue", "-1" = "deepskyblue3", "0" = "grey50", "1" = "red", "2" = "darkred")
  gp <- ggplot(df, aes(x = loc, y = copynum)) + 
    geom_point(size = 0.1, aes(color = factor(call))) +
    theme_bw() +
    facet_grid(.~ chr, scales = "free", space = "free") +
    scale_color_manual(values = pointCols) +
    ylab("Log Ratio") +
    labs(title = title) +
    ylim(ylim) +
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          strip.text.x = element_text(size = 20),
          title = element_text(size = 20),
          axis.text.y.left = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          panel.grid = element_blank(),
          axis.ticks.x = element_blank(), 
          legend.text = element_blank(),
          legend.title = element_blank(),
          legend.position = 'none')
  if (plotLine) gp <- gp + geom_line(aes(x = loc, y = segs), size = 1)
  if (print) print(gp)
  else if (return) return(gp)
  else message("specify to print or return")
}
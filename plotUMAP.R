library(ggplot2)
library(ggrepel)

# plot a UMAP
plotUMAP <- function(umap, 
                     group = NULL, 
                     label = F, 
                     density = F, 
                     point.size = 1, 
                     geom.density.bins = 10){
  df <- data.frame(x = umap$layout[, 1], 
                   y = umap$layout[, 2])
  if (is.null(group)) group <- rep("None", nrow(df))
  df$group <- group
  if (label) {
    gp <- ggplot(df, aes(x = x, y = y, color = group, label = rownames(df))) +
      geom_point(size = point.size) +
      theme_bw() +
      geom_text_repel()
  }
  else {
    gp <- ggplot(df, aes(x = x, y = y, color = group)) +
      geom_point(size = point.size) +
      theme_bw()
  }
  if (density){
    gp <- gp + geom_density_2d(inherit.aes = T, linemitre = 2, bins = geom.density.bins)
  }
  return(gp)
}

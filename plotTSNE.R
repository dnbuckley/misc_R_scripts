library(ggplot2)
library(ggrepel)

# plot a TSNE
plotTSNE <- function(ts, group = NULL, label = F, density = F){
  df <- data.frame(x = ts[, 1], 
                   y = ts[, 2])
  if (is.null(group)) group <- rep("None", nrow(df))
  df$group <- group
  
  if (label) {
    gp <- ggplot(df, aes(x = x, y = y, color = group, label = rownames(df))) +
      geom_point() +
      theme_bw() +
      geom_text_repel()
  }
  else {
    gp <- ggplot(df, aes(x = x, y = y, color = group)) +
      geom_point() +
      theme_bw()
  }
  if (density){
    gp <- gp + geom_density_2d(inherit.aes = T, linemitre = 2)
  }
  return(gp)
}

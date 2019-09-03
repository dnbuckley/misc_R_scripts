library(ggplot2)
library(ggrepel)

# plot a UMAP
plotUMAP <- function(umap, group = NULL, label = F){
  df <- data.frame(x = umap$layout[, 1], 
                   y = umap$layout[, 2])
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
  return(gp)
}

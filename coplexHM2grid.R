library(ComplexHeatmap)
library(gridExtra)

# functions copied from stack overflow
.grab_grob <- function(){
  grid.grab()
}
.drawGridHeatmap  <- function(hm) {
  draw(hm)
  .grab_grob()
}

# should return object that can be added to a
# gridExtra style grob, it will also plot the hm
# idk how it works so I can't stop that part
complexHM2grid <- function(heatmap){
  heatmap <- lapply(list(heatmap), .drawGridHeatmap)[[1]]
  return(heatmap)
}
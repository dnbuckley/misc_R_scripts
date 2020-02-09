library(rtracklayer)

bigWigToBedGraph <- function(fileName){
  gr <- import(fileName)
  name <- gsub("\\.bw$", ".bg")
  export.bedGraph(gr, name)
}
library(rtracklayer)

bigWigToBedGraph <- function(fileName){
  gr <- import(fileName)
  name <- gsub("\\.bw$", ".bg", fileName)
  export.bedGraph(gr, name)
}

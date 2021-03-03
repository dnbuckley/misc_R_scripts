#combine all GRanges within a given range
combineInRange <- function(gr, range){
  end(gr) <- end(gr) + range
  gr <- reduce(gr)
  end(gr) <- end(gr) - range
  return(gr)
}

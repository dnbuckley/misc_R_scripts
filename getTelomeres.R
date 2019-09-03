library(BSgenome)
library(rtracklayer)

.chrAsNum <- function(tbl){
  tbl$chrom <- gsub("chr", "", tbl$chrom)
  tbl$chrom[tbl$chrom=="X"] <- 23
  tbl$chrom[tbl$chrom=="Y"] <- 24
  tbl$chrom <- as.numeric(tbl$chrom)
  tbl[order(tbl$chrom),]
}

getTelomeres <- function( genome="hg19" ){
  mySession <- try(browserSession("UCSC"), silent=TRUE)
  # In case of failure, try another mirror
  if(inherits(mySession, "try-error"))
    mySession <- browserSession("UCSC",
                                url="http://genome-euro.ucsc.edu/cgi-bin/")
  genome(mySession) <- genome
  obj <- ucscTableQuery(mySession, table="gap")
  tbl <- getTable(obj)
  tbl <- tbl[tbl$type=="telomere", c("chrom", "chromStart", "chromEnd")]
  colnames(tbl)[2:3] <- c("telomereStart", "telomereEnd")
  x <- .chrAsNum(tbl)
  x$chrom <- paste0("chr", x$chrom)
  x$chrom[grepl("chr23", x$chrom)] <- "chrX"
  x$chrom[grepl("chr24", x$chrom)] <- "chrY"
  return(GRanges(x))
}
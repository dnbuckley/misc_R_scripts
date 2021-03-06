# read in bsrate from Andrew's bsrate program
# in methPipe, summary = T will only read the first
# 3 lines, otherwise it will read the whole table
readBSrate <- function(file, summary = F) {
  if (summary) {
    s <- read.table(file, nrows = 3, sep = "\n", stringsAsFactors = F)
    s <- strsplit(s$V1, " = ")
    pos <- unlist(strsplit(s[[2]][2], "\t"))
    neg <- unlist(strsplit(s[[3]][2], "\t"))
    df <- data.frame(ovarall.conv.rate = as.numeric(s[[1]][2]),
                     pos.conv.rate = as.numeric(pos[1]),
                     neg.conv.rate = as.numeric(neg[1]),
                     pos.conv.count = as.numeric(pos[2]),
                     neg.conv.count = as.numeric(neg[2]),
                     stringsAsFactors = F)
  } else {
    df <- read.table(file, skip = 3, header = T)
  }
  return(df)
}
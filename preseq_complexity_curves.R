library(ggplot2)
library(reshape2)

# plots complexity curve from preseq output
preseqComplexityCurve <- function(files, sampleNames = NULL){
  crvs <- lapply(files, read.table, header = T)
  domains <- lapply(crvs, function(x) x$total_reads)
  domain <- max(unlist(domains))
  step <- unique(sapply(domains, function(x){(max(x)-min(x))/(length(x)-1)}))
  stopifnot(length(step) == 1)
  mat <- matrix(nrow = domain/step + 1, ncol = length(crvs) + 1)
  if (is.null (sampleNames)){
    colnames(mat) <- c("total_reads", gsub("(\\.txt|.*\\/)", "", files))
  } else {
    colnames(mat) <- c("total_reads", sampleNames)
  }
  mat[, 1] <- seq(0, domain, step)
  for (i in 1:length(crvs)){
    crv <- crvs[[i]]
    max <- max(crv$total_reads)
    v <- c(crv$distinct_reads, rep(NA, (domain-max)/step))
    mat[, i+1] <- v
  }
  df <- melt(as.data.frame(mat), id.vars = "total_reads")
  ggplot(df, aes(x = total_reads, y = value, color = variable)) +
    geom_line(size = 1) +
    ylab("distinct_reads") +
    theme_bw() +
    geom_abline(slope = 1, 
                intercept = 0,
                color = "dodgerblue", 
                linetype = "dashed",
                size = 1)
}



library(ggplot2)
library(reshape2)
library(parallel)
library(pspline)
library(gridExtra)

# series of functions designed to take data generated
# by preseq lc_extrap and analyze it

# I think this aproximates F'(x)
.getInstSlope <- function(x, y){
  pspl <- smooth.Pspline(x, y)
  f1 <- predict(pspl, x, nderiv=1)
  return(f1)
}

.getUniqueAlignments <- function(multiqcFolder, sampleNames){
  qcFiles <- list.files(multiqcFolder, recursive = T, full.names = T)
  bisFile <- qcFiles[grep("bismark_alignment.txt", qcFiles)]
  bisRep <- read.table(bisFile, header = T)
  if (nrow(bisRep) != length(sampleNames)) {
    return(NULL)
  } else {
    idx <- vector()
    for (i in 1:length(sampleNames)){
      grp <- grep(sampleNames[i], bisRep$Sample)
      if (length(grp) == 0) {
        return(NULL)
      } else {
        idx[i] <- grp
      }
    }
    if (length(unique(idx)) != length(idx) || length(idx) != length(sampleNames)) {
      return(NULL)
    } else {
      return(bisRep$aligned_reads[idx])
    }
  }
}

.getFittedPlot <- function(df, nReads, slopeThreshold = 0.15,
                           dropoff = 0.5, name = "") {
  df <- df[df$f1 >= slopeThreshold, ]
  col <- ifelse(df$f1 >= dropoff, "blue", "red")
  if (!is.na(nReads)) {
    col[df$TOTAL_READS <= nReads] <- "green"
    str <- paste0("Current # Reads = ", nReads,
                  "\nRecommended # Reads = ", max(df$TOTAL_READS[df$f1 >= dropoff]),
                  "\nIncrease of ", round(((max(df$TOTAL_READS[df$f1 >= dropoff])-nReads)/nReads)*100), "%")
  } else {
    str <- paste0("Current # Reads = ", nReads,
                  "\nRecommended # Reads = ", max(df$TOTAL_READS[df$f1 >= dropoff]))
  } 
  gp1 <- ggplot(df, aes(x = TOTAL_READS, y = EXPECTED_DISTINCT)) +
    geom_line(color = col, size = 1) +
    theme_bw() +
    ylim(0, max(df$TOTAL_READS)) +
    geom_abline(intercept = 0, 
                slope = 1, 
                linetype = "dashed", 
                color = "firebrick", 
                size = 0.75) +
    labs(title = name) +
    theme(axis.text = element_blank())
  gp2 <- ggplot(df, aes(x = TOTAL_READS, y = f1)) +
    geom_line(color = col, size = 1) +
    theme_bw() +
    ylab("d/dx(EXPECTED_DISTINCT)") +
    labs(title = "f'(x)") +
    theme(axis.text = element_blank())
  grb <- arrangeGrob(gp1, gp2, ncol = 2, bottom = str)
  return(grb)
}

extrapAnalyze <- function(files, 
                          sampleNames = NULL, 
                          multiqcFolder = NULL,
                          cores = 10){
  ext <- mclapply(files, read.table, header = T, mc.cores = cores)
  if (is.null(sampleNames)){
    sampleNames <- files
  }
  names(ext) <- sampleNames
  
  ext <- mapply(function(x, n){
    x$origin <- n
    return(x)
  }, ext, sampleNames, SIMPLIFY = F)
  
  # get actual # of reads if you can
  if (!is.null(multiqcFolder)) {
    numReads <- .getUniqueAlignments(multiqcFolder, sampleNames)
    if (is.null(numReads)) {
      message("WARNING: Unable to map sampleNames to MultiQC Data")
      numReads <- rep(NA, length(sampleNames))
    }
  } else {
    numReads <- rep(NA, length(sampleNames))
  }
  
  # calculate slope at each point
  ext <- lapply(ext, function(x){
    x$f1 <- .getInstSlope(x = x$TOTAL_READS, y = x$EXPECTED_DISTINCT)
    return(x)
  })
  grbs <- mcmapply(.getFittedPlot, 
                   df = ext, 
                   nReads = numReads, 
                   name = sampleNames,
                   mc.cores = cores,
                   SIMPLIFY = F)
  return(grbs)
}

if (0){
  # Demonstration of pspline
  library(pspline)
  library(ggplot2)
  df <- data.frame(x = seq(-5, 5, 0.01))
  df$y <- df$x^2
  pspl <- smooth.Pspline(df$x, df$y)
  df$f1 <- predict(pspl, df$x, nderiv = 1)
  df$f2 <- predict(pspl, df$x, nderiv = 2)
  diff <- as.numeric(predict(pspl, df$x, nderiv = 0)) - df$y
  summary(diff)
  hist(diff)
  df <- melt(df, id.vars = "x")
  ggplot(df, aes(x = x, y = value, color = variable)) +
    geom_line(size = 1)
}
  
  
  
  
  
  
  
  
  



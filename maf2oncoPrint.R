library(maftools)

# I don't like oncoprint, so this should... hopefully...
# Convert a maf to a matrix that can be used by oncoPlot.
# The VCF files should be annotated by VEP before converting
# to maf, if you did it right maf@data should have 133 columns

# requires unique names in the "Tumor_Sample_Barcode" field

# debug
if (0){
  setwd("~/Desktop/salhia_lab/poetic_paper/vcf_annotated/mafAnalysis/")
  z <- read.delim("maf_files/P01-015-T-1-FHC.maf", sep = '\t', skip = 1, header = T)
  z <- z[1:50000, ]
  dat <- z
}

.getWorstMutation <- function(muts, ref){
  if (length(unique(muts)) == 1){
    return(unique(muts))
  } else {
    muts <- unique(muts)
  }
  ref <- ref[ref$class %in% muts, ]
  # if there are 2 or more mutations more deleterious than 
  # a missense mutation return "multi-hit"
  if (nrow(ref[ref$badness >= 9, ]) >= 2){
    return("multi-hit")
  } else{
    return(ref$class[which(ref$badness %in% max(ref$badness))])
  }
}

maf2oncoPrint <- function(maf){
  if (class(maf) == "MAF"){
    dat <- maf@data
  } else if (class(maf) == "data.frame"){
    dat <- maf
  } else {
    message("ERROR: invalid class")
    stop()
  }
  # ranking badness of variant classifications, open to interpretation
  ref <- openxlsx::read.xlsx("~/misc_R_scripts/script_files/vep_classification_ranking.xlsx")
  ref <- ref[order(ref$class), ]
  dat <- dat[!is.na(dat$Tumor_Sample_Barcode), ]
  dat <- dat[order(dat$Tumor_Sample_Barcode), ]
  mat <- matrix(nrow = length(unique(dat$Hugo_Symbol)),
                ncol = length(unique(dat$Tumor_Sample_Barcode)))
  colnames(mat) <- unique(dat$Tumor_Sample_Barcode)
  rownames(mat) <- unique(dat$Hugo_Symbol)
  dat <- split(dat, dat$Tumor_Sample_Barcode)
  stopifnot(identical(names(dat), colnames(mat)))
  mat <- as.matrix(mat[order(rownames(mat)), ])
  # lmao good luck figuring out how this works
  for (i in 1:ncol(mat)){
    df <- dat[[i]]
    df$Hugo_Symbol <- as.character(df$Hugo_Symbol)
    df$Variant_Classification <- as.character(df$Variant_Classification)
    df <- df[order(df$Hugo_Symbol), ]
    mat[, i] <- rownames(mat) %in% df$Hugo_Symbol
    mutLst <- split(df$Variant_Classification, df$Hugo_Symbol)
    con <- lapply(mutLst, .getWorstMutation, ref)
    stopifnot(identical(rownames(mat)[mat[, i]], names(con)))
    mat[mat[, i], i] <- unlist(con)
  }
  mat[is.na(mat)] <- ""
  return(mat)
}












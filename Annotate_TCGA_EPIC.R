library(GenomicRanges)
# source("~/misc_R_scripts/annotateInfinium.R")

# provide filename of TCGA file, will annotate and either return or send
annotateTCGA_EPIC <- function(filename, return = F, save = F, 
                              genome = NULL, refAnnotations = NULL, 
                              annoDir = "~/misc_R_scripts/script_files/", 
                              saveDir = "", skipLines = 0){
  # idiot proofing
  if (!return & !save){
    message("select return or save")
    return(NULL)
  }
  if (is.null(genome) & is.null(refAnnotations)){
    message("provide either a genome or a refAnnotations granges object (Wanding's RDS objects)")
    return(NULL)
  }
  if (!is.null(genome) & !is.null(refAnnotations)){
    warning("you provided both a refAnnotation object and genome, using refAnnotation obj for annotation")
  }
  # load annotations if not provided
  if (is.null(refAnnotations) & !is.null(genome)){
    # attempt to guess the chip type
    chip <- ifelse(grepl("HumanMethylation450", filename), "hm450",
                   ifelse(grepl("HumanMethylation27", filename), "hm27", "EPIC"))
    refAnnotations <- readRDS(paste0(annoDir, chip, ".", genome, ".manifest.rds"))
  }
  df <- read.table(filename, sep = "\t", header = T, skip = skipLines)
  df <- df[, !grepl("(Chromosome|Start|End|Genomic_Coordinate)", names(df))]
  names(df)[1] <- "CpG_ID"
  df$CpG_ID <- as.character(df$CpG_ID)
  refAnnotations <- refAnnotations[order(names(refAnnotations))]
  df <- df[order(df$CpG_ID), ]
  commonCpGs <- intersect(names(refAnnotations), df$CpG_ID)
  refAnnotations <- refAnnotations[names(refAnnotations) %in% commonCpGs]
  df <- df[df$CpG_ID %in% commonCpGs, ]
  stopifnot(identical(df$CpG_ID, names(refAnnotations)))
  df <- cbind(as.data.frame(granges(refAnnotations)),
              df,
              mcols(refAnnotations))
  gr <- GRanges(df)
  gr <- gr[order(gr)]
  if (save){
    if (!is.null(genome)) outfile <- gsub("_hg.*", paste0(".annotated.", genome, ".rds"), filename)
    else outfile <- gsub("_hg.*", ".annotated.rds", filename)
    outfile <- gsub(".txt$", ".rds", outfile)
    outfile <- gsub(".*\\/", "", outfile)
    outfile <- paste0(saveDir, outfile)
    message("saving to ", outfile)
    saveRDS(gr, outfile)
    rm(gr, df, refAnnotations)
  }
  if (return){
    return(gr)
  }
  return(NULL)
}

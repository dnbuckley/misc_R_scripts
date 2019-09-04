library(GenomicRanges)
source("~/misc_R_scripts/annotateInfinium.R")

# provide filename of TCGA file, will annotate and either return or send
annotateTCGA_EPIC <- function(filename, return = F, save = F, 
                              genome = NULL, refAnnotations = NULL, 
                              annoDir = "~/misc_R_scripts/script_files/"){
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
  df <- read.table(filename, sep = "\t", header = T)
  df <- df[, !grepl("(Chromosome|Start|End)", names(df))]
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
    saveRDS(gr, gsub("_hg.*", ".annotated.", genome, ".rds"))
  }
  if (return){
    return(gr)
  }
}

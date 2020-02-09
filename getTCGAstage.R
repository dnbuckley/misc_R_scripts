# get TCGA stage from file name
name <- names[1]
getTCGAstage <- function(name, clinDat = NULL, simplifyStage = T){
  if (is.null(clinDat)){
    clinDat <- read.table("~/misc_R_scripts/script_files/TCGA_450k_clinical_data.tsv", 
                          sep = "\t", header = T)
  }
  name <- gsub(".*lvl-[0-9]\\.", "", name)
  name <- gsub("-[0-9]{2}[A-Z]-.*", "", name)
  row <- clinDat[clinDat$submitter_id %in% name, ]
  row <- unique(row)
  stopifnot(nrow(row) == 1)
  stage <- row$primary_diagnosis
}
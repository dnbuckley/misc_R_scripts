library(openxlsx)

getTCGAtype <- function(barcode, xl = NULL, abriv = F){
  if (is.null(xl)){
    xl <- read.xlsx("~/misc_R_scripts/script_files/tcga_sample_type_codes.xlsx")
  }
  xl$definition <- gsub(" ", "_", xl$definition)
  barcode <- gsub("-[0-9]{2}D-.*", "", barcode)
  barcode <- gsub(".*-", "", barcode)
  num <- as.numeric(gsub("[A-Z]", "", barcode))
  row <- xl[xl$code %in% num, ]
  if (nrow(row) != 1) {
    return("Undefined")
  } else if (abriv) {
    return(row$short_letter)
  } else {
    return(row$definition)
  }
}

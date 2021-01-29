library(openxlsx)

xlsx2txt <- function(file) {
  xl <- read.xlsx(file)
  write.table(xl, row.names = F, col.names = T, quote = F, sep = "\t",
              file = gsub(".xlsx$", ".txt", file))
}

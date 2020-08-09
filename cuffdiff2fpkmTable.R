cuffdiff2fpkmTable <- function(df) {
  samples <- df[, grepl("^sample", names(df))]
  fpkms <- df[, grepl("^value", names(df))]
  stopifnot(identical(dim(samples), dim(fpkms)))
  IDs <- as.character(samples[1, ])
  names(fpkms) <- paste0("value_", IDs)
  return(fpkms)
}
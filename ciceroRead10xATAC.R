library(cicero)

# create ATAC CDS from 10x output
ciceroRead10xATAC <- function(countsDir, doPreProc = T) {
  # find files
  message("reading ATAC files in ", countsDir, "...")
  countsDir <- paste0(countsDir, "/")
  mtx.file <- list.files(countsDir, "matrix.mtx", recursive = T)
  mtx.file <- mtx.file[grepl("filtered", mtx.file) & !grepl("analysis", mtx.file)]
  bc.file <- list.files(countsDir, "barcodes.tsv", recursive = T)
  bc.file <- bc.file[grepl("filtered", bc.file) & !grepl("analysis", bc.file)]
  peaks.file <- list.files(countsDir, "peaks.bed", recursive = T)
  feat.file <- list.files(countsDir, "features.tsv", recursive = T)
  feat.file <- feat.file[grepl("filtered", feat.file) & !grepl("analysis", feat.file)]
  stopifnot(length(c(mtx.file, bc.file, peaks.file, feat.file)) == 4)
  
  features <- read.delim(paste0(countsDir, feat.file), header = F)
  
  indata <- Matrix::readMM(paste0(countsDir, mtx.file)) 
  indata <- indata[features$V3 == "Peaks", ]
  indata@x[indata@x > 0] <- 1
  
  cellinfo <- read.table(paste0(countsDir, bc.file))
  row.names(cellinfo) <- cellinfo$V1
  names(cellinfo) <- "cells"
  
  peakinfo <- read.table(paste0(countsDir, peaks.file))
  names(peakinfo) <- c("chr", "bp1", "bp2")
  peakinfo$site_name <- paste(peakinfo$chr, peakinfo$bp1, peakinfo$bp2, sep="_")
  row.names(peakinfo) <- peakinfo$site_name
  
  row.names(indata) <- row.names(peakinfo)
  colnames(indata) <- row.names(cellinfo)
  
  # make CDS
  input_cds <-  suppressWarnings(new_cell_data_set(indata,
                                                   cell_metadata = cellinfo,
                                                   gene_metadata = peakinfo))
  
  input_cds <- monocle3::detect_genes(input_cds)
  #Ensure there are no peaks included with zero reads
  input_cds <- input_cds[Matrix::rowSums(exprs(input_cds)) != 0,]
  
  if (doPreProc) {
    message("Running preprocessing...")
    input_cds <- detect_genes(input_cds)
    input_cds <- estimate_size_factors(input_cds)
    input_cds <- preprocess_cds(input_cds, method = "LSI", verbose = T)
    input_cds <- reduce_dimension(input_cds, reduction_method = 'UMAP',
                                  preprocess_method = "LSI", verbose = T)
  }
  
  return(input_cds)
}


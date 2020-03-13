library(Rsamtools)

setwd("~/Desktop/salhia_lab/misc_R_scripts/")
bam <- "../reprod_ovarian/umi_dup_rates/all_bams/run1/0-6-2_S18_R1_001_bismark_bt2_pe.hg38.sort.bam"
# this file is not what chris thinks it is,
# no reads lie within these regions, all are too long
gr <- readRDS("../stitch_ovarian_together/ovarian_formatted_regions_2.rds")
gr <- gr[gr$plex == 2]

filterBamInRange <- function(bam, gr){
  if (!file.exists(paste0(bam, ".bai"))){
    message("creating index for ", bam)
    indexBam(bam)
  }
  bai <- paste0(bam, ".bai")
  # params <- ScanBamParam(what = c("qname", "rname", "pos", "qwidth"))
  bam <- scanBam(bam, bai)
  readNames <- lapply(bam, function(x){x$qname})
  starts <- lapply(bam, function(x){x$pos})
  seqNames <- lapply(bam, function(x){x$rname})
  widths <- lapply(bam, function(x){x$qwidth})
  df <- data.frame(chromosome = lapply(seqNames, function(x){x}),
                   start = lapply(starts, function(x){x}),
                   width = lapply(widths, function(x){x}),
                   name = lapply(readNames, function(x){x}))
  names(df) <- c("chromosome", "start", "width", "name")
  df$end <- df$start+df$width
  reads <- GRanges(df)
  ovlp <- findOverlaps(gr, reads, type = "within")
  reads <- reads[to(findOverlaps(gr, reads))]
}
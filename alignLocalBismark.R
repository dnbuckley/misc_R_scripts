library(GenomicRanges)
library(BSgenome)
library(seqinr)
library(Rsamtools)
library(chromstaR)
source("~/misc_R_scripts/getBSgenomeObj.R")

# uses bismark to find where a primer is in a DMR region
# because chris is too lazy to do it as he makes the primers

# using bismark for this purpose: https://i.imgur.com/WQ6K6iE.jpg

# gr <- GRanges("chr8:37683362-37683591")
# query.seq <- "TTTTTAGGCGGGTTTGAGGTAAGGGAG"
genome <- "hg19"
buffer <- 100

.makeFQ <- function(seq, tmpFile) {
  # dummy read name to trick bismark
  fq <- c("@M05368:125:000000000-D8BW6:1:1101:13163:1837", 
          seq, "+", paste0(rep("I", nchar(seq)), collapse = ""))
  write.table(fq, tmpFile, sep = "\n", col.names = F, row.names = F, quote = F)
}

alignLocalBismark <- function(gr, query.seq, genome, buffer = NULL,
                              bispath = "/Users/dbuckley/Desktop/bin/Bismark/",
                              bt2path = "/Users/dbuckley/Desktop/bin/") {
  bsg <- getBSgenomeObj(genome = genome)
  bis <- paste0(bispath, "bismark")
  bis.build <- paste0(bispath, "bismark_genome_preparation")
  temp.fq <- "temp.fastq"
  tempdir <- paste0("tempbisdir_", query.seq)
  
  if (!is.null(buffer)) {
    start(gr) <- start(gr)-buffer
    end(gr) <- end(gr)+buffer
  }
  ref.seq <- getSeq(bsg, gr)
  dir.create(tempdir, showWarnings = F)
  setwd(tempdir)
  dir.create("genome", showWarnings = F)
  # export(ref.seq, "genome/temp.fa", format = "fasta")
  write.fasta(ref.seq, "chrZ", "genome/temp.fa", nbchar = 80)
  cmd <- paste0(bis.build, " --path_to_aligner ", bt2path, " genome/")
  system(cmd)
  .makeFQ(query.seq, temp.fq)
  cmd <- paste0(bis, " --path_to_bowtie2 ", bt2path, " --se ", temp.fq, " genome/")
  system(cmd, show.output.on.console = "invisible")
  bam <- list.files(".", "bam$")
  mr <- as.numeric(system(paste0("samtools view ", bam, " | wc -l"), intern = T))
  if (mr == 1) {
    bam <- readBamFileAsGRanges(bam)
    primer.loc <- GRanges(seqnames = seqnames(gr), 
                          ranges = IRanges(start = start(gr)+start(bam)-1,
                                           end = start(gr)+end(bam)-1),
                          strand = strand(bam))
    primer.loc$mapq <- primer.loc$mapq
  } else {
    primer.loc <- NA
  }
  setwd("../")
  system(paste0("rm -r ", tempdir))
  return(primer.loc)
}







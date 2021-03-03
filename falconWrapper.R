library(Canopy)
library(falcon)

# multithreaded wrapper functions for processing canopy ASCN;
# preproc pipeline: bt2 -> picard MarkDuplicates -> GATK HaplotypeCaller -> GATK ASEReadCounter

# called primary & relapse for simplicity, can slso be mets
falconWrapper <- function(primaryFile, relapseFile, normalFile, 
                          cores = 10, outDir = "falcon_output/", doQC = T) {
  outDir <- paste0(outDir, "/")
  dir.create(outDir)
  primary <- read.delim(primaryFile)
  relapse <- read.delim(relapseFile)
  normal <- read.delim(normalFile)
  primary <- .preprocessRC(primary, normal)
  relapse <- .preprocessRC(relapse, normal)
  # calculate depth ratio (total read counts of tumor versus normal)
  rdep.relapse=sum(relapse$Tumor_ReadCount_Total)/sum(relapse$Normal_ReadCount_Total)
  rdep.primary=sum(primary$Tumor_ReadCount_Total)/sum(primary$Normal_ReadCount_Total)
  primary <- split(primary, primary$Chromosome)
  relapse <- split(relapse, relapse$Chromosome)
  # for GRCh and hg
  possibleChrs <- c(paste0("chr", 1:22), 1:22)
  primary <- primary[names(primary) %in% possibleChrs]
  relapse <- relapse[names(relapse) %in% possibleChrs]
  stopifnot(identical(names(primary), names(relapse)))
  mcmapply(.runFalcon, primary, relapse, names(primary), 
           MoreArgs = list(rdep.primary = rdep.primary,
                           rdep.relapse = rdep.relapse,
                           outDir = outDir, doQC = doQC),
           mc.cores = cores)
}

.preprocessRC <- function(t, n) {
  rownames(t) <- paste0(t$contig, ":", t$position)
  rownames(n) <- paste0(n$contig, ":", n$position)
  common <- intersect(rownames(t), rownames(n))
  t <- t[rownames(t) %in% common, ]
  n <- n[rownames(n) %in% common, ]
  stopifnot(identical(rownames(t), rownames(n)))
  stopifnot(identical(t$refAllele, n$refAllele))
  ret <- data.frame(Chromosome = t$contig,
                    Start_position = t$position,
                    End_position = t$position,
                    Reference_Allele = t$refAllele,
                    TumorSeq_Allele1 = t$refAllele,
                    TumorSeq_Allele2 = t$altAllele,
                    Tumor_ReadCount_Total = t$totalCount, 
                    Tumor_ReadCount_Ref = t$refCount,
                    Tumor_ReadCount_Alt = t$altCount,
                    Match_Norm_Seq_Allele1 = n$refAllele,
                    Match_Norm_Seq_Allele2 = n$altAllele,
                    Normal_ReadCount_Total = n$totalCount,
                    Normal_ReadCount_Ref = n$refCount,
                    Normal_ReadCount_Alt = n$altCount)
  return(ret)
}

.runFalcon <- function(primary.chr, relapse.chr, chr, rdep.primary, rdep.relapse, outDir, doQC) {
  source('~/misc_R_scripts/falcon.qc.R')
  relapse.chr=relapse.chr[relapse.chr[,'Match_Norm_Seq_Allele1']!=' ',]
  relapse.chr=relapse.chr[relapse.chr[,'Match_Norm_Seq_Allele2']!=' ',]
  relapse.chr=relapse.chr[relapse.chr[,'Reference_Allele']!=' ',]
  relapse.chr=relapse.chr[relapse.chr[,'TumorSeq_Allele1']!=' ',]
  relapse.chr=relapse.chr[relapse.chr[,'TumorSeq_Allele2']!=' ',]
  # get germline heterozygous loci (normal allele1 != normal allele2)
  relapse.chr=relapse.chr[(as.matrix(relapse.chr[,'Match_Norm_Seq_Allele1'])!=as.matrix(relapse.chr[,'Match_Norm_Seq_Allele2'])),]
  ############################################################
  # QC procedures to remove false neg and false pos variants.
  # The thresholds can be adjusted.
  ############################################################
  # remove indels (this can be relaxed but we think indels are harder to call than SNPs)
  indel.filter1=nchar(as.matrix(relapse.chr[,'Reference_Allele']))<=1
  indel.filter2=nchar(as.matrix(relapse.chr[,'Match_Norm_Seq_Allele1']))<=1
  indel.filter3=nchar(as.matrix(relapse.chr[,'Match_Norm_Seq_Allele2']))<=1
  indel.filter4=nchar(as.matrix(relapse.chr[,'TumorSeq_Allele1']))<=1
  indel.filter5=nchar(as.matrix(relapse.chr[,'TumorSeq_Allele2']))<=1
  relapse.chr=relapse.chr[indel.filter1 & indel.filter2 & indel.filter3 & indel.filter4 & indel.filter5,]
  # total number of reads greater than 30 in both tumor and normal
  depth.filter1=(relapse.chr[,"Normal_ReadCount_Ref"]+relapse.chr[,"Normal_ReadCount_Alt"])>=30
  depth.filter2=(relapse.chr[,"Tumor_ReadCount_Ref"]+relapse.chr[,"Tumor_ReadCount_Alt"])>=30
  relapse.chr=relapse.chr[depth.filter1 & depth.filter2,]
  #########################
  # Generate FALCON input.
  #########################
  # Data frame with four columns: tumor ref, tumor alt, normal ref, normal alt.
  readMatrix.relapse=as.data.frame(relapse.chr[,c('Tumor_ReadCount_Ref',
                                                  'Tumor_ReadCount_Alt',
                                                  'Normal_ReadCount_Ref',
                                                  'Normal_ReadCount_Alt')])
  colnames(readMatrix.relapse)=c('AT','BT','AN','BN')
  ###############################
  # Run FALCON and view results.
  ###############################
  tauhat.relapse=getChangepoints(readMatrix.relapse)
  cn.relapse = getASCN(readMatrix.relapse, tauhat=tauhat.relapse, rdep = rdep.relapse, threshold = 0.3)
  # Chromosomal view of segmentation results.
  pdf(file=paste(outDir, 'falcon.relapse.',chr,'.pdf',sep=''),width=5,height=8)
  view(cn.relapse,pos=relapse.chr[,'Start_position'], rdep = rdep.relapse)
  dev.off()
  # save image file.
  save.image(file=paste(outDir, 'falcon_relapse_',chr,'.rda',sep=''))
  ########################################
  # Further curate FALCON's segmentation.
  ########################################
  if(length(tauhat.relapse)>0){
    length.thres=10^6  # Threshold for length of segments, in base pair.
    delta.cn.thres=0.3  # Threshold of absolute copy number difference between consecutive segments.
    # source('~/misc_R_scripts/falcon.qc.R') # Can be downloaded from
    # https://github.com/yuchaojiang/Canopy/blob/master/instruction/falcon.qc.R
    falcon.qc.list = falcon.qc(readMatrix = readMatrix.relapse,
                               tauhat = tauhat.relapse,
                               cn = cn.relapse,
                               st_bp = relapse.chr[,"Start_position"],
                               end_bp = relapse.chr[,"End_position"],
                               rdep = rdep.relapse,
                               length.thres = length.thres,
                               delta.cn.thres = delta.cn.thres)
    
    tauhat.relapse=falcon.qc.list$tauhat
    cn.relapse=falcon.qc.list$cn
  }
  # Chromosomal view of QC'ed segmentation results.
  if (doQC) {
    pdf(file=paste(outDir, 'falcon.relapse.qc.',chr,'.pdf',sep=''),width=5,height=8)
    view(cn.relapse,pos=relapse.chr[,'Start_position'], rdep = rdep.relapse)
    dev.off()
  }
  #################################################
  # Generate Canopy's input with s.d. measurement.
  #################################################
  # This is to generate table output including genomic locations for 
  # segment boudaries.
  # For Canopy's input, we use Bootstrap-based method to estimate the
  # standard deviations for the allele-specific copy numbers.
  source('~/misc_R_scripts/falcon.output.R')
  falcon.output=falcon.output(readMatrix = readMatrix.relapse,
                              tauhat = tauhat.relapse,
                              cn = cn.relapse,
                              st_bp = relapse.chr[,"Start_position"],
                              end_bp = relapse.chr[,"End_position"],
                              nboot = 5000)
  falcon.output = cbind(chr=rep(chr,nrow(falcon.output)), falcon.output)
  write.table(falcon.output, file=paste(outDir, 'falcon.relapse.output.',chr,'.txt',sep=''), 
              col.names =T, row.names = F, sep='\t', quote = F)
  
  ###########################################
  ###########################################
  #
  #        Primary tumor
  #
  ###########################################
  ###########################################  
  
  ###########################################
  # Focus on germline heterozygous variants.
  ###########################################
  
  # remove variants with missing genotype
  primary.chr=primary.chr[primary.chr[,'Match_Norm_Seq_Allele1']!=' ',]
  primary.chr=primary.chr[primary.chr[,'Match_Norm_Seq_Allele2']!=' ',]
  primary.chr=primary.chr[primary.chr[,'Reference_Allele']!=' ',]
  primary.chr=primary.chr[primary.chr[,'TumorSeq_Allele1']!=' ',]
  primary.chr=primary.chr[primary.chr[,'TumorSeq_Allele2']!=' ',]
  # get germline heterozygous loci (normal allele1 != normal allele2)
  primary.chr=primary.chr[(as.matrix(primary.chr[,'Match_Norm_Seq_Allele1'])!=as.matrix(primary.chr[,'Match_Norm_Seq_Allele2'])),]
  ############################################################
  # QC procedures to remove false neg and false pos variants.
  # The thresholds can be adjusted.
  ############################################################
  # remove indels (this can be relaxed but we think indels are harder to call than SNPs)
  indel.filter1=nchar(as.matrix(primary.chr[,'Reference_Allele']))<=1
  indel.filter2=nchar(as.matrix(primary.chr[,'Match_Norm_Seq_Allele1']))<=1
  indel.filter3=nchar(as.matrix(primary.chr[,'Match_Norm_Seq_Allele2']))<=1
  indel.filter4=nchar(as.matrix(primary.chr[,'TumorSeq_Allele1']))<=1
  indel.filter5=nchar(as.matrix(primary.chr[,'TumorSeq_Allele2']))<=1
  primary.chr=primary.chr[indel.filter1 & indel.filter2 & indel.filter3 & indel.filter4 & indel.filter5,]
  # total number of reads greater than 30 in both tumor and normal
  depth.filter1=(primary.chr[,"Normal_ReadCount_Ref"]+primary.chr[,"Normal_ReadCount_Alt"])>=30
  depth.filter2=(primary.chr[,"Tumor_ReadCount_Ref"]+primary.chr[,"Tumor_ReadCount_Alt"])>=30
  primary.chr=primary.chr[depth.filter1 & depth.filter2,]
  #########################
  # Generate FALCON input.
  #########################
  # Data frame with four columns: tumor ref, tumor alt, normal ref, normal alt.
  readMatrix.primary=as.data.frame(primary.chr[,c('Tumor_ReadCount_Ref',
                                                  'Tumor_ReadCount_Alt',
                                                  'Normal_ReadCount_Ref',
                                                  'Normal_ReadCount_Alt')])
  colnames(readMatrix.primary)=c('AT','BT','AN','BN')
  dim(readMatrix.primary); dim(primary.chr)
  ###############################
  # Run FALCON and view results.
  ###############################
  tauhat.primary=getChangepoints(readMatrix.primary)
  cn.primary = getASCN(readMatrix.primary, tauhat=tauhat.primary, rdep = rdep.primary, threshold = 0.3)
  
  # Chromosomal view of segmentation results.
  pdf(file=paste(outDir, 'falcon.primary.',chr,'.pdf',sep=''),width=5,height=8)
  view(cn.primary,pos=primary.chr[,'Start_position'], rdep = rdep.primary)
  dev.off()
  # save image file.
  save.image(file=paste(outDir, 'falcon_primary_',chr,'.rda',sep=''))
  ########################################
  # Further curate FALCON's segmentation.
  ########################################
  if(length(tauhat.primary)>0){
    length.thres=10^6  # Threshold for length of segments, in base pair.
    delta.cn.thres=0.3  # Threshold of absolute copy number difference between consecutive segments.
    # source('falcon_demo/falcon.qc.R') # Can be downloaded from
    # https://github.com/yuchaojiang/Canopy/blob/master/instruction/falcon.qc.R
    falcon.qc.list = falcon.qc(readMatrix = readMatrix.primary,
                               tauhat = tauhat.primary,
                               cn = cn.primary,
                               st_bp = primary.chr[,"Start_position"],
                               end_bp = primary.chr[,"End_position"],
                               rdep = rdep.primary,
                               length.thres = length.thres,
                               delta.cn.thres = delta.cn.thres)
    
    tauhat.primary=falcon.qc.list$tauhat
    cn.primary=falcon.qc.list$cn
  }
  # Chromosomal view of QC'ed segmentation results.
  if (doQC) {
    pdf(file=paste(outDir, 'falcon.primary.qc.',chr,'.pdf',sep=''),width=5,height=8)
    view(cn.primary,pos=primary.chr[,'Start_position'], rdep = rdep.primary)
    dev.off()
  }
  #################################################
  # Generate Canopy's input with s.d. measurement.
  #################################################
  # This is to generate table output including genomic locations for 
  # segment boudaries.
  # For Canopy's input, we use Bootstrap-based method to estimate the
  # standard deviations for the allele-specific copy numbers.
  # source('falcon_demo/falcon.output.R') # Can be downloaded from
  # https://github.com/yuchaojiang/Canopy/blob/master/instruction/falcon.output.R
  source('~/misc_R_scripts/falcon.output.R')
  falcon.output=falcon.output(readMatrix = readMatrix.primary,
                              tauhat = tauhat.primary,
                              cn = cn.primary,
                              st_bp = primary.chr[,"Start_position"],
                              end_bp = primary.chr[,"End_position"],
                              nboot = 5000)
  falcon.output = cbind(chr=rep(chr,nrow(falcon.output)), falcon.output)
  write.table(falcon.output, file=paste(outDir ,'falcon.primary.output.',chr,'.txt',sep=''), 
              col.names = T, row.names = F, sep='\t', quote = F)  
}









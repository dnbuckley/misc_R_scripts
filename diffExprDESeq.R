library(DESeq2)

# get expression data from desseq
diffExprDESeq <- function(grpA, grpB, grpAlab = "grpA", grpBlab = "grpB",
                          lfcShrink = T, 
                          saveObjName = NULL) {
  allFiles <- data.frame(path = c(grpA, grpB))
  allFiles$sample <- gsub("\\/.*", "", allFiles$path)
  # allFiles$sample <- gsub("_.*", "", allFiles$path)
  allFiles$group <- ifelse(allFiles$path %in% grpA, grpAlab, 
                           ifelse(allFiles$path %in% grpB, grpBlab, NA))
  allFiles <- allFiles[!is.na(allFiles$group), ]
  message("Testing ", nrow(allFiles), " samples.")
  sampleTable <- data.frame(sampleName = allFiles$sample,
                            fileName = allFiles$path,
                            condition = allFiles$group)
  sampleTable$condition <- factor(sampleTable$condition)
  ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                         directory = directory,
                                         design = ~ condition)
  dds <- DESeq(ddsHTSeq, parallel = T)
  res <- results(dds, contrast = c('condition', grpAlab, grpBlab), parallel = T)
  if (lfcShrink) {
    res <- lfcShrink(dds, contrast = c('condition', grpAlab, grpBlab), res=res, type = 'ashr', parallel = T)
  }
  if (!is.null(saveObjName)) {
    saveRDS(dds, saveObjName)
  }
  res <- as.data.frame(res)
  res$gene <- rownames(res)
  return(res)
}

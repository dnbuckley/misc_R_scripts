library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(ggplot2)
library(chromVAR)
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)
library(BSgenome.Hsapiens.UCSC.hg38)

# scripts here are drawn from the Seurat web page, basic steps:
#### 1) Read in data
#### 2) Filter cells
#### 3) Reduce dimensions
#### 4) Evaluate motifs

readSeuratMultiome <- function(dir, 
                               showQC = T, 
                               alignment_group = NULL,
                               annotations = NULL) {
  
  if (!grepl("outs", dir)) {
    dir <- paste0(dir, "/outs/")
  }
  message("reading from ", dir)
  
  h5.10x <- Read10X_h5(paste0(dir ,"/filtered_feature_bc_matrix.h5"))
  rna_counts <- h5.10x$`Gene Expression`
  atac_counts <- h5.10x$Peaks
  
  # Create Seurat object
  seurat.obj <- CreateSeuratObject(counts = rna_counts)
  seurat.obj[["percent.mt"]] <- PercentageFeatureSet(seurat.obj, pattern = "^MT-")
  
  # Now add in the ATAC-seq data
  # we'll only use peaks in standard chromosomes
  grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
  grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
  atac_counts <- atac_counts[as.vector(grange.use), ]
  if(is.null(annotations)) {
    annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
  }
  seqlevelsStyle(annotations) <- 'UCSC'
  genome(annotations) <- "hg38"
  frag.file <- paste0(dir, "/atac_fragments.tsv.gz")
  chrom_assay <- CreateChromatinAssay(counts = atac_counts, 
                                      sep = c(":", "-"), 
                                      genome = 'hg38',
                                      fragments = frag.file,
                                      min.cells = 10,
                                      annotation = annotations)
  seurat.obj[["ATAC"]] <- chrom_assay
  if (showQC) {
    message("Plotting qc...")
    vp <- VlnPlot(seurat.obj, features = c("nCount_ATAC", "nCount_RNA","percent.mt"), ncol = 3,
                  log = TRUE, pt.size = 0) + 
      NoLegend() +
      labs(title = dir)
    plot(vp)
  }
  if (!is.null(alignment_group)) {
    seurat.obj$alignment_group <- alignment_group
  }
  return(seurat.obj)
}

filterSeurat <- function(seurat.obj, 
                         max.atac = 7e4,
                         min.atac = 1e3,
                         max.rna = 2.5e4,
                         min.rna = 1e3,
                         max.pct.mt = 20,
                         showQC = T) {
  seurat.obj <- subset(
    x = seurat.obj,
    subset = nCount_ATAC < max.atac &
      nCount_ATAC > min.atac &
      nCount_RNA < max.rna &
      nCount_RNA > min.rna &
      percent.mt < max.pct.mt
  )
  if (showQC) {
    message("Plotting qc...")
    vp <- VlnPlot(seurat.obj, features = c("nCount_ATAC", "nCount_RNA","percent.mt"), ncol = 3,
                  log = TRUE, pt.size = 0) + 
      NoLegend()
    plot(vp)
  }
  return(seurat.obj)
}

reduceDimsSeruat <- function(seurat.obj, 
                             plotUMAP = T, 
                             color.by = "seurat_clusters",
                             conserve.memory = F) {
  # RNA analysis
  DefaultAssay(seurat.obj) <- "RNA"
  seurat.obj <- SCTransform(seurat.obj, 
                            verbose = T, 
                            conserve.memory = conserve.memory) %>% 
    RunPCA() %>% 
    RunUMAP(dims = 1:50, 
            reduction.name = 'umap.rna', 
            reduction.key = 'rnaUMAP_')
  # ATAC analysis
  # We exclude the first dimension as this is typically correlated with sequencing depth
  DefaultAssay(seurat.obj) <- "ATAC"
  seurat.obj <- RunTFIDF(seurat.obj)
  seurat.obj <- FindTopFeatures(seurat.obj, min.cutoff = 'q0')
  seurat.obj <- RunSVD(seurat.obj, verbose = T)
  seurat.obj <- RunUMAP(seurat.obj, reduction = 'lsi', 
                        dims = 2:50, 
                        reduction.name = "umap.atac", 
                        reduction.key = "atacUMAP_")
  
  seurat.obj <- FindMultiModalNeighbors(seurat.obj, 
                                        reduction.list = list("pca", "lsi"), 
                                        dims.list = list(1:50, 2:50), 
                                        verbose = T)
  seurat.obj <- RunUMAP(seurat.obj, 
                        nn.name = "weighted.nn", 
                        reduction.name = "wnn.umap", 
                        reduction.key = "wnnUMAP_", 
                        verbose = T)
  seurat.obj <- FindClusters(seurat.obj, 
                             graph.name = "wsnn", 
                             algorithm = 3, 
                             verbose = T)
  if (plotUMAP) {
    p <- plotSeuratUmaps(seurat.obj, color.by = color.by)
  }
  return(seurat.obj)
}

evalMotifsSeurat <- function(seurat.obj) {
  # Scan the DNA sequence of each peak for the presence of each motif, 
  # and create a Motif object
  DefaultAssay(seurat.obj) <- "ATAC"
  pwm_set <- getMatrixSet(x = JASPAR2020, 
                          opts = list(species = 9606, all_versions = FALSE))
  motif.matrix <- CreateMotifMatrix(features = granges(seurat.obj), 
                                    pwm = pwm_set, 
                                    genome = 'hg38', 
                                    use.counts = FALSE)
  motif.object <- CreateMotifObject(data = motif.matrix, 
                                    pwm = pwm_set)
  seurat.obj <- SetAssayData(seurat.obj, assay = 'ATAC', 
                             slot = 'motifs', 
                             new.data = motif.object)
  # Note that this step can take 30-60 minutes 
  seurat.obj <- RunChromVAR(object = seurat.obj,
                            genome = BSgenome.Hsapiens.UCSC.hg38)
  return(seurat.obj)
}

plotSeuratUmaps <- function(seurat.obj, color.by = "seurat_clusters") {
  p1 <- DimPlot(seurat.obj, 
                reduction = "umap.rna", 
                label = TRUE, 
                label.size = 2.5, 
                repel = TRUE, 
                group.by = color.by) + 
    ggtitle("RNA")
  p2 <- DimPlot(seurat.obj, 
                reduction = "umap.atac", 
                label = TRUE, 
                label.size = 2.5, 
                repel = TRUE, 
                group.by = color.by) + 
    ggtitle("ATAC")
  p3 <- DimPlot(seurat.obj, 
                reduction = "wnn.umap", 
                label = TRUE, 
                label.size = 2.5, 
                repel = TRUE,
                group.by = color.by) + 
    ggtitle("WNN")
  p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
}

preFetchAnnotations <- function() {
  annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
  return(annotations)
}







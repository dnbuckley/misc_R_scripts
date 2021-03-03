library(monocle3)

# standard preprocessing as recommended by Trapnell lab
# http://cole-trapnell-lab.github.io/monocle-release/
preProcCDS <- function(cds, num_dim = 100, alignment_group = NULL, cores = 1) {
  cds <- preprocess_cds(cds, num_dim = num_dim, verbose = T)
  if (!is.null(alignment_group)) {
    message("align_cds by ", alignment_group, "...")
    cds = align_cds(cds, num_dim = num_dim, 
                    alignment_group = alignment_group, 
                    verbose = T)
  }
  message("reducing dimensions...")
  cds <- reduce_dimension(cds, verbose = T, cores = cores)
  return(cds)
}

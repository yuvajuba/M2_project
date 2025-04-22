#' Normalize and extract log-transformed counts for selected genes
#'
#' This function performs size-factor normalization followed by log10 transformation
#' of raw count data. It returns both scaled (z-score) and unscaled log-normalized
#' values for a selected subset of genes.
#'
#' @param RCounts A numeric matrix or data.frame of raw counts (genes x samples)
#' @param Sel_genes A named vector of genes (names must match rownames of RCounts)
#'
#' @return A list with two elements:
#' \describe{
#'   \item{scaled}{Z-score scaled log-normalized counts}
#'   \item{unscaled}{Raw log-normalized counts without scaling}
#' }
#'
#' @importFrom DESeq2 estimateSizeFactorsForMatrix
#' @export
#'
#' @examples
#' # Assuming 'counts' is your matrix and 'sel' a named vector of selected genes
#' data_out <- nlog_scale_data(counts, sel)
#' scaled_data <- data_out$scaled
#' unscaled_data <- data_out$unscaled
#'

nlog_scale_data <- function(RCounts, Sel_genes){
  
  # -- Checks ----
  if (!is.data.frame(RCounts) & !is.matrix(RCounts)) {
    stop("RCounts must be a data.frame or matrix.")
  }
  if(all(apply(RCounts,2,class) != "numeric")){
    stop("All your variables must be numeric values !")
  }
  
  if(is.null(names(Sel_genes)) || any(names(Sel_genes) == "")){
    stop("Sel_genes must be a named vector with the rownames of your RCounts")
  }
  if(all(!(names(Sel_genes) %in% rownames(RCounts)))){
    stop(paste("There are some unmatched genes between RCounts and Sel_genes",
               "rownames(RCounts) need to match names(Sel_genes)",
               sep = "\n"))
  }
  
  # -- Proceed ----
  ## 1. with scaling :
  n_fact <- estimateSizeFactorsForMatrix(RCounts)
  n_count <- sweep(RCounts, 2, n_fact, FUN = "/")
  n_count <- log10(n_count+1)
  n_count <- n_count[names(Sel_genes),]
  n_count <- t(apply(n_count,1,scale))
  colnames(n_count) <- colnames(RCounts)
  
  ## 2. without scaling
  log_counts <- sweep(RCounts, 2, n_fact, FUN = "/")
  log_counts <- log10(log_counts+1)
  log_counts <- log_counts[names(Sel_genes),]
  
  # -- return ----
  return(list(scaled = n_count,
              unscaled = log_counts))
}





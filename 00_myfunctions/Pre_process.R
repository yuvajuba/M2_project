#' Pre_process RNA-seq data
#'
#' @param RCounts Raw count matrix (genes x samples).
#' @param Exp.Conditions Data frame containing sample metadata.
#' @param condition_col Name of the column indicating experimental conditions.
#' @param replicate_col Name of the column indicating biological replicates.
#' @param keep_above Raw count threshold for filtering genes (default = 40).
#' @param cpm_count_threshold CPM threshold for filtering (default = 0.5).
#' @param apply_filt_samp Minimum number of samples meeting CPM threshold (default = 4).
#' @param pca_p.size Point size in PCA plot (default = 5).
#'
#' @return A list with:
#'   \item{raw_filt}{Filtered count matrix}
#'   \item{dge}{DGEList object with TMM normalization}
#'   \item{MDS_plt}{MDS plot object}
#'   \item{PCA_plt}{PCA ggplot object}
#'   \item{PCA_obj}{Raw PCA object}
#'
#' @export

Pre_process <- function(RCounts,
                        Exp.Conditions,
                        condition_col = "Conditions",
                        replicate_col = "Replicate",
                        keep_above = 40,
                        cpm_count_threshold = 0.5,
                        apply_filt_samp = 4,
                        pca_p.size = 5) {
  
  # ---- Checks ----
  if (!is.data.frame(RCounts) & !is.matrix(RCounts)) {
    stop("RCounts must be a data.frame or matrix.")
  }
  if (!all(c(condition_col, replicate_col) %in% colnames(Exp.Conditions))) {
    stop(paste("Exp.Conditions must contain columns:",
               condition_col, "and", replicate_col))
  }
  
  
  # ---- Filtering raw counts ----
  raw_filt <- RCounts[rowMaxs(as.matrix(RCounts)) >= keep_above, ]
  cpm_counts <- cpm(raw_filt)
  raw_filt <- raw_filt[rowSums(cpm_counts > cpm_count_threshold) >= apply_filt_samp, ]
  
  
  # ---- Normalization ----
  dge <- DGEList(raw_filt)
  dge <- calcNormFactors(dge, method = "TMM")
  log_cpm <- cpm(dge, log = TRUE, prior.count = 1) %>% as.data.frame()
  
  
  # ---- MDS ----
  group <- factor(Exp.Conditions[[condition_col]])
  mds <- plotMDS(dge,
                 labels = colnames(log_cpm),
                 col = as.numeric(group))
  
  
  # ---- PCA ----
  pca <- prcomp(t(log_cpm), scale. = TRUE)
  pca_plt <- autoplot(pca,
                      data = Exp.Conditions,
                      colour = replicate_col,
                      shape = condition_col,
                      size = pca_p.size) +
    labs(title = "PCA") +
    My_theme() +
    geom_hline(yintercept = 0, linetype = 2, color = "gray", linewidth = 0.5) +
    geom_vline(xintercept = 0, linetype = 2, color = "gray", linewidth = 0.5) +
    scale_colour_manual(values = MyPalette)
  
  
  # ---- Return ----
  list(
    raw_filt = raw_filt,
    dge = dge,
    MDS_plt = mds,
    PCA_plt = pca_plt,
    PCA_obj = pca
  )
}

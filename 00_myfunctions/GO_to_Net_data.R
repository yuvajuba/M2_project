#' Prepare GO & DEA Results for Network Plot
#'
#' This function merges a GO enrichment object with DEA results and prepares a tidy dataframe
#' that can be used for downstream visualization (e.g. network plots like cnetplot).
#'
#' @param go.object A `clusterProfiler` enrichResult object (from enrichGO, gseGO, etc.)
#' @param dea.results A data.frame of DEA results (typically the result of DESeq2 or limma)
#' @param condition.name Optional. A string to label the cluster/condition (e.g., "siC8"). Default is "Condition1".
#' @param go.term.colname Column name in the GO results corresponding to the GO term (default: "Description")
#' @param go.gene.colname Column name in the GO results with genes (default: "geneID")
#' @param go.pval.colname Column name for GO adjusted p-values (default: "p.adjust")
#' @param dea.gene.colname Column name for gene names in DEA table (default: "Gene_name")
#' @param dea.logFC.colname Column name for log2 fold change in DEA table (default: "log2FoldChange")
#' @param dea.pval.colname Column name for adjusted p-values in DEA table (default: "padj")
#'
#' @return A tidy data.frame with columns: Term, Gene, Cluster, and logFC
#' @export
#'
#' @examples
#' df <- Go_to_Net.data(go.object = myGO,
#'                      dea.results = myDEGs,
#'                      condition.name = "siC8")


Go_to_Net.data <- function(go.object,
                           dea.results,
                           condition.name = NULL,
                           go.term.colname = "Description",
                           go.gene.colname = "geneID",
                           go.pval.colname = "p.adjust",
                           dea.gene.colname = "Gene_name",
                           dea.logFC.colname = "log2FoldChange",
                           dea.pval.colname = "padj"){
  
  if (!inherits(go.object, "enrichResult")) {
    stop("You need to provide a GO object from enrichGO or similar.")
  }
  
  if(!(is.data.frame(dea.results))){
    stop("dea.results needs to be a dataframe")
  }
  
  GO.res <- go.object@result %>% 
    dplyr::filter(.data[[go.pval.colname]] < 0.05) %>% 
    dplyr::select(Term = go.term.colname, 
                  Gene = go.gene.colname)
  
  DEA.res <- dea.results %>% 
    dplyr::filter(.data[[dea.pval.colname]] < 0.05) %>% 
    dplyr::select(Gene = dea.gene.colname, 
                  logFC = dea.logFC.colname)
  
  
  dat <- GO.res %>% 
    dplyr::mutate(Cluster = ifelse(is.null(condition.name), "Condition1", condition.name))  %>% 
    tidyr::separate_rows(Gene, sep = "/") %>% 
    left_join(DEA.res, by = "Gene")
  
  
  
  return(dat)
  
}
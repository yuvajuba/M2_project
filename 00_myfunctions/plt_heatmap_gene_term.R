#' @title Heatmap of Genes across Enriched Terms
#'
#' @description
#' Visualize the presence/absence of selected genes across enriched GO terms as a heatmap.
#' Optionally, display the associated log2 fold-change values if available.
#'
#' @param obj A dataframe of enrichment results containing at least a GO term column and a gene list column (e.g., result from `enrichGO`).
#' @param terms_col The name of the column containing the GO terms (default: `"Description"`).
#' @param genelist_col The name of the column containing the gene IDs (default: `"geneID"`).
#' @param sel_genes A character vector of gene names to include. If `NULL`, the top `n_genes` are selected from the object (default: `NULL`).
#' @param n_genes Number of genes to show (used only if `sel_genes` is `NULL`, default: 50).
#' @param fc Logical. If `TRUE`, colors tiles by log2 fold change. If `FALSE`, shows presence/absence only (default: `FALSE`).
#' @param FCdata A dataframe with fold change values to be merged with the matrix (only used if `fc = TRUE`).
#' @param FC_gene_col Column name of the gene names in `FCdata` (default: `"Gene_name"`).
#' @param FC_col Column name of the log2FC values in `FCdata` (default: `"log2FoldChange"`).
#'
#' @return A ggplot2 heatmap object.
#'
#' @examples
#' \dontrun{
#' # Presence/absence heatmap
#' Hm_term.genes(go_result)
#'
#' # With log2FC
#' Hm_term.genes(go_result, fc = TRUE, FCdata = degs_df)
#' }
#'
#' @export



Hm_term.genes <- function(obj, 
                          terms_col = "Description", 
                          genelist_col = "geneID",
                          sel_genes = NULL,
                          n_genes = 50,
                          fc = F,
                          FCdata,
                          FC_gene_col = "Gene_name",
                          FC_col = "log2FoldChange"){
  
  library(tidyr)
  
  dat <- obj %>% 
    dplyr::select(terms_col, genelist_col) %>% 
    tidyr::separate_rows(geneID, sep = "/")
  
  if(is.null(sel_genes)){
    selected_genes <- unique(dat[[genelist_col]])[1:n_genes]
    
  } else {
    selected_genes <- sel_genes
  }
  
  ## set binary table for geom_tile
  dat <- dat %>% 
    dplyr::mutate(Present = 1) %>% 
    pivot_wider(names_from = genelist_col, values_from = Present, values_fill = 0) %>% 
    pivot_longer(-terms_col, names_to = "Genes", values_to = "Present")
  
  if(isFALSE(fc)){
    hm <- dat %>% 
      dplyr::filter(Genes %in% selected_genes) %>% 
      ggplot(aes(x= Genes, y= Description, fill= factor(Present)))+
      geom_tile(color = "white")+
      scale_fill_manual(values = c("0" = "white", "1" = "darkred"))+
      theme_minimal()+
      labs(y="", x="", fill="")+
      guides(fill = "none")+
      theme(axis.text.x = element_text(angle = 90, face = "bold", size = 10, colour = "black"),
            axis.text.y = element_text(face = "bold", size = 10, colour = "black"),
            panel.background = element_rect(colour = "black", linewidth = 0.4))
    
  } else {
    
    ## add fold change
    dat <- dat %>% 
      left_join(FCdata %>% dplyr::select(Genes = FC_gene_col, logFC = FC_col), by = "Genes") %>% 
      dplyr::mutate(logFC = case_when(Present == 1 ~ logFC,
                                      Present == 0 ~ 0,
                                      TRUE ~ 0))
    hm <- dat %>% 
      dplyr::filter(Genes %in% selected_genes) %>% 
      ggplot(aes(x= Genes, y= Description, fill= logFC))+
      geom_tile(color = "white")+
      scale_fill_gradient2(low = "midnightblue",
                           mid = "white", midpoint = 0,
                           high = "darkred", 
                           name = "log2FC")+
      theme_minimal()+
      labs(y="", x="")+
      theme(axis.text.x = element_text(angle = 90, face = "bold", size = 10, colour = "black"),
            axis.text.y = element_text(face = "bold", size = 10, colour = "black"),
            panel.background = element_rect(colour = "black", linewidth = 0.4))
    
  }
  
  
  return(hm)
}
#' My_volcano: Custom volcano plot for DEA results
#'
#' @param df A data frame with differential expression results.
#' @param x Column name for log2 fold change (default: "log2FoldChange").
#' @param y Column name for adjusted p-value (default: "padj").
#' @param lab Column name for gene labels (default: "Gene_name").
#' @param point_size Size of scatter points.
#' @param label_genes Whether to label selected genes (default: FALSE).
#' @param lab_size Font size for gene labels.
#' @param lab_max_overlap Max number of overlaps for label repel.
#' @param lab_box.padding Padding for labels.
#' @param lab_colour Color of label text.
#' @param lab_fill Background fill color for labels.
#' @param lab_force_pull Force pull for label repel.
#' @param selected_genes Optional vector of gene names to label.
#' @param logFC_thresh Threshold for log2 fold change.
#' @param padj_thresh Threshold for adjusted p-value.
#' @param col_palette Named color vector (optional).
#' @param extend_xlim Extra margin on x-axis.
#'
#' @return A ggplot2 volcano plot
#' @export

My_volcano <- function(df, 
                       x= "log2FoldChange", 
                       y= "padj", 
                       lab= "Gene_name",
                       point_size= 2,
                       label_genes= F,
                       lab_size= 3,
                       lab_max_overlap= 20,
                       lab_box.padding= 0.5,
                       lab_colour = "navy", 
                       lab_fill = "white", 
                       lab_force_pull = .5,
                       selected_genes= NULL,
                       logFC_thresh= 1.5,
                       padj_thresh= 1e-5,
                       col_palette= NULL,
                       extend_xlim= 0.5){
  
  if(!(is.data.frame(df))) stop("your object isn't a dataframe !")
  if(!all(c(x, y) %in% colnames(df))){
    stop(paste0("<",x,"> or <",y,"> isn't found in your dataframe !"))
  }
  
  ## Setup thresholds :
  ## ---------------- :
  hline <- -log10(padj_thresh)
  vline <- logFC_thresh
  
  ## Default gene selection :
  ## ---------------------- :
  if(is.null(selected_genes)){
    selected_genes <- df %>% 
      dplyr::filter(.data[[y]] < 1e-10, abs(.data[[x]]) > 2) %>% 
      dplyr::pull(.data[[lab]])
  }
  
  ## Default color palette :
  ## --------------------- :
  default_palette <- c("pval & logFC" = "darkred",
                       "logFC" = "darkgreen",
                       "p-value" = "midnightblue",
                       "NS" = "lightgray")
  
  col_palette <- if (is.null(col_palette) || length(col_palette) < 4) {
    default_palette
  } else {
    setNames(col_palette[1:4], names(default_palette))
  }
  
  
  ## Add color & label columns :
  ## ------------------------- :
  df_plot <- df %>%
    dplyr::mutate(
      Col = case_when(
        abs(.data[[x]]) > logFC_thresh & .data[[y]] < padj_thresh ~ "pval & logFC",
        abs(.data[[x]]) > logFC_thresh & .data[[y]] > padj_thresh ~ "logFC",
        abs(.data[[x]]) < logFC_thresh & .data[[y]] < padj_thresh ~ "p-value",
        TRUE ~ "NS"
      ),
      Sel_genes = ifelse(.data[[lab]] %in% selected_genes, .data[[lab]], NA)
    )
  
  
  
  ## Base plot :
  ## --------- :
  plt <- ggplot(df_plot, aes(x = .data[[x]], y = -log10(.data[[y]]))) +
    geom_vline(xintercept = c(-vline, vline), linetype = 2, color = "black", linewidth = 0.4) +
    geom_hline(yintercept = hline, linetype = 2, color = "black", linewidth = 0.4) +
    geom_point(aes(colour = Col), size = point_size, alpha = 0.6) +
    labs(title = "Volcano Plot",
         caption = paste0("Number of genes : ", nrow(df))) +
    scale_color_manual(values = col_palette) +
    xlim(c(min(df[[x]]) - extend_xlim, max(df[[x]]) + extend_xlim)) +
    guides(colour = guide_legend(title = NULL, override.aes = list(size = 5)))+
    theme(legend.position = "top",
          legend.margin = margin(b=.05, unit = "in"),
          legend.key = element_rect(colour = "white"),
          panel.background = element_rect(fill = "white",
                                          colour = "darkgray",
                                          linewidth = .5),
          legend.text = element_text(size = 13,
                                     color = "black",
                                     face = "bold"),
          axis.title = element_text(size = 14,
                                    face = "bold",
                                    colour = "darkred"),
          axis.title.x = element_text(margin = margin(t=12)),
          axis.title.y = element_text(margin = margin(r=12)),
          axis.text = element_text(size = 12,
                                   colour = "black"),
          plot.title = element_text(size = 16,
                                    face = "bold",
                                    colour = "darkred",
                                    margin = margin(b=12)),
          plot.caption = element_text(size = 13,
                                      colour = "black"))
  
  
  ## Add labels if needed :
  ## -------------------- :
  if (label_genes) {
    plt <- plt +
      ggrepel::geom_label_repel(aes(label = Sel_genes),
                                size = lab_size,
                                max.overlaps = lab_max_overlap,
                                box.padding = lab_box.padding,
                                label.r = 0.4,
                                label.size = 0.25,
                                force_pull = lab_force_pull,
                                colour = lab_colour,
                                fill = lab_fill)
  }
  
  
  return(plt)
}
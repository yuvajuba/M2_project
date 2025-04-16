# Clean start :
# ----------- :
rm(list = ls())
graphics.off()
cat("\014")

# SETUP ENVIRONNEMENT :
# =================== :
source("00_setup.R")

# Preprocessing parameters :
# ------------------------ :
min_count <- 40
cpm_thresh <- 0.5
samples_apply <- 4

# DEA filtering parameters :
# ------------------------ :
logfc_thresh <- 1
pval_thresh <- 0.05

# FUNCTIONS : ####

My_theme <- function(){
  theme_minimal()+
    theme(
      # customize the plot title
      plot.title = element_text(size = 15,
                                colour = "darkred",
                                face = "bold",
                                hjust = 0.5,
                                margin = margin(b = 12)),
      plot.subtitle = element_text(size = 13,
                                   colour = "black",
                                   hjust = 0.5),
      plot.background = element_rect(fill = "white",
                                     colour = "white",
                                     linewidth = 1),
      
      # customize the panels
      panel.background = element_rect(colour = "black",
                                      fill = "white",
                                      linewidth = 1),
      panel.grid = element_blank(),
      
      # customize axis
      axis.title = element_text(size = 13,
                                colour = "darkred",
                                face = "bold"),
      axis.title.x = element_text(margin = margin(t = 10)),
      axis.title.y = element_text(margin = margin(r = 10)),
      axis.text = element_text(size = 12,
                               colour = "black"),
      
      # customize legend
      legend.title = element_text(size = 13,
                                  colour = "darkred",
                                  face = "bold",
                                  margin = margin(b = 10)),
      legend.text = element_text(size = 12,
                                 colour = "black"),
    )
}

Pre_process <- function(RCounts,
                        Conditions,
                        keep_above = 40,
                        cpm_count_threshold = 0.5,
                        apply_filt_samp = 4,
                        pca_p.size = 5){
  
  rawcounts <- RCounts
  Experimental <- Conditions
  
  # filtering the rawcounts :
  raw_filt <- rawcounts %>% 
    dplyr::mutate(M = apply(rawcounts, 1, max)) %>% 
    dplyr::filter(M >= keep_above) %>% 
    dplyr::select(-M)
  
  cpm_counts <- cpm(raw_filt)
  raw_filt <- raw_filt[rowSums(cpm_counts > cpm_count_threshold) >= apply_filt_samp,]
  
  # normalizing the data :
  dge <- DGEList(raw_filt)
  dge <- calcNormFactors(dge, method = "TMM")
  log_cpm <- cpm(dge, log = T, prior.count = 1) %>% as.data.frame()
  
  # MDS plot :
  group <- factor(Experimental$Conditions)
  mds <- plotMDS(dge,
                 labels= colnames(log_cpm),
                 col= as.numeric(group))
  
  # PCA plot :
  pca <- prcomp(t(log_cpm), scale. = TRUE)
  pca_plt <- autoplot(pca, 
                      data = Experimental, 
                      colour = "Replicate",
                      shape = "Conditions",
                      size = pca_p.size)+
    labs(title = "PCA")+
    My_theme()+
    geom_hline(yintercept = 0, linetype = 2, color = "gray", linewidth = .5)+
    geom_vline(xintercept = 0, linetype = 2, color = "gray", linewidth = .5)+
    scale_colour_manual(values = MyPalette)
  
  
  ## Regrouping the outputs :
  Outs <- list(raw_filt = raw_filt,
               dge = dge,
               MDS_plt = mds,
               PCA_plt = pca_plt,
               PCA_obj = pca)
  
  return(Outs)
}

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
  
  if(!(is.data.frame(df))){
    stop("your object isn't a dataframe !")
  }
  
  if(!all(c(x, y) %in% colnames(df))){
    stop(paste0("<",x,"> or <",y,"> isn't found in your dataframe !"))
  }
  
  hline <- -log10(padj_thresh)
  vline <- logFC_thresh
  
  ## Select genes :
  if(is.null(selected_genes)){
    selected_genes <- df %>% 
      dplyr::filter(.data[[y]] < 1e-10,
                    abs(.data[[x]]) > 2) %>% 
      dplyr::pull(.data[[lab]])
  }
  
  ## Color palette :
  if(is.null(col_palette)){
    col_palette <- c("pval & logFC"="darkred",
                     "logFC"="darkgreen",
                     "p-value"="midnightblue",
                     "NS"="lightgray")
  } else {
    if(length(col_palette) < 4){
      warning("The provided palette is of length < 4 ; therefor we're using default palette")
      col_palette <- c("pval & logFC"="darkred",
                       "logFC"="darkgreen",
                       "p-value"="midnightblue",
                       "NS"="lightgray")
      
    } else {
      col_palette <- c("pval & logFC" = col_palette[1],
                       "logFC" = col_palette[2],
                       "p-value" = col_palette[3],
                       "NS" = col_palette[4])
    }
  }
  
  
  ## plot :
  df %>% 
    dplyr::mutate(Col = case_when(abs(.data[[x]]) > logFC_thresh & 
                                    .data[[y]] < padj_thresh ~ "pval & logFC",
                                  abs(.data[[x]]) > logFC_thresh & 
                                    .data[[y]] > padj_thresh ~ "logFC",
                                  abs(.data[[x]]) < logFC_thresh & 
                                    .data[[y]] < padj_thresh ~ "p-value",
                                  abs(.data[[x]]) < logFC_thresh & 
                                    .data[[y]] > padj_thresh ~ "NS",
                                  TRUE ~ "NS"),
                  Sel_genes = ifelse(.data[[lab]] %in% selected_genes,
                                     .data[[lab]],
                                     NA)) -> df.plot
  
  
  
  if(isFALSE(label_genes)){
    
    plt <- df.plot %>% 
      ggplot(aes(x= .data[[x]],
                 y= -log10(.data[[y]])))+
      geom_vline(xintercept = -(vline), linetype = 2, color = "black", linewidth = .4)+
      geom_vline(xintercept = vline, linetype = 2, color = "black", linewidth = .4)+
      geom_hline(yintercept = hline, linetype = 2, color = "black", linewidth = .4)+
      geom_point(aes(colour = Col),
                 size = point_size,
                 alpha = 0.6)+
      labs(title = "Volcano Plot",
           caption = paste0("Number of genes : ", nrow(df)))+
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
                                        colour = "black"))+
      guides(colour = guide_legend(title = NULL,
                                   override.aes = list(size = 5)))+
      scale_color_manual(values = col_palette)+
      xlim(c((min(df[[x]]) - extend_xlim),(max(df[[x]]) + extend_xlim)))
    
  } else {
    
    plt <- df.plot %>% 
      ggplot(aes(x= .data[[x]],
                 y= -log10(.data[[y]])))+
      geom_vline(xintercept = -(vline), linetype = 2, color = "black", linewidth = .4)+
      geom_vline(xintercept = vline, linetype = 2, color = "black", linewidth = .4)+
      geom_hline(yintercept = hline, linetype = 2, color = "black", linewidth = .4)+
      geom_point(aes(colour = Col),
                 size = point_size,
                 alpha = 0.6)+
      labs(title = "Volcano Plot",
           caption = paste0("Number of genes : ", nrow(df)))+
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
                                        colour = "black"))+
      guides(colour = guide_legend(title = NULL,
                                   override.aes = list(size = 5)))+
      scale_color_manual(values = col_palette)+
      xlim(c((min(df[[x]]) - extend_xlim),(max(df[[x]]) + extend_xlim)))+
      geom_label_repel(aes(label = .data[["Sel_genes"]]),
                       size = lab_size,
                       max.overlaps = lab_max_overlap, 
                       box.padding = lab_box.padding,
                       label.r = 0.4,
                       label.size = .25, 
                       force_pull = lab_force_pull,
                       colour = lab_colour, 
                       fill = lab_fill)
    
  }
  
  
  return(plt)
}


# IMPORT DATA : ####

rawcounts <- read_xlsx(paste0(data_dir,"C8_rawcounts.xlsx")) %>% column_to_rownames("ID")
metadata <- read_xlsx(paste0(data_dir,"C8_metadata.xlsx"))
Experimental <- data.frame(row.names = colnames(rawcounts),
                           Conditions = str_split_i(colnames(rawcounts), "_", 1),
                           Replicate = str_split_i(colnames(rawcounts), "_", 2))

# PROCESSING : ####

## Counts distribution : ####
## ----------------------- ##
rawcounts %>% 
  pivot_longer(cols = colnames(rawcounts),
               names_to = "Samples",
               values_to = "Counts") %>% 
  mutate(Conditions = str_split_i(Samples, "_", 1)) %>% 
  ggplot()+
  geom_bar(aes(x= Samples,
               y= Counts/10^6,
               fill= Conditions),
           stat = "identity")+
  scale_fill_manual(values = c("gold","brown"))+
  My_theme()+
  theme(axis.text.x = element_text(size = 13, angle = 25, margin = margin(t=10)),
        axis.title.x = element_blank(),
        legend.title = element_text(size = 12,
                                    face = "bold",
                                    colour = "darkred"),
        legend.text = element_text(size = 10,
                                   colour = "black"))+
  labs(title = "Read Counts",
       y = "Counts in million",
       x = "") 


## Processing : ####
## -------------- ##

## 1. Dim reduction :
## ---------------- :
PreProcess_res <- Pre_process(RCounts = rawcounts,
                              Conditions = Experimental,
                              keep_above = min_count,
                              cpm_count_threshold = cpm_thresh,
                              apply_filt_samp = samples_apply)

## 2. Filtering :
## ------------ :
Experimental <- Experimental[which(!(Experimental$Replicate %in% c("n2","n5"))),]
rawcounts <- rawcounts[,rownames(Experimental)]

## 3. Dim reduction again :
## ---------------------- :
PreProcess_res <- Pre_process(RCounts = rawcounts,
                              Conditions = Experimental,
                              keep_above = min_count,
                              cpm_count_threshold = cpm_thresh,
                              apply_filt_samp = samples_apply,
                              pca_p.size = 5)

## 4. Final filtered counts :
## ------------------------ :
raw_filt <- PreProcess_res$raw_filt
dge <- PreProcess_res$dge


# DEA : ####

## Get results : 
## ----------- :
DESeqObj <- DESeqDataSetFromMatrix(countData = dge, 
                                   colData = Experimental,
                                   design = ~ Conditions)
DESeqObj$Conditions <- relevel(DESeqObj$Conditions, ref = "siLuc")
dds <- DESeq(DESeqObj)
res <- results(dds, 
               contrast =  c("Conditions","siC8","siLuc"),
               alpha = 0.05) %>%
  as.data.frame() %>% 
  rownames_to_column(var = "ID") %>% 
  inner_join(metadata[,c("ID","Gene_name")],by = "ID") %>% 
  dplyr::filter(!(is.na(padj)))

raw_filt <- raw_filt[res$ID,] # removing NULL p values from rawcounts

## Filtering results :
## ----------------- :
res_filt <- res %>% 
  dplyr::filter(padj < pval_thresh,
                abs(log2FoldChange) > logfc_thresh) %>% 
  dplyr::select(ID, Gene_name, baseMean, log2FoldChange, padj) %>% 
  dplyr::arrange(desc(log2FoldChange))

## Merging results :
## --------------- :
AllResults <- res  %>%  
  dplyr::mutate(Log_baseMean = log10(baseMean),
                Rank = case_when(padj < 1e-100 ~ 
                                   round(log2FoldChange*(-log10(1e-50))/20, 4),
                                 padj >= 1e-100 & padj < 1e-30 ~ 
                                   round(log2FoldChange*(-log10(1e-30))/20, 4),
                                 padj >= 1e-30 & padj < 1e-15 ~ 
                                   round(log2FoldChange*(-log10(1e-15))/20, 4),
                                 padj >= 1e-15 & padj < 1e-10 ~ 
                                   round(log2FoldChange*(-log10(1e-10))/20, 4),
                                 padj >= 1e-10 & padj < 1e-05 ~ 
                                   round(log2FoldChange*(-log10(1e-05))/20, 4),
                                 padj >= 1e-05 & padj < 5e-02 ~ 
                                   round(log2FoldChange*(-log10(padj))/20, 4),
                                 padj >= 5e-02 & padj < 0.5 ~ 
                                   round(log2FoldChange*(-log10(0.25))/20, 4),
                                 padj >= 0.25 ~ 0),
                DEGs = case_when(padj < pval_thresh & 
                                   abs(log2FoldChange) > logfc_thresh ~ T,
                                 TRUE ~ F)) %>% 
  dplyr::arrange(desc(Rank)) %>%  
  dplyr::select(ID,Gene_name,log2FoldChange,padj,Rank,DEGs,baseMean,Log_baseMean) %>% 
  column_to_rownames(var = "ID")


AllResults <- merge(AllResults, raw_filt, by = 0) %>% 
  dplyr::rename(ID = Row.names) %>% 
  dplyr::arrange(desc(Rank))

## Saving results :
## -------------- :
write_xlsx(AllResults, path = paste0(obj_dir,"C8_AllResults.xlsx"))


# VISUALIZATION : ####

## Volcano plot :
## ------------ :
My_volcano(df = res, 
           x= "log2FoldChange", 
           y= "padj", 
           lab= "Gene_name",
           label_genes= T,
           point_size = 2,
           lab_size= 3,
           lab_max_overlap= 20,
           lab_box.padding= 0.5,
           lab_colour = "navy", 
           lab_fill = "white", 
           lab_force_pull = .5,
           selected_genes= NULL,
           logFC_thresh= logfc_thresh,
           padj_thresh= pval_thresh,
           col_palette= NULL,
           extend_xlim= 0.5)+
  labs(caption = paste(paste0("Total number of genes : ", nrow(res)),
                       paste0("Number of up-regulated genes : ", nrow(res_filt[res_filt$log2FoldChange > 0,])),
                       paste0("Number of down-regulated genes : ", nrow(res_filt[res_filt$log2FoldChange < 0,])), 
                       sep = "\n"))+
  theme(plot.caption = element_text(size = 13, face = "bold.italic", colour = "#333"))



# PEA : ####

## Running GO :
## ---------- :
GOobj <- enrichGO(gene = res_filt$Gene_name,
                  OrgDb = "org.Hs.eg.db",
                  keyType = "SYMBOL",
                  minGSSize = 10,
                  maxGSSize = 1500,
                  readable = T,
                  ont = "BP")

## Extract results :
## --------------- :
GO.res <- GOobj@result %>% 
  dplyr::filter(p.adjust < 0.05) %>% 
  dplyr::mutate(RichFactor = round(Count / as.numeric(sub("/\\d+","",BgRatio)),5),
                BgGenes = as.numeric(sub("/\\d+","",BgRatio))) %>% 
  dplyr::select(Description, RichFactor, p.adjust, GeneRatio, BgRatio, Count, BgGenes, geneID) %>% 
  dplyr::arrange(p.adjust)

## Dot plot visualization :
## ---------------------- :
GO.res %>% 
  head(20) %>% 
  dplyr::mutate(Description = ifelse(nchar(Description) <= 60,
                                     Description,
                                     paste0(substr(Description,1,56), "...."))) %>% 
  ggplot(aes(x= RichFactor, y= fct_reorder(Description, RichFactor)))+
  geom_segment(aes(xend= 0, yend= Description))+
  geom_point(aes(color= p.adjust, size= Count))+
  scale_color_viridis_c(guide = guide_colorbar(reverse = T))+
  scale_size_continuous(range = c(3,10))+
  theme_linedraw()+
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5, colour = "darkred", 
                                  margin = margin(b=0.1, unit = "in")),
        plot.subtitle = element_text(size = 14, face = "bold", hjust = 0.5, colour = "brown", 
                                     margin = margin(b=0.2, unit = "in")),
        plot.caption = element_text(size = 14, face = "bold.italic", colour = "#333",
                                    margin = margin(t=0.2, unit = "in")),
        plot.margin = unit(c(0.1,0.1,0.1,0.1), "in"),
        axis.title.x = element_text(size = 16, face = "bold", colour = "darkred", 
                                    margin = margin(t=0.2, unit = "in")),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 15, face = "bold", colour = "darkred",
                                    margin = margin(b=0.2, unit = "in")),
        legend.text = element_text(size = 13),
        legend.box.margin = margin(l=0.2, unit = "in"))+
  labs(title = "GO Enriched Terms",
       subtitle = "Top 20 most significant terms",
       caption = paste(paste0("Total number of enriched terms : ", nrow(GOobj@result)),
                       paste0("Number of significant terms : ", nrow(GO.res)), 
                       sep = "\n"))

## Saving results :
## -------------- :
saveRDS(GOobj, file = paste0(obj_dir, "C8_GOobject.rds"))






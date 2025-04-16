# Clean start :
# ----------- :
rm(list = ls())
graphics.off()
cat("\014")

# SETUP ENVIRONNEMENT :
# =================== :
source("00_setup.R")

# Importing data ----

data <- list()
for(i in Conditions){
  
  rawcounts <- read_xlsx(paste0(data_dir,i,"_rawcounts.xlsx")) %>% column_to_rownames(var = "ID")
  metadata <- read_xlsx(paste0(data_dir,i,"_metadata.xlsx"))
  Experimental <- data.frame(row.names = colnames(rawcounts),
                             Conditions = str_split_i(colnames(rawcounts), "_", 1),
                             Replicate = str_split_i(colnames(rawcounts), "_", 2))
  
  if(i == "C8"){
    Experimental <- Experimental[which(!(Experimental$Replicate %in% c("n2","n5"))),]
    rawcounts <- rawcounts[,rownames(Experimental)]
  } else {
    Experimental <- Experimental[which(!(Experimental$Replicate %in% c("n1"))),]
    rawcounts <- rawcounts[,rownames(Experimental)]
  }
  
  PreProcess_res <- Pre_process(RCounts = rawcounts,
                                Conditions = Experimental,
                                keep_above = 40,
                                cpm_count_threshold = 0.5,
                                apply_filt_samp = ifelse(i == "C8", 4, 3),
                                pca_p.size = 5)
  
  raw_filt <- PreProcess_res$raw_filt
  dge <- PreProcess_res$dge
  
  DESeqObj <- DESeqDataSetFromMatrix(countData = dge, 
                                     colData = Experimental,
                                     design = ~ Conditions)
  DESeqObj$Conditions <- relevel(DESeqObj$Conditions, ref = "siLuc")
  dds <- DESeq(DESeqObj)
  if(i == "C8"){
    res <- results(dds, 
                   contrast =  c("Conditions","siC8","siLuc"),
                   alpha = 0.05) %>%
      as.data.frame() 
  } else {
    res <- results(dds, 
                   contrast =  c("Conditions","siCTSB","siLuc"),
                   alpha = 0.05) %>%
      as.data.frame() 
  }
  
  res <- res %>% 
    rownames_to_column(var = "ID") %>% 
    inner_join(metadata[,c("ID","Gene_name")],by = "ID") %>% 
    dplyr::filter(!(is.na(padj)))
  
  raw_filt <- raw_filt[res$ID,]
  
  res_filt <- res %>% 
    dplyr::filter(padj < 0.05,
                  abs(log2FoldChange) > 1) %>% 
    dplyr::select(ID, Gene_name, baseMean, log2FoldChange, padj) %>% 
    dplyr::arrange(desc(log2FoldChange))
  
  if(i == "C8"){
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
                    DEGs = case_when(padj < 0.05 & 
                                       abs(log2FoldChange) > 1 ~ T,
                                     TRUE ~ F)) %>% 
      dplyr::arrange(desc(Rank)) %>%  
      dplyr::select(ID,Gene_name,log2FoldChange,padj,Rank,DEGs,baseMean,Log_baseMean) %>% 
      column_to_rownames(var = "ID")
    
    AllResults <- merge(AllResults, raw_filt, by = 0) %>% 
      dplyr::rename(ID = Row.names) %>% 
      dplyr::arrange(desc(Rank))
    
  } else {
    AllResults <- res  %>%  
      dplyr::mutate(Log_baseMean = log10(baseMean),
                    Rank = case_when(padj < 1e-20 ~ 
                                       round(log2FoldChange*(-log10(1e-20))/20, 4),
                                     padj >= 1e-20 & padj < 1e-05 ~ 
                                       round(log2FoldChange*(-log10(1e-10))/20, 4),
                                     padj >= 1e-05 & padj < 5e-02 ~ 
                                       round(log2FoldChange*(-log10(padj))/20, 4),
                                     padj >= 5e-02 & padj < 0.5 ~ 
                                       round(log2FoldChange*(-log10(0.25))/20, 4),
                                     padj >= 0.25 ~ 0),
                    DEGs = case_when(padj < 0.05 & 
                                       abs(log2FoldChange) > 1 ~ T,
                                     TRUE ~ F)) %>% 
      dplyr::arrange(desc(Rank)) %>%  
      dplyr::select(ID,Gene_name,log2FoldChange,padj,Rank,DEGs,baseMean,Log_baseMean) %>% 
      column_to_rownames(var = "ID")
    
    AllResults <- merge(AllResults, raw_filt, by = 0) %>% 
      dplyr::rename(ID = Row.names) %>% 
      dplyr::arrange(desc(Rank))
  }
  
  data[[paste0(i,"_AllResults")]] <- AllResults
  
  GOobj <- enrichGO(gene = res_filt$Gene_name,
                    OrgDb = "org.Hs.eg.db",
                    keyType = "SYMBOL",
                    minGSSize = 10,
                    maxGSSize = 1500,
                    readable = T,
                    ont = "BP")
  
  data[[paste0(i,"_GOobject")]] <- GOobj
  
}

## --- Remove unnecessary variables :
rm(AllResults, dds, DESeqObj, dge, GOobj, Experimental, metadata, PreProcess_res, 
   raw_filt, rawcounts, res, res_filt, i, List_packages, function_files, cond, 
   load_or_install, Pre_process)

## --- Get the data to use
res.c8 <- data$C8_AllResults
res.ctsb <- data$CTSB_AllResults

go.c8 <- data$C8_GOobject@result[data$C8_GOobject@result$p.adjust < 0.05,]
go.ctsb <- data$CTSB_GOobject@result[data$CTSB_GOobject@result$p.adjust < 0.05,]



# Selecting terms ----

GOterms_comparison <- Display_Venn(Markers = list(C8 = go.c8$Description, 
                                                  CTSB = go.ctsb$Description), 
                                   colpalette = NULL, 
                                   set.names = NULL, 
                                   set.name.size = 5,
                                   text.size = 5,
                                   Padding = 0.06)

GOterms_comparison$plot





## We selected these terms below among the shared pathways between the 2 analysis
Selected_terms <- readLines("selected_terms.txt")
all(Selected_terms %in% GOterms_comparison$intersections$`C8 ∩ CTSB`)


netdat.c8 <- Go_to_Net.data(go.object = data$C8_GOobject, 
                         dea.results =res.c8,
                         condition.name = "siC8")
netdat.ctsb <- Go_to_Net.data(go.object = data$CTSB_GOobject, 
                           dea.results =res.ctsb,
                           condition.name = "siCTSB",
                           dea.gene.colname = "Gene_name")

netdata <- rbind(netdat.c8, netdat.ctsb) %>% 
  dplyr::filter(Term %in% Selected_terms)
  

# Global overview : 
# --------------- :

b.dat1 <- netdata %>% distinct(Gene, Cluster) %>% dplyr::group_by(Gene) %>% 
  dplyr::summarise(group = ifelse(n()>1, "Shared", Cluster[1]), .groups = "drop")

b.data <- netdata %>% 
  dplyr::left_join(b.dat1, by = "Gene") %>% 
  dplyr::group_by(Term, group) %>% 
  dplyr::summarise(n.genes = n(), .groups = "drop")

rm(b.dat1)

b.data %>% 
  ggplot()+
  geom_bar(aes(x = n.genes,
               y = Term,
               fill = group),
           stat = "identity",
           position = "stack")+
  labs(x = "Number of genes",
       y = "",
       fill = "Genes")+
  theme(panel.background = element_rect(fill = "white", colour = "#ddd"),
        axis.title.x = element_text(size = 13, colour = "darkred", face = "bold", 
                                    margin = margin(t=0.2, unit = "in")),
        axis.text = element_text(size = 12, colour = "black", face = "bold"),
        legend.title = element_text(size = 13, color = "darkred", face = "bold",
                                    margin = margin(b=0.2, unit = "in")))+
  scale_fill_manual(values = c("#600", "#880", "#088"))



## Network
Net_plot2(object = netdata, 
          select_terms = Selected_terms[c(3,5,7,13)],
          net.layout = "fr",
          edge.alpha = 0.2,
          edge.colors = NULL,
          edge.width = .5,
          node.term.size1 = 9,
          node.term.size2 = 5,
          node.term.color1 = "black",
          node.term.color2 = "black",
          node.gene.size = 3,
          node.gene.alpha = .8,
          label.term.show = F,
          label.term.size = 3,
          label.term.color = "navy",
          label.term.fill = "white",
          label.term.padding = 0.2,
          label.gene.size = 3,
          label.gene.colors = c("black","blue","red")) 




# Selecting genes ----

## 1. Pre-selection : 
##    Genes shared between these 4 pathways (excluding distinct genes)
term_list <- list()
for(i in Selected_terms[c(3,5,7,13)]){
  term_list[[i]] <- unique(netdata$Gene[which(netdata$Term == i)])
}

Genes_comparison <- Display_Venn(Markers = term_list, 
                                 colpalette = NULL, 
                                 set.names = NULL, 
                                 set.name.size = 5,
                                 text.size = 5,
                                 Padding = 0.03)

Genes_comparison$plot

gene_to_remove <- unique(c(Genes_comparison$group_specific$`unique each`$`inflammatory response`,
                           Genes_comparison$group_specific$`unique each`$`response to cytokine`,
                           Genes_comparison$group_specific$`unique each`$`leukocyte activation`,
                           Genes_comparison$group_specific$`unique each`$`response to bacterium`))

to_keep <- setdiff(c(term_list$`inflammatory response`,
                     term_list$`response to cytokine`,
                     term_list$`leukocyte activation`,
                     term_list$`response to bacterium`), gene_to_remove) %>% unique()


## 2. Final selection
##    Filtering the pre-selected list by log2FC (1.5)

hm.data <- go.c8[go.c8$Description %in% Selected_terms,]
tmp <- res.c8[res.c8$Gene_name %in% to_keep, c("Gene_name","log2FoldChange")]

selected_genes <- tmp$Gene_name[which(abs(tmp$log2FoldChange) > 1.5)]



## Heatmap
########## #
Hm_term.genes(obj = hm.data, 
              sel_genes = selected_genes,
              fc = T,
              FCdata = res.c8,
              FC_gene_col = "Gene_name",
              FC_col = "log2FoldChange")


## Volcano 
########## #
My_volcano(df = res.c8, 
           x= "log2FoldChange", 
           y= "padj", 
           lab= "Gene_name",
           point_size= 2,
           label_genes= T,
           lab_size= 3,
           lab_max_overlap= 20,
           lab_box.padding= 0.5,
           lab_colour = "navy", 
           lab_fill = "white", 
           lab_force_pull = .5,
           selected_genes= selected_genes,
           logFC_thresh= 1,
           padj_thresh= 5e-2,
           col_palette= NULL,
           extend_xlim= 0.5)





# ## Hm : Exp-Gene ###
# ## --------------- -
# 
# c8.counts <- res.c8[,c(1,10:17)] %>% column_to_rownames(var = "ENSEMBL") 
# list_genes <- selected_genes %>% 
#   setNames(res.c8$ENSEMBL[which(res.c8$Gene_name %in% selected_genes)])
# 
# 
# n_fact <- estimateSizeFactorsForMatrix(c8.counts)
# n_count <- sweep(c8.counts, 2, n_fact, FUN="/")
# n_count <- log10(n_count+1)
# n_count <- n_count[names(list_genes),]
# n_count <- t(apply(n_count,1,scale))
# colnames(n_count) <- colnames(c8.counts)
# 
# log_counts <- sweep(c8.counts, 2, n_fact, FUN="/")
# log_counts <- log10(log_counts+1)
# log_counts <- log_counts[names(list_genes),]
# 
# Cond <- as.factor(str_split_i(colnames(n_count),"_",1))
# ColConditions <- setNames(c("purple","khaki4"), levels(Cond))
# 
# HAnnot <- HeatmapAnnotation(Condition = Cond,
#                             col = list(Condition = ColConditions),
#                             annotation_name_side = "left",
#                             show_annotation_name = F,
#                             show_legend = F,
#                             annotation_name_gp = list(fontsize = 10,
#                                                       col = "navy",
#                                                       fontface = "bold"))
# 
# HAnnot2 <- HeatmapAnnotation(LogCounts = anno_boxplot(log_counts, 
#                                                       axis = TRUE, 
#                                                       gp = gpar(fill = "orange", col = "black"),
#                                                       height = unit(2, "cm")),
#                              annotation_name_side = "left",
#                              annotation_name_gp = gpar(fontsize = 9,
#                                                        col = "darkred",
#                                                        fontface = "bold"),
#                              annotation_name_rot = 90)
# 
# 
# 
# h1 <- Heatmap(n_count,
#               cluster_rows = T, 
#               bottom_annotation = HAnnot,
#               top_annotation = HAnnot2,
#               column_title = "Heatmap", 
#               column_title_gp = gpar(fontsize = 16, 
#                                      col = "darkred",
#                                      fontface = "bold"),
#               show_row_dend = F,
#               show_column_dend = F,
#               column_names_gp = gpar(col = "darkred", 
#                                      fontface = "bold",
#                                      fontsize = 11),
#               column_names_rot = 60,
#               column_labels = colnames(n_count),
#               row_labels = res.c8$Gene_name[which(res.c8$ENSEMBL %in% names(list_genes))],
#               name = "Z-score",
#               cluster_columns = F,
#               col = colorRamp2(c(-2,0,2), 
#                                c("darkgreen","white","darkred")))
# 
# h1



# Saving instance ----

sink("session_info.txt")
devtools::session_info()
sink()








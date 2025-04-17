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
samples_apply <- 3

# DEA filtering parameters :
# ------------------------ :
logfc_thresh <- 1
pval_thresh <- 0.05



# IMPORT DATA : ####

rawcounts <- read_xlsx(paste0(data_dir,"CTSB_rawcounts.xlsx")) %>% column_to_rownames(var = "ID")
metadata <- read_xlsx(paste0(data_dir,"CTSB_metadata.xlsx"))
Experimental <- data.frame(row.names = colnames(rawcounts),
                           Conditions = str_split_i(colnames(rawcounts), "_", 1),
                           Replicate = str_split_i(colnames(rawcounts), "_", 2))



# Counts distribution : ####
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

# Processing : ####
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
Experimental <- Experimental[which(!(Experimental$Replicate %in% c("n1"))),]
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
               contrast =  c("Conditions","siCTSB","siLuc"),
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
                Rank = case_when(padj < 1e-20 ~ 
                                   round(log2FoldChange*(-log10(1e-20))/20, 4),
                                 padj >= 1e-20 & padj < 1e-05 ~ 
                                   round(log2FoldChange*(-log10(1e-10))/20, 4),
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
write_xlsx(AllResults, path = paste0(obj_dir,"CTSB_AllResults.xlsx"))




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
  theme(plot.caption = element_text(size = 13, face = "bold.italic", colour = "#333"))+
  ylim(c(0,12))


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
saveRDS(GOobj, file = paste0(obj_dir, "CTSB_GOobject.rds"))

















































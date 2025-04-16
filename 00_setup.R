# ================================================
#   Project Setup - RNAseq Analysis
# ================================================

## 1. Load or Install Packages
load_or_install <- function(pkgs){
  for(pkg in pkgs){
    if (!require(pkg, character.only = TRUE)) {
      if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
      library(pkg, character.only = TRUE)
    }
  }
}

# -- List of packages for this project
List_packages <- c(
  "DESeq2", "edgeR", "dplyr", "tibble", "stringr", "ggplot2", 
  "writexl", "readxl", "tidyr", "ggfortify", "ggrepel",
  "ComplexHeatmap", "circlize", "forcats", "clusterProfiler",
  "enrichplot", "AnnotationDbi", "org.Hs.eg.db", "purrr"
)

# -- Load all required packages
load_or_install(List_packages)


## 2. Load custom project functions
function_files <- list.files("00_myfunctions", full.names = TRUE, pattern = "\\.R$")
sapply(function_files, source)


## 3. Define output folders (optional but recommended)
fig_dir <- "~/Bureau/Projects/RNAseq/M2_project/out_fig/"
obj_dir <- "~/Bureau/Projects/RNAseq/M2_project/out_obj/"
data_dir <- "~/Bureau/Projects/RNAseq/M2_project/data/"


## 4. Define variables :

# --- Color palette :
MyPalette <- c("#9933aa","#ffdd22","#aa4400","#ff0000","#337722","#00ff66","#005566","#002277",
               "#441144","#aa0077","#00bbff","#003333","#4422cc","#116611","#330077","#111111",
               "#667700","#ddaa00","#33ffff","#ff22ff","#ffff33","#00ff00","#0000ff","#444444")


# --- Conditions vector :
Conditions <- c("C8", "CTSB")


## 5. Define some functions :

# --- ggplot theme function
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

# --- pre-processing data
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

# --- Venn diagram for comparison
Display_Venn <- function(Markers, 
                         colpalette = NULL, 
                         set.names = NULL, 
                         set.name.size = 5,
                         text.size = 5,
                         Padding = 0.03) {
  
  # package install checking
  pkgs <- c("ggvenn","purrr")
  for(pkg in pkgs){
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg)
    }
  }
  
  library(ggvenn)
  library(purrr)
  
  # Arguments validation
  if (is.null(Markers)) {
    stop("Argument 'Markers' is NULL. Please provide a list of sets.")
  }
  
  if (!is.list(Markers)) {
    stop("Argument 'Markers' must be a list of vectors.")
  }
  
  if (!all(sapply(Markers, is.vector))) {
    stop("All elements in the 'Markers' list must be vectors.")
  }
  
  
  n <- ifelse(length(Markers) >= 4, 4, length(Markers))
  
  if (n < 2) {
    stop("You need at least 2 sets to make a comparison.")
  }
  
  if (length(Markers) > 4) {
    warning("The function supports up to 4 sets only. Taking the first 4 vectors as arguments")
  }
  
  # Handling colours
  default_palette <- c("navy", "red", "darkgreen", "violet")
  if (is.null(colpalette) || length(colpalette) < n) {
    warning("The provided colour palette is too short or NULL. Using default colours.")
    colpalette <- default_palette[1:n]
  } else if(length(colpalette) >= n){
    colpalette <- colpalette[1:n]
  }
  
  # Handling sets
  sets <- Markers[1:n]
  if (!(is.null(set.names))) {
    if (length(set.names) != n) {
      stop("The length of 'set.names' must match the number of sets in 'Markers'")
    }
    names(sets) <- set.names
  } else if (is.null(names(sets))) {
    names(sets) <- paste0("Set", 1:n)
  }
  
  # Calculate intersections
  common_genes <- list()
  for(i in 2:n){
    x <- combn(sets, i, simplify = F)
    for(s in 1:length(x)){
      if(i == 2){
        common_genes[[paste(names(x[[s]]), collapse = " ∩ ")]] <- intersect(x[[s]][[1]], x[[s]][[2]])
      }
      
      if(i == 3){
        common_genes[[paste(names(x[[s]]), collapse = " ∩ ")]] <- purrr::reduce(list(x[[s]][[1]],
                                                                                     x[[s]][[2]],
                                                                                     x[[s]][[3]]), 
                                                                                intersect)
      }
      
      if(i == 4){
        common_genes[["Common all"]] <- purrr::reduce(x[[s]], intersect)
      }
    }
  }
  
  # Calculate group specific genes
  specific_genes <- list()
  for(i in 2:n){
    x <- combn(sets, i, simplify = F)
    for(s in 1:length(x)){
      if(i == 2){
        specific_genes[[paste(names(x[[s]]), collapse = " | ")]] <- setNames(list(setdiff(x[[s]][[1]], 
                                                                                          intersect(x[[s]][[1]], 
                                                                                                    x[[s]][[2]])),
                                                                                  setdiff(x[[s]][[2]], 
                                                                                          intersect(x[[s]][[1]], 
                                                                                                    x[[s]][[2]]))),
                                                                             c(names(x[[s]][1]), 
                                                                               names(x[[s]][2])))
      }
      
      if(i == 3){
        specific_genes[[paste(names(x[[s]]), collapse = " | ")]] <- setNames(list(setdiff(x[[s]][[1]], 
                                                                                          union(x[[s]][[2]], 
                                                                                                x[[s]][[3]])),
                                                                                  setdiff(x[[s]][[2]], 
                                                                                          union(x[[s]][[1]], 
                                                                                                x[[s]][[3]])),
                                                                                  setdiff(x[[s]][[3]], 
                                                                                          union(x[[s]][[2]], 
                                                                                                x[[s]][[1]]))),
                                                                             c(names(x[[s]][1]), 
                                                                               names(x[[s]][2]), 
                                                                               names(x[[s]][3])))
      }
      
      if(i == 4){
        specific_genes[["unique each"]] <- setNames(list(setdiff(x[[s]][[1]],
                                                                 purrr::reduce(list(x[[s]][[2]], 
                                                                                    x[[s]][[3]], 
                                                                                    x[[s]][[4]]), 
                                                                               union)),
                                                         setdiff(x[[s]][[2]],
                                                                 purrr::reduce(list(x[[s]][[1]], 
                                                                                    x[[s]][[3]], 
                                                                                    x[[s]][[4]]), 
                                                                               union)),
                                                         setdiff(x[[s]][[3]],
                                                                 purrr::reduce(list(x[[s]][[2]], 
                                                                                    x[[s]][[1]], 
                                                                                    x[[s]][[4]]), 
                                                                               union)),
                                                         setdiff(x[[s]][[4]],
                                                                 purrr::reduce(list(x[[s]][[2]], 
                                                                                    x[[s]][[3]], 
                                                                                    x[[s]][[1]]), 
                                                                               union))),
                                                    c(names(x[[s]][1]),
                                                      names(x[[s]][2]),
                                                      names(x[[s]][3]),
                                                      names(x[[s]][4])))
      }
    }
  }
  
  # Plot
  p <- ggvenn(
    sets,
    fill_color = colpalette,
    stroke_size = .8,
    show_elements = F,
    stroke_linetype = "solid",
    set_name_color = colpalette,
    set_name_size = set.name.size,
    text_color = "black",
    text_size = text.size,
    padding = Padding, 
    show_stats = "c", 
    fill_alpha = 0.4
  )
  
  return(list(plot = p, intersections = common_genes, group_specific = specific_genes))
}

# --- GOobject to Net data function
Go_to_Net.data <- function(go.object,
                           dea.results,
                           condition.name = NULL,
                           go.term.colname = "Description",
                           go.gene.colname = "geneID",
                           go.pval.colname = "p.adjust",
                           dea.gene.colname = "Gene_name",
                           dea.logFC.colname = "log2FoldChange",
                           dea.pval.colname = "padj"){
  
  ## go.object :      provide a GOobject 
  ## dea.results :    provide a DEA result dataframe
  
  if(!(class(go.object) == "enrichResult")){
    stop("You need to provide a GO object")
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


# --- Personalized network plots
Net_plot <- function(object,
                     select_terms = NULL,
                     net.layout = "fr",
                     edge.alpha = 0.3,
                     edge.colors = NULL,
                     edge.width = .7,
                     node.term.size1 = 10,
                     node.term.size2 = 7,
                     node.term.color = "black",
                     node.gene.size = 4,
                     node.gene.colors = NULL,
                     node.gene.alpha = .8,
                     label.term.show = F,
                     label.term.size = 3,
                     label.term.color = "white",
                     label.term.fill = "purple",
                     label.term.padding = 0.2,
                     label.gene.size = 3){
  
  ## object:        A dataframe containing at least a "Term" and "Gene" column.
  ##                In addition it can have a "Cluster" column for the experimental condition
  ##                and a "logFC" column as well. As an example, the df could be the result 
  ##                of the tidyr::separate_rows(geneID, sep = "/") function.
  ## select_terms:  A vector containing a list of terms you want to filter with your df.
  ##                Default is NULL, so no filtering would be applied.
  
  
  # --------------------- #
  # Checking for packages #
  # --------------------- #
  required_pkgs <- c("igraph", "ggraph", "tidygraph", "dplyr", "tidyr", "ggplot2")
  installed <- rownames(installed.packages())
  for (pkg in required_pkgs) {
    if (!pkg %in% installed) install.packages(pkg)
  }
  lapply(required_pkgs, library, character.only = TRUE)
  
  
  # ---------------- #
  # Debuging section #
  # ---------------- #
  if(!(is.data.frame(object))){
    stop("object must be a dataframe !!")
  }
  
  if(!all(c("Term","Gene") %in% colnames(object))){
    stop("object must have at least 'Term' and 'Gene' column !!")
  }
  
  if(is.null(select_terms)){
    selected_terms <- unique(object[["Term"]])
  } else {
    selected_terms <- select_terms
  }
  
  advanced <- ifelse(all(c("Cluster","logFC") %in% colnames(object)), T, F)
  
  
  
  # ----------------- #
  # The core function #
  # ----------------- #
  
  ## 1. Initiate a df with the selected terms :
  ########################################### #
  dat <- object[object$Term %in% selected_terms,]
  
  ## 2. In case we have the Cluster and logFC column :
  ################################################## #
  if(isTRUE(advanced)){
    d1 <- dat %>% 
      dplyr::distinct(Gene, Cluster) %>%
      dplyr::group_by(Gene) %>%
      dplyr::summarise(group = ifelse(n() > 1, 
                                      "Shared", 
                                      paste0(Cluster[1])), 
                       .groups = "drop")
    
    d2 <- dat %>% 
      dplyr::select(Gene, Cluster, logFC) %>%
      dplyr::group_by(Gene) %>%
      dplyr::summarise(logFC = if (n() == 1) {
        logFC[1]
      } else if (length(unique(sign(logFC))) == 1) {
        mean(logFC) 
      } else {
        NA  
      },
      .groups = "drop")
    
    Genes_info <- dplyr::left_join(d1, d2, by = "Gene")
    
    ### 2.1. Create the graph :
    ######################### #
    graph_data <- igraph::graph_from_data_frame(dat[, c("Gene", "Term")], directed = FALSE)
    
    ### 2.2. Add node type (Gene or Term) :
    ##################################### #
    V(graph_data)$type <- ifelse(V(graph_data)$name %in% dat$Gene, "Gene", "Term")
    
    ### 2.3. Add the attributes (group & logFC) :
    ########################################### #
    V(graph_data)$group <- Genes_info$group[match(V(graph_data)$name, Genes_info$Gene)]
    V(graph_data)$logFC <- Genes_info$logFC[match(V(graph_data)$name, Genes_info$Gene)]
    V(graph_data)$reg <- case_when(V(graph_data)$logFC > 0 ~ "Gene up",
                                   V(graph_data)$logFC < 0 ~ "Gene down",
                                   TRUE ~ NA)
    V(graph_data)$pathways <- case_when(V(graph_data)$type == "Term" ~ V(graph_data)$name,
                                        TRUE ~ NA)
    
    ### 2.4. Add Term names to the edges :
    #################################### #
    E(graph_data)$term <- dat$Term[match(
      paste0(dat$Gene, dat$Term),
      paste0(ends(graph_data, es = E(graph_data))[, 1],
             ends(graph_data, es = E(graph_data))[, 2])
    )]
    
  } else {
    
    ## 3. In case of a simple network :
    ################################# #
    
    graph_data <- igraph::graph_from_data_frame(dat[, c("Gene", "Term")], directed = FALSE)
    V(graph_data)$type <- ifelse(V(graph_data)$name %in% dat$Gene, "Gene", "Term")
    V(graph_data)$group <- ifelse(V(graph_data)$type == "Gene", "Gene", "Term")
    V(graph_data)$logFC <- NA
    
    E(graph_data)$term <- dat$Term[match(
      paste0(dat$Gene, dat$Term),
      paste0(ends(graph_data, es = E(graph_data))[, 1],
             ends(graph_data, es = E(graph_data))[, 2])
    )]
  }
  
  
  
  # ---------------------------- #
  # Graph prep and visualization #
  # ---------------------------- #
  
  graph_tbl <- tidygraph::as_tbl_graph(graph_data)
  
  nterm <- length(selected_terms)
  MyPalette <- c("#9933aa","#aa4400","#ff0000","#337722","#002277",
                 "#441144","#005566","#aa0077","#00bbff","#003333","#4422cc","#116611",
                 "#667700","#ddaa00","#ff22ff","#00ff00","#330077","#0000ff")
  
  if(is.null(edge.colors)){
    e.colors <- MyPalette[1:nterm]
  } else {
    e.colors <- edge.colors
  }
  
  if(is.null(node.gene.colors)){
    n.colors <- setdiff(MyPalette, e.colors)
  } else {
    n.colors <- node.gene.colors
  }
  
  g <- ggraph(graph_tbl, layout = net.layout)+
    theme_void()+
    
    ## Edges :
    geom_edge_link(aes(color = term), alpha = edge.alpha, edge_width = edge.width)+
    scale_edge_color_manual(values = e.colors)+
    guides(edge_color = guide_legend(title = "Selected terms",
                                     override.aes = list(edge_width = 1.5),
                                     position = "right",
                                     direction = "vertical"))+
    
    ## Nodes :
    geom_node_point(colour = ifelse(V(graph_tbl)$type == "Term", "black", "khaki"),
                    alpha = ifelse(V(graph_tbl)$type == "Term", .3, 0),
                    size = ifelse(V(graph_tbl)$type == "Term", node.term.size1, 3))+
    geom_node_point(colour = ifelse(V(graph_tbl)$type == "Term", node.term.color, "khaki"),
                    alpha = ifelse(V(graph_tbl)$type == "Term", 1, 0),
                    size = ifelse(V(graph_tbl)$type == "Term", node.term.size2, 3))+
    
    geom_node_point(aes(color = group,
                        shape = reg), 
                    alpha = node.gene.alpha,
                    size = ifelse(V(graph_tbl)$type == "Term", 
                                  node.term.size2, 
                                  node.gene.size))+
    
    scale_color_manual(values = n.colors, na.translate = F)+
    scale_shape_manual(values = c(16,18), na.translate = F)+
    guides(color = guide_legend(title = "Genes",
                                override.aes = list(size = 4)),
           shape = guide_legend(title = "logFC",
                                override.aes = list(size = 4)))+
    
    ## Label genes :
    geom_node_text(aes(label = name,
                       colour = group),
                   repel = T, 
                   show.legend = F,
                   fontface = "bold.italic",
                   size = label.gene.size)+
    
    ## theme :
    theme(legend.title = element_text(size = 13, 
                                      face = "bold", 
                                      color = "darkred",
                                      margin = margin(t=0.3, b=0.1, unit = "in")),
          legend.text = element_text(size = 10,
                                     face = "bold",
                                     color = "black"),
          legend.margin = margin(l=0.2, unit = "in"))
  
  
  
  if(isTRUE(label.term.show)){
    g <- g+
      geom_node_label(aes(label = pathways), 
                      show.legend = F,
                      repel = T,
                      label.padding = label.term.padding,
                      label.size = .5,
                      size = label.term.size,
                      fontface = "bold.italic",
                      colour = label.term.color,
                      fill = label.term.fill)
  }
  
  
  
  return(g)
  
}

Net_plot2 <- function(object,
                      select_terms = NULL,
                      net.layout = "fr",
                      edge.alpha = 0.3,
                      edge.colors = NULL,
                      edge.width = .7,
                      node.term.size1 = 10,
                      node.term.size2 = 7,
                      node.term.color1 = "black",
                      node.term.color2 = "black",
                      node.gene.size = 4,
                      node.gene.alpha = .8,
                      label.term.show = F,
                      label.term.size = 3,
                      label.term.color = "white",
                      label.term.fill = "purple",
                      label.term.padding = 0.2,
                      label.gene.size = 3,
                      label.gene.colors = NULL){
  
  ## object:        A dataframe containing at least a "Term" and "Gene" column.
  ##                In addition it can have a "Cluster" column for the experimental condition
  ##                and a "logFC" column as well. As an example, the df could be the result 
  ##                of the tidyr::separate_rows(geneID, sep = "/") function.
  ## select_terms:  A vector containing a list of terms you want to filter with your df.
  ##                Default is NULL, so no filtering would be applied.
  
  
  # --------------------- #
  # Checking for packages #
  # --------------------- #
  required_pkgs <- c("igraph", "ggraph", "tidygraph", "dplyr", "tidyr", "ggplot2")
  installed <- rownames(installed.packages())
  for (pkg in required_pkgs) {
    if (!pkg %in% installed) install.packages(pkg)
  }
  lapply(required_pkgs, library, character.only = TRUE)
  
  
  # ---------------- #
  # Debuging section #
  # ---------------- #
  if(!(is.data.frame(object))){
    stop("object must be a dataframe !!")
  }
  
  if(!all(c("Term","Gene") %in% colnames(object))){
    stop("object must have at least 'Term' and 'Gene' column !!")
  }
  
  if(is.null(select_terms)){
    selected_terms <- unique(object[["Term"]])
  } else {
    selected_terms <- select_terms
  }
  
  advanced <- ifelse(all(c("Cluster","logFC") %in% colnames(object)), T, F)
  
  
  
  # ----------------- #
  # The core function #
  # ----------------- #
  
  ## 1. Initiate a df with the selected terms :
  ########################################### #
  dat <- object[object$Term %in% selected_terms,]
  
  ## 2. In case we have the Cluster and logFC column :
  ################################################## #
  if(isTRUE(advanced)){
    d1 <- dat %>% 
      dplyr::distinct(Gene, Cluster) %>%
      dplyr::group_by(Gene) %>%
      dplyr::summarise(group = ifelse(n() > 1, 
                                      "Shared", 
                                      paste0(Cluster[1])), 
                       .groups = "drop")
    
    d2 <- dat %>% 
      dplyr::select(Gene, Cluster, logFC) %>%
      dplyr::group_by(Gene) %>%
      dplyr::summarise(logFC = if (n() == 1) {
        logFC[1]
      } else if (length(unique(sign(logFC))) == 1) {
        mean(logFC) 
      } else {
        NA  
      },
      .groups = "drop")
    
    Genes_info <- dplyr::left_join(d1, d2, by = "Gene")
    
    ### 2.1. Create the graph :
    ######################### #
    graph_data <- igraph::graph_from_data_frame(dat[, c("Gene", "Term")], directed = FALSE)
    
    ### 2.2. Add node type (Gene or Term) :
    ##################################### #
    V(graph_data)$type <- ifelse(V(graph_data)$name %in% dat$Gene, "Gene", "Term")
    
    ### 2.3. Add the attributes (group & logFC) :
    ########################################### #
    V(graph_data)$genes <- case_when(V(graph_data)$type == "Gene" ~ V(graph_data)$name,
                                     TRUE ~ NA)
    V(graph_data)$terms <- case_when(V(graph_data)$type == "Term" ~ V(graph_data)$name,
                                     TRUE ~ NA)
    V(graph_data)$group <- Genes_info$group[match(V(graph_data)$name, Genes_info$Gene)]
    V(graph_data)$logFC <- Genes_info$logFC[match(V(graph_data)$name, Genes_info$Gene)]
    V(graph_data)$reg <- case_when(V(graph_data)$logFC > 0 ~ "Gene up",
                                   V(graph_data)$logFC < 0 ~ "Gene down",
                                   TRUE ~ NA)
    
    
    ### 2.4. Add Term names to the edges :
    #################################### #
    E(graph_data)$term <- dat$Term[match(
      paste0(dat$Gene, dat$Term),
      paste0(ends(graph_data, es = E(graph_data))[, 1],
             ends(graph_data, es = E(graph_data))[, 2])
    )]
    
  } else {
    
    ## 3. In case of a simple network :
    ################################# #
    
    graph_data <- igraph::graph_from_data_frame(dat[, c("Gene", "Term")], directed = FALSE)
    V(graph_data)$type <- ifelse(V(graph_data)$name %in% dat$Gene, "Gene", "Term")
    V(graph_data)$group <- ifelse(V(graph_data)$type == "Gene", "Gene", "Term")
    V(graph_data)$logFC <- NA
    
    E(graph_data)$term <- dat$Term[match(
      paste0(dat$Gene, dat$Term),
      paste0(ends(graph_data, es = E(graph_data))[, 1],
             ends(graph_data, es = E(graph_data))[, 2])
    )]
  }
  
  
  
  # ---------------------------- #
  # Graph prep and visualization #
  # ---------------------------- #
  
  graph_tbl <- tidygraph::as_tbl_graph(graph_data)
  
  nterm <- length(selected_terms)
  MyPalette <- c("#9933aa","#aa4400","#ff0000","#337722","#002277",
                 "#441144","#005566","#aa0077","#00bbff","#003333","#4422cc","#116611",
                 "#667700","#ddaa00","#ff22ff","#00ff00","#330077","#0000ff")
  
  if(is.null(edge.colors)){
    e.colors <- MyPalette[1:nterm]
  } else {
    e.colors <- edge.colors
  }
  
  if(is.null(label.gene.colors)){
    l.colors <- setdiff(MyPalette, e.colors)
  } else {
    l.colors <- label.gene.colors
  }
  
  g <- ggraph(graph_tbl, layout = net.layout)+
    theme_void()+
    
    ## Edges :
    geom_edge_link(aes(color = term), alpha = edge.alpha, edge_width = edge.width)+
    scale_edge_color_manual(values = e.colors)+
    guides(edge_color = guide_legend(title = "Selected terms",
                                     override.aes = list(edge_width = 1.5),
                                     position = "right",
                                     direction = "vertical"))+
    
    ## Nodes :
    geom_node_point(colour = ifelse(!(is.na(V(graph_tbl)$terms)), node.term.color1, "white"),
                    alpha = ifelse(!(is.na(V(graph_tbl)$terms)), .3, 0),
                    size = ifelse(!(is.na(V(graph_tbl)$terms)), node.term.size1, 0))+
    geom_node_point(colour = ifelse(!(is.na(V(graph_tbl)$terms)), node.term.color2, "white"),
                    alpha = ifelse(!(is.na(V(graph_tbl)$terms)), 1, 0),
                    size = ifelse(!(is.na(V(graph_tbl)$terms)), node.term.size2, 0))+
    
    geom_node_point(aes(colour = logFC,
                        shape = group), 
                    alpha = node.gene.alpha,
                    size = ifelse(V(graph_tbl)$type == "Term", 
                                  node.term.size2, 
                                  node.gene.size))+
    
    scale_colour_gradientn(
      colours = c("#007700", "#77ff77", "#ffffff", "#ff7777", "#770000"),
      values = scales::rescale(c(-3, -1.5, 0, 1.5, 3.5)),
      limits = c(-3, 3.5)
    )+
    
    scale_shape_manual(values = c(17,19,18),
                       na.translate = F)+
    
    guides(colour = guide_colourbar(title = "log2FC",
                                    barwidth = unit(.6, "cm"),
                                    barheight = unit(3, "cm")),
           shape = guide_legend(title = "Genes",
                                override.aes = list(size = 4)))+
    
    
    
    ## Label genes :
    geom_node_text(aes(label = genes),
                   colour = ifelse(V(graph_tbl)$group == unique(V(graph_tbl)$group)[1],
                                   l.colors[1],
                                   ifelse(V(graph_tbl)$group == unique(V(graph_tbl)$group)[2],
                                          l.colors[2],
                                          l.colors[3])),
                   repel = T, 
                   show.legend = F,
                   fontface = "bold.italic",
                   size = label.gene.size)+
    
    ## theme :
    theme(legend.title = element_text(size = 13, 
                                      face = "bold", 
                                      color = "darkred",
                                      margin = margin(t=0.3, b=0.1, unit = "in")),
          legend.text = element_text(size = 10,
                                     face = "bold",
                                     color = "black"),
          legend.margin = margin(l=0.2, unit = "in"),
          legend.justification = "top")
  
  
  if(isTRUE(label.term.show)){
    g <- g+
      geom_node_label(aes(label = terms), 
                      show.legend = F,
                      repel = T,
                      label.padding = label.term.padding,
                      label.size = .5,
                      size = label.term.size,
                      fontface = "bold.italic",
                      colour = label.term.color,
                      fill = label.term.fill)
  }
  
  
  
  return(g)
  
}

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




#' @title Custom Cnetplot-like Network Visualization
#' 
#' @description
#' Create an interactive-style gene-term network plot based on a dataframe containing Gene Ontology terms and associated genes.
#' Supports custom coloring by condition, logFC, and shared genes.
#'
#' @param object A dataframe containing at least the columns `"Term"` and `"Gene"`. It can also contain `"Cluster"` and `"logFC"` for enhanced plotting.
#' @param select_terms A character vector of selected GO terms to display. If NULL (default), all terms are included.
#' @param net.layout Layout method used for plotting (default `"fr"` – Fruchterman-Reingold).
#' @param edge.alpha Opacity of edges (default 0.3).
#' @param edge.colors Manual color vector for edges (one color per GO term). If NULL, uses default palette.
#' @param edge.width Width of the edges (default 0.7).
#' @param node.term.size1 Outer node size for terms (default 10).
#' @param node.term.size2 Inner node size for terms (default 7).
#' @param node.term.color Outer color for term nodes (default `"black"`).
#' @param node.gene.size Size of the gene nodes (default 4).
#' @param node.gene.colors Named color vector for gene node coloring (groups like `"siC8"`, `"Shared"`...). If NULL, uses default palette.
#' @param node.gene.alpha Transparency of gene nodes (default 0.8).
#' @param label.term.show Logical. Whether to label GO terms (default `FALSE`).
#' @param label.term.size Size of the GO term labels (default 3).
#' @param label.term.color Color of term label text (default `"white"`).
#' @param label.term.fill Background fill of term labels (default `"purple"`).
#' @param label.term.padding Padding around the term label text (default 0.2).
#' @param label.gene.size Size of the gene labels (default 3).
#'
#' @return A ggplot2 object with the gene-term network plotted via ggraph.
#'
#' @examples
#' \dontrun{
#' Net_plot2(df_go_genes)
#' Net_plot2(df_go_genes, select_terms = c("immune response", "T cell activation"))
#' }
#'
#' @export



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
  
  # ---------------- #
  # Debuging section #
  # ---------------- #
  if(!(is.data.frame(object))) stop("object must be a dataframe !!")
  
  if(!all(c("Term","Gene") %in% colnames(object))){
    stop("object must have at least 'Term' and 'Gene' column !!")
  }
  
  selected_terms <- if(is.null(select_terms)) unique(object[["Term"]]) else select_terms
  
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
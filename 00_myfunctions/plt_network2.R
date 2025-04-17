#' Custom Cnetplot-like Network Visualization
#'
#' Builds a gene-pathway network from a dataframe containing pathway terms and gene IDs, 
#' optionally enhanced with clusters and logFC values. Allows custom styling for edge colors, 
#' node shapes, sizes, labels, and color scales.
#'
#' @param object Dataframe with at least "Term" and "Gene" columns. Optionally includes "Cluster" and "logFC".
#' @param select_terms Character vector of pathway terms to filter and display. If NULL, all terms are used.
#' @param net.layout Graph layout (default: "fr"). Other options: "kk", "circle", etc.
#' @param edge.alpha Transparency for edges.
#' @param edge.colors Color palette for pathway terms.
#' @param edge.width Width of the edge links.
#' @param node.term.size1 Outer glow size for term nodes.
#' @param node.term.size2 Inner point size for term nodes.
#' @param node.term.color1 Outer color for term nodes.
#' @param node.term.color2 Inner color for term nodes.
#' @param node.gene.size Size of gene nodes.
#' @param node.gene.alpha Alpha for gene nodes.
#' @param label.term.show Whether to show labels for the GO terms.
#' @param label.term.size Font size for term labels.
#' @param label.term.color Color of the text labels for terms.
#' @param label.term.fill Fill color for the term label boxes.
#' @param label.term.padding Padding for term labels.
#' @param label.gene.size Size for gene labels.
#' @param label.gene.colors Color vector for each gene group (Shared, siC8, etc.)
#'
#' @return A `ggplot2` object of the network graph.
#' @export
#'
#' @examples
#' Net_plot(df, select_terms = c("immune response", "cell cycle"))


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
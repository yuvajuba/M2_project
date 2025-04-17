#' My_theme: Custom ggplot2 theme for RNA-seq visualizations
#'
#' A clean and publication-ready theme based on `theme_minimal()`, with
#' adjustments to title, axes, legend, background, and panel styles.
#'
#' @return A ggplot2 theme object
#' @export

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
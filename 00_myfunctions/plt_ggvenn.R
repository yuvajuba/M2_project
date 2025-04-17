#' Display_Venn
#'
#' Generate a Venn diagram (up to 4 sets) and return intersections and specific sets.
#'
#' @param Markers A named list of vectors, each representing a group of elements (e.g., genes).
#' @param colpalette A vector of colors to use for the different sets (length 2–4). Default palette available.
#' @param set.names A custom character vector to rename the sets. Must match the number of sets (optional).
#' @param set.name.size Size of the set labels in the plot. Default is 5.
#' @param text.size Size of the intersection count labels. Default is 5.
#' @param Padding Padding between elements and the plot border. Default is 0.04.
#' @param show_elements Logical, whether to display set elements in the plot. Default is FALSE.
#'
#' @return A list with:
#' \item{plot}{The Venn diagram (ggplot object).}
#' \item{intersections}{A list of intersected elements between 2 or more sets.}
#' \item{group_specific}{A list of elements unique to each group.}
#'
#' @examples
#' my_list <- list(A = c("g1","g2","g3"), B = c("g2","g4"), C = c("g3","g5"))
#' Display_Venn(my_list)
#'
#' @export

Display_Venn <- function(Markers, 
                         colpalette = NULL, 
                         set.names = NULL, 
                         set.name.size = 5,
                         text.size = 5,
                         Padding = 0.04,
                         show_elements = FALSE) {
  
  
  
  # 1. Input checks :
  # --------------- :
  if (!is.list(Markers)) stop("Argument 'Markers' must be a list of vectors.")
  if (!all(sapply(Markers, is.vector))) stop("All elements in the 'Markers' list must be vectors.")
  
  
  n <- min(length(Markers), 4)
  if (n < 2) stop("You need at least 2 sets to make a Venn diagram.")
  if (length(Markers) > 4) warning("Only the first 4 sets will be used.")
  
  
  # 2. Handling sets :
  # ---------------- :
  sets <- Markers[1:n]
  if (!(is.null(set.names))) {
    if (length(set.names) != n) {
      stop("The length of 'set.names' must match the number of sets in 'Markers'")
    }
    names(sets) <- set.names
  } else if (is.null(names(sets))) {
    names(sets) <- paste0("Set", 1:n)
  }
  
  
  # 3. Handling colours :
  # ------------------- :
  default_palette <- c("navy", "red", "darkgreen", "violet")
  colpalette <- if (is.null(colpalette) || length(colpalette) < n) {
    warning("Invalid color palette. Using default.")
    default_palette[1:n]
  } else {
    colpalette[1:n]
  }
  
  
  
  # 4. Calculate intersections :
  # -------------------------- :
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
  
  
  # 5. Calculate group specific genes :
  # --------------------------------- :
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
  p <- ggvenn::ggvenn(
    sets,
    fill_color = colpalette,
    stroke_size = 0.8,
    set_name_color = colpalette,
    set_name_size = set.name.size,
    text_size = text.size,
    text_color = "black",
    stroke_linetype = "solid",
    padding = Padding,
    fill_alpha = 0.4,
    show_stats = "c",
    show_elements = show_elements
  )
  
  return(list(plot = p, intersections = common_genes, group_specific = specific_genes))
}
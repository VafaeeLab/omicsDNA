
#' Plot Filtered ABACUS Community Network
#'
#' This function performs community detection using the ABACUS algorithm and plots
#' a selected set of layers from a multilayer network. It supports both "multiforce"
#' and "circular" layouts and automatically adjusts the grid layout based on the
#' number of layers selected.
#'
#' @param net A multilayer network object of class `ml.network`.
#' @param layer_names A character vector of layer names to include in the plot (e.g., `c("E1", "E2", "M1")`).
#' @param min.actors Minimum number of actors per community. Default: 15.
#' @param min.layers Minimum number of layers a community must span. Default: 2.
#' @param layout The layout to use for plotting. One of `"multiforce"` (default) or `"circular"`.
#'
#' @return A filtered data frame of communities returned by ABACUS.
#'
#' @examples
#' \dontrun{
#' plot_filtered_ABACUS_network(net, layer_names = c("E1", "E2", "M1"))
#' }
#'
#' @importFrom multinet abacus_ml layout_multiforce_ml layout_circular_ml layers_ml
#' @importFrom dplyr filter
#' @export
plot_filtered_ABACUS_network <- function(net,
                                         layer_names,
                                         min.actors = 15,
                                         min.layers = 2,
                                         layout = "multiforce") {

  layout_multiforce <- multinet::layout_multiforce_ml(net, gravity = 0.3)
  layout_circular <- multinet::layout_circular_ml(net)

  selected_layout <- switch(layout,
                            "multiforce" = layout_multiforce,
                            "circular" = layout_circular,
                            stop("Error: Invalid layout type. Choose 'multiforce' or 'circular'."))

  message("Note: Running ABACUS community detection...")
  ABACUS_communities_net <- multinet::abacus_ml(net, min.actors = min.actors, min.layers = min.layers)

  flt_ABACUS <- dplyr::filter(ABACUS_communities_net, layer %in% layer_names)

  num_layers <- length(layer_names)
  grid_rows <- floor(sqrt(num_layers))
  grid_cols <- ceiling(num_layers / grid_rows)

  message(sprintf("Note: Plotting layout grid: %dx%d", grid_rows, grid_cols))

  plot(net, layout = selected_layout,
       grid = c(grid_rows, grid_cols),
       layers = layer_names,
       com = flt_ABACUS,
       vertex.cex = 1.5, vertex.size = 5,
       vertex.labels = NA)

  return(flt_ABACUS)
}

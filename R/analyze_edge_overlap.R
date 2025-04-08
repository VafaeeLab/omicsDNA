
#' Analyze Edge Overlap Across Network Layers
#'
#' This function calculates the Jaccard similarity index for gene-gene relationships (edges)
#' between all pairs of layers in a multilayer network and visualizes the overlap as a heatmap.
#'
#' @param net A multilayer network object of class `ml.network`, typically
#'        created using create_multilayer_network().
#'
#' @return A symmetric matrix containing Jaccard similarity scores for edge overlap
#'         between each pair of layers.
#'
#' @details
#' - Computes Jaccard index for edges using `layer_comparison_ml()` with `method = "jaccard.edges"`.
#' - Converts the matrix to long format and creates a heatmap.
#'
#' @examples
#' \dontrun{
#' edge_jaccard_matrix <- analyze_edge_overlap(net)
#' }
#'
#' @importFrom multinet layer_comparison_ml layers_ml
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_gradient theme_minimal
#'             labs ggtitle element_text
#' @export
analyze_edge_overlap <- function(net) {
  message("Note: Computing Jaccard similarity for edges (gene-gene relationships) between layers...")

  comparison_matrix <- as.matrix(
    multinet::layer_comparison_ml(net,
                                  layers = multinet::layers_ml(net),
                                  method = "jaccard.edges",
                                  mode = "all", K = 0)
  )

  df <- as.data.frame(as.table(comparison_matrix))

  plot <- ggplot2::ggplot(df, ggplot2::aes(Var1, Var2, fill = Freq)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::scale_fill_gradient(low = "white", high = "blue") +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = "Layers", y = "Layers", fill = "Jaccard Index") +
    ggplot2::ggtitle("Edge Overlap Between Layers - Jaccard") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))

  print(plot)
  return(comparison_matrix)
}

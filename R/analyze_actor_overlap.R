
#' Analyze Actor (Gene) Overlap Across Network Layers
#'
#' This function calculates the Jaccard similarity index for actors (genes)
#' between all pairs of layers in a multilayer network using the
#' `multinet` package. The results are visualized as a heatmap.
#'
#' @param net A multilayer network object of class `ml.network`, typically
#'        created using create_multilayer_network().
#'
#' @return A symmetric matrix containing Jaccard similarity scores for actor
#'         (gene) overlap between each pair of layers.
#'
#' @details
#' - Computes pairwise Jaccard indices of shared actors across layers.
#' - Converts the similarity matrix to long format.
#' - Plots the result as a heatmap using `ggplot2`.
#'
#' @examples
#' \dontrun{
#' actor_jaccard_matrix <- analyze_actor_overlap(net)
#' }
#'
#' @importFrom multinet layer_comparison_ml layers_ml
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_gradient theme_minimal
#'             labs ggtitle element_text
#' @export
analyze_actor_overlap <- function(net) {
  message("Note: Computing Jaccard similarity for actors (genes) between layers...")

  comparison_matrix <- as.matrix(
    multinet::layer_comparison_ml(net,
                                  layers = multinet::layers_ml(net),
                                  method = "jaccard.actors",
                                  mode = "all", K = 0)
  )

  df <- as.data.frame(as.table(comparison_matrix))

  plot <- ggplot2::ggplot(df, ggplot2::aes(Var1, Var2, fill = Freq)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::scale_fill_gradient(low = "white", high = "blue") +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = "Layers", y = "Layers", fill = "Jaccard Index") +
    ggplot2::ggtitle("Actors (Genes) Overlap Between Layers - Jaccard") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))

  print(plot)
  return(comparison_matrix)
}

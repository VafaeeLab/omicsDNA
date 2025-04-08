
#' Analyze Community Counts vs. Sample Sizes Across Layers
#'
#' This function analyzes the relationship between the number of detected
#' communities (CIDs) in each network layer and the number of samples
#' per group. It visualizes this relationship using a scatter plot with
#' a linear regression line.
#'
#' @param ABACUS_communities_net A data frame returned from `abacus_ml()`
#'        containing community detection results, with at least the columns
#'        `layer` and `cid`.
#' @param num_samples A named numeric vector where names are layer identifiers
#'        (e.g., "E1", "M1") and values represent the number of samples per group.
#'
#' @return A named vector where each element is the number of unique community
#'         IDs (CIDs) detected in the corresponding layer.
#'
#' @examples
#' \dontrun{
#' num_samples <- c(E1 = 46, E2 = 22, M1 = 16, M2 = 27)
#' analyze_layer_communities(ABACUS_communities_net, num_samples)
#' }
#'
#' @importFrom dplyr group_by summarise
#' @importFrom ggplot2 ggplot aes geom_point geom_smooth labs theme_minimal
#' @export
analyze_layer_communities <- function(ABACUS_communities_net, num_samples) {
  layer_cid_counts <- ABACUS_communities_net %>%
    dplyr::group_by(layer) %>%
    dplyr::summarise(unique_cid_count = dplyr::n_distinct(cid), .groups = "drop")

  message("Note: Number of distinct communities detected in each layer:")
  print(layer_cid_counts)

  layer_cid_vector <- setNames(layer_cid_counts$unique_cid_count, layer_cid_counts$layer)

  if (!is.null(names(num_samples)) && identical(sort(names(num_samples)), sort(names(layer_cid_vector)))) {
    data <- data.frame(
      cids = layer_cid_vector,
      samples = num_samples[names(layer_cid_vector)]
    )

    print(
      ggplot2::ggplot(data, ggplot2::aes(x = cids, y = samples)) +
        ggplot2::geom_point(size = 3, color = "blue") +
        ggplot2::geom_smooth(method = "lm", se = FALSE, color = "red") +
        ggplot2::labs(x = "Number of CIDs",
                      y = "Number of Samples per Group",
                      title = "Relationship between CIDs and Samples per Group") +
        ggplot2::theme_minimal()
    )

    return(layer_cid_vector)
  } else {
    stop("Error: Layer names in `num_samples` do not match the ones in ABACUS result.")
  }
}

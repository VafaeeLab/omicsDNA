
#' Select Common Edges Within Each Group Based on Repetition Threshold
#'
#' Identifies and retains edges that appear consistently across multiple repetitions
#' within each group and computes their average weights.
#'
#' @param edge_lists A nested list of edge lists for each group (e.g., age group).
#'        The structure should be: `list[[group]][[repetition]]`, where each element
#'        is a data frame with columns `from`, `to`, and `weight`.
#' @param threshold_ratio A numeric value (between 0 and 1) indicating the minimum
#'        proportion of repetitions in which an edge must appear to be retained.
#'        Default is `0.7`, meaning an edge must appear in at least 70% of repetitions.
#'
#' @return A named list (`selected_edges`) of data frames. Each data frame contains
#'         the filtered and averaged edges for a specific group, with columns:
#'         `from`, `to`, and `weight` (mean weight across selected repetitions).
#'
#' @details This function is typically used after constructing multiple networks per group
#' (e.g., via bootstrapping or subsampling) and allows filtering for reproducible interactions.
#'
#' @examples
#' \dontrun{
#' selected_edges <- select_edges(edge_lists, threshold_ratio = 0.7)
#' head(selected_edges[["Group1"]])
#' }
#'
#' @importFrom stats aggregate
#' @export
select_edges <- function(edge_lists, threshold_ratio = 0.7) {
  num_repetitions <- length(edge_lists[[1]])
  threshold <- num_repetitions * threshold_ratio

  selected_edges <- vector("list", length(names(edge_lists)))
  names(selected_edges) <- names(edge_lists)

  for (group in names(selected_edges)) {
    edges <- do.call("rbind", edge_lists[[group]])

    edge_counts <- table(paste(edges$from, edges$to, sep = "_"))
    repeated_edges <- edge_counts[edge_counts >= threshold]

    filtered_edges <- edges[paste(edges$from, edges$to, sep = "_") %in% names(repeated_edges), ]
    selected_edges[[group]] <- aggregate(weight ~ from + to, data = filtered_edges, FUN = mean)
  }

  return(selected_edges)
}

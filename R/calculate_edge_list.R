
#' Convert Weighted Adjacency Matrix to Edge List with Preserved Signs
#'
#' This function transforms a weighted adjacency matrix into an edge list
#' format, preserving the sign (positive or negative) of the original correlations.
#'
#' @param weighted_adj_mat A square numeric matrix representing the weighted adjacency matrix
#'                         of a graph. Typically generated from a correlation matrix with
#'                         possible negative values.
#'
#' @return A data frame with three columns:
#' \describe{
#'   \item{from}{The source node (gene name).}
#'   \item{to}{The target node (gene name).}
#'   \item{weight}{The edge weight (correlation), preserving the original sign.}
#' }
#'
#' @details The function first creates an `igraph` object using the absolute values of the
#' weights, then restores the original sign based on the original matrix. This allows
#' correct representation of both positive and negative correlations in downstream analysis.
#'
#' @examples
#' adj_matrix <- matrix(c(0, 0.8, -0.5,
#'                        0.8, 0, 0.6,
#'                       -0.5, 0.6, 0),
#'                      nrow = 3, byrow = TRUE)
#' rownames(adj_matrix) <- colnames(adj_matrix) <- c("Gene1", "Gene2", "Gene3")
#' edge_list <- calculate_edge_list(adj_matrix)
#'
#' @importFrom igraph graph_from_adjacency_matrix as_edgelist E
#' @export
calculate_edge_list <- function(weighted_adj_mat) {
  negative_edges <- weighted_adj_mat < 0
  adjusted_matrix <- abs(weighted_adj_mat)

  graph <- igraph::graph_from_adjacency_matrix(adjusted_matrix, mode = "undirected", weighted = TRUE)
  edges <- igraph::as_edgelist(graph)
  weights <- igraph::E(graph)$weight

  for (i in seq_len(nrow(edges))) {
    node1 <- edges[i, 1]
    node2 <- edges[i, 2]

    if (negative_edges[node1, node2]) {
      weights[i] <- -weights[i]
    }
  }

  edgelist_df <- data.frame(from = edges[, 1], to = edges[, 2], weight = weights)
  return(edgelist_df)
}

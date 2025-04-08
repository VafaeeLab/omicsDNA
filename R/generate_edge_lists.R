
#' Generate Edge Lists from Adjacency Matrices Across Repetitions
#'
#' Converts a nested list of adjacency matrices (grouped by repetition and group)
#' into corresponding edge lists using a specified edge list generation function.
#'
#' @param adjacency_matrices A nested list where each element represents a repetition,
#'        and within each repetition, a list of weighted adjacency matrices for each group.
#'        The structure should be: `list[[repetition]][[Group]]`.
#' @param num_repetitions An integer indicating how many repetitions were used to generate
#'        the adjacency matrices. Default is 50.
#' @param func_edgelist A function that takes a weighted adjacency matrix and returns an
#'        edge list data frame. Default is `calculate_edge_list`.
#'
#' @return A nested list of edge lists. The outer list is indexed by group, and each
#'         element is a list of edge lists for each repetition. The structure is:
#'         `list[[Group]][[repetition]]`.
#'
#' @examples
#' \dontrun{
#' # Assuming adjacency_matrices has already been generated:
#' edge_lists <- generate_edge_lists(adjacency_matrices)
#' head(edge_lists[["Group1"]][[1]])
#' }
#'
#' @export
generate_edge_lists <- function(adjacency_matrices, num_repetitions = 50, func_edgelist = calculate_edge_list) {

  edge_lists <- list()
  Groups <- names(adjacency_matrices[[1]])

  for (i in 1:num_repetitions) {
    for (Group in Groups) {

      adjacency_matrix <- adjacency_matrices[[i]][[Group]]
      edge_list <- func_edgelist(adjacency_matrix)
      edge_list <- data.frame(edge_list, Group = Group)

      if (!(Group %in% names(edge_lists))) {
        edge_lists[[Group]] <- list()
      }

      edge_lists[[Group]][[i]] <- edge_list
    }
  }

  return(edge_lists)
}



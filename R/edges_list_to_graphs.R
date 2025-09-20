
# ------------------------------------------------------------------------------
# 6 - Convert a list of edge data frames to a list of igraphs
# ------------------------------------------------------------------------------

#' Convert per‑layer edge tables to per‑layer `igraph` graphs
#'
#' @description
#' Takes a **named list** of edge data frames (e.g., the output of
#' `consensusEdges(..., as_list = TRUE)`) and returns a **named list** of
#' `igraph` graphs—one per layer. Column names must include `from` and `to`;
#' `weight` is optional (default `1`).
#'
#' @param edges_list Named list where each element is a data frame with columns
#'   `from`, `to`, and optional `weight`.
#' @param directed Logical; build directed graphs? Default `FALSE`.
#' @return Named list of `igraph` objects.
#'
#' @examples
#' \dontrun{
#' list_graphs <- edges_list_to_graphs(cons_list, directed = FALSE)
#' }
#' @importFrom igraph graph_from_data_frame
#' @export
edges_list_to_graphs <- function(edges_list, directed = FALSE) {
  stopifnot(is.list(edges_list), length(edges_list) >= 1)
  if (is.null(names(edges_list)) || any(!nzchar(names(edges_list)))) {
    names(edges_list) <- paste0("layer", seq_along(edges_list))
  }

  lapply(edges_list, function(df) {
    if (!is.data.frame(df) || !all(c("from","to") %in% names(df))) {
      stop("Each element must be a data frame containing columns `from` and `to`.")
    }
    if (!("weight" %in% names(df))) df$weight <- 1
    df <- data.frame(
      from   = as.character(df$from),
      to     = as.character(df$to),
      weight = suppressWarnings(as.numeric(df$weight)),
      stringsAsFactors = FALSE
    )
    if (all(is.na(df$weight))) df$weight <- 1
    igraph::graph_from_data_frame(df[, c("from","to","weight")], directed = directed)
  })
}


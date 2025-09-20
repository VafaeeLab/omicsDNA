
# ------------------------------------------------------------------------------
# 7 - Convert a list of igraphs to a list of edge data frames
# ------------------------------------------------------------------------------

#' Convert per‑layer `igraph` graphs to per‑layer edge tables
#'
#' @description
#' Takes a **named list** of `igraph` graphs (one per layer) and returns a
#' **named list** of data frames with columns `from`, `to`, and `weight`. This is
#' the inverse of `edges_list_to_graphs()`.
#'
#' - If a graph has no `weight` edge attribute, weights are set to `1`.
#' - If vertex `name` is missing, numeric vertex IDs are assigned and coerced to
#'   character.
#' - For **undirected** graphs, you can optionally canonicalize endpoints
#'   (`from = min(u,v)`, `to = max(u,v)`) to make the orientation deterministic.
#'
#' @param graphs_list Named list of `igraph` objects (one per layer).
#' @param weight_attr Edge attribute to map into the `weight` column. Default `"weight"`.
#' @param default_weight Numeric value used when `weight_attr` is missing or `NA`. Default `1`.
#' @param canonicalize_undirected Logical; if `TRUE`, canonicalize endpoints for undirected graphs.
#'   Default `TRUE`.
#' @param keep_edge_attrs Optional character vector of additional edge attributes to include
#'   in the output tables (columns are appended after `weight`). Attributes not present in a
#'   layer are silently skipped (set `verbose = TRUE` to be notified).
#' @param verbose Logical; print helpful messages (e.g., when vertex names are created or
#'   requested extra attributes are missing). Default `TRUE`.
#'
#' @return A **named list** of data frames (one per layer) with columns
#'   `from`, `to`, `weight`, and any requested extras.
#'
#' @examples
#' # Two small undirected layer graphs
#' g1 <- igraph::graph_from_data_frame(
#'   data.frame(from=c("g1","g2","g3"), to=c("g2","g3","g4"), weight=c(0.8,0.6,0.7)),
#'   directed = FALSE
#' )
#' g2 <- igraph::graph_from_data_frame(
#'   data.frame(from=c("g1","g2"), to=c("g3","g4"), weight=c(0.9,0.65)),
#'   directed = FALSE
#' )
#'
#' edges_list <- graphs_to_edges_list(list(E1 = g1, E2 = g2))
#' str(edges_list$E1)
#'
#' # Round‑trip (for undirected graphs, use directed = FALSE)
#' list_graphs2 <- edges_list_to_graphs(edges_list, directed = FALSE)
#'
#' @importFrom igraph as_data_frame
#' @export
graphs_to_edges_list <- function(
    graphs_list,
    weight_attr             = "weight",
    default_weight          = 1,
    canonicalize_undirected = TRUE,
    keep_edge_attrs         = NULL,
    verbose                 = TRUE
) {
  stopifnot(is.list(graphs_list), length(graphs_list) >= 1)
  if (is.null(names(graphs_list)) || any(!nzchar(names(graphs_list)))) {
    names(graphs_list) <- paste0("layer", seq_along(graphs_list))
  }

  out <- setNames(vector("list", length(graphs_list)), names(graphs_list))

  for (ln in names(graphs_list)) {
    g <- graphs_list[[ln]]
    if (!inherits(g, "igraph")) stop("Element '", ln, "' is not an igraph object.")

    # Ensure vertex names exist
    if (is.null(igraph::V(g)$name)) {
      igraph::V(g)$name <- as.character(seq_len(igraph::vcount(g)))
      if (isTRUE(verbose)) message("Layer ", ln, ": vertex `name` was missing; numeric IDs assigned.")
    }

    # Empty graph → empty edge table
    if (igraph::ecount(g) == 0L) {
      df <- data.frame(from = character(),
                       to   = character(),
                       weight = numeric(),
                       stringsAsFactors = FALSE)
      # (Optionally add requested extra columns as 0‑length vectors)
      if (length(keep_edge_attrs)) {
        for (att in keep_edge_attrs) df[[att]] <- df$from # creates character(0)
      }
      out[[ln]] <- df
      next
    }

    # Pull edges (uses vertex names if present)
    edf <- igraph::as_data_frame(g, what = "edges")
    # Ensure character endpoints
    edf$from <- as.character(edf$from)
    edf$to   <- as.character(edf$to)

    # Prepare weights
    if (!(weight_attr %in% names(edf))) {
      wt <- rep(default_weight, nrow(edf))
    } else {
      wt <- suppressWarnings(as.numeric(edf[[weight_attr]]))
      wt[is.na(wt)] <- default_weight
    }

    # Base output
    df <- data.frame(from = edf$from, to = edf$to, weight = wt, stringsAsFactors = FALSE)

    # Canonicalize for undirected graphs (deterministic orientation)
    if (!igraph::is_directed(g) && isTRUE(canonicalize_undirected) && nrow(df) > 0L) {
      a <- pmin(df$from, df$to); b <- pmax(df$from, df$to)
      df$from <- a; df$to <- b
    }

    # Append any requested extra edge attributes present in this layer
    if (length(keep_edge_attrs)) {
      extras <- intersect(keep_edge_attrs, setdiff(names(edf), c("from", "to", weight_attr)))
      for (att in extras) df[[att]] <- edf[[att]]
      missing <- setdiff(keep_edge_attrs, extras)
      if (isTRUE(verbose) && length(missing)) {
        message("Layer ", ln, ": requested edge attribute(s) not found: ",
                paste(missing, collapse = ", "))
      }
    }

    out[[ln]] <- df
  }

  out
}

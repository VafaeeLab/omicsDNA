

# ------------------------------------------------------------------------------
# 5 - Build a multilayer network from per-layer igraph objects (aligned with build_multiNet)
# ------------------------------------------------------------------------------

#' Build a multilayer network from per‑layer **igraph** graphs
#'
#' @description
#' Construct a multilayer network—**one layer per graph**—from a named list of
#' `igraph` objects and return a `multinet::ml.network`. This sibling of
#' `build_multiNet()` is tuned so that, given equivalent inputs, the resulting
#' multilayer object is the **same** as what `build_multiNet()` would produce.
#'
#' **Alignment with `build_multiNet()`**
#' - **Layer ordering:** if `layerOrder = NULL`, layers are ordered **alphabetically**.
#' - **Vertex attributes:** when `nodesMetadata` is supplied, attributes from
#'   metadata **overwrite** existing graph vertex attributes of the same name
#'   (default behavior in `build_multiNet()`).
#' - **No automatic simplification or actor normalization** by default.
#'
#' **Input expectations**
#' - `graphsPerLayer` is a **named list** where each element is an `igraph`.
#'   If names are missing, synthetic names (`"layer1"`, `"layer2"`, …) are added.
#' - Each graph should have vertex attribute `name` for actor IDs; if missing,
#'   numeric IDs are assigned and coerced to character.
#' - If an edge attribute `weight` is absent, it is added and set to `1`.
#'
#' **Persistence**
#' When `save_to_rds = TRUE`, the multilayer network is saved under
#' `results_dir` (default: `getOption("mlnet.results_dir","omicsDNA_results")`)
#' with an auto‑stamped filename unless `rds_file` is supplied.
#'
#' @param graphsPerLayer Named list of `igraph` objects (one per layer).
#' @param layerOrder Optional character vector specifying which layers to keep
#'   and in what order. If `NULL`, layers are **alphabetically** ordered.
#' @param nodesMetadata Optional data frame of node (actor) attributes.
#' @param featureID_col Column in `nodesMetadata` holding actor IDs to match
#'   against vertex `name`. Default `"feature_id"`.
#' @param nodeAttrCols Character vector of columns from `nodesMetadata` to attach
#'   as vertex attributes. Default `NULL` = attach all except `featureID_col`.
#' @param attr_conflict Conflict policy when a vertex attribute exists both in a
#'   graph and in `nodesMetadata`. One of `"overwrite_with_meta"` (default),
#'   `"keep_graph"`, or `"coalesce"` (fill graph NAs with metadata).
#' @param actor_normalize Optional character vector of normalization steps
#'   applied to vertex names within each graph. Any subset of
#'   `c("strip_version","trim","tolower")`. Default `NULL` = disabled (matches
#'   `build_multiNet()` behavior).
#' @param simplify Logical; if `TRUE`, call `igraph::simplify()` per layer to
#'   merge parallel edges and drop loops. Default `FALSE`.
#' @param simplify_weight How to combine the `weight` attribute when simplifying.
#'   One of `"sum"`, `"mean"`, `"median"`, `"min"`, `"max"`, `"first"`, `"last"`.
#'   Default `"sum"`.
#' @param results_dir Directory for outputs. Default
#'   `getOption("mlnet.results_dir","omicsDNA_results")`.
#' @param save_to_rds Logical; save the `ml.network` object to RDS. Default `TRUE`.
#' @param rds_file Optional filename for the RDS; if relative, placed under
#'   `results_dir`. When `NULL`, a timestamped name is generated.
#' @param save_layers_graphml Logical; also export each layer as GraphML. Default `FALSE`.
#' @param graphml_prefix Basename prefix for GraphML files (layer name and a
#'   timestamp are appended). Default `"layer"`.
#' @param verbose Logical; print per‑layer summaries and file paths. Default `TRUE`.
#'
#' @return A `multinet::ml.network` object with one layer per input graph.
#'
#' @seealso \code{\link{build_multiNet}} (edge tables or adjacencies);
#'   \code{\link{consensusEdges}}, \code{\link{edgesFromAdjacency}}.
#'
#' @examples
#' \dontrun{
#' g1 <- igraph::graph_from_data_frame(
#'   data.frame(from=c("a","b"), to=c("b","c"), weight=c(0.8, 0.6)), directed = FALSE)
#' g2 <- igraph::graph_from_data_frame(
#'   data.frame(from=c("a","c"), to=c("c","d"), weight=c(0.7, 0.9)), directed = FALSE)
#'
#' net <- build_multinet_from_graphs(
#'   graphsPerLayer = list(E1 = g1, E2 = g2),
#'   layerOrder     = NULL,            # alphabetical by default
#'   nodesMetadata  = NULL,
#'   save_to_rds    = FALSE
#' )
#' }
#'
#' @importFrom igraph V E vertex_attr_names set_vertex_attr simplify write_graph
#' @importFrom multinet ml_empty add_igraph_layer_ml
#' @export
build_multinet_from_graphs <- function(
    graphsPerLayer,
    layerOrder          = NULL,
    nodesMetadata       = NULL,
    featureID_col       = "feature_id",
    nodeAttrCols        = NULL,
    attr_conflict       = c("overwrite_with_meta","keep_graph","coalesce"),
    actor_normalize     = NULL,
    simplify            = FALSE,
    simplify_weight     = c("sum","mean","median","min","max","first","last"),
    results_dir         = getOption("mlnet.results_dir","omicsDNA_results"),
    save_to_rds         = TRUE,
    rds_file            = NULL,
    save_layers_graphml = FALSE,
    graphml_prefix      = "layer",
    verbose             = TRUE
) {
  attr_conflict   <- match.arg(attr_conflict)
  simplify_weight <- match.arg(simplify_weight)
  if (is.null(actor_normalize)) actor_normalize <- character(0)
  stopifnot(is.list(graphsPerLayer), length(graphsPerLayer) >= 1)

  .ensure_dir <- function(d) if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  .is_abs     <- function(p) grepl("^(/|[A-Za-z]:[\\/])", p)
  .norm <- function(x, steps) {
    x <- as.character(x)
    for (s in steps) {
      x <- switch(s,
                  "strip_version" = sub("\\.\\d+$", "", x),
                  "trim"          = trimws(x),
                  "tolower"       = tolower(x),
                  x)
    }
    x
  }

  # Validate & name layers
  lay_names <- names(graphsPerLayer)
  if (is.null(lay_names) || any(!nzchar(lay_names))) {
    lay_names <- paste0("layer", seq_along(graphsPerLayer))
    names(graphsPerLayer) <- lay_names
  }

  # Order layers: alphabetical when not specified (matches build_multiNet)
  if (is.null(layerOrder)) {
    lay_names <- sort(unique(lay_names))
    graphsPerLayer <- graphsPerLayer[lay_names]
  } else {
    lay_names <- layerOrder[layerOrder %in% lay_names]
    if (!length(lay_names)) stop("None of the requested `layerOrder` names match the provided graphs.")
    graphsPerLayer <- graphsPerLayer[lay_names]
  }

  # Clean graphs: ensure names & weight; optional normalization/simplification
  cleaned <- vector("list", length(graphsPerLayer)); names(cleaned) <- lay_names
  for (ln in lay_names) {
    g <- graphsPerLayer[[ln]]
    if (!inherits(g, "igraph")) stop("Layer '", ln, "' is not an igraph object.")

    # Ensure vertex names
    if (is.null(igraph::V(g)$name)) {
      igraph::V(g)$name <- as.character(seq_len(igraph::vcount(g)))
      if (verbose) message("Layer ", ln, ": vertex `name` was missing; numeric IDs assigned.")
    }

    # Normalize (to match typical upstream handling in some pipelines; disabled by default)
    if (length(actor_normalize)) {
      vn2 <- .norm(igraph::V(g)$name, actor_normalize)
      if (anyDuplicated(vn2)) {
        if (verbose) {
          dups <- unique(vn2[duplicated(vn2)])
          message("Layer ", ln, ": actor_normalize produced duplicate names; disambiguating with make.unique(). ",
                  "Examples: ", paste(utils::head(dups, 5), collapse = ", "),
                  if (length(dups) > 5) " ..." else "")
        }
        vn2 <- make.unique(vn2, sep = "_")
      }
      igraph::V(g)$name <- vn2
    }

    # Ensure weight exists and is numeric
    if (!("weight" %in% igraph::edge_attr_names(g))) {
      igraph::E(g)$weight <- 1
    } else {
      igraph::E(g)$weight <- suppressWarnings(as.numeric(igraph::E(g)$weight))
      if (all(is.na(igraph::E(g)$weight))) igraph::E(g)$weight <- 1
    }

    # Simplify if requested
    if (isTRUE(simplify)) {
      comb <- list(weight = switch(simplify_weight,
                                   sum="sum", mean="mean", median="median",
                                   min="min", max="max", first="first", last="last"),
                   "ignore")
      g <- igraph::simplify(g, remove.multiple = TRUE, remove.loops = TRUE, edge.attr.comb = comb)
    }

    cleaned[[ln]] <- g
  }

  # Attach vertex attributes from metadata (overwrite by default, to match build_multiNet)
  if (!is.null(nodesMetadata)) {
    if (!(featureID_col %in% names(nodesMetadata))) {
      stop("`nodesMetadata` must contain the ID column named in `featureID_col` (got '",
           featureID_col, "').")
    }
    attrs <- if (is.null(nodeAttrCols)) setdiff(names(nodesMetadata), featureID_col) else nodeAttrCols
    miss  <- setdiff(attrs, names(nodesMetadata))
    if (length(miss)) stop("nodeAttrCols not found in nodesMetadata: ", paste(miss, collapse = ", "))
    key <- as.character(nodesMetadata[[featureID_col]])

    for (ln in lay_names) {
      g  <- cleaned[[ln]]
      vn <- igraph::V(g)$name
      idx <- match(vn, key)

      for (attr in attrs) {
        meta_vals <- nodesMetadata[[attr]][idx]
        if (attr_conflict == "keep_graph" && attr %in% igraph::vertex_attr_names(g)) {
          # keep as-is
        } else if (attr_conflict == "coalesce" && attr %in% igraph::vertex_attr_names(g)) {
          cur <- igraph::vertex_attr(g, attr)
          fill <- ifelse(is.na(cur) | cur == "", meta_vals, cur)
          g <- igraph::set_vertex_attr(g, name = attr, value = fill)
        } else { # overwrite_with_meta OR new attribute
          g <- igraph::set_vertex_attr(g, name = attr, value = meta_vals)
        }
      }
      cleaned[[ln]] <- g
      if (verbose) {
        cov <- mean(!is.na(idx))
        message(sprintf("Layer %-8s | vertices: %5d | matched in metadata: %5.1f%%",
                        ln, length(vn), 100 * cov))
      }
    }
  } else if (verbose) {
    for (ln in lay_names) {
      message(sprintf("Layer %-8s | vertices: %5d | edges: %5d | directed: %s",
                      ln,
                      igraph::vcount(cleaned[[ln]]),
                      igraph::ecount(cleaned[[ln]]),
                      if (igraph::is_directed(cleaned[[ln]])) "yes" else "no"))
    }
  }

  # Assemble multilayer network
  net <- multinet::ml_empty()
  for (ln in lay_names) multinet::add_igraph_layer_ml(net, cleaned[[ln]], ln)

  # Save artifacts
  .ensure_dir(results_dir)
  stamp <- format(Sys.time(), "%Y-%m-%d_%H%M%S")

  if (isTRUE(save_to_rds)) {
    if (is.null(rds_file)) rds_file <- sprintf("multiNet_%s.rds", stamp)  # match build_multiNet default basename
    if (!.is_abs(rds_file)) rds_file <- file.path(results_dir, rds_file)
    saveRDS(net, rds_file)
    if (verbose) message("Saved multilayer network RDS: ", normalizePath(rds_file, FALSE))
  }

  if (isTRUE(save_layers_graphml)) {
    for (ln in lay_names) {
      fn <- file.path(results_dir, sprintf("%s_%s_%s.graphml", graphml_prefix, ln, stamp))
      igraph::write_graph(cleaned[[ln]], fn, format = "graphml")
    }
    if (verbose) message("Saved per-layer GraphML files under: ", normalizePath(results_dir, FALSE))
  }

  net
}

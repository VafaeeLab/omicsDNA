
# ------------------------------------------------------------------------------
# 4 - Build a multi-layer network from per-layer edge tables or adjacency matrices
# ------------------------------------------------------------------------------

#' Build a multilayer network from per‑layer edge tables or adjacency matrices
#'
#' @description
#' This routine constructs a multilayer network—one layer per group/condition—
#' and returns a `multinet::ml.network`. It is designed to accept the most common
#' representations used in practice and to normalise them into the structure that
#' `multinet` expects.
#'
#' **Accepted input formats**
#' - **Per‑layer edge tables:** a **named list** where each element is a data
#'   frame containing `from`, `to`, and (optionally) `weight`.
#' - **Flattened edge table:** a single data frame with columns `from`, `to`,
#'   optional `weight`, a layer column (auto‑detected), and optional repetition
#'   column (also auto‑detected).
#' - **Adjacency matrices:** a single square matrix, a list of matrices (one per
#'   layer), or a nested list (repeat × layer or layer × repeat). These are
#'   first converted to edge tables via `edgesFromAdjacency()` (which must be
#'   available when using adjacency inputs).
#'
#' @details
#' **Workflow and internal logic**
#' 1. **Normalisation to a layer‑indexed list.** Regardless of the input form,
#'    the function produces a **named list keyed by layer**, with each element a
#'    tidy data frame of `from`, `to`, `weight`.
#'    - For **flattened tables**, the function detects the layer column from
#'      `c("layer","Layer","group","Group")`. The repetition column (if present)
#'      is detected from `c("rep","repeat","repetition","bootstrap","iter",
#'      "iteration","Rep","Repeat")`. If not found, repetitions default to a
#'      single level `"1"`.
#'    - For **nested lists**, the option `list_order` controls whether the top
#'      level indexes **layers** or **repetitions**. With `list_order = "auto"`,
#'      the function inspects names to infer a sensible layout. Lists deeper than
#'      two levels are not supported for adjacency inputs.
#' 2. **Handling repetitions (optional).** If a repetition column is present in
#'    a flattened table, you may either:
#'    - **Filter repetitions** with `rep_filter` (named labels or numeric indices
#'      over the unique repetition labels), or
#'    - **Collapse repetitions within each layer** using `rep_collapse`
#'      (`"mean"`, `"median"`, `"sum"`, or `"none"`). When `"none"`, rows from
#'      different repetitions are kept as‑is; duplicates may therefore persist to
#'      the next step.
#' 3. **Directedness and duplicate aggregation.** For **undirected** layers
#'    (`directed = FALSE`), endpoints are canonicalised prior to deduplication
#'    (`from = pmin(from,to)`, `to = pmax(from,to)`), ensuring `A–B` and `B–A`
#'    are treated as the same edge. Duplicate `from`–`to` rows **within a layer**
#'    are combined according to `aggregate_duplicates` (`"none"`, `"mean"`,
#'    `"median"`, `"sum"`). If `"none"`, multiple parallel edges may be created
#'    in the resulting `igraph` layer.
#' 4. **Per‑layer graph construction.** Each layer is converted to an `igraph`
#'    object via `igraph::graph_from_data_frame()`. Edge weights default to `1`
#'    when missing and are coerced to numeric.
#' 5. **Node attributes (optional).** If `nodesMetadata` is supplied, vertex
#'    attributes are attached by matching `nodesMetadata[[featureID_col]]` to the
#'    vertex names. When `nodeAttrCols = NULL`, all columns except `featureID_col`
#'    are attached. Unmatched vertices receive `NA` for those attributes.
#' 6. **Assembly of the multilayer object.** The per‑layer `igraph`s are inserted
#'    into a new `multinet::ml.network` using `multinet::add_igraph_layer_ml()`.
#'
#' **Layer ordering and naming**
#' - If `layerOrder` is **not** provided, layers are ordered alphabetically by
#'   their names. If `layerOrder` **is** provided, only the listed layers are
#'   kept and they appear in the specified order.
#'
#' **Notes and limitations**
#' - Vertices appear only if referenced by at least one edge in that layer. If
#'   you require isolates, add them upstream (e.g., by inserting zero‑weight
#'   self‑loops or extending this function to accept an explicit vertex set).
#' - When using **adjacency** inputs, this function relies on
#'   `edgesFromAdjacency()`. For sparse matrices, zeros are not stored; hence a
#'   zero entry will not generate an edge.
#'
#' **Persistence and exports**
#' - When `save_to_rds = TRUE`, the multilayer network is saved as an RDS file
#'   under `results_dir` (default:
#'   `getOption("mlnet.results_dir", "omicsDNA_results")`). The filename is
#'   auto‑generated unless `rds_file` is provided (relative paths are resolved
#'   under `results_dir`; absolute paths are respected).
#' - Setting `save_layers_graphml = TRUE` additionally writes each layer as a
#'   GraphML file using `igraph::write_graph()`, with filenames prefixed by
#'   `graphml_prefix`, suffixed by the layer name and a timestamp.
#'
#' @param edgeListPerLayer A named list of per‑layer edge data frames; or a
#'   flattened edge data frame; or a single adjacency matrix; or a (nested) list
#'   of adjacency matrices. Edge data frames must contain `from`, `to`, and
#'   optional `weight` (defaulted to `1` if missing). Adjacency inputs require
#'   `edgesFromAdjacency()`.
#' @param nodesMetadata Optional data frame of node attributes. Its identifier
#'   column must be named in `featureID_col` and must match the vertex names used
#'   in the edges.
#' @param featureID_col Name of the identifier column in `nodesMetadata` used to
#'   match attributes to vertices. Default `"feature_id"`.
#' @param nodeAttrCols Character vector naming columns of `nodesMetadata` to
#'   attach as vertex attributes. If `NULL` and `nodesMetadata` is provided, all
#'   columns except `featureID_col` are attached.
#' @param layerOrder Optional character vector of layer names. If supplied, only
#'   these layers are retained and they are ordered accordingly; otherwise layers
#'   are ordered alphabetically.
#' @param directed Logical; whether per‑layer graphs are directed. Default `FALSE`.
#' @param verbose Logical; print per‑layer summaries (vertex counts and attribute
#'   match coverage) and file‑saving messages. Default `TRUE`.
#' @param layer_col,rep_col For flattened data‑frame input only: explicit column
#'   names for layer and repetition identifiers. If `NULL`, they are auto‑detected
#'   as described above.
#' @param rep_filter Optional vector (labels or numeric indices over the unique
#'   repetition labels) selecting which repetitions to keep for flattened input.
#' @param list_order For nested lists: `"auto"`, `"layer_rep"` (top level indexes
#'   layers), or `"rep_layer"` (top level indexes repetitions). Default `"auto"`.
#' @param rep_collapse How to combine multiple repetitions within each layer when
#'   present and not filtered. One of `"mean"`, `"median"`, `"sum"`, `"none"`.
#'   Default `"mean"`.
#' @param aggregate_duplicates How to combine duplicate `from`–`to` rows within a
#'   layer **after** repetition handling. One of `"none"`, `"mean"`, `"median"`,
#'   `"sum"`. Default `"none"`.
#' @param results_dir Directory to write outputs. Default
#'   `getOption("mlnet.results_dir", "omicsDNA_results")`.
#' @param save_to_rds Logical; save the `ml.network` object as an RDS file under
#'   `results_dir`. Default `TRUE`.
#' @param rds_file Optional filename for the RDS. If relative, it is created
#'   under `results_dir`; absolute paths are respected. Default `NULL` (auto‑named
#'   with a timestamp).
#' @param save_layers_graphml Logical; also export each layer as a GraphML file.
#'   Default `FALSE`.
#' @param graphml_prefix Basename prefix for per‑layer GraphML files (layer name
#'   and a timestamp are appended).
#'
#' @return A `multinet::ml.network` object containing one layer per group/condition.
#'
#' @section Practical notes:
#' - Edge weights are coerced to numeric; non‑numeric values will produce `NA`
#'   weights, which are passed through to `igraph` unchanged.
#' - If your input list lacks names, synthetic layer names (`"layer1"`, `"layer2"`, …)
#'   are generated.
#'
#' @examples
#' \dontrun{
#' ## From consensusEdges(..., as_list = TRUE):
#' net <- build_multiNet(
#'   cons_list,
#'   nodesMetadata   = genes_info,
#'   featureID_col   = "GeneName",
#'   nodeAttrCols    = "GeneType",
#'   directed        = FALSE,
#'   aggregate_duplicates = "mean"
#' )
#'
#' ## From flattened edges (e.g., edgesFromAdjacency(..., flatten = TRUE)):
#' net <- build_multiNet(
#'   edges_df,
#'   rep_collapse    = "median"
#' )
#'
#' ## From adjacency matrices (list[layer] or nested repeat × layer):
#' net <- build_multiNet(
#'   adjacency,
#'   list_order      = "rep_layer"
#' )
#' }
#'
#' @seealso
#'   \code{\link{edgesFromAdjacency}} for converting adjacencies to edge tables;
#'   \code{\link{consensusEdges}} for obtaining stable, per‑layer edge sets.
#'
#' @importFrom igraph graph_from_data_frame set_vertex_attr V write_graph
#' @importFrom multinet ml_empty add_igraph_layer_ml
#' @export
build_multiNet <- function(
    edgeListPerLayer,
    nodesMetadata         = NULL,
    featureID_col         = "feature_id",
    nodeAttrCols          = NULL,
    layerOrder            = NULL,
    directed              = FALSE,
    verbose               = TRUE,
    layer_col             = NULL,
    rep_col               = NULL,
    rep_filter            = NULL,
    list_order            = c("auto","layer_rep","rep_layer"),
    rep_collapse          = c("mean","median","sum","none"),
    aggregate_duplicates  = c("none","mean","median","sum"),
    results_dir           = getOption("mlnet.results_dir","omicsDNA_results"),
    save_to_rds           = TRUE,
    rds_file              = NULL,
    save_layers_graphml   = FALSE,
    graphml_prefix        = "layer"
) {
  list_order           <- match.arg(list_order)
  rep_collapse         <- match.arg(rep_collapse)
  aggregate_duplicates <- match.arg(aggregate_duplicates)

  # ---- tiny helpers ----
  .ensure_dir <- function(d) if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  .is_abs     <- function(p) grepl("^(/|[A-Za-z]:[\\/])", p)
  .collapse_fun <- function(which) switch(which,
                                          mean   = function(z) mean(z, na.rm = TRUE),
                                          median = function(z) stats::median(z, na.rm = TRUE),
                                          sum    = function(z) sum(z, na.rm = TRUE),
                                          none   = function(z) z[1]
  )
  .aggregate_edges <- function(df, how = "none", directed = FALSE) {
    if (!nrow(df)) return(df)
    # canonicalize endpoints for undirected graphs
    if (!directed) {
      a <- pmin(df$from, df$to); b <- pmax(df$from, df$to)
      df$from <- a; df$to <- b
    }
    if (how == "none") return(df)
    FUN <- .collapse_fun(how)
    stats::aggregate(weight ~ from + to, data = df, FUN = FUN)
  }
  .is_adj_matrix <- function(x) (is.matrix(x) || inherits(x, "Matrix")) && is.numeric(x) && nrow(x) == ncol(x)
  .is_list_node  <- function(x) is.list(x) && !is.data.frame(x)
  .list_depth    <- function(x) if (!.is_list_node(x) || !length(x)) 0L else 1L + max(vapply(x, .list_depth, integer(1L)))
  .get_leaf      <- function(x) { while (.is_list_node(x) && length(x)) x <- x[[1]]; x }

  # ---- normalize input to named list-by-layer of (from,to,weight) ----
  to_list_by_layer <- function(x) {
    # A) flattened edge data frame
    if (is.data.frame(x)) {
      stopifnot(all(c("from","to") %in% names(x)))
      if (!("weight" %in% names(x))) x$weight <- 1
      # detect layer/rep
      lc <- layer_col
      if (is.null(lc)) {
        cands <- c("layer","Layer","group","Group")
        hit <- cands[cands %in% names(x)][1]
        if (is.na(hit)) stop("Could not detect a layer column. Use `layer_col`.")
        lc <- hit
      }
      rc <- rep_col
      if (is.null(rc)) {
        rcands <- c("rep","repeat","repetition","bootstrap","iter","iteration","Rep","Repeat")
        hit <- rcands[rcands %in% names(x)][1]
        rc  <- if (is.na(hit)) NULL else hit
      }
      # standardize
      x <- data.frame(
        layer  = as.character(x[[lc]]),
        from   = as.character(x$from),
        to     = as.character(x$to),
        weight = as.numeric(x$weight),
        rep    = if (is.null(rc)) "1" else as.character(x[[rc]]),
        stringsAsFactors = FALSE
      )
      x <- x[!is.na(x$from) & !is.na(x$to) & !is.na(x$weight), , drop = FALSE]

      # reps: filter or collapse
      if (!is.null(rc)) {
        reps <- unique(x$rep)
        if (!is.null(rep_filter)) {
          # labels or indices
          keep_reps <- if (is.numeric(rep_filter)) reps[rep_filter] else rep_filter
          x <- x[x$rep %in% keep_reps, , drop = FALSE]
          x$rep <- NULL
        } else {
          if (length(reps) > 1 && rep_collapse != "none") {
            FUN <- .collapse_fun(rep_collapse)
            x <- stats::aggregate(weight ~ layer + from + to, data = x, FUN = FUN)
          } else {
            x$rep <- NULL
          }
        }
      } else {
        x$rep <- NULL
      }
      # split by layer
      return(split(x[, c("from","to","weight")], x$layer))
    }

    # B) single adjacency matrix
    if (.is_adj_matrix(x)) {
      df <- edgesFromAdjacency(x, flatten = FALSE, save_to_rds = FALSE)
      if (is.list(df)) df <- df[[1]]
      if (!nrow(df)) df <- data.frame(from=character(), to=character(), weight=numeric())
      return(list(`1` = df))
    }

    # C) lists (of matrices or edge data frames)
    if (.is_list_node(x)) {
      depth <- .list_depth(x)
      leaf  <- .get_leaf(x)
      if (.is_adj_matrix(leaf)) {
        # list of matrices → flatten with edgesFromAdjacency then reuse df path
        if (depth == 1L) {
          tmp <- edgesFromAdjacency(x, flatten = TRUE, id_cols = "layer", save_to_rds = FALSE)
          return(to_list_by_layer(tmp))
        }
        if (depth == 2L) {
          if (list_order == "auto") {
            inner_name_sets <- lapply(x, function(xx) if (.is_list_node(xx)) names(xx) else character())
            common_inners   <- Reduce(intersect, inner_name_sets)
            list_order      <- if (length(common_inners) >= 1) "rep_layer" else "layer_rep"
          }
          id_cols <- if (list_order == "rep_layer") c("rep","layer") else c("layer","rep")
          tmp <- edgesFromAdjacency(x, flatten = TRUE, id_cols = id_cols, save_to_rds = FALSE)
          return(to_list_by_layer(tmp))
        }
        stop("Nested adjacency lists deeper than 2 are not supported.")
      }
      if (is.data.frame(leaf)) {
        # depth 1: named list of per-layer edge data frames
        if (depth == 1L) {
          if (is.null(names(x)) || any(names(x) == "")) names(x) <- paste0("layer", seq_along(x))
          for (nm in names(x)) {
            df <- x[[nm]]
            if (!all(c("from","to") %in% names(df))) stop("Layer '", nm, "' must have columns `from` and `to`.")
            if (!("weight" %in% names(df))) df$weight <- 1
            x[[nm]] <- data.frame(from   = as.character(df$from),
                                  to     = as.character(df$to),
                                  weight = as.numeric(df$weight),
                                  stringsAsFactors = FALSE)
          }
          return(x)
        }
        # depth 2: list-of-lists -> build a flat df (layer + rep) then reuse df path
        if (depth == 2L) {
          if (list_order == "auto") {
            inner_name_sets <- lapply(x, function(xx) if (.is_list_node(xx)) names(xx) else character())
            common_inners   <- Reduce(intersect, inner_name_sets)
            list_order      <- if (length(common_inners) >= 1) "rep_layer" else "layer_rep"
          }
          rows <- list()
          for (i in seq_along(x)) {
            top_name <- names(x)[i]; if (is.null(top_name) || top_name == "") top_name <- as.character(i)
            inner    <- x[[i]]
            for (j in seq_along(inner)) {
              inner_name <- names(inner)[j]; if (is.null(inner_name) || inner_name == "") inner_name <- as.character(j)
              df <- inner[[j]]
              if (!("weight" %in% names(df))) df$weight <- 1
              df <- data.frame(from=as.character(df$from), to=as.character(df$to), weight=as.numeric(df$weight), stringsAsFactors = FALSE)
              if (list_order == "rep_layer") {
                rows[[length(rows)+1L]] <- cbind(data.frame(rep=top_name, layer=inner_name, stringsAsFactors = FALSE), df)
              } else {
                rows[[length(rows)+1L]] <- cbind(data.frame(layer=top_name, rep=inner_name, stringsAsFactors = FALSE), df)
              }
            }
          }
          flat <- do.call(rbind, rows)
          return(to_list_by_layer(flat))
        }
      }
      stop("Unsupported list content for `edgeListPerLayer`.")
    }

    stop("Unsupported input type for `edgeListPerLayer`.")
  }

  # ---- normalize ----
  layers <- to_list_by_layer(edgeListPerLayer)
  lay_names <- names(layers)
  if (is.null(lay_names) || any(!nzchar(lay_names))) {
    lay_names <- paste0("L", seq_along(layers))
    names(layers) <- lay_names
  }

  # order layers
  if (is.null(layerOrder)) {
    lay_names <- sort(unique(lay_names))
  } else {
    lay_names <- layerOrder[layerOrder %in% lay_names]
    if (!length(lay_names)) stop("No valid layer names after applying `layerOrder`.")
  }

  # ---- build per-layer igraphs ----
  graphs_list <- vector("list", length(lay_names)); names(graphs_list) <- lay_names
  for (ln in lay_names) {
    df <- layers[[ln]]
    # clean & default
    if (!nrow(df)) {
      graphs_list[[ln]] <- igraph::graph_from_data_frame(
        data.frame(from=character(), to=character(), weight=numeric()),
        directed = directed
      )
      next
    }
    if (!("weight" %in% names(df))) df$weight <- 1
    df <- data.frame(from=as.character(df$from), to=as.character(df$to), weight=as.numeric(df$weight), stringsAsFactors = FALSE)

    # aggregate duplicates (also canonicalize for undirected)
    df <- .aggregate_edges(df, how = aggregate_duplicates, directed = directed)

    g <- igraph::graph_from_data_frame(df[, c("from","to","weight")], directed = directed)

    # attach node attributes
    if (!is.null(nodesMetadata)) {
      if (!(featureID_col %in% names(nodesMetadata))) {
        stop("`nodesMetadata` must contain the ID column: ", featureID_col)
      }
      vn  <- igraph::V(g)$name
      key <- as.character(nodesMetadata[[featureID_col]])
      idx <- match(vn, key)

      attrs <- if (is.null(nodeAttrCols)) setdiff(names(nodesMetadata), featureID_col) else nodeAttrCols
      miss  <- setdiff(attrs, names(nodesMetadata))
      if (length(miss)) stop("nodeAttrCols not found in nodesMetadata: ", paste(miss, collapse = ", "))

      for (attr in attrs) {
        vals <- ifelse(is.na(idx), NA, nodesMetadata[[attr]][idx])
        g <- igraph::set_vertex_attr(g, name = attr, value = vals)
      }

      if (verbose) {
        cov <- mean(!is.na(idx))
        message(sprintf("Layer %-6s | vertices: %5d | attribute match: %5.1f%%",
                        ln, length(vn), 100 * cov))
      }
    } else if (verbose) {
      message(sprintf("Layer %-6s | vertices: %5d", ln, length(igraph::V(g))))
    }

    graphs_list[[ln]] <- g
  }

  # ---- assemble multilayer network ----
  net <- multinet::ml_empty()
  for (ln in lay_names) multinet::add_igraph_layer_ml(net, graphs_list[[ln]], ln)

  if (verbose) {
    cat("Layer summary (unique actors per layer):\n")
    actor_counts <- vapply(lay_names, function(ln) length(igraph::V(graphs_list[[ln]])$name), integer(1))
    print(data.frame(Layer = lay_names, ActorCount = actor_counts, row.names = NULL))
  }

  # ---- save artifacts ----
  .ensure_dir(results_dir)
  stamp <- format(Sys.time(), "%Y-%m-%d_%H%M%S")

  if (isTRUE(save_to_rds)) {
    if (is.null(rds_file)) rds_file <- sprintf("multiNet_%s.rds", stamp)
    if (!.is_abs(rds_file)) rds_file <- file.path(results_dir, rds_file)
    saveRDS(net, rds_file)
    if (verbose) message("Saved multilayer network RDS: ", normalizePath(rds_file, FALSE))
  }

  if (isTRUE(save_layers_graphml)) {
    for (ln in lay_names) {
      fn <- file.path(results_dir, sprintf("%s_%s_%s.graphml", graphml_prefix, ln, stamp))
      igraph::write_graph(graphs_list[[ln]], fn, format = "graphml")
    }
    if (verbose) message("Saved per-layer GraphML files under: ", normalizePath(results_dir, FALSE))
  }

  net
}


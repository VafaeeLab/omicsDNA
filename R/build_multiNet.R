# ------------------------------------------------------------------------------
# 4 - Build a multi-layer network from per-layer edge tables or adjacency matrices
# ------------------------------------------------------------------------------

#' Build a multilayer network from edge lists or adjacency matrices
#'
#' @description
#' Constructs a \pkg{multinet} multilayer network in which each layer represents a
#' group/condition (e.g., time point, omics layer, experimental group). The function
#' accepts common network input formats (edge tables or adjacency matrices), converts
#' them into per-layer \pkg{igraph} objects, and then assembles them into a single
#' \code{multinet::ml.network}.
#'
#' @details
#' \section{What is the input to this function?}{
#' The core input is \code{edgeListPerLayer}, which can be provided in several practical formats:
#' \itemize{
#'   \item \strong{Named list of per-layer edge tables:} each element is a data frame
#'         with columns \code{from}, \code{to}, and optionally \code{weight}.
#'   \item \strong{Single “flattened” edge table:} one data frame containing \code{from}, \code{to},
#'         optional \code{weight}, and a layer identifier column (automatically detected or specified
#'         by \code{layer_col}). A repetition column may also be present (detected or specified by
#'         \code{rep_col}).
#'   \item \strong{Adjacency matrices:} either a single square matrix, a list of matrices (one per layer),
#'         or a nested list representing repetitions and layers. Adjacency inputs are converted to edge
#'         tables using \code{edgesFromAdjacency()}.
#' }
#'
#' Optionally, \code{nodesMetadata} can be supplied to attach node (actor) attributes to each layer.
#' Vertices are matched using \code{nodesMetadata[[featureID_col]]} against vertex names.
#' }
#'
#' \section{How the network is built (high-level logic)}{
#' \enumerate{
#'   \item \strong{Standardise inputs:} all supported formats are converted into a named list
#'         \code{layer -> data.frame(from,to,weight)}.
#'   \item \strong{Handle repetitions (if present):} repetitions can be filtered (\code{rep_filter})
#'         and/or collapsed within each layer (\code{rep_collapse}).
#'   \item \strong{Resolve duplicates:} duplicate edges within a layer can be combined using
#'         \code{aggregate_duplicates}. For undirected networks, endpoints are canonicalised so that
#'         \code{A–B} and \code{B–A} are treated as the same edge.
#'   \item \strong{Build layer graphs:} each layer is converted to an \pkg{igraph} object.
#'   \item \strong{Attach node attributes (optional):} selected columns from \code{nodesMetadata} are
#'         assigned as vertex attributes.
#'   \item \strong{Assemble multilayer network:} layer graphs are inserted into a new
#'         \code{multinet::ml.network}.
#' }
#' }
#'
#' \section{What files does this function write?}{
#' Writing files is optional and controlled by flags:
#' \itemize{
#'   \item If \code{save_to_rds = TRUE}, the multilayer network object is saved as an \code{.rds}
#'         file. If \code{rds_file} is \code{NULL}, the function auto-generates a timestamped filename
#'         (e.g., \code{multiNet_<timestamp>.rds}) under \code{results_dir}.
#'   \item If \code{save_layers_graphml = TRUE}, each layer graph is exported as a \code{.graphml} file
#'         under \code{results_dir} with names of the form
#'         \code{<graphml_prefix>_<layer>_<timestamp>.graphml}.
#' }
#' If both flags are \code{FALSE}, no files are written and the network exists only in memory.
#' }
#'
#' @param edgeListPerLayer Network definition provided as either (i) a named list of per-layer edge
#'   tables, (ii) a single flattened edge table containing a layer column, (iii) an adjacency matrix,
#'   or (iv) a list/nested list of adjacency matrices. Edge tables must contain \code{from}, \code{to},
#'   and optionally \code{weight} (default weight is 1 when missing). Adjacency inputs require
#'   \code{edgesFromAdjacency()} to be available.
#' @param nodesMetadata Optional data frame of node (actor) metadata to attach as vertex attributes.
#'   Row identifiers are matched to vertex names via \code{featureID_col}.
#' @param featureID_col Column name in \code{nodesMetadata} that contains the vertex identifiers used
#'   for matching. Default is \code{"feature_id"}.
#' @param nodeAttrCols Character vector naming which columns of \code{nodesMetadata} to attach as
#'   vertex attributes. If \code{NULL}, all columns except \code{featureID_col} are attached.
#' @param layerOrder Optional vector defining the layer order. If provided, only listed layers are used
#'   and they appear in this order; otherwise layers are ordered alphabetically.
#' @param directed Logical; if \code{TRUE}, treat edges as directed within each layer. Default \code{FALSE}.
#' @param verbose Logical; if \code{TRUE}, print per-layer summaries and file-saving messages.
#' @param layer_col,rep_col For flattened edge-table input only: explicit column names identifying the
#'   layer and repetition variables. If \code{NULL}, common names are auto-detected.
#' @param rep_filter Optional repetition filter for flattened input: can be repetition labels or numeric
#'   indices over the unique repetition labels.
#' @param list_order For nested list inputs: specifies whether the top level indexes layers or repetitions.
#'   One of \code{"auto"}, \code{"layer_rep"}, or \code{"rep_layer"}.
#' @param rep_collapse How to combine repetitions within each layer when multiple repetitions exist.
#'   One of \code{"mean"}, \code{"median"}, \code{"sum"}, or \code{"none"}.
#' @param aggregate_duplicates How to combine duplicated edges within a layer after repetition handling.
#'   One of \code{"none"}, \code{"mean"}, \code{"median"}, or \code{"sum"}.
#' @param results_dir Directory used for output files (RDS/GraphML) when saving is enabled.
#' @param save_to_rds Logical; if \code{TRUE}, save the resulting \code{ml.network} as an \code{.rds}.
#' @param rds_file Optional filename for the RDS output. Relative paths are created under \code{results_dir};
#'   absolute paths are respected. If \code{NULL}, a timestamped name is used.
#' @param save_layers_graphml Logical; if \code{TRUE}, also export each layer to GraphML format.
#' @param graphml_prefix Prefix used in GraphML filenames.
#'
#' @return A \code{multinet::ml.network} object containing one layer per condition/group.
#'
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


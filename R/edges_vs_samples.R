
#------------------------------------------------------
# 13 - Check relationship between #samples & #edges | Scatter plot
#------------------------------------------------------

#' Edges vs. sample size across layers
#'
#' @title Plot the relationship between per-layer edge counts and sample size
#' @description
#' For each network layer, compute the number of edges and draw a scatter plot
#' (y = edges, x = samples) with an OLS trend line. Returns a named integer vector
#' (one value per layer) and attaches the ggplot object under attr(x, "plot").
#'
#' @param communities Optional data.frame with a 'layer' column (used only as a
#'   hint to restrict/reorder layers for compatibility with older code).
#' @param num_samples Named numeric vector of sample counts per layer (non-negative).
#' @param net Optional multilayer object that can be coerced to a list of igraphs
#'   via as.list(net); edges are counted with igraph::ecount().
#' @param edges Optional data.frame of edges with endpoints and layer information.
#'   Endpoints auto-detected among common pairs (from/to, source/target, etc.).
#'   If from_layer/to_layer exist, only intra-layer rows are kept.
#' @param num_edges Optional named vector of precomputed edge counts per layer.
#'   If supplied, it takes precedence over `net`/`edges`.
#' @param layers Optional character vector to restrict/reorder layers.
#' @param directed Logical; if FALSE (default), undirected duplicates are
#'   canonicalised to min(u,v)–max(u,v) before counting.
#'
#' @return Named integer vector of edge counts per layer with attribute "plot".
#' @importFrom ggplot2 ggplot aes geom_point geom_smooth labs theme_minimal
#' @importFrom igraph ecount
#' @export
edges_vs_samples <- function(communities,
                             num_samples,
                             net       = NULL,
                             edges     = NULL,
                             num_edges = NULL,
                             layers    = NULL,
                             directed  = FALSE) {

  ## ---- validate num_samples (preserve names!) ----
  if (is.null(num_samples) || !length(num_samples))
    stop("`num_samples` must be a non-empty named numeric vector.")
  if (is.null(names(num_samples)) || any(!nzchar(names(num_samples))))
    stop("`num_samples` must have non-empty names (layer IDs).")

  ns_names <- as.character(names(num_samples))
  ns_vals  <- suppressWarnings(as.numeric(num_samples))
  if (any(!is.finite(ns_vals) | ns_vals < 0))
    stop("`num_samples` must be finite and non-negative.")
  names(ns_vals) <- ns_names

  ## ---- print number of samples per layer to the console ----
  cat("Samples per layer (n):\n")
  for (i in seq_along(ns_vals)) {
    cat(sprintf("  • %s: %s\n", ns_names[i], format(ns_vals[i], trim = TRUE, scientific = FALSE)))
  }

  ## ---- optional layer hints from communities ----
  layer_hint <- NULL
  if (!is.null(communities)) {
    if (!is.data.frame(communities))
      stop("`communities` must be a data.frame when provided.")
    if (!"layer" %in% names(communities))
      stop("`communities` must contain a 'layer' column when provided.")
    layer_hint <- unique(as.character(communities$layer))
  }

  ## ---- compute per-layer edge counts: num_edges | net | edges ----
  edge_counts <- NULL

  # (A) explicit counts win
  if (!is.null(num_edges)) {
    if (is.null(names(num_edges)) || any(!nzchar(names(num_edges))))
      stop("`num_edges` must have names (layer IDs).")
    edge_counts <- as.integer(round(num_edges))
    names(edge_counts) <- as.character(names(num_edges))
  }

  # (B) otherwise from multinet/igraph list
  if (is.null(edge_counts) && !is.null(net)) {
    glist <- try(as.list(net), silent = TRUE)
    if (inherits(glist, "try-error") || !is.list(glist))
      stop("Could not coerce `net` to a list of igraphs via as.list(net).")

    # drop non-igraph entries and reserved '_flat_' if present
    if (!is.null(names(glist))) {
      glist <- glist[names(glist) != "_flat_"]
    }
    glist <- Filter(function(g) inherits(g, "igraph"), glist)
    if (!length(glist))
      stop("No igraph layers found in `net` after filtering.")

    # FIX: ecount() returns double; coerce to integer to satisfy vapply's type
    edge_counts <- vapply(
      glist,
      function(g) as.integer(igraph::ecount(g)),
      integer(1)
    )
    if (is.null(names(edge_counts)) && !is.null(names(glist)))
      names(edge_counts) <- names(glist)
  }

  # (C) otherwise from an edges data frame
  if (is.null(edge_counts) && !is.null(edges)) {
    if (!is.data.frame(edges))
      stop("`edges` must be a data.frame when provided.")

    nm <- names(edges)
    # heuristics for endpoint columns
    endpoint_pairs <- list(
      c("from","to"), c("source","target"),
      c("actor1","actor2"), c("i","j"), c("v1","v2")
    )
    a_col <- b_col <- NA_character_
    for (p in endpoint_pairs) if (all(p %in% nm)) { a_col <- p[1]; b_col <- p[2]; break }
    if (is.na(a_col)) {
      layer_like <- c("layer","Layer","from_layer","to_layer","group","Group","l1","l2")
      char_cols  <- which(vapply(edges, function(x) is.character(x) || is.factor(x), logical(1)))
      char_cols  <- setdiff(char_cols, match(layer_like, nm, nomatch = 0))
      if (length(char_cols) < 2)
        stop("Could not identify two endpoint columns in `edges`.")
      a_col <- nm[char_cols[1]]; b_col <- nm[char_cols[2]]
    }

    # normalize / derive layer column
    layer_col <- if ("layer" %in% nm) {
      "layer"
    } else if (all(c("from_layer","to_layer") %in% nm)) {
      keep <- as.character(edges$from_layer) == as.character(edges$to_layer)
      edges <- edges[keep, , drop = FALSE]
      "from_layer"
    } else if ("Layer" %in% nm) {
      "Layer"
    } else if ("group" %in% nm) {
      "group"
    } else if ("Group" %in% nm) {
      "Group"
    } else {
      edges$layer <- "L1"; "layer"
    }

    if (!nrow(edges)) stop("No intra-layer edges remain after filtering.")

    # canonicalize endpoints if undirected
    if (!isTRUE(directed)) {
      u <- pmin(as.character(edges[[a_col]]), as.character(edges[[b_col]]))
      v <- pmax(as.character(edges[[a_col]]), as.character(edges[[b_col]]))
    } else {
      u <- as.character(edges[[a_col]])
      v <- as.character(edges[[b_col]])
    }
    key <- paste(u, v, sep = "\t")
    lay <- as.character(edges[[layer_col]])

    df <- unique(data.frame(layer = lay, key = key, stringsAsFactors = FALSE))
    tab <- table(df$layer)
    edge_counts <- setNames(as.integer(tab), names(tab))
  }

  if (is.null(edge_counts))
    stop("Provide one of: `num_edges`, `net`, or `edges` to compute per-layer edge counts.")

  ## ---- restrict/reorder layers if requested/hinted ----
  if (length(layers)) {
    keep <- intersect(as.character(layers), names(edge_counts))
    if (!length(keep))
      stop("None of the requested `layers` are present in the edge source.")
    edge_counts <- edge_counts[keep]
  } else if (length(layer_hint)) {
    keep <- intersect(as.character(layer_hint), names(edge_counts))
    if (length(keep)) edge_counts <- edge_counts[keep]
  }

  ## ---- strict name matching against num_samples ----
  if (!setequal(names(edge_counts), names(ns_vals))) {
    stop(
      "Layer name mismatch between edge counts and `num_samples`.\n",
      "  edges:      {", paste(sort(names(edge_counts)), collapse = ", "), "}\n",
      "  num_samples:{", paste(sort(ns_names),        collapse = ", "), "}"
    )
  }

  # align order
  edge_counts <- edge_counts[ns_names]

  ## ---- plot: x=samples, y=edges ----
  df_plot <- data.frame(
    layer   = ns_names,
    samples = as.numeric(ns_vals),
    edges   = as.integer(edge_counts),
    stringsAsFactors = FALSE
  )

  plt <- ggplot2::ggplot(df_plot, ggplot2::aes(x = samples, y = edges)) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_smooth(method = "lm", se = FALSE) +
    ggplot2::labs(
      x = "Number of samples (per group)",
      y = "Number of edges (per layer network)",
      title = "Edges vs. sample size across layers"
    ) +
    ggplot2::theme_minimal()

  print(plt)

  attr(edge_counts, "plot") <- plt
  invisible(edge_counts)
}


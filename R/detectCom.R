
# ---------------------------------------
# 11 - detectCom
# ---------------------------------------

#' Detect communities on a multilayer network via a layer‑collapsed supra‑graph
#'
#' @description
#' This function performs community detection on a **single, undirected supra‑graph**
#' obtained by collapsing the selected layers of a multilayer network. Parallel
#' edges across layers are aggregated either by **presence count** (in how many
#' layers an edge appears) or by **weight sum** (sum of edge weights across
#' layers). A community partition is then computed once on the supra‑graph and
#' mapped back to each layer. Communities are labelled `C1`, `C2`, … with `C1`
#' corresponding to the largest community by default.
#'
#' Supported methods: `"louvain"` (default), `"infomap"`, a simple `"clique"`
#' (k‑clique union followed by connected components), and—when available—
#' `"abacus"` via `multinet::abacus_ml()`.
#'
#' @details
#' **Input normalisation and scope**
#' - Edges are retrieved from `net` using `multinet::edges_ml()`. If both
#'   `from_layer` and `to_layer` columns are present, only **intra‑layer**
#'   edges are retained (`from_layer == to_layer`). Otherwise, a single `layer`
#'   column is used (or created as `"L1"` if absent).
#' - The `layers` argument limits the analysis to a subset of layers; by default
#'   all available layers are used (after intra‑layer restriction).
#' - Edge weights are taken from the first column found among
#'   `c("weight","Weight","w","value","score")`. If none is found, weights
#'   default to `1` per edge.
#'
#' **Supra‑graph construction**
#' - For each selected layer, edges are reduced to `(from, to, weight)`, then
#'   endpoint order is canonicalised (`from = pmin(from, to)`,
#'   `to = pmax(from, to)`), so the supra‑graph is undirected.
#' - If `edgeWeight = "count"`, an edge’s weight equals the **number of distinct
#'   layers** in which that undirected pair appears. If `edgeWeight = "sum"`,
#'   weights are **summed across layers**. The supra‑graph is built with
#'   `igraph::graph_from_data_frame()` and simplified with `weight = "sum"`.
#'
#' **Community detection methods**
#' - `"louvain"`: modularity maximisation via `igraph::cluster_louvain()` using
#'   supra‑graph edge weights.
#' - `"infomap"`: flow‑based clustering via `igraph::cluster_infomap()` using
#'   supra‑graph edge weights.
#' - `"clique"`: compute all cliques of size ≥ `clique.k` on the supra‑graph,
#'   connect every pair of actors co‑occurring in any such clique (k‑clique
#'   union), and take connected components as communities. This is a pragmatic
#'   fallback and may be costly for dense graphs and large `k`.
#' - `"abacus"`: if `multinet::abacus_ml()` is exported in your `multinet` build,
#'   community assignments are obtained directly from ABACUS, then **filtered**
#'   by `min.actors` and `min.layers` (after any layer subsetting). The same
#'   relabelling and file‑output logic is applied.
#'
#' **Mapping back to layers and filtering**
#' - After detecting communities on the supra‑graph, assignments are mapped back
#'   to each selected layer by keeping only actors that actually appear in that
#'   layer’s edges. Communities are then **filtered** to retain only those
#'   meeting both criteria: at least `min.actors` unique actors (pooled across
#'   all selected layers) and presence in at least `min.layers` **distinct**
#'   layers.
#'
#' **Reproducibility and labels**
#' - If `seed` is provided, the random seed is set before running Louvain or
#'   Infomap. (Note that some implementations may still exhibit mild
#'   non‑determinism.) When `relabel_by_size = TRUE`, community IDs are re‑coded
#'   so that `C1` is largest, `C2` next, and so on.
#'
#' **Persistence**
#' - When `save_to_rds = TRUE`, the assignments are written to
#'   `file.path(results_dir, sprintf("communities_%s_%s_<timestamp>.rds",
#'   tolower(method), tolower(edgeWeight)))`. If `write_csv = TRUE`, a CSV with
#'   the same basename is also produced. If `write_summary_csv = TRUE`, an
#'   additional CSV summarises each community’s size and layer span.
#'
#' @param net        A multilayer network compatible with `multinet::edges_ml()`.
#' @param method     Community detection method: one of
#'   `c("louvain","infomap","clique","abacus")`. Default `"louvain"`.
#' @param layers     Character vector naming layers to include. Default: all
#'   layers present after intra‑layer filtering.
#' @param edgeWeight How to aggregate parallel edges across layers:
#'   `"count"` (number of layers where an edge appears) or `"sum"` (sum of
#'   edge weights across layers). Default `"count"`.
#' @param min.actors Minimum number of unique actors required for a community to
#'   be kept (after mapping back to layers). Default `15`.
#' @param min.layers Minimum number of **distinct layers** in which a community
#'   must appear to be kept. Default `2`.
#' @param clique.k   For `method = "clique"`, minimum clique size \(k\). Default `3`.
#' @param seed       Optional integer seed for Louvain/Infomap runs. Default `NULL`.
#' @param relabel_by_size Logical; if `TRUE`, relabel communities by decreasing
#'   size so `C1` is largest. Default `TRUE`.
#' @param results_dir Output directory for saved artefacts. Default
#'   `getOption("mlnet.results_dir", "omicsDNA_results")`.
#' @param save_to_rds Logical; write the assignments as an RDS file. Default `TRUE`.
#' @param rds_file   Optional file name for the RDS output; relative paths are
#'   resolved under `results_dir`, absolute paths are respected. Default `NULL`
#'   (auto‑named).
#' @param write_csv  Logical; additionally write a CSV of assignments. Default `FALSE`.
#' @param csv_prefix Basename prefix for CSV files. Default `"communities"`.
#' @param write_summary_csv Logical; write a summary CSV (community size and
#'   layer span). Default `TRUE`.
#' @param verbose    Logical; print progress and file locations. Default `TRUE`.
#'
#' @return A data frame with columns:
#'   - `actor`: actor identifier,
#'   - `com`: community label (`C1`, `C2`, …),
#'   - `layer`: layer name,
#'   - `method`: the method used.
#'
#'   The return value carries attributes:
#'   - `"community_sizes"` — named integer vector of community sizes (ordered),
#'   - `"layer_span"`      — named integer vector of number of layers per community,
#'   - `"files"`           — named list of output paths (if artefacts were written).
#'
#' @section Practical notes
#' - Cross‑layer edges (if present) are ignored; only intra‑layer edges contribute
#'   to the supra‑graph and to per‑layer mappings.
#' - If no weight column is detected, all edges are treated as weight `1`.
#' - The `"clique"` method may be computationally expensive for dense graphs or
#'   larger `clique.k` due to clique enumeration.
#'
#' @examples
#' \dontrun{
#' # Louvain on all layers; edge weight = number of layers in which an edge appears
#' comm <- detectCom(
#'   net,
#'   method      = "louvain",
#'   edgeWeight  = "count",
#'   min.actors  = 15,
#'   min.layers  = 2,
#'   seed        = 1,
#'   save_to_rds = TRUE,
#'   write_csv   = FALSE,
#'   verbose     = TRUE
#' )
#'
#' # Clique-based communities using k = 4, focusing on selected layers
#' comm_clq <- detectCom(
#'   net,
#'   method     = "clique",
#'   layers     = c("GroupA","GroupB"),
#'   clique.k   = 4,
#'   edgeWeight = "sum",
#'   write_summary_csv = TRUE
#' )
#' }
#'
#' @seealso
#'   \code{\link{build_multiNet}} to assemble multilayer graphs,
#'   \code{\link{consensusEdges}} for robust per‑layer edges,
#'   \code{\link{add_network_attributes}} to enrich networks prior to clustering.
#'
#' @importFrom multinet edges_ml
#' @importFrom igraph graph_from_data_frame simplify cluster_louvain cluster_infomap cliques components E V membership ecount
#' @export
detectCom <- function(net,
                      method            = c("louvain","infomap","clique","abacus"),
                      layers            = NULL,
                      edgeWeight        = c("count","sum"),
                      min.actors        = 15,
                      min.layers        = 2,
                      clique.k          = 3,
                      seed              = NULL,
                      relabel_by_size   = TRUE,
                      results_dir       = getOption("mlnet.results_dir","omicsDNA_results"),
                      save_to_rds       = TRUE,
                      rds_file          = NULL,
                      write_csv         = FALSE,
                      csv_prefix        = "communities",
                      write_summary_csv = TRUE,
                      verbose           = TRUE) {

  method     <- match.arg(method)
  edgeWeight <- match.arg(edgeWeight)
  stopifnot(is.numeric(min.actors), min.actors >= 1,
            is.numeric(min.layers), min.layers >= 1)
  if (method == "clique") stopifnot(is.numeric(clique.k), clique.k >= 2)

  .ensure_dir <- function(d) if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  .is_abs     <- function(p) grepl("^(/|[A-Za-z]:[\\/])", p)
  .sanitize   <- function(s) gsub("[^[:alnum:]_.-]+", "_", s)
  .pick <- function(cands, nms) { z <- cands[cands %in% nms]; if (length(z)) z[1] else NA_character_ }

  # ---- 0) Pull all multilayer edges once ----
  Eraw <- try(multinet::edges_ml(net), silent = TRUE)
  if (inherits(Eraw, "try-error") || is.null(Eraw)) {
    stop("Could not retrieve edges from `net` via multinet::edges_ml().")
  }
  Eall <- if (is.data.frame(Eraw)) Eraw else {
    tmp <- try(as.data.frame(Eraw, stringsAsFactors = FALSE), silent = TRUE)
    if (inherits(tmp, "try-error") || is.null(tmp) || !nrow(tmp)) {
      stop("`edges_ml(net)` did not return a coercible table of edges.")
    }
    tmp
  }
  nm <- names(Eall)

  # Endpoint columns (robust heuristics)
  pairs <- list(
    c("from_actor","to_actor"), c("from","to"), c("source","target"),
    c("actor1","actor2"), c("i","j"), c("v1","v2")
  )
  a_col <- b_col <- NA_character_
  for (p in pairs) if (all(p %in% nm)) { a_col <- p[1]; b_col <- p[2]; break }
  if (is.na(a_col)) {
    layer_like <- c("layer","Layer","from_layer","to_layer","l1","l2")
    char_cols  <- which(vapply(Eall, function(x) is.character(x) || is.factor(x), logical(1)))
    char_cols  <- setdiff(char_cols, match(layer_like, nm, nomatch = 0))
    if (length(char_cols) < 2)
      stop("Could not identify two endpoint columns in the global edge table.")
    a_col <- nm[char_cols[1]]; b_col <- nm[char_cols[2]]
  }

  # Layer column(s); keep intra-layer
  if ("from_layer" %in% nm && "to_layer" %in% nm) {
    Eall <- Eall[Eall$from_layer == Eall$to_layer, , drop = FALSE]
    Eall$layer <- as.character(Eall$from_layer)
  } else if ("layer" %in% nm) {
    Eall$layer <- as.character(Eall$layer)
  } else {
    Eall$layer <- "L1"
  }
  if (!nrow(Eall)) stop("No intra-layer edges found in the network.")

  # Select layers
  if (is.null(layers)) {
    layers <- sort(unique(Eall$layer))
  } else {
    layers <- intersect(as.character(layers), unique(Eall$layer))
  }
  if (!length(layers)) stop("No layers available after intersection with requested `layers`.")

  # ---- 1) Minimal per-layer tables ----
  w_col_global <- .pick(c("weight","Weight","w","value","score"), nm)
  edge_list <- lapply(layers, function(ly) {
    ed <- Eall[Eall$layer == ly, , drop = FALSE]
    if (!nrow(ed)) return(ed[FALSE, , drop = FALSE])
    from <- as.character(ed[[a_col]])
    to   <- as.character(ed[[b_col]])
    w    <- if (!is.na(w_col_global) && w_col_global %in% names(ed)) as.numeric(ed[[w_col_global]]) else rep(1, length(from))
    data.frame(layer = ly, from = from, to = to, weight = w, stringsAsFactors = FALSE)
  })
  names(edge_list) <- layers

  # ---- 2) Supra-graph (undirected) ----
  edges_all <- do.call(rbind, edge_list)
  if (is.null(edges_all) || !nrow(edges_all))
    stop("No edges found in the selected layers.")

  # canonicalize pairs (undirected)
  aa <- pmin(edges_all$from, edges_all$to)
  bb <- pmax(edges_all$from, edges_all$to)
  edges_all$from <- aa; edges_all$to <- bb

  if (edgeWeight == "count") {
    per_layer_unique <- unique(edges_all[, c("layer","from","to")])
    agg <- stats::aggregate(rep(1L, nrow(per_layer_unique)),
                            by = per_layer_unique[c("from","to")], FUN = sum)
    names(agg)[3] <- "weight"
  } else { # "sum"
    agg <- stats::aggregate(weight ~ from + to, data = edges_all, FUN = sum, na.rm = TRUE)
  }

  g <- igraph::graph_from_data_frame(agg, directed = FALSE)
  if (igraph::ecount(g) == 0) stop("Supra-graph has no edges after aggregation.")
  g <- igraph::simplify(g, edge.attr.comb = list(weight = "sum", "ignore"))

  if (!is.null(seed)) set.seed(as.integer(seed))
  if (verbose) message("• Built supra-graph: |V|=", length(igraph::V(g)), " |E|=", igraph::ecount(g),
                       "  (edgeWeight=", edgeWeight, ", method=", method, ")")

  # ---- 3) Community detection ----
  memb <- NULL
  if (method == "louvain") {
    co <- igraph::cluster_louvain(g, weights = igraph::E(g)$weight)
    memb <- igraph::membership(co)
  } else if (method == "infomap") {
    co <- igraph::cluster_infomap(g, e.weights = igraph::E(g)$weight)
    memb <- igraph::membership(co)
  } else if (method == "clique") {
    cls <- igraph::cliques(g, min = clique.k)
    if (!length(cls)) stop("No cliques of size >= ", clique.k, " found.")
    # Build clique-union graph: connect all pairs that co-occur in any k-clique
    pairs_df <- do.call(rbind, lapply(cls, function(cl) {
      vs <- igraph::V(g)$name[cl]
      if (length(vs) < 2) return(NULL)
      m <- utils::combn(vs, 2)
      data.frame(from = m[1,], to = m[2,], stringsAsFactors = FALSE)
    }))
    pairs_df <- unique(pairs_df)
    h <- igraph::graph_from_data_frame(pairs_df, directed = FALSE, vertices = igraph::V(g)$name)
    comp <- igraph::components(h)
    memb <- comp$membership; names(memb) <- igraph::V(h)$name
  } else if (method == "abacus") {
    if (!"abacus_ml" %in% getNamespaceExports("multinet"))
      stop("`abacus_ml` is not available in this multinet build.")
    ab <- multinet::abacus_ml(net, min.actors = min.actors, min.layers = min.layers)
    if (!all(c("actor","com","layer") %in% names(ab)))
      stop("Unexpected ABACUS output format.")
    # If user selected a subset of layers, filter and re-check span/size.
    if (!is.null(layers)) ab <- ab[ab$layer %in% layers, , drop = FALSE]
    size_by <- tapply(ab$actor, ab$com, function(v) length(unique(v)))
    span_by <- tapply(ab$layer, ab$com, function(v) length(unique(v)))
    keep    <- names(size_by)[size_by >= min.actors & span_by >= min.layers]
    res <- ab[ab$com %in% keep, c("actor","com","layer")]
    res$method <- "abacus"
    # relabel by size if requested
    if (relabel_by_size && length(keep)) {
      ord <- order(size_by[keep], decreasing = TRUE)
      map <- setNames(paste0("C", seq_along(ord)), names(size_by[keep])[ord])
      res$com <- unname(map[res$com])
      size_by <- setNames(as.integer(size_by[keep][ord]), names(map))   # temporarily old names
      names(size_by) <- unname(map[names(size_by)])                     # rename to C1, C2, ...
      span_by <- span_by[keep][ord]; names(span_by) <- names(size_by)
    }
    attr(res, "community_sizes") <- size_by[unique(names(size_by))]
    attr(res, "layer_span")      <- span_by[unique(names(span_by))]
    # save artifacts
    files <- list()
    .ensure_dir(results_dir); stamp <- format(Sys.time(), "%Y-%m-%d_%H%M%S")
    method_tag <- tolower(method); weight_tag <- tolower(edgeWeight)
    if (isTRUE(save_to_rds)) {
      if (is.null(rds_file)) rds_file <- sprintf("communities_%s_%s_%s.rds", method_tag, weight_tag, stamp)
      if (!.is_abs(rds_file)) rds_file <- file.path(results_dir, rds_file)
      saveRDS(res, rds_file); files$rds <- rds_file
      if (verbose) message("Saved RDS: ", normalizePath(rds_file, FALSE))
    }
    if (isTRUE(write_csv)) {
      csv_file <- file.path(results_dir, sprintf("%s_%s_%s_%s.csv", csv_prefix, method_tag, weight_tag, stamp))
      utils::write.csv(res, csv_file, row.names = FALSE); files$csv <- csv_file
      if (verbose) message("Saved CSV: ", normalizePath(csv_file, FALSE))
    }
    if (isTRUE(write_summary_csv) && length(attr(res, "community_sizes"))) {
      sm <- data.frame(
        com   = names(attr(res, "community_sizes")),
        size  = as.integer(attr(res, "community_sizes")),
        span  = as.integer(attr(res, "layer_span")[names(attr(res, "community_sizes"))]),
        stringsAsFactors = FALSE
      )
      sum_file <- file.path(results_dir, sprintf("%s_summary_%s_%s_%s.csv", csv_prefix, method_tag, weight_tag, stamp))
      utils::write.csv(sm, sum_file, row.names = FALSE); files$summary_csv <- sum_file
      if (verbose) message("Saved summary CSV: ", normalizePath(sum_file, FALSE))
    }
    attr(res, "files") <- files
    return(res)
  }
  if (is.null(memb) || !length(memb)) stop("Community detection failed.")

  # ---- 4) Relabel by size (C1 largest, …) ----
  memb_vec <- as.integer(memb)
  names(memb_vec) <- names(memb)
  if (relabel_by_size) {
    sizes <- sort(table(memb_vec), decreasing = TRUE)
    map   <- setNames(paste0("C", seq_along(sizes)), names(sizes))
    com_str <- unname(map[as.character(memb_vec)])
  } else {
    com_str <- paste0("C", memb_vec)
  }
  memb_df <- data.frame(actor = names(memb_vec), com = com_str, stringsAsFactors = FALSE)

  # ---- 5) Expand to each layer & filter ----
  actors_by_layer <- lapply(edge_list, function(ed) unique(c(ed$from, ed$to)))
  names(actors_by_layer) <- layers
  pieces <- lapply(names(actors_by_layer), function(ly) {
    a <- actors_by_layer[[ly]]
    if (!length(a)) return(NULL)
    df <- memb_df[memb_df$actor %in% a, , drop = FALSE]
    if (!nrow(df)) return(NULL)
    transform(df, layer = ly, stringsAsFactors = FALSE)
  })
  communities <- do.call(rbind, Filter(Negate(is.null), pieces))
  if (is.null(communities) || !nrow(communities))
    stop("No community assignments found for selected layers (after mapping actors).")
  communities$method <- method

  size_by <- tapply(communities$actor, communities$com, function(v) length(unique(v)))
  span_by <- tapply(communities$layer, communities$com, function(v) length(unique(v)))
  keep    <- names(size_by)[size_by >= min.actors & span_by >= min.layers]
  communities <- communities[communities$com %in% keep, , drop = FALSE]

  # Keep attributes (sorted by size desc for readability)
  if (length(keep)) {
    ord <- order(size_by[keep], decreasing = TRUE)
    size_by <- size_by[keep][ord]
    span_by <- span_by[keep][ord]
  } else {
    size_by <- integer(0); span_by <- integer(0)
  }
  attr(communities, "community_sizes") <- size_by
  attr(communities, "layer_span")      <- span_by

  if (verbose) {
    message("• Communities kept: ", length(size_by),
            " (min.actors=", min.actors, ", min.layers=", min.layers, ")")
  }

  # ---- 6) Save artifacts under results/ ----
  files <- list()
  .ensure_dir(results_dir)
  stamp <- format(Sys.time(), "%Y-%m-%d_%H%M%S")
  method_tag <- tolower(method); weight_tag <- tolower(edgeWeight)

  if (isTRUE(save_to_rds)) {
    if (is.null(rds_file)) rds_file <- sprintf("communities_%s_%s_%s.rds", method_tag, weight_tag, stamp)
    if (!.is_abs(rds_file)) rds_file <- file.path(results_dir, rds_file)
    saveRDS(communities, rds_file); files$rds <- rds_file
    if (verbose) message("Saved RDS: ", normalizePath(rds_file, FALSE))
  }
  if (isTRUE(write_csv)) {
    csv_file <- file.path(results_dir, sprintf("%s_%s_%s_%s.csv", csv_prefix, method_tag, weight_tag, stamp))
    utils::write.csv(communities, csv_file, row.names = FALSE); files$csv <- csv_file
    if (verbose) message("Saved CSV: ", normalizePath(csv_file, FALSE))
  }
  if (isTRUE(write_summary_csv) && length(size_by)) {
    sm <- data.frame(
      com  = names(size_by),
      size = as.integer(size_by),
      span = as.integer(span_by[names(size_by)]),
      stringsAsFactors = FALSE
    )
    sum_file <- file.path(results_dir, sprintf("%s_summary_%s_%s_%s.csv", csv_prefix, method_tag, weight_tag, stamp))
    utils::write.csv(sm, sum_file, row.names = FALSE); files$summary_csv <- sum_file
    if (verbose) message("Saved summary CSV: ", normalizePath(sum_file, FALSE))
  }

  attr(communities, "files") <- files
  communities
}


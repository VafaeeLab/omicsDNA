# ---------------------------------------
# 11 - detectCom
# ---------------------------------------

#' Detect communities on a multilayer network via an undirected, layer‑collapsed supra‑graph
#'
#' @description
#' This function performs community detection on a **single, layer‑collapsed supra‑graph**
#' obtained by aggregating the selected layers of a multilayer network. **By construction,
#' the supra‑graph is undirected**: endpoint order is canonicalised (A–B ≡ B–A) and any
#' direction in the input is **ignored**. Parallel edges across layers are aggregated by
#' either **presence count** (in how many layers an edge appears) or **weight sum** (sum of
#' edge weights over layers). A community partition is computed once on the supra‑graph
#' and then mapped back to each layer. Communities are labelled `C1`, `C2`, … with `C1`
#' corresponding to the largest community (when `relabel_by_size = TRUE`).
#'
#' Supported methods:
#' - `"glouvain"`: Louvain‑style modularity optimisation on the undirected supra‑graph
#'   (implemented via `igraph::cluster_louvain()`).
#' - `"louvain"`: **deprecated alias** of `"glouvain"`; accepted for backward compatibility
#'   with a warning.
#' - `"infomap"`: flow‑based clustering via `igraph::cluster_infomap()`.
#' - `"clique"`: k‑clique union (cliques of size ≥ `clique.k`) followed by connected
#'   components.
#' - `"abacus"`: when available, communities are read directly from
#'   `multinet::abacus_ml()` and then filtered by `min.actors`/`min.layers`.
#'
#' @details
#' **Input normalisation and scope**
#' - Edges are retrieved from `net` using `multinet::edges_ml()`. If both `from_layer` and
#'   `to_layer` exist, only *intra‑layer* edges are retained (`from_layer == to_layer`);
#'   otherwise a single `layer` column is used (or created as `"L1"` if absent).
#' - `layers` restricts analysis to a subset of layers; by default, all available layers
#'   after the intra‑layer restriction are used.
#' - Edge weights are taken from the first present of
#'   `c("weight","Weight","w","value","score")`; if none is found, weights default to `1`.
#'
#' **Supra‑graph construction (always undirected)**
#' - For each selected layer, edges are reduced to `(from, to, weight)` and endpoint order
#'   is canonicalised (`from = pmin(from, to)`, `to = pmax(from, to)`), so the supra‑graph
#'   is **undirected** even if the original multilayer network was directed.
#' - If `edgeWeight = "count"`, an edge’s weight is the number of **distinct layers** in
#'   which that undirected pair occurs. If `edgeWeight = "sum"`, weights are summed across
#'   layers. The supra‑graph is built with `igraph::graph_from_data_frame()` and simplified
#'   with `weight = "sum"`.
#'
#' **Community detection methods**
#' - `"glouvain"`: calls `igraph::cluster_louvain()` with supra‑graph edge weights.
#' - `"infomap"`: calls `igraph::cluster_infomap()` with supra‑graph edge weights.
#' - `"clique"`: enumerates all cliques of size ≥ `clique.k`; connects every pair of actors
#'   that co‑occur in any such clique (k‑clique union); returns connected components as
#'   communities. This can be expensive for dense graphs and larger `clique.k`.
#' - `"abacus"`: if available in your `multinet` build, community assignments are obtained
#'   from `abacus_ml()` and then **filtered** by `min.actors` and `min.layers` after any
#'   layer subsetting; the same relabelling and output logic applies.
#'
#' **Mapping back to layers & filtering**
#' - After detection on the supra‑graph, assignments are mapped back to each layer by
#'   retaining only actors that appear in that layer’s edges. Communities are then
#'   **filtered** to keep those with at least `min.actors` unique actors (pooled across the
#'   selected layers) and present in at least `min.layers` **distinct** layers.
#' - For `"clique"`, note the distinction: `clique.k` controls *how* communities are formed
#'   (cohesion threshold), whereas `min.actors`/`min.layers` control *which* detected
#'   communities are retained.
#'
#' **Reproducibility & labels**
#' - If `seed` is provided, the random seed is set before running glouvain/infomap (mild
#'   non‑determinism may persist depending on BLAS/threading). When `relabel_by_size = TRUE`,
#'   community IDs are recoded so that `C1` is largest, `C2` next, etc.
#'
#' **Persistence**
#' - When `save_to_rds = TRUE`, an RDS is written under `results_dir` with a timestamped
#'   name derived from the method and edge‑weight aggregation. If `write_csv = TRUE`, a CSV
#'   of assignments is also produced. If `write_summary_csv = TRUE`, a per‑community summary
#'   (size and layer span) is saved.
#'
#' @param net        A multilayer network compatible with `multinet::edges_ml()`.
#' @param method     Community detection method. One of
#'   `c("glouvain","louvain","infomap","clique","abacus")`. Default `"glouvain"`.
#'   (Note: `"louvain"` is a deprecated alias for `"glouvain"`.)
#' @param layers     Character vector naming layers to include. Default: all layers present
#'   after intra‑layer filtering.
#' @param edgeWeight How to aggregate parallel edges across layers:
#'   `"count"` (number of layers in which an edge appears) or `"sum"` (sum of edge weights
#'   across layers). Default `"count"`.
#' @param min.actors Minimum number of unique actors required for a community to be kept
#'   (after mapping back to layers). Default `15`.
#' @param min.layers Minimum number of **distinct layers** in which a community must appear
#'   to be kept. Default `2`.
#' @param clique.k   For `method = "clique"`, minimum clique size k (k‑clique union).
#'   Default `3`.
#' @param seed       Optional integer seed for glouvain/infomap runs. Default `NULL`.
#' @param relabel_by_size Logical; if `TRUE`, relabel communities by decreasing size so `C1`
#'   is largest. Default `TRUE`.
#' @param results_dir Output directory for saved artefacts. Default
#'   `getOption("mlnet.results_dir", "omicsDNA_results")`.
#' @param save_to_rds Logical; write the assignments as an RDS file. Default `TRUE`.
#' @param rds_file   Optional file name for the RDS output; relative paths are resolved
#'   under `results_dir`, absolute paths are respected. Default `NULL` (auto‑named).
#' @param write_csv  Logical; additionally write a CSV of assignments. Default `FALSE`.
#' @param csv_prefix Basename prefix for CSV files. Default `"communities"`.
#' @param write_summary_csv Logical; write a summary CSV (community size and layer span).
#'   Default `TRUE`.
#' @param verbose    Logical; print progress and file locations. Default `TRUE`.
#'
#' @return
#' A data frame with columns:
#' \itemize{
#'   \item `actor`: actor identifier,
#'   \item `com`: community label (`C1`, `C2`, …),
#'   \item `layer`: layer name,
#'   \item `method`: the method used (e.g., `"glouvain"`).
#' }
#' The return value carries attributes:
#' \itemize{
#'   \item `"community_sizes"` — named integer vector of community sizes (ordered),
#'   \item `"layer_span"`      — named integer vector of number of layers per community,
#'   \item `"files"`           — named list of output paths (if artefacts were written).
#' }
#'
#' @section Practical notes:
#' - Cross‑layer edges (if present) are ignored; only intra‑layer edges contribute to both
#'   the supra‑graph and per‑layer mappings.
#' - If no weight column is detected, all edges are treated as weight `1`.
#' - `"clique"` may be computationally expensive for dense graphs or larger `clique.k`.
#'
#' @examples
#' \dontrun{
#' # Glouvain (Louvain‑style) on all layers; edge weight = count of layers
#' comm <- detectCom(
#'   net,
#'   method      = "glouvain",
#'   edgeWeight  = "count",
#'   min.actors  = 15,
#'   min.layers  = 2,
#'   seed        = 1,
#'   save_to_rds = TRUE,
#'   write_csv   = FALSE,
#'   verbose     = TRUE
#' )
#'
#' # Backward‑compatible alias (deprecated): will warn and run as glouvain
#' comm2 <- detectCom(net, method = "louvain")
#'
#' # Clique‑based communities using k = 4; retain communities with >= 20 actors and in >= 2 layers
#' comm_clq <- detectCom(
#'   net,
#'   method      = "clique",
#'   clique.k    = 4,
#'   min.actors  = 20,
#'   min.layers  = 2,
#'   edgeWeight  = "sum",
#'   write_summary_csv = TRUE
#' )
#'
#' # Infomap with selected layers
#' sel <- c("Endothelium","Immune","N.LoH")
#' comm_infomap <- detectCom(net, method = "infomap", layers = sel, edgeWeight = "count")
#' }
#'
#' @seealso
#'   \code{\link{build_multiNet}} to assemble multilayer graphs,
#'   \code{\link{consensusEdges}} for robust per‑layer edges,
#'   \code{\link{add_network_attributes}} to enrich networks prior to clustering.
#'
#' @importFrom multinet edges_ml
#' @importFrom igraph graph_from_data_frame simplify cluster_louvain cluster_infomap cliques components E V membership ecount
#' @importFrom utils write.csv
#' @importFrom stats aggregate
#' @export
detectCom <- function(net,
                      method            = c("glouvain","louvain","infomap","clique","abacus"),
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

  # ---- method/args normalization ----
  method     <- match.arg(method)
  if (identical(method, "louvain")) {
    warning("`method = \"louvain\"` is deprecated; use `method = \"glouvain\"`.", call. = FALSE)
    method <- "glouvain"
  }
  edgeWeight <- match.arg(edgeWeight)
  stopifnot(is.numeric(min.actors), min.actors >= 1,
            is.numeric(min.layers), min.layers >= 1)
  if (method == "clique") stopifnot(is.numeric(clique.k), clique.k >= 2)

  # ---- helpers ----
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

  # ---- 2) Supra-graph (always undirected) ----
  if (isTRUE(verbose)) {
    message("Building a layer-collapsed supra-graph as UNDIRECTED: ",
            "input edge directions (if any) are ignored.")
  }
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
  if (verbose) message("• Built supra-graph: |V|=", length(igraph::V(g)),
                       " |E|=", igraph::ecount(g),
                       "  (edgeWeight=", edgeWeight, ", method=", method, ")")

  # ---- 3) Community detection ----
  memb <- NULL
  if (method == "glouvain") {
    co <- igraph::cluster_louvain(g, weights = igraph::E(g)$weight)
    memb <- igraph::membership(co)

  } else if (method == "infomap") {
    co <- igraph::cluster_infomap(g, e.weights = igraph::E(g)$weight)
    memb <- igraph::membership(co)

  } else if (method == "clique") {
    if (isTRUE(verbose))
      message("Clique-based detection: using clique.k = ", clique.k,
              " (increase for more stringent, denser communities).")
    cls <- igraph::cliques(g, min = clique.k)
    if (!length(cls)) stop("No cliques of size >= ", clique.k, " found.")

    # Build k-clique union graph: connect pairs co-occurring in any k-clique
    pairs_df <- do.call(rbind, lapply(cls, function(cl) {
      vs <- igraph::V(g)$name[cl]
      if (length(vs) < 2) return(NULL)
      m <- utils::combn(vs, 2)
      data.frame(from = m[1,], to = m[2,], stringsAsFactors = FALSE)
    }))
    pairs_df <- unique(pairs_df)
    h <- igraph::graph_from_data_frame(pairs_df, directed = FALSE, vertices = igraph::V(g)$name)
    comp <- igraph::components(h)
    if (isTRUE(verbose)) {
      message("k-clique union produced ", comp$no, " components before size/layer filtering.")
    }
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

    # save artifacts
    files <- list()
    .ensure_dir(results_dir); stamp <- format(Sys.time(), "%Y-%m-%d_%H%M%S")
    method_tag <- "glouvain" %in% tolower("abacus") # dummy to avoid R CMD check note
    method_tag <- "abacus"  # label consistently
    weight_tag <- tolower(edgeWeight)

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
    if (isTRUE(write_summary_csv) && length(size_by)) {
      sm <- data.frame(
        com   = names(size_by),
        size  = as.integer(size_by),
        span  = as.integer(span_by[names(size_by)]),
        stringsAsFactors = FALSE
      )
      sum_file <- file.path(results_dir, sprintf("%s_summary_%s_%s_%s.csv", csv_prefix, method_tag, weight_tag, stamp))
      utils::write.csv(sm, sum_file, row.names = FALSE); files$summary_csv <- sum_file
      if (isTRUE(verbose)) message("Saved summary CSV: ", normalizePath(sum_file, FALSE))
    }

    attr(res, "community_sizes") <- if (length(keep)) as.integer(size_by) else integer(0)
    attr(res, "layer_span")      <- if (length(keep)) as.integer(span_by[names(size_by)]) else integer(0)
    attr(res, "files")           <- files
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

  # ---- 6) Save artifacts under results_dir ----
  files <- list()
  .ensure_dir(results_dir)
  stamp <- format(Sys.time(), "%Y-%m-%d_%H%M%S")
  method_tag <- tolower(method)
  method_tag <- sub("^louvain$", "glouvain", method_tag)  # normalize alias in filenames
  weight_tag <- tolower(edgeWeight)

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

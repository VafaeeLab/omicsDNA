#' Community detection in multilayer networks (multinet)
#'
#' @md
#' @description
#' Detect communities in a multilayer network of class [multinet::ml.network].
#'
#' This function supports **two analysis modes**:
#'
#' 1. **Multilayer-native mode** (`supra_graph = FALSE`, default): calls a
#'    multilayer community algorithm from the **`multinet`** package directly:
#'    - [multinet::glouvain_ml()]
#'    - [multinet::infomap_ml()]
#'    - [multinet::clique_percolation_ml()]
#'    - [multinet::abacus_ml()]
#'
#' 2. **Supra-graph (layer-collapsed) mode** (`supra_graph = TRUE`): collapses the
#'    chosen layers into a single undirected graph (a *supra-graph*), runs a
#'    community method from **`igraph`**, then maps community membership back to
#'    each layer.
#'    - `method = "glouvain"` -> [igraph::cluster_louvain()]
#'    - `method = "infomap"`  -> [igraph::cluster_infomap()]
#'    - `method = "clique"`   -> clique-union components (based on
#'      [igraph::cliques()] + [igraph::components()])
#'
#' In both modes, the result is standardized to a **tidy, edge-independent**
#' membership table (one row per *actor–layer* membership) and optionally written
#' to disk.
#'
#' @details
#' ## Output format (always the same)
#' The returned object is a data frame with columns:
#' - `actor`: actor identifier (e.g., gene name)
#' - `com`: community id (e.g., `"C1"`, `"C2"`, ... if relabeled; otherwise the
#'   algorithm’s raw ids coerced to character)
#' - `layer`: layer name
#' - `method`: the method used (`"glouvain"`, `"infomap"`, `"clique"`, `"abacus"`)
#'
#' The object also carries attributes:
#' - `community_sizes`: named integer vector = number of **unique actors** per community
#' - `layer_span`: named integer vector = number of **unique layers** spanned by each community
#' - `files`: list with file paths that were written (`rds`, `csv`, `summary_csv`)
#'
#' ## Multilayer-native mode (`supra_graph = FALSE`)
#' In this mode, `detectCom()` delegates community detection to the corresponding
#' `multinet::*_ml()` function using your parameters (e.g., `gamma`, `omega` for
#' multilayer generalized Louvain).
#'
#' **Important practical note:** even in multilayer-native mode, `detectCom()`
#' may **change** the final output compared to the raw `multinet::*_ml()` result
#' because it can apply **post-filtering** (`min.actors`, `min.layers`) and
#' **optional relabeling** (`relabel_by_size`).
#'
#' If you want the wrapper output to be as close as possible to the direct
#' multilayer algorithm output, use:
#' - `min.actors = 1`, `min.layers = 1`, and
#' - `relabel_by_size = FALSE`.
#'
#' ## Supra-graph mode (`supra_graph = TRUE`)
#' In supra-graph mode, the function:
#' 1. Extracts intra-layer edges from the multilayer network.
#' 2. Builds a single undirected graph where vertices are actors and edges are
#'    aggregated across layers.
#' 3. Runs an `igraph` community method on this collapsed graph.
#' 4. Maps each actor’s community membership back to the layers in which the
#'    actor appears.
#'
#' Edge aggregation across layers is controlled by `edgeWeight`:
#' - `"count"`: supra-edge weight = number of layers in which the edge appears
#'   at least once.
#' - `"sum"`: supra-edge weight = sum of per-layer edge weights across layers
#'   (when available; otherwise edges default to weight 1).
#'
#' **Important:** supra-graph mode is a different analysis pipeline than
#' multilayer-native mode. You should not expect identical partitions between
#' these modes, even with similar settings.
#'
#' ## Filtering and relabeling (applies to all modes)
#' After community detection, communities are retained only if:
#' - number of **unique actors** in the community `>= min.actors`, and
#' - number of **distinct layers** represented in the community `>= min.layers`.
#'
#' If `relabel_by_size = TRUE`, communities are renamed to `"C1"`, `"C2"`, ...
#' in decreasing order of community size (ties broken arbitrarily).
#'
#' @param net A multilayer network of class [multinet::ml.network].
#'
#' @param method Community detection method. One of:
#' - `"glouvain"` (default): generalized Louvain
#' - `"infomap"`: Infomap community detection
#' - `"clique"`: clique-based communities
#'   - multilayer-native mode: k-clique percolation via `multinet`
#'   - supra-graph mode: clique-union components via `igraph`
#' - `"abacus"`: ABACUS multilayer community detection (multilayer-native only)
#' - `"louvain"` is accepted as an alias of `"glouvain"` (deprecated; converted
#'   internally with a warning).
#'
#' @param supra_graph Logical.
#' - `FALSE` (default): **multilayer-native** algorithms from `multinet` are used.
#' - `TRUE`: layers are collapsed into a single supra-graph and an `igraph`
#'   method is used (only for `method = "glouvain"`, `"infomap"`, `"clique"`).
#'
#' @param layers Optional character vector of layer names to include.
#' - If `NULL`, all layers are used.
#' - If provided, only those layers (intersected with available layers) are used.
#'
#' @param min.actors Integer (>= 1). Minimum number of **unique actors** required
#' for a community to be kept. Applied as a post-filter to all methods.
#'
#' @param min.layers Integer (>= 1). Minimum number of **distinct layers**
#' required for a community to be kept. Applied as a post-filter to all methods.
#'
#' @param relabel_by_size Logical.
#' - `TRUE` (default): relabel communities as `"C1"`, `"C2"`, ... ordered by
#'   decreasing community size.
#' - `FALSE`: keep the algorithm’s original community ids (converted to character).
#'
#' @param seed Optional integer. If supplied, `set.seed(seed)` is called before
#' running algorithms that may be stochastic.
#' - In multilayer-native mode: used for `"glouvain"` and `"infomap"`.
#' - In supra-graph mode: used before the `igraph` clustering step.
#'
#' @param glouvain_gamma Numeric (> 0). Resolution parameter passed to
#' [multinet::glouvain_ml()] **only when** `method = "glouvain"` and
#' `supra_graph = FALSE`. Ignored otherwise.
#'
#' @param glouvain_omega Numeric (>= 0). Inter-layer coupling parameter passed to
#' [multinet::glouvain_ml()] **only when** `method = "glouvain"` and
#' `supra_graph = FALSE`. Ignored otherwise.
#'
#' @param clique.k Integer (>= 3). Clique size threshold.
#' - Multilayer-native mode (`supra_graph = FALSE`): forwarded as `k` to
#'   [multinet::clique_percolation_ml()].
#' - Supra-graph mode (`supra_graph = TRUE`): used as the minimum clique size in
#'   [igraph::cliques()].
#'
#' @param clique.m Integer (>= 1). Multilayer-only clique parameter forwarded as
#' `m` to [multinet::clique_percolation_ml()] **only when**
#' `method = "clique"` and `supra_graph = FALSE`. Ignored in supra-graph mode.
#'
#' @param infomap_overlapping Logical. Passed to [multinet::infomap_ml()] as
#' `overlapping` **only when** `method = "infomap"` and `supra_graph = FALSE`.
#' Ignored in supra-graph mode.
#'
#' @param infomap_directed Logical. Passed to [multinet::infomap_ml()] as
#' `directed` **only when** `method = "infomap"` and `supra_graph = FALSE`.
#' Ignored in supra-graph mode.
#'
#' @param infomap_self_links Logical. Passed to [multinet::infomap_ml()] as
#' `self.links` **only when** `method = "infomap"` and `supra_graph = FALSE`.
#' Ignored in supra-graph mode.
#'
#' @param edgeWeight How to aggregate edge weights across layers in supra-graph mode.
#' One of `"count"` or `"sum"`.
#' - `"count"`: weight = number of layers in which the edge appears at least once
#' - `"sum"`: weight = sum of per-layer weights across layers (when available)
#' Ignored when `supra_graph = FALSE`.
#'
#' @param results_dir Output directory for optional files (RDS/CSV/summary CSV).
#' If it does not exist, it is created. Defaults to
#' `getOption("mlnet.results_dir", "omicsDNA_results")`.
#'
#' @param save_to_rds Logical. If `TRUE`, saves the community membership table as
#' an `.rds` file in `results_dir` (unless `rds_file` is absolute).
#'
#' @param rds_file Optional character filename for the RDS output.
#' - If `NULL`, a timestamped filename is generated automatically.
#' - If relative, it is created inside `results_dir`.
#' - If absolute, it is used as-is.
#'
#' @param write_csv Logical. If `TRUE`, writes the membership table to a CSV file.
#'
#' @param csv_prefix Character prefix used to build output CSV filenames.
#'
#' @param write_summary_csv Logical. If `TRUE`, writes a small summary CSV with:
#' `com`, `size` (unique actors), `span` (layers).
#'
#' @param verbose Logical. If `TRUE`, prints progress messages and file paths.
#'
#' @return A data frame with columns `actor`, `com`, `layer`, `method`, plus
#' attributes `community_sizes`, `layer_span`, and `files`.
#'
#' @examples
#' # Build a tiny 2-layer multiplex from igraph graphs
#' g1 <- igraph::graph_from_data_frame(
#'   data.frame(from = c("a","b","c"), to = c("b","c","d"), w_ = c(1,1,1)),
#'   directed = FALSE
#' )
#' g2 <- igraph::graph_from_data_frame(
#'   data.frame(from = c("a","a","d"), to = c("c","d","e"), w_ = c(1,1,1)),
#'   directed = FALSE
#' )
#'
#' net <- multinet::ml_empty()
#' multinet::add_igraph_layer_ml(net, g1, name = "E1")
#' multinet::add_igraph_layer_ml(net, g2, name = "E2")
#'
#' # Multilayer-native GLouvain (supra_graph = FALSE)
#' comm_ml <- detectCom(
#'   net,
#'   method          = "glouvain",
#'   supra_graph     = FALSE,
#'   glouvain_gamma  = 1,
#'   glouvain_omega  = 1,
#'   min.actors      = 1,
#'   min.layers      = 1,
#'   relabel_by_size = FALSE,
#'   seed            = 1,
#'   save_to_rds     = FALSE,
#'   write_summary_csv = FALSE,
#'   verbose         = FALSE
#' )
#'
#' # Supra-graph GLouvain (supra_graph = TRUE)
#' comm_supra <- detectCom(
#'   net,
#'   method          = "glouvain",
#'   supra_graph     = TRUE,
#'   edgeWeight      = "count",
#'   min.actors      = 1,
#'   min.layers      = 1,
#'   relabel_by_size = FALSE,
#'   seed            = 1,
#'   save_to_rds     = FALSE,
#'   write_summary_csv = FALSE,
#'   verbose         = FALSE
#' )
#'
#' @importFrom multinet layers_ml ml_empty add_igraph_layer_ml
#' @importFrom multinet edges_ml glouvain_ml infomap_ml clique_percolation_ml abacus_ml
#' @importFrom igraph graph_from_data_frame simplify cliques components cluster_louvain cluster_infomap V E
#' @importFrom stats aggregate
#' @export
detectCom <- function(
    net,
    method              = c("glouvain","louvain","infomap","clique","abacus"),
    supra_graph         = FALSE,
    layers              = NULL,
    min.actors          = 15,
    min.layers          = 2,
    relabel_by_size     = TRUE,
    seed                = NULL,
    glouvain_gamma      = 1,
    glouvain_omega      = 1,
    clique.k            = 3,
    clique.m            = 1,
    infomap_overlapping = FALSE,
    infomap_directed    = FALSE,
    infomap_self_links  = TRUE,
    edgeWeight          = c("count","sum"),
    results_dir         = getOption("mlnet.results_dir","omicsDNA_results"),
    save_to_rds         = TRUE,
    rds_file            = NULL,
    write_csv           = FALSE,
    csv_prefix          = "communities",
    write_summary_csv   = TRUE,
    verbose             = TRUE
) {

  # ---- normalize and validate ----
  method <- match.arg(method)
  if (identical(method, "louvain")) {
    warning("`method = \"louvain\"` is deprecated; use `method = \"glouvain\"`.", call. = FALSE)
    method <- "glouvain"
  }
  edgeWeight <- match.arg(edgeWeight)

  stopifnot(is.numeric(min.actors), length(min.actors) == 1L, min.actors >= 1,
            is.numeric(min.layers), length(min.layers) == 1L, min.layers >= 1)

  stopifnot(is.logical(relabel_by_size), length(relabel_by_size) == 1L,
            is.logical(supra_graph), length(supra_graph) == 1L,
            is.logical(verbose), length(verbose) == 1L)

  stopifnot(is.numeric(glouvain_gamma), length(glouvain_gamma) == 1L, glouvain_gamma > 0,
            is.numeric(glouvain_omega), length(glouvain_omega) == 1L, glouvain_omega >= 0)

  if (method == "clique") {
    stopifnot(is.numeric(clique.k), length(clique.k) == 1L, clique.k >= 3,
              is.numeric(clique.m), length(clique.m) == 1L, clique.m >= 1)
  }

  stopifnot(is.logical(infomap_overlapping), length(infomap_overlapping) == 1L,
            is.logical(infomap_directed), length(infomap_directed) == 1L,
            is.logical(infomap_self_links), length(infomap_self_links) == 1L)

  # ---- helpers ----
  .ensure_dir <- function(d) if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  .is_abs     <- function(p) grepl("^(/|[A-Za-z]:[\\/])", p)
  .pick       <- function(cands, nms) { z <- cands[cands %in% nms]; if (length(z)) z[1] else NA_character_ }

  # Subset layers by rebuilding a temporary ml.network (so the detection itself
  # only sees the requested layers).
  .subset_layers_ml <- function(n, layers_vec) {
    if (is.null(layers_vec)) return(n)
    avail <- try(multinet::layers_ml(n), silent = TRUE)
    if (inherits(avail, "try-error") || is.null(avail)) return(n)

    layers_vec <- intersect(as.character(layers_vec), as.character(avail))
    if (!length(layers_vec)) stop("No layers available after intersection with requested `layers`.")

    # no-op if already all
    if (length(layers_vec) == length(avail) && all(sort(layers_vec) == sort(avail))) return(n)

    n2 <- multinet::ml_empty()
    for (ly in layers_vec) {
      g <- igraph::as.igraph(n, layers = ly, merge.actors = TRUE, all.actors = FALSE)
      multinet::add_igraph_layer_ml(n2, g, name = ly)
    }
    n2
  }

  # Build and return the standardized result + write optional files.
  .finalize_result <- function(res_df, size_by, span_by, method_used, mode_tag, weight_tag = NULL) {
    res_df$method <- method_used

    files <- list()
    .ensure_dir(results_dir)
    stamp <- format(Sys.time(), "%Y-%m-%d_%H%M%S")

    # Filenames include mode to avoid accidental overwrites when comparing
    # multilayer-native vs supra-graph runs.
    if (is.null(weight_tag)) {
      base <- sprintf("%s_%s_%s_%s", csv_prefix, tolower(method_used), mode_tag, stamp)
    } else {
      base <- sprintf("%s_%s_%s_%s_%s", csv_prefix, tolower(method_used), mode_tag, tolower(weight_tag), stamp)
    }

    if (isTRUE(save_to_rds)) {
      if (is.null(rds_file)) {
        rds_file_use <- paste0(base, ".rds")
      } else {
        rds_file_use <- rds_file
      }
      if (!.is_abs(rds_file_use)) rds_file_use <- file.path(results_dir, rds_file_use)
      saveRDS(res_df, rds_file_use)
      files$rds <- rds_file_use
      if (verbose) message("Saved RDS: ", normalizePath(rds_file_use, winslash = "/", mustWork = FALSE))
    }

    if (isTRUE(write_csv)) {
      csv_file <- file.path(results_dir, paste0(base, ".csv"))
      utils::write.csv(res_df, csv_file, row.names = FALSE)
      files$csv <- csv_file
      if (verbose) message("Saved CSV: ", normalizePath(csv_file, winslash = "/", mustWork = FALSE))
    }

    if (isTRUE(write_summary_csv) && length(size_by)) {
      sm <- data.frame(
        com  = names(size_by),
        size = as.integer(size_by),
        span = as.integer(span_by[names(size_by)]),
        stringsAsFactors = FALSE
      )
      sum_file <- file.path(results_dir, paste0(base, "_summary.csv"))
      utils::write.csv(sm, sum_file, row.names = FALSE)
      files$summary_csv <- sum_file
      if (verbose) message("Saved summary CSV: ", normalizePath(sum_file, winslash = "/", mustWork = FALSE))
    }

    attr(res_df, "community_sizes") <- if (length(size_by)) as.integer(size_by) else integer(0)
    attr(res_df, "layer_span")      <- if (length(size_by)) as.integer(span_by[names(size_by)]) else integer(0)
    attr(res_df, "files")           <- files

    res_df
  }

  # ---- decide mode ----
  # abacus has no supra-graph variant here; always run multilayer-native.
  use_ml_native <- (!isTRUE(supra_graph)) || method == "abacus"

  # ========================================================================
  # 1) Multilayer-native mode: delegate to multinet::*_ml()
  # ========================================================================
  if (use_ml_native) {
    if (verbose) {
      msg <- paste0("Running multilayer-native community detection (multinet), method = \"", method, "\".")
      if (isTRUE(supra_graph) && method == "abacus") {
        msg <- paste0(msg, " Note: `supra_graph = TRUE` is ignored for method = \"abacus\".")
      }
      if (!isTRUE(supra_graph)) {
        msg <- paste0(msg, " Note: `edgeWeight` is ignored when supra_graph = FALSE.")
      }
      message(msg)
    }

    net_use <- .subset_layers_ml(net, layers)

    if (!is.null(seed) && method %in% c("glouvain","infomap")) {
      set.seed(as.integer(seed))
    }

    cm_raw <- switch(method,
                     glouvain = multinet::glouvain_ml(net_use, gamma = glouvain_gamma, omega = glouvain_omega),
                     infomap  = multinet::infomap_ml(net_use,
                                                     overlapping = infomap_overlapping,
                                                     directed    = infomap_directed,
                                                     self.links  = infomap_self_links),
                     clique   = multinet::clique_percolation_ml(net_use, k = clique.k, m = clique.m),
                     abacus   = multinet::abacus_ml(net_use, min.actors = min.actors, min.layers = min.layers),
                     stop("Unsupported method: ", method)
    )

    if (!is.data.frame(cm_raw)) {
      cm_raw <- try(as.data.frame(cm_raw, stringsAsFactors = FALSE), silent = TRUE)
      if (inherits(cm_raw, "try-error") || is.null(cm_raw))
        stop("Community detection did not return a data frame.")
    }

    if (!nrow(cm_raw)) {
      res     <- data.frame(actor = character(), com = character(), layer = character(), stringsAsFactors = FALSE)
      size_by <- integer(0)
      span_by <- integer(0)
      return(.finalize_result(res, size_by, span_by, method, mode_tag = "ml"))
    }

    nm <- names(cm_raw)
    actor_col <- .pick(c("actor","Actor","node","Node","vertex","Vertex"), nm)
    layer_col <- .pick(c("layer","Layer"), nm)
    com_col   <- .pick(c("com","community","cluster","cid","community_id","Community"), nm)
    if (any(is.na(c(actor_col, layer_col, com_col)))) {
      stop("Unexpected output format from multinet community algorithm: could not detect actor/layer/community columns.")
    }

    res <- data.frame(
      actor = as.character(cm_raw[[actor_col]]),
      com   = as.character(cm_raw[[com_col]]),
      layer = as.character(cm_raw[[layer_col]]),
      stringsAsFactors = FALSE
    )

    # Safety: filter to requested layers if the upstream network could not be subset
    if (!is.null(layers)) {
      layers_use <- intersect(as.character(layers), unique(res$layer))
      res <- res[res$layer %in% layers_use, , drop = FALSE]
    }

    if (!nrow(res)) {
      if (verbose) message("No community assignments remain after applying `layers` filter.")
      res     <- data.frame(actor = character(), com = character(), layer = character(), stringsAsFactors = FALSE)
      size_by <- integer(0)
      span_by <- integer(0)
      return(.finalize_result(res, size_by, span_by, method, mode_tag = "ml"))
    }

    # Post-filter by size/span (even if abacus already filtered)
    size_by <- tapply(res$actor, res$com, function(v) length(unique(v)))
    span_by <- tapply(res$layer, res$com, function(v) length(unique(v)))
    keep    <- names(size_by)[size_by >= min.actors & span_by >= min.layers]
    res     <- res[res$com %in% keep, , drop = FALSE]

    if (!length(keep) || !nrow(res)) {
      if (verbose) message("No communities satisfy min.actors/min.layers after filtering.")
      res     <- data.frame(actor = character(), com = character(), layer = character(), stringsAsFactors = FALSE)
      size_by <- integer(0)
      span_by <- integer(0)
      return(.finalize_result(res, size_by, span_by, method, mode_tag = "ml"))
    }

    if (isTRUE(relabel_by_size)) {
      ord <- order(size_by[keep], decreasing = TRUE)
      map <- setNames(paste0("C", seq_along(ord)), names(size_by[keep])[ord])
      res$com <- unname(map[res$com])
      size_by <- size_by[keep][ord]; names(size_by) <- unname(map[names(size_by)])
      span_by <- span_by[keep][ord]; names(span_by) <- names(size_by)
    } else {
      size_by <- size_by[keep]
      span_by <- span_by[keep]
      ord <- order(size_by, decreasing = TRUE)
      size_by <- size_by[ord]
      span_by <- span_by[ord]
    }

    return(.finalize_result(res, size_by, span_by, method, mode_tag = "ml"))
  }

  # ========================================================================
  # 2) Supra-graph mode: collapse layers and run igraph algorithms
  # ========================================================================

  if (!method %in% c("glouvain","infomap","clique")) {
    stop("`supra_graph = TRUE` is only supported for methods 'glouvain', 'infomap', and 'clique'.")
  }

  if (verbose) {
    message("Running supra-graph community detection (layer-collapsed), method = \"", method, "\".")
    if (method == "glouvain") message("Note: glouvain_gamma/glouvain_omega are ignored in supra-graph mode.")
    if (method == "infomap")  message("Note: infomap_* parameters are ignored in supra-graph mode.")
    if (method == "clique")   message("Note: clique.m is ignored in supra-graph mode.")
  }

  # Pull multilayer edge table
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

  # Identify endpoint columns
  pairs <- list(
    c("from_actor","to_actor"), c("from","to"), c("source","target"),
    c("actor1","actor2"), c("i","j"), c("v1","v2")
  )
  a_col <- b_col <- NA_character_
  for (p in pairs) {
    if (all(p %in% nm)) { a_col <- p[1]; b_col <- p[2]; break }
  }
  if (is.na(a_col)) {
    layer_like <- c("layer","Layer","from_layer","to_layer","l1","l2")
    char_cols  <- which(vapply(Eall, function(x) is.character(x) || is.factor(x), logical(1)))
    char_cols  <- setdiff(char_cols, match(layer_like, nm, nomatch = 0))
    if (length(char_cols) < 2) stop("Could not identify two endpoint columns in the global edge table.")
    a_col <- nm[char_cols[1]]
    b_col <- nm[char_cols[2]]
  }

  # Layer columns; keep intra-layer edges
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
    layers_use <- sort(unique(Eall$layer))
  } else {
    layers_use <- intersect(as.character(layers), unique(Eall$layer))
  }
  if (!length(layers_use)) stop("No layers available after intersection with requested `layers`.")

  # Choose weight column if present
  w_col_global <- .pick(c("weight","Weight","w","w_","value","score","correlation"), nm)

  edge_list <- lapply(layers_use, function(ly) {
    ed <- Eall[Eall$layer == ly, , drop = FALSE]
    if (!nrow(ed)) return(ed[FALSE, , drop = FALSE])
    from <- as.character(ed[[a_col]])
    to   <- as.character(ed[[b_col]])
    w    <- if (!is.na(w_col_global) && w_col_global %in% names(ed)) as.numeric(ed[[w_col_global]]) else rep(1, length(from))
    data.frame(layer = ly, from = from, to = to, weight = w, stringsAsFactors = FALSE)
  })
  names(edge_list) <- layers_use

  edges_all <- do.call(rbind, edge_list)
  if (is.null(edges_all) || !nrow(edges_all)) stop("No edges found in the selected layers.")

  # Canonicalize for undirected supra-graph
  aa <- pmin(edges_all$from, edges_all$to)
  bb <- pmax(edges_all$from, edges_all$to)
  edges_all$from <- aa
  edges_all$to   <- bb

  # Aggregate weights across layers
  if (edgeWeight == "count") {
    per_layer_unique <- unique(edges_all[, c("layer","from","to")])
    agg <- stats::aggregate(rep(1L, nrow(per_layer_unique)),
                            by = per_layer_unique[c("from","to")], FUN = sum)
    names(agg)[3] <- "weight"
  } else {
    agg <- stats::aggregate(weight ~ from + to, data = edges_all, FUN = sum, na.rm = TRUE)
  }

  g <- igraph::graph_from_data_frame(agg, directed = FALSE)
  if (igraph::ecount(g) == 0) stop("Supra-graph has no edges after aggregation.")
  g <- igraph::simplify(g, edge.attr.comb = list(weight = "sum", "ignore"))

  if (!is.null(seed)) set.seed(as.integer(seed))
  if (verbose) message("Built supra-graph: |V|=", length(igraph::V(g)), " |E|=", igraph::ecount(g),
                       " (edgeWeight=", edgeWeight, ")")

  memb <- NULL
  if (method == "glouvain") {
    co <- igraph::cluster_louvain(g, weights = igraph::E(g)$weight)
    memb <- igraph::membership(co)
  } else if (method == "infomap") {
    co <- igraph::cluster_infomap(g, e.weights = igraph::E(g)$weight)
    memb <- igraph::membership(co)
  } else if (method == "clique") {
    cls <- igraph::cliques(g, min = clique.k)
    if (!length(cls)) stop("No cliques of size >= ", clique.k, " found.")

    pairs_df <- do.call(rbind, lapply(cls, function(cl) {
      vs <- igraph::V(g)$name[cl]
      if (length(vs) < 2) return(NULL)
      m <- utils::combn(vs, 2)
      data.frame(from = m[1, ], to = m[2, ], stringsAsFactors = FALSE)
    }))
    pairs_df <- unique(pairs_df)

    h <- igraph::graph_from_data_frame(pairs_df, directed = FALSE, vertices = igraph::V(g)$name)
    comp <- igraph::components(h)
    memb <- comp$membership
    names(memb) <- igraph::V(h)$name
  }

  if (is.null(memb) || !length(memb)) stop("Community detection failed.")

  # Standardize membership to a character community id (numeric labels).
  memb_vec <- as.integer(memb)
  names(memb_vec) <- names(memb)
  memb_df <- data.frame(actor = names(memb_vec), com = as.character(memb_vec), stringsAsFactors = FALSE)

  # Map back to each layer
  actors_by_layer <- lapply(edge_list, function(ed) unique(c(ed$from, ed$to)))
  pieces <- lapply(names(actors_by_layer), function(ly) {
    a <- actors_by_layer[[ly]]
    if (!length(a)) return(NULL)
    df <- memb_df[memb_df$actor %in% a, , drop = FALSE]
    if (!nrow(df)) return(NULL)
    transform(df, layer = ly, stringsAsFactors = FALSE)
  })
  res <- do.call(rbind, Filter(Negate(is.null), pieces))
  if (is.null(res) || !nrow(res)) stop("No community assignments found after mapping back to layers.")

  # Filter by size/span across layers
  size_by <- tapply(res$actor, res$com, function(v) length(unique(v)))
  span_by <- tapply(res$layer, res$com, function(v) length(unique(v)))
  keep    <- names(size_by)[size_by >= min.actors & span_by >= min.layers]
  res     <- res[res$com %in% keep, , drop = FALSE]

  if (!length(keep) || !nrow(res)) {
    if (verbose) message("No communities satisfy min.actors/min.layers after filtering.")
    res     <- data.frame(actor = character(), com = character(), layer = character(), stringsAsFactors = FALSE)
    size_by <- integer(0)
    span_by <- integer(0)
    return(.finalize_result(res, size_by, span_by, method, mode_tag = "supra", weight_tag = edgeWeight))
  }

  # Keep summary stats only for retained communities and order by size
  size_by <- size_by[keep]
  span_by <- span_by[keep]
  ord     <- order(size_by, decreasing = TRUE)
  size_by <- size_by[ord]
  span_by <- span_by[ord]

  # Optionally relabel kept communities to C1, C2, ... (consecutive)
  if (isTRUE(relabel_by_size)) {
    map <- setNames(paste0("C", seq_along(size_by)), names(size_by))
    res$com  <- unname(map[res$com])
    names(size_by) <- unname(map[names(size_by)])
    names(span_by) <- names(size_by)
  }

  .finalize_result(res, size_by, span_by, method, mode_tag = "supra", weight_tag = edgeWeight)
}

# ---------------------------------------------------------------------------
# detectCom.R
# Strict community detection wrapper for multilayer networks (multinet)
# ---------------------------------------------------------------------------

#' Strict community detection wrapper for multilayer networks (multinet)
#'
#' @md
#' @description
#' `detectCom()` supports **two analysis modes** for a multilayer network of class
#' [multinet::ml.network].
#'
#' ## 1) Strict multilayer-native passthrough (`supra_graph = FALSE`, default)
#' In this mode, `detectCom()` is a **thin wrapper** around `multinet::*_ml()` and
#' returns the **raw multinet output object unchanged**:
#' - **no post-filtering**
#' - **no relabeling**
#' - **no column renaming**
#' - **no added `method` column**
#' - **no attributes attached**
#'
#' Mappings:
#' - `method = "glouvain"` -> [multinet::glouvain_ml()]
#' - `method = "infomap"`  -> [multinet::infomap_ml()]
#' - `method = "clique"`   -> [multinet::clique_percolation_ml()]
#' - `method = "abacus"`   -> [multinet::abacus_ml()]
#'
#' Notes for strict mode:
#' - `layers` is **ignored** (to guarantee identity with calling `multinet::*_ml(net, ...)`).
#' - `min.actors`, `min.layers`, and `relabel_by_size` are **not applied** as wrapper
#'   post-processing to the returned object.
#' - For `method="abacus"`, `min.actors` and `min.layers` are passed through to
#'   [multinet::abacus_ml()] because they are part of that algorithm’s inputs.
#'
#' ## 2) Supra-graph mode (`supra_graph = TRUE`)
#' In this mode, selected layers are collapsed into an undirected actor–actor
#' supra-graph, an **igraph** community method is applied, and memberships are
#' mapped back to layers where each actor appears.
#'
#' Mappings:
#' - `method = "glouvain"` -> [igraph::cluster_louvain()]
#' - `method = "infomap"`  -> [igraph::cluster_infomap()]
#' - `method = "clique"`   -> clique-union components derived from
#'   [igraph::cliques()] + [igraph::components()]
#'
#' Post-filtering (`min.actors`, `min.layers`) and optional relabeling
#' (`relabel_by_size`) are applied **only in supra-graph mode**.
#'
#' @details
#' ## Summary CSV output (matches your attached format)
#' If `write_summary_csv = TRUE`, a summary CSV is written with columns:
#' - `com`: community id
#' - `size`: number of **unique actors** in the community
#' - `span`: number of **unique layers** represented by that community
#'
#' The summary is ordered by decreasing `size`.
#'
#' - In **strict passthrough mode** (`supra_graph=FALSE`), the summary is computed
#'   from a **coerced copy** of the raw multinet output (when possible), but the
#'   returned object remains unchanged.
#' - If `relabel_by_size = TRUE`, the **summary file only** uses `C1, C2, ...`
#'   ordered by decreasing `size` (this does not modify the returned object).
#'
#' ## Edge weights in supra-graph mode
#' When `supra_graph = TRUE`, inter-layer aggregation uses `edgeWeight`:
#' - `"count"`: supra-edge weight = number of layers where the edge appears
#' - `"sum"`: supra-edge weight = sum of per-layer weights (when available)
#'
#' ## Backward-compatible aliases (optional convenience)
#' This implementation also accepts legacy names passed via `...`:
#' - `glouvain_gamma` -> `gamma`, `glouvain_omega` -> `omega`
#' - `clique.k` -> `k`, `clique.m` -> `m`
#' - `infomap_overlapping` -> `overlapping`, `infomap_directed` -> `directed`,
#'   `infomap_self_links` -> `self.links`
#'
#' @param net A multilayer network of class [multinet::ml.network].
#'
#' @param method Community detection method. One of `"glouvain"`, `"louvain"` (alias),
#' `"infomap"`, `"clique"`, `"abacus"`.
#'
#' @param supra_graph Logical. Default `FALSE`.
#' - `FALSE`: strict passthrough to `multinet::*_ml()` (returned object unchanged)
#' - `TRUE`: supra-graph mode (igraph) with optional post-processing
#'
#' @param layers Optional character vector of layer names to include **only in supra-graph mode**.
#' Ignored in strict passthrough mode.
#'
#' @param gamma Numeric (>0). Resolution parameter for [multinet::glouvain_ml()]
#' used when `method="glouvain"` and `supra_graph=FALSE`.
#'
#' @param omega Numeric (>=0). Inter-layer coupling for [multinet::glouvain_ml()]
#' used when `method="glouvain"` and `supra_graph=FALSE`.
#'
#' @param overlapping Logical. Passed to [multinet::infomap_ml()] when
#' `method="infomap"` and `supra_graph=FALSE`.
#'
#' @param directed Logical. Passed to [multinet::infomap_ml()] when
#' `method="infomap"` and `supra_graph=FALSE`.
#'
#' @param self.links Logical. Passed to [multinet::infomap_ml()] when
#' `method="infomap"` and `supra_graph=FALSE`.
#'
#' @param k Integer (>=3). Clique size threshold. Passed to
#' [multinet::clique_percolation_ml()] in strict mode; used as minimum clique
#' size for [igraph::cliques()] in supra-graph mode.
#'
#' @param m Integer (>=1). Multilayer clique parameter passed to
#' [multinet::clique_percolation_ml()] in strict mode. Ignored in supra-graph mode.
#'
#' @param edgeWeight How to aggregate edges across layers in supra-graph mode:
#' one of `"count"` or `"sum"`. Ignored when `supra_graph=FALSE`.
#'
#' @param min.actors Integer (>=1). Post-filtering threshold for supra-graph mode.
#' For `method="abacus"` in strict mode, passed through to [multinet::abacus_ml()].
#'
#' @param min.layers Integer (>=1). Post-filtering threshold for supra-graph mode.
#' For `method="abacus"` in strict mode, passed through to [multinet::abacus_ml()].
#'
#' @param relabel_by_size Logical. In supra-graph mode, optionally relabels kept
#' communities to `"C1","C2",...` in decreasing size order.
#' In strict passthrough mode, does **not** change the returned object, but controls
#' relabeling inside the **summary CSV** (if written).
#'
#' @param seed Optional integer. If provided, `set.seed(seed)` is called before
#' running methods that may be stochastic.
#'
#' @param results_dir Output directory for optional files. Created if missing.
#'
#' @param save_to_rds Logical. If `TRUE`, saves the returned object (`supra_graph=FALSE`)
#' or membership table (`supra_graph=TRUE`) as an `.rds`.
#'
#' @param rds_file Optional filename for the RDS output. If `NULL`, a timestamped
#' name is generated. Relative paths are placed under `results_dir`.
#'
#' @param write_csv Logical. If `TRUE`, writes a CSV (coercing to data frame when needed).
#'
#' @param csv_prefix Character prefix used to build output filenames.
#'
#' @param write_summary_csv Logical. If `TRUE`, writes the summary CSV with columns
#' `com,size,span` (ordered by decreasing size).
#'
#' @param verbose Logical. If `TRUE`, prints progress messages and written file paths.
#'
#' @param ... Additional arguments passed to the underlying `multinet::*_ml()`
#' function when in strict mode (`supra_graph=FALSE`). Also used to accept the
#' backward-compatible aliases listed above.
#'
#' @return
#' - If `supra_graph = FALSE`: returns the **raw object** returned by `multinet::*_ml()`
#'   unchanged.
#' - If `supra_graph = TRUE`: returns a data frame with columns `actor`, `com`, `layer`, `method`
#'   and attributes `community_sizes`, `layer_span`, and `files`.
#'
#' @examples
#' # ---- Build a tiny 2-layer multiplex ----
#' g1 <- igraph::graph_from_data_frame(
#'   data.frame(from = c("a","a","b","c"),
#'              to   = c("b","c","c","d"),
#'              w_   = 1),
#'   directed = FALSE
#' )
#' g2 <- igraph::graph_from_data_frame(
#'   data.frame(from = c("a","a","b","c"),
#'              to   = c("b","c","c","e"),
#'              w_   = 1),
#'   directed = FALSE
#' )
#'
#' net_A <- multinet::ml_empty()
#' multinet::add_igraph_layer_ml(net_A, g1, name = "L1")
#' multinet::add_igraph_layer_ml(net_A, g2, name = "L2")
#'
#' # ============================================================
#' # Example 1: STRICT passthrough (glouvain, supra_graph = FALSE)
#' # ============================================================
#' set.seed(1)
#' cm_mn <- multinet::glouvain_ml(net_A, gamma = 1, omega = 1)
#'
#' cm_dc <- detectCom(
#'   net_A,
#'   method      = "glouvain",
#'   supra_graph = FALSE,
#'   gamma       = 1,
#'   omega       = 1,
#'   seed        = 1,
#'   verbose     = FALSE
#' )
#'
#' all.equal(cm_mn, cm_dc)
#'
#' # ============================================================
#' # Example 2: STRICT passthrough (clique, supra_graph = FALSE)
#' #           + write summary CSV (com,size,span)
#' # ============================================================
#' cm_mn2 <- multinet::clique_percolation_ml(net_A, k = 3, m = 2)
#'
#' tmpdir <- file.path(tempdir(), "detectCom_example")
#' dir.create(tmpdir, showWarnings = FALSE, recursive = TRUE)
#'
#' cm_dc2 <- detectCom(
#'   net_A,
#'   method            = "clique",
#'   supra_graph       = FALSE,
#'   k                 = 3,
#'   m                 = 2,
#'   results_dir       = tmpdir,
#'   write_summary_csv = TRUE,
#'   relabel_by_size   = TRUE,  # summary file uses C1,C2,...; return object unchanged
#'   verbose           = FALSE
#' )
#'
#' all.equal(cm_mn2, cm_dc2)
#' list.files(tmpdir, pattern = "_summary\\.csv$", full.names = TRUE)
#'
#' # ============================================================
#' # Example 3: STRICT passthrough (infomap, supra_graph = FALSE)
#' # ============================================================
#' set.seed(7)
#' cm_mn3 <- multinet::infomap_ml(net_A, overlapping = FALSE, directed = FALSE, self.links = TRUE)
#'
#' cm_dc3 <- detectCom(
#'   net_A,
#'   method      = "infomap",
#'   supra_graph = FALSE,
#'   overlapping = FALSE,
#'   directed    = FALSE,
#'   self.links  = TRUE,
#'   seed        = 7,
#'   verbose     = FALSE
#' )
#'
#' all.equal(cm_mn3, cm_dc3)
#'
#' @importFrom multinet ml_empty add_igraph_layer_ml edges_ml
#' @importFrom multinet glouvain_ml infomap_ml clique_percolation_ml abacus_ml
#' @importFrom igraph graph_from_data_frame simplify cliques components cluster_louvain cluster_infomap membership V E vcount ecount
#' @importFrom stats aggregate
#' @importFrom utils write.csv
#' @export
detectCom <- function(
    net,
    method = c("glouvain", "louvain", "infomap", "clique", "abacus"),
    supra_graph = FALSE,
    layers = NULL,

    # ---- multinet algorithm parameters (names match multinet) ----
    gamma       = 1,
    omega       = 1,
    overlapping = FALSE,
    directed    = FALSE,
    self.links  = TRUE,
    k           = 3,
    m           = 1,

    # ---- supra-graph options ----
    edgeWeight = c("count", "sum"),

    # ---- post-processing (ONLY used when supra_graph = TRUE) ----
    min.actors      = 15,
    min.layers      = 2,
    relabel_by_size = TRUE,

    # ---- reproducibility ----
    seed = NULL,

    # ---- optional outputs ----
    results_dir       = getOption("mlnet.results_dir", "omicsDNA_results"),
    save_to_rds       = FALSE,
    rds_file          = NULL,
    write_csv         = FALSE,
    csv_prefix        = "communities",
    write_summary_csv = FALSE,
    verbose           = TRUE,

    ...
) {
  # ---- dependencies (safe even inside package; helpful when sourced standalone) ----
  if (!requireNamespace("multinet", quietly = TRUE)) {
    stop("Package 'multinet' is required. Install it with install.packages('multinet').", call. = FALSE)
  }
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package 'igraph' is required. Install it with install.packages('igraph').", call. = FALSE)
  }

  # ---- normalize args ----
  method <- match.arg(method)
  if (identical(method, "louvain")) {
    warning("`method = \"louvain\"` is deprecated; using `method = \"glouvain\"`.", call. = FALSE)
    method <- "glouvain"
  }
  edgeWeight <- match.arg(edgeWeight)

  # ---- capture ... and map legacy aliases (if supplied) ----
  dots <- list(...)

  # legacy glouvain names
  if (!is.null(dots$glouvain_gamma)) {
    if (missing(gamma)) gamma <- dots$glouvain_gamma
    dots$glouvain_gamma <- NULL
  }
  if (!is.null(dots$glouvain_omega)) {
    if (missing(omega)) omega <- dots$glouvain_omega
    dots$glouvain_omega <- NULL
  }

  # legacy clique names
  if (!is.null(dots$clique.k)) {
    if (missing(k)) k <- dots$clique.k
    dots$clique.k <- NULL
  }
  if (!is.null(dots$clique.m)) {
    if (missing(m)) m <- dots$clique.m
    dots$clique.m <- NULL
  }

  # legacy infomap names
  if (!is.null(dots$infomap_overlapping)) {
    if (missing(overlapping)) overlapping <- dots$infomap_overlapping
    dots$infomap_overlapping <- NULL
  }
  if (!is.null(dots$infomap_directed)) {
    if (missing(directed)) directed <- dots$infomap_directed
    dots$infomap_directed <- NULL
  }
  if (!is.null(dots$infomap_self_links)) {
    if (missing(self.links)) self.links <- dots$infomap_self_links
    dots$infomap_self_links <- NULL
  }

  # ---- validate ----
  stopifnot(is.logical(supra_graph), length(supra_graph) == 1L)
  stopifnot(is.logical(verbose), length(verbose) == 1L)

  stopifnot(is.numeric(gamma), length(gamma) == 1L, gamma > 0)
  stopifnot(is.numeric(omega), length(omega) == 1L, omega >= 0)

  stopifnot(is.logical(overlapping), length(overlapping) == 1L)
  stopifnot(is.logical(directed), length(directed) == 1L)
  stopifnot(is.logical(self.links), length(self.links) == 1L)

  stopifnot(is.numeric(k), length(k) == 1L, k >= 3)
  stopifnot(is.numeric(m), length(m) == 1L, m >= 1)

  stopifnot(is.numeric(min.actors), length(min.actors) == 1L, min.actors >= 1)
  stopifnot(is.numeric(min.layers), length(min.layers) == 1L, min.layers >= 1)
  stopifnot(is.logical(relabel_by_size), length(relabel_by_size) == 1L)

  # ---- helpers ----
  .ensure_dir <- function(d) {
    if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  }
  .is_abs <- function(p) {
    is.character(p) && length(p) == 1L && grepl("^(/|[A-Za-z]:[\\/])", p)
  }
  .pick <- function(cands, nms) {
    z <- cands[cands %in% nms]
    if (length(z)) z[1] else NA_character_
  }
  .stamp <- function() format(Sys.time(), "%Y-%m-%d_%H%M%S")

  .write_summary_like_attachment <- function(df, base) {
    nm <- names(df)
    actor_col <- .pick(c("actor","Actor","node","Node","vertex","Vertex"), nm)
    layer_col <- .pick(c("layer","Layer"), nm)
    com_col   <- .pick(c("com","community","cluster","cid","community_id","Community"), nm)

    if (any(is.na(c(actor_col, layer_col, com_col)))) {
      if (isTRUE(verbose)) message("Skipping summary CSV: could not infer actor/layer/community columns.")
      return(invisible(NULL))
    }

    # dedupe rows to avoid double counting
    tmp <- unique(data.frame(
      actor = as.character(df[[actor_col]]),
      layer = as.character(df[[layer_col]]),
      com   = as.character(df[[com_col]]),
      stringsAsFactors = FALSE
    ))

    size_by <- tapply(tmp$actor, tmp$com, function(v) length(unique(v)))
    span_by <- tapply(tmp$layer, tmp$com, function(v) length(unique(v)))

    size_by <- as.integer(size_by)
    span_by <- as.integer(span_by[names(size_by)])

    ord <- order(size_by, decreasing = TRUE)
    ids <- names(size_by)[ord]

    if (isTRUE(relabel_by_size)) {
      com_out <- paste0("C", seq_along(ids))
    } else {
      com_out <- as.character(ids)
    }

    sm <- data.frame(
      com  = com_out,
      size = as.integer(size_by[ids]),
      span = as.integer(span_by[ids]),
      stringsAsFactors = FALSE
    )

    sum_file <- file.path(results_dir, paste0(base, "_summary.csv"))
    utils::write.csv(sm, sum_file, row.names = FALSE)
    if (isTRUE(verbose)) message("Saved summary CSV: ", normalizePath(sum_file, winslash = "/", mustWork = FALSE))

    invisible(sum_file)
  }

  # Save outputs WITHOUT modifying the returned object (strict mode)
  .save_raw_outputs <- function(obj, method_used, mode_tag, weight_tag = NULL) {
    if (!isTRUE(save_to_rds) && !isTRUE(write_csv) && !isTRUE(write_summary_csv)) return(invisible(NULL))

    .ensure_dir(results_dir)
    st <- .stamp()

    base <- if (is.null(weight_tag)) {
      sprintf("%s_%s_%s_%s", csv_prefix, tolower(method_used), mode_tag, st)
    } else {
      sprintf("%s_%s_%s_%s_%s", csv_prefix, tolower(method_used), mode_tag, tolower(weight_tag), st)
    }

    # RDS
    if (isTRUE(save_to_rds)) {
      rds_use <- if (is.null(rds_file)) paste0(base, ".rds") else rds_file
      if (!.is_abs(rds_use)) rds_use <- file.path(results_dir, rds_use)
      saveRDS(obj, rds_use)
      if (isTRUE(verbose)) message("Saved RDS: ", normalizePath(rds_use, winslash = "/", mustWork = FALSE))
    }

    # CSV (coerce copy only)
    if (isTRUE(write_csv)) {
      csv_file <- file.path(results_dir, paste0(base, ".csv"))
      obj_df <- if (is.data.frame(obj)) obj else try(as.data.frame(obj, stringsAsFactors = FALSE), silent = TRUE)
      if (inherits(obj_df, "try-error") || is.null(obj_df)) {
        warning("write_csv=TRUE but raw output could not be coerced to data.frame; skipping CSV.", call. = FALSE)
      } else {
        utils::write.csv(obj_df, csv_file, row.names = FALSE)
        if (isTRUE(verbose)) message("Saved CSV: ", normalizePath(csv_file, winslash = "/", mustWork = FALSE))
      }
    }

    # Summary CSV (coerce copy only)
    if (isTRUE(write_summary_csv)) {
      obj_df <- if (is.data.frame(obj)) obj else try(as.data.frame(obj, stringsAsFactors = FALSE), silent = TRUE)
      if (inherits(obj_df, "try-error") || is.null(obj_df) || !nrow(obj_df)) {
        if (isTRUE(verbose)) message("Skipping summary CSV: raw output not coercible or empty.")
      } else {
        .write_summary_like_attachment(obj_df, base)
      }
    }

    invisible(NULL)
  }

  # Finalize supra-mode output (attach attrs + save)
  .finalize_supra <- function(res_df, method_used, weight_tag = NULL) {
    .ensure_dir(results_dir)
    st <- .stamp()

    base <- if (is.null(weight_tag)) {
      sprintf("%s_%s_supra_%s", csv_prefix, tolower(method_used), st)
    } else {
      sprintf("%s_%s_supra_%s_%s", csv_prefix, tolower(method_used), tolower(weight_tag), st)
    }

    files <- list()

    if (isTRUE(save_to_rds)) {
      rds_use <- if (is.null(rds_file)) paste0(base, ".rds") else rds_file
      if (!.is_abs(rds_use)) rds_use <- file.path(results_dir, rds_use)
      saveRDS(res_df, rds_use)
      files$rds <- rds_use
      if (isTRUE(verbose)) message("Saved RDS: ", normalizePath(rds_use, winslash = "/", mustWork = FALSE))
    }

    if (isTRUE(write_csv)) {
      csv_file <- file.path(results_dir, paste0(base, ".csv"))
      utils::write.csv(res_df, csv_file, row.names = FALSE)
      files$csv <- csv_file
      if (isTRUE(verbose)) message("Saved CSV: ", normalizePath(csv_file, winslash = "/", mustWork = FALSE))
    }

    if (isTRUE(write_summary_csv) && nrow(res_df)) {
      sum_file <- .write_summary_like_attachment(res_df, base)
      if (!is.null(sum_file)) files$summary_csv <- sum_file
    }

    # Always attach sizes/spans for supra output
    if (nrow(res_df)) {
      size_by <- tapply(res_df$actor, res_df$com, function(v) length(unique(v)))
      span_by <- tapply(res_df$layer, res_df$com, function(v) length(unique(v)))
      ord <- order(as.integer(size_by), decreasing = TRUE)
      size_by <- as.integer(size_by[ord])
      span_by <- as.integer(span_by[names(size_by)])
    } else {
      size_by <- integer(0)
      span_by <- integer(0)
    }

    attr(res_df, "community_sizes") <- size_by
    attr(res_df, "layer_span")      <- span_by
    attr(res_df, "files")           <- files

    res_df
  }

  # abacus has no supra-graph variant here
  use_ml_native <- (!isTRUE(supra_graph)) || method == "abacus"

  # =======================================================================
  # 1) STRICT multilayer-native passthrough (supra_graph = FALSE)
  # =======================================================================
  if (use_ml_native) {
    if (!is.null(seed)) set.seed(as.integer(seed))

    if (isTRUE(verbose)) {
      msg <- paste0("Running STRICT passthrough multinet::*_ml(), method = \"", method, "\".")
      if (isTRUE(supra_graph) && method == "abacus") {
        msg <- paste0(msg, " Note: supra_graph=TRUE is ignored for method=\"abacus\".")
      }
      if (!is.null(layers) && length(layers)) {
        msg <- paste0(msg, " Note: `layers` is ignored in strict passthrough mode (subset upstream if needed).")
      }
      msg <- paste0(msg, " Note: no wrapper post-filtering/relabeling is applied when supra_graph=FALSE.")
      message(msg)
    }

    cm_raw <- switch(
      method,
      glouvain = do.call(multinet::glouvain_ml,
                         c(list(n = net, gamma = gamma, omega = omega), dots)),
      infomap  = do.call(multinet::infomap_ml,
                         c(list(n = net,
                                overlapping = overlapping,
                                directed    = directed,
                                self.links  = self.links), dots)),
      clique   = do.call(multinet::clique_percolation_ml,
                         c(list(n = net, k = as.integer(k), m = as.integer(m)), dots)),
      abacus   = do.call(multinet::abacus_ml,
                         c(list(n = net, min.actors = min.actors, min.layers = min.layers), dots)),
      stop("Unsupported method: ", method, call. = FALSE)
    )

    # optional saving (does not modify cm_raw)
    .save_raw_outputs(cm_raw, method_used = method, mode_tag = "ml")

    # strict guarantee: return raw object unchanged
    return(cm_raw)
  }

  # =======================================================================
  # 2) Supra-graph mode (supra_graph = TRUE; method in glouvain/infomap/clique)
  # =======================================================================
  if (!method %in% c("glouvain", "infomap", "clique")) {
    stop("supra_graph=TRUE is only supported for method = 'glouvain', 'infomap', or 'clique'.", call. = FALSE)
  }
  if (!is.null(seed)) set.seed(as.integer(seed))

  if (isTRUE(verbose)) {
    message("Running supra-graph mode (layer-collapsed), method = \"", method, "\".")
    message("Post-filtering/relabeling may be applied here (min.actors/min.layers/relabel_by_size).")
  }

  # ---- pull edge table from multinet ----
  Eraw <- try(multinet::edges_ml(net), silent = TRUE)
  if (inherits(Eraw, "try-error") || is.null(Eraw)) {
    stop("Could not retrieve edges from `net` via multinet::edges_ml().", call. = FALSE)
  }
  Eall <- if (is.data.frame(Eraw)) Eraw else {
    tmp <- try(as.data.frame(Eraw, stringsAsFactors = FALSE), silent = TRUE)
    if (inherits(tmp, "try-error") || is.null(tmp)) {
      stop("edges_ml(net) did not return a coercible table of edges.", call. = FALSE)
    }
    tmp
  }
  if (!nrow(Eall)) stop("No edges found in `net`.", call. = FALSE)

  nm <- names(Eall)

  # ---- identify endpoint columns ----
  pairs <- list(
    c("from_actor","to_actor"),
    c("from","to"),
    c("source","target"),
    c("actor1","actor2"),
    c("i","j"),
    c("v1","v2")
  )
  a_col <- b_col <- NA_character_
  for (p in pairs) {
    if (all(p %in% nm)) { a_col <- p[1]; b_col <- p[2]; break }
  }
  if (is.na(a_col)) {
    layer_like <- c("layer","Layer","from_layer","to_layer","l1","l2")
    char_cols <- which(vapply(Eall, function(x) is.character(x) || is.factor(x), logical(1)))
    char_cols <- setdiff(char_cols, match(layer_like, nm, nomatch = 0))
    if (length(char_cols) < 2) stop("Could not identify two endpoint columns in edges table.", call. = FALSE)
    a_col <- nm[char_cols[1]]
    b_col <- nm[char_cols[2]]
  }

  # ---- determine layer and keep intra-layer edges ----
  if ("from_layer" %in% nm && "to_layer" %in% nm) {
    Eall <- Eall[Eall$from_layer == Eall$to_layer, , drop = FALSE]
    Eall$layer <- as.character(Eall$from_layer)
  } else if ("layer" %in% nm) {
    Eall$layer <- as.character(Eall$layer)
  } else if ("Layer" %in% nm) {
    Eall$layer <- as.character(Eall$Layer)
  } else {
    Eall$layer <- "L1"
  }
  if (!nrow(Eall)) stop("No intra-layer edges found after filtering.", call. = FALSE)

  # ---- choose layers ----
  if (is.null(layers)) {
    layers_use <- sort(unique(Eall$layer))
  } else {
    layers_use <- intersect(as.character(layers), unique(Eall$layer))
  }
  if (!length(layers_use)) stop("No layers available after intersecting with `layers`.", call. = FALSE)

  Eall <- Eall[Eall$layer %in% layers_use, , drop = FALSE]
  if (!nrow(Eall)) stop("No edges left after applying `layers` selection.", call. = FALSE)

  # ---- choose weight column if present ----
  w_col <- .pick(c("weight","Weight","w","w_","value","score","correlation"), nm)

  edges_all <- data.frame(
    layer  = as.character(Eall$layer),
    from   = as.character(Eall[[a_col]]),
    to     = as.character(Eall[[b_col]]),
    weight = if (!is.na(w_col) && w_col %in% names(Eall)) as.numeric(Eall[[w_col]]) else 1,
    stringsAsFactors = FALSE
  )

  # canonicalize for undirected supra-graph
  aa <- pmin(edges_all$from, edges_all$to)
  bb <- pmax(edges_all$from, edges_all$to)
  edges_all$from <- aa
  edges_all$to   <- bb

  # ---- aggregate edge weights across layers ----
  if (edgeWeight == "count") {
    per_layer_unique <- unique(edges_all[, c("layer","from","to")])
    agg <- stats::aggregate(rep(1L, nrow(per_layer_unique)),
                            by = per_layer_unique[c("from","to")], FUN = sum)
    names(agg) <- c("from","to","weight")
  } else {
    agg <- stats::aggregate(weight ~ from + to, data = edges_all, FUN = sum, na.rm = TRUE)
  }

  g <- igraph::graph_from_data_frame(agg, directed = FALSE)
  if (igraph::ecount(g) == 0) stop("Supra-graph has no edges after aggregation.", call. = FALSE)
  g <- igraph::simplify(g, edge.attr.comb = list(weight = "sum", "ignore"))

  if (isTRUE(verbose)) {
    message("Built supra-graph: |V|=", igraph::vcount(g), " |E|=", igraph::ecount(g),
            " (edgeWeight=", edgeWeight, ")")
  }

  # ---- run supra-graph community detection ----
  memb <- NULL
  if (method == "glouvain") {
    co <- igraph::cluster_louvain(g, weights = igraph::E(g)$weight)
    memb <- igraph::membership(co)
  } else if (method == "infomap") {
    co <- igraph::cluster_infomap(g, e.weights = igraph::E(g)$weight)
    memb <- igraph::membership(co)
  } else if (method == "clique") {
    cls <- igraph::cliques(g, min = as.integer(k))
    if (!length(cls)) stop("No cliques of size >= ", as.integer(k), " found in supra-graph.", call. = FALSE)

    pairs_df <- do.call(rbind, lapply(cls, function(cl) {
      vs <- igraph::V(g)$name[cl]
      if (length(vs) < 2) return(NULL)
      mm <- utils::combn(vs, 2)
      data.frame(from = mm[1, ], to = mm[2, ], stringsAsFactors = FALSE)
    }))
    pairs_df <- unique(pairs_df)

    h <- igraph::graph_from_data_frame(pairs_df, directed = FALSE, vertices = igraph::V(g)$name)
    comp <- igraph::components(h)
    memb <- comp$membership
    names(memb) <- igraph::V(h)$name
  }

  if (is.null(memb) || !length(memb)) stop("Supra-graph community detection failed.", call. = FALSE)

  memb_vec <- as.integer(memb)
  names(memb_vec) <- names(memb)

  memb_df <- data.frame(
    actor = names(memb_vec),
    com   = as.character(memb_vec),
    stringsAsFactors = FALSE
  )

  # ---- map communities back to each layer ----
  actors_by_layer <- lapply(layers_use, function(ly) {
    ed <- edges_all[edges_all$layer == ly, , drop = FALSE]
    unique(c(ed$from, ed$to))
  })
  names(actors_by_layer) <- layers_use

  res <- do.call(rbind, lapply(layers_use, function(ly) {
    a <- actors_by_layer[[ly]]
    df <- memb_df[memb_df$actor %in% a, , drop = FALSE]
    if (!nrow(df)) return(NULL)
    df$layer  <- ly
    df$method <- method
    df
  }))
  if (is.null(res) || !nrow(res)) stop("No supra-graph community assignments after mapping back to layers.", call. = FALSE)

  # ---- post-filtering (supra_graph=TRUE only) ----
  size_by <- tapply(res$actor, res$com, function(v) length(unique(v)))
  span_by <- tapply(res$layer, res$com, function(v) length(unique(v)))
  keep <- names(size_by)[as.integer(size_by) >= min.actors & as.integer(span_by) >= min.layers]
  res <- res[res$com %in% keep, , drop = FALSE]

  if (!length(keep) || !nrow(res)) {
    if (isTRUE(verbose)) message("No supra-graph communities satisfy min.actors/min.layers after filtering.")
    empty <- data.frame(actor = character(), com = character(), layer = character(), method = character(),
                        stringsAsFactors = FALSE)
    return(.finalize_supra(empty, method_used = method, weight_tag = edgeWeight))
  }

  # ---- relabel by decreasing size (optional; supra_graph=TRUE only) ----
  if (isTRUE(relabel_by_size)) {
    size2 <- tapply(res$actor, res$com, function(v) length(unique(v)))
    ord <- order(as.integer(size2), decreasing = TRUE)
    old_ids <- names(size2)[ord]
    new_ids <- paste0("C", seq_along(old_ids))
    map <- stats::setNames(new_ids, old_ids)
    res$com <- unname(map[res$com])
  }

  .finalize_supra(res, method_used = method, weight_tag = edgeWeight)
}

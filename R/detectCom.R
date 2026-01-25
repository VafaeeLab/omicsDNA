# ---------------------------------------------------------------------------
# detectCom.R
# Community detection wrapper for multilayer networks (multinet)
# ---------------------------------------------------------------------------

#' Community detection wrapper for multilayer networks (multinet)
#'
#' @md
#' @description
#' `detectCom()` supports **two analysis modes** for a multilayer network of class
#' [multinet::ml.network].
#'
#' ## 1) Multilayer-native passthrough (`supra_graph = FALSE`, default)
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

    # --- legacy / package-style names (kept for compatibility) ---
    glouvain_gamma = 1,
    glouvain_omega = 1,
    clique.k = 3,
    clique.m = 1,
    infomap_overlapping = FALSE,
    infomap_directed    = FALSE,
    infomap_self_links  = TRUE,

    # --- strict-mode friendly aliases (your examples use these) ---
    gamma       = NULL,
    omega       = NULL,
    k           = NULL,
    m           = NULL,
    overlapping = NULL,
    directed    = NULL,
    self.links  = NULL,

    # --- supra-graph options ---
    edgeWeight = c("count", "sum"),

    # --- post-processing (ONLY applied when supra_graph = TRUE) ---
    min.actors      = 15,
    min.layers      = 2,
    relabel_by_size = TRUE,

    # --- reproducibility ---
    seed = NULL,

    # --- outputs ---
    results_dir       = getOption("mlnet.results_dir", "omicsDNA_results"),
    save_to_rds       = FALSE,
    rds_file          = NULL,
    write_csv         = FALSE,
    csv_prefix        = "communities",
    write_summary_csv = TRUE,   # <- default TRUE so omicsDNA_results is created + summary written
    verbose           = TRUE,

    ...
) {
  # ---------------------------
  # normalize / validate
  # ---------------------------
  method <- match.arg(method)
  if (identical(method, "louvain")) {
    warning("`method = \"louvain\"` is deprecated; using `method = \"glouvain\"`.", call. = FALSE)
    method <- "glouvain"
  }
  edgeWeight <- match.arg(edgeWeight)

  stopifnot(is.logical(supra_graph), length(supra_graph) == 1L)
  stopifnot(is.logical(verbose), length(verbose) == 1L)
  stopifnot(is.numeric(min.actors), length(min.actors) == 1L, min.actors >= 1)
  stopifnot(is.numeric(min.layers), length(min.layers) == 1L, min.layers >= 1)
  stopifnot(is.logical(relabel_by_size), length(relabel_by_size) == 1L)

  # resolve aliases (if user provided gamma/k/etc, those override legacy names)
  gamma_use <- if (!is.null(gamma)) gamma else glouvain_gamma
  omega_use <- if (!is.null(omega)) omega else glouvain_omega

  k_use <- if (!is.null(k)) k else clique.k
  m_use <- if (!is.null(m)) m else clique.m

  overlapping_use <- if (!is.null(overlapping)) overlapping else infomap_overlapping
  directed_use    <- if (!is.null(directed))    directed    else infomap_directed
  self_links_use  <- if (!is.null(self.links))  self.links  else infomap_self_links

  stopifnot(is.numeric(gamma_use), length(gamma_use) == 1L, gamma_use > 0)
  stopifnot(is.numeric(omega_use), length(omega_use) == 1L, omega_use >= 0)
  stopifnot(is.numeric(k_use), length(k_use) == 1L, k_use >= 3)
  stopifnot(is.numeric(m_use), length(m_use) == 1L, m_use >= 1)

  stopifnot(is.logical(overlapping_use), length(overlapping_use) == 1L)
  stopifnot(is.logical(directed_use), length(directed_use) == 1L)
  stopifnot(is.logical(self_links_use), length(self_links_use) == 1L)

  # abacus has no supra-graph implementation here
  use_ml_native <- (!isTRUE(supra_graph)) || method == "abacus"

  # ---------------------------
  # helpers
  # ---------------------------
  .ensure_dir <- function(d) {
    if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  }
  .is_abs <- function(p) {
    is.character(p) && length(p) == 1L && grepl("^(/|[A-Za-z]:[\\/])", p)
  }
  .stamp <- function() format(Sys.time(), "%Y-%m-%d_%H%M%S")
  .pick <- function(cands, nms) {
    z <- cands[cands %in% nms]
    if (length(z)) z[1] else NA_character_
  }

  .write_empty_summary <- function(base) {
    .ensure_dir(results_dir)
    sum_file <- file.path(results_dir, paste0(base, "_summary.csv"))
    sm <- data.frame(com = character(), size = integer(), span = integer(), stringsAsFactors = FALSE)
    utils::write.csv(sm, sum_file, row.names = FALSE)
    if (isTRUE(verbose)) message("Saved summary CSV (empty): ", normalizePath(sum_file, winslash = "/", mustWork = FALSE))
    sum_file
  }

  .write_summary_like_attachment <- function(df, base) {
    # Always create dir and summary file path
    .ensure_dir(results_dir)
    sum_file <- file.path(results_dir, paste0(base, "_summary.csv"))

    # Empty input -> write empty summary
    if (is.null(df) || !nrow(df)) {
      sm <- data.frame(com = character(), size = integer(), span = integer(), stringsAsFactors = FALSE)
      utils::write.csv(sm, sum_file, row.names = FALSE)
      if (isTRUE(verbose)) message("Saved summary CSV (empty): ", normalizePath(sum_file, winslash = "/", mustWork = FALSE))
      return(sum_file)
    }

    nm <- names(df)
    actor_col <- .pick(c("actor","Actor","node","Node","vertex","Vertex"), nm)
    layer_col <- .pick(c("layer","Layer"), nm)
    com_col   <- .pick(c("com","community","cluster","cid","community_id","Community"), nm)

    # If we cannot infer columns, still write empty summary (do NOT error)
    if (any(is.na(c(actor_col, layer_col, com_col)))) {
      if (isTRUE(verbose)) message("Summary CSV: could not infer actor/layer/com columns; writing empty summary.")
      sm <- data.frame(com = character(), size = integer(), span = integer(), stringsAsFactors = FALSE)
      utils::write.csv(sm, sum_file, row.names = FALSE)
      return(sum_file)
    }

    tmp <- unique(data.frame(
      actor = as.character(df[[actor_col]]),
      layer = as.character(df[[layer_col]]),
      com   = as.character(df[[com_col]]),
      stringsAsFactors = FALSE
    ))

    if (!nrow(tmp)) {
      sm <- data.frame(com = character(), size = integer(), span = integer(), stringsAsFactors = FALSE)
      utils::write.csv(sm, sum_file, row.names = FALSE)
      if (isTRUE(verbose)) message("Saved summary CSV (empty): ", normalizePath(sum_file, winslash = "/", mustWork = FALSE))
      return(sum_file)
    }

    size_by <- tapply(tmp$actor, tmp$com, function(v) length(unique(v)))
    span_by <- tapply(tmp$layer, tmp$com, function(v) length(unique(v)))

    # No communities -> empty summary
    if (length(size_by) == 0L) {
      sm <- data.frame(com = character(), size = integer(), span = integer(), stringsAsFactors = FALSE)
      utils::write.csv(sm, sum_file, row.names = FALSE)
      if (isTRUE(verbose)) message("Saved summary CSV (empty): ", normalizePath(sum_file, winslash = "/", mustWork = FALSE))
      return(sum_file)
    }

    size_by <- as.integer(size_by)
    span_by <- as.integer(span_by[names(size_by)])

    ord <- order(size_by, decreasing = TRUE, na.last = TRUE)
    ids <- names(size_by)[ord]
    ids <- ids[!is.na(ids)]

    # if after removing NA we have nothing, still empty summary
    if (!length(ids)) {
      sm <- data.frame(com = character(), size = integer(), span = integer(), stringsAsFactors = FALSE)
      utils::write.csv(sm, sum_file, row.names = FALSE)
      if (isTRUE(verbose)) message("Saved summary CSV (empty): ", normalizePath(sum_file, winslash = "/", mustWork = FALSE))
      return(sum_file)
    }

    size_vec <- as.integer(size_by[ids])
    span_vec <- as.integer(span_by[ids])

    com_out <- if (isTRUE(relabel_by_size)) paste0("C", seq_along(ids)) else as.character(ids)

    sm <- data.frame(
      com  = com_out,
      size = as.integer(size_vec),
      span = as.integer(span_vec),
      stringsAsFactors = FALSE
    )

    utils::write.csv(sm, sum_file, row.names = FALSE)
    if (isTRUE(verbose)) message("Saved summary CSV: ", normalizePath(sum_file, winslash = "/", mustWork = FALSE))
    sum_file
  }

  .save_outputs_raw <- function(obj, method_used, mode_tag, weight_tag = NULL) {
    if (!isTRUE(save_to_rds) && !isTRUE(write_csv) && !isTRUE(write_summary_csv)) return(invisible(NULL))

    .ensure_dir(results_dir)
    base <- if (is.null(weight_tag)) {
      sprintf("%s_%s_%s_%s", csv_prefix, tolower(method_used), mode_tag, .stamp())
    } else {
      sprintf("%s_%s_%s_%s_%s", csv_prefix, tolower(method_used), mode_tag, tolower(weight_tag), .stamp())
    }

    # RDS (save raw object unchanged)
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

    # Summary CSV (coerce copy only; robust to empties)
    if (isTRUE(write_summary_csv)) {
      obj_df <- if (is.data.frame(obj)) obj else try(as.data.frame(obj, stringsAsFactors = FALSE), silent = TRUE)
      if (inherits(obj_df, "try-error") || is.null(obj_df)) {
        if (isTRUE(verbose)) message("Summary CSV: raw output not coercible; writing empty summary.")
        .write_empty_summary(base)
      } else {
        .write_summary_like_attachment(obj_df, base)
      }
    }

    invisible(NULL)
  }

  .finalize_supra <- function(res_df, method_used, weight_tag = NULL) {
    if (!isTRUE(save_to_rds) && !isTRUE(write_csv) && !isTRUE(write_summary_csv)) {
      # still attach attrs even if no writing
      files <- list()
    } else {
      .ensure_dir(results_dir)
      files <- list()
    }

    base <- if (is.null(weight_tag)) {
      sprintf("%s_%s_supra_%s", csv_prefix, tolower(method_used), .stamp())
    } else {
      sprintf("%s_%s_supra_%s_%s", csv_prefix, tolower(method_used), tolower(weight_tag), .stamp())
    }

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

    if (isTRUE(write_summary_csv)) {
      sum_file <- .write_summary_like_attachment(res_df, base)
      files$summary_csv <- sum_file
    }

    # attributes (always)
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

  # ---------------------------
  # RNG seed
  # ---------------------------
  if (!is.null(seed)) set.seed(as.integer(seed))

  # =======================================================================
  # 1) STRICT multilayer-native passthrough (supra_graph = FALSE OR abacus)
  # =======================================================================
  if (use_ml_native) {
    if (isTRUE(verbose)) {
      msg <- paste0("Running STRICT passthrough multinet::*_ml(), method = \"", method, "\".")
      if (!is.null(layers) && length(layers)) {
        msg <- paste0(msg, " Note: `layers` is ignored in strict mode (subset upstream if needed).")
      }
      message(msg)
    }

    cm_raw <- switch(
      method,
      glouvain = multinet::glouvain_ml(net, gamma = gamma_use, omega = omega_use, ...),
      infomap  = multinet::infomap_ml(net,
                                      overlapping = overlapping_use,
                                      directed    = directed_use,
                                      self.links  = self_links_use, ...),
      clique   = multinet::clique_percolation_ml(net, k = as.integer(k_use), m = as.integer(m_use), ...),
      abacus   = multinet::abacus_ml(net, min.actors = min.actors, min.layers = min.layers, ...),
      stop("Unsupported method: ", method, call. = FALSE)
    )

    # side-effect outputs only; return object unchanged
    .save_outputs_raw(cm_raw, method_used = method, mode_tag = "ml", weight_tag = NULL)

    return(cm_raw)
  }

  # =======================================================================
  # 2) Supra-graph mode (supra_graph = TRUE; glouvain/infomap/clique)
  # =======================================================================
  if (!method %in% c("glouvain", "infomap", "clique")) {
    stop("supra_graph=TRUE is only supported for method = 'glouvain', 'infomap', or 'clique'.", call. = FALSE)
  }

  if (isTRUE(verbose)) {
    message("Running supra-graph mode (layer-collapsed), method = \"", method, "\".")
  }

  # ---- pull multilayer edges ----
  Eraw <- try(multinet::edges_ml(net), silent = TRUE)
  if (inherits(Eraw, "try-error") || is.null(Eraw)) {
    stop("Could not retrieve edges from `net` via multinet::edges_ml().", call. = FALSE)
  }
  Eall <- if (is.data.frame(Eraw)) Eraw else {
    tmp <- try(as.data.frame(Eraw, stringsAsFactors = FALSE), silent = TRUE)
    if (inherits(tmp, "try-error") || is.null(tmp)) stop("edges_ml(net) did not return a coercible table.", call. = FALSE)
    tmp
  }
  if (!nrow(Eall)) stop("No edges found in `net`.", call. = FALSE)

  nm <- names(Eall)

  # endpoints
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
    char_cols <- which(vapply(Eall, function(x) is.character(x) || is.factor(x), logical(1)))
    char_cols <- setdiff(char_cols, match(layer_like, nm, nomatch = 0))
    if (length(char_cols) < 2) stop("Could not identify two endpoint columns in edges table.", call. = FALSE)
    a_col <- nm[char_cols[1]]
    b_col <- nm[char_cols[2]]
  }

  # layer + intra-layer filtering
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
  if (!nrow(Eall)) stop("No intra-layer edges found.", call. = FALSE)

  # choose layers
  if (is.null(layers)) {
    layers_use <- sort(unique(Eall$layer))
  } else {
    layers_use <- intersect(as.character(layers), unique(Eall$layer))
  }
  if (!length(layers_use)) stop("No layers available after intersecting with `layers`.", call. = FALSE)
  Eall <- Eall[Eall$layer %in% layers_use, , drop = FALSE]
  if (!nrow(Eall)) stop("No edges left after applying `layers` selection.", call. = FALSE)

  # weight col (optional)
  w_col <- .pick(c("weight","Weight","w","w_","value","score","correlation"), nm)

  edges_all <- data.frame(
    layer  = as.character(Eall$layer),
    from   = as.character(Eall[[a_col]]),
    to     = as.character(Eall[[b_col]]),
    weight = if (!is.na(w_col) && w_col %in% names(Eall)) as.numeric(Eall[[w_col]]) else 1,
    stringsAsFactors = FALSE
  )

  # canonicalize undirected
  aa <- pmin(edges_all$from, edges_all$to)
  bb <- pmax(edges_all$from, edges_all$to)
  edges_all$from <- aa
  edges_all$to   <- bb

  # aggregate
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

  # detect communities on supra graph
  memb <- NULL
  if (method == "glouvain") {
    co <- igraph::cluster_louvain(g, weights = igraph::E(g)$weight)
    memb <- igraph::membership(co)
  } else if (method == "infomap") {
    co <- igraph::cluster_infomap(g, e.weights = igraph::E(g)$weight)
    memb <- igraph::membership(co)
  } else if (method == "clique") {
    cls <- igraph::cliques(g, min = as.integer(k_use))
    if (!length(cls)) stop("No cliques of size >= ", as.integer(k_use), " found in supra-graph.", call. = FALSE)

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
  memb_df <- data.frame(actor = names(memb_vec), com = as.character(memb_vec), stringsAsFactors = FALSE)

  # map back to layers
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

  if (is.null(res) || !nrow(res)) {
    empty <- data.frame(actor = character(), com = character(), layer = character(), method = character(),
                        stringsAsFactors = FALSE)
    return(.finalize_supra(empty, method_used = method, weight_tag = edgeWeight))
  }

  # post-filter (supra only)
  size_by <- tapply(res$actor, res$com, function(v) length(unique(v)))
  span_by <- tapply(res$layer, res$com, function(v) length(unique(v)))
  keep <- names(size_by)[as.integer(size_by) >= min.actors & as.integer(span_by) >= min.layers]
  res <- res[res$com %in% keep, , drop = FALSE]

  if (!length(keep) || !nrow(res)) {
    empty <- data.frame(actor = character(), com = character(), layer = character(), method = character(),
                        stringsAsFactors = FALSE)
    return(.finalize_supra(empty, method_used = method, weight_tag = edgeWeight))
  }

  # relabel (supra only)
  if (isTRUE(relabel_by_size)) {
    size2 <- tapply(res$actor, res$com, function(v) length(unique(v)))
    ord <- order(as.integer(size2), decreasing = TRUE)
    old_ids <- names(size2)[ord]
    new_ids <- paste0("C", seq_along(old_ids))
    map <- stats::setNames(new_ids, old_ids)
    res$com <- unname(map[res$com])
  }

  # standardize column order/types
  res <- res[, c("actor","com","layer","method"), drop = FALSE]
  res$actor  <- as.character(res$actor)
  res$com    <- as.character(res$com)
  res$layer  <- as.character(res$layer)
  res$method <- as.character(res$method)

  .finalize_supra(res, method_used = method, weight_tag = edgeWeight)
}

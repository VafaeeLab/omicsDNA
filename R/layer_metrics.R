
#' Compute per‑layer network and node metrics, and save results to a run folder:
#'
#' @description
#' Builds an \pkg{igraph} graph for each layer of a \pkg{multinet} object,
#' computes a standard set of **node‑level** and **layer‑level** metrics, and
#' writes all artefacts into a timestamped run folder. This is designed to be a
#' drop‑in analysis step you can re‑run reproducibly.
#'
#' @details
#' **Input normalisation and scope**
#' - Edges are retrieved once via `multinet::edges_ml(net)` and coerced to a data
#'   frame. If both `from_layer` and `to_layer` are present, **only intra‑layer**
#'   edges are kept (`from_layer == to_layer`) and a single `layer` column is set
#'   to `from_layer`. Otherwise, an existing `layer`/`Layer` column is used; if
#'   neither exists, a single layer `"L1"` is assumed.
#' - Endpoint columns are identified robustly by name (tries common pairs such as
#'   `from`/`to`, `source`/`target`, `actor1`/`actor2`). If these are not found,
#'   the **first two character columns** (excluding obvious layer columns) are
#'   used as endpoints.
#'
#' **Graph construction and directedness**
#' - For each selected layer, an \pkg{igraph} object is created with vertices
#'   defined by the endpoints appearing in that layer’s edges. If `directed = TRUE`,
#'   a directed graph is built and passed to all directed‑aware metrics; otherwise
#'   an undirected graph is used.
#'
#' **Node‑level metrics (per layer)**
#' - `degree` (mode `"all"`; for directed graphs this is in + out),
#' - `betweenness` (normalised; respects `directed`),
#' - `eigenvector` centrality (scaled; respects `directed`; errors are caught and
#'   return `NA`),
#' - `closeness` (mode `"all"`; non‑finite values coerced to `0`),
#' - `pagerank` (respects `directed`; returns the PageRank vector),
#' - local `clustering` coefficient (`transitivity(type = "local", isolates = "zero")`),
#' - `coreness` (k‑core index; mode `"all"`).
#'
#' **Layer‑level summaries**
#' - `n_nodes`, `n_edges`, `density` (`edge_density(loops = FALSE)`),
#' - number of weak components (`components(mode = "weak")`),
#' - `avg_path_len` and `diameter` with `unconnected = TRUE` (guarded by `tryCatch`;
#'   may return `NA` on degenerate graphs),
#' - degree distribution summaries: `mean`, `sd`, `median`, `min`, `max`,
#'   coefficient of variation (`sd/mean`, `NA` if mean = 0), and **skewness**
#'   (third central moment over `sd^3`; requires `n ≥ 3`, otherwise `NA`).
#'
#' **Saving strategy (run‑folder style)**
#' - Uses the global results directory `getOption("mlnet.results_dir", "omicsDNA_results")`.
#' - Creates a **run‑specific subfolder** named
#'   `"<file_prefix>_<YYYY-mm-dd_HHMMSS>"` inside that directory.
#' - Writes:
#'   - `layer_metrics_summary.csv` (one row per layer);
#'   - node metrics either as a single workbook
#'     `layer_metrics_nodes.xlsx` (one sheet per layer; requires **writexl**),
#'     **or** as `nodes_per_layer_csv/` containing one CSV per layer.
#'
#' **Notes**
#' - All metrics are computed on **unweighted** graphs. If weighted metrics are
#'   needed, extend the \pkg{igraph} calls accordingly.
#' - Layers with no nodes/edges yield empty node tables and `NA` layer summaries
#'   where appropriate.
#'
#' @param net       A multilayer network compatible with
#'                  `multinet::edges_ml()` (and typically `multinet::layers_ml()`).
#' @param layers    Optional character vector of layer names to analyse. Default
#'                  `NULL` (analyse **all** layers found after intra‑layer filtering).
#' @param directed  Logical; build directed graphs (`TRUE`) or undirected graphs
#'                  (`FALSE`). Default `FALSE`.
#' @param file_prefix Basename used for the run folder and output files inside it.
#'                    Default `"layer_metrics"`.
#' @param verbose   Logical; print per‑layer progress and file paths. Default `TRUE`.
#'
#' @return Invisibly returns a list with:
#' \itemize{
#'   \item `summary` — data frame with one row per layer and the layer‑level metrics;
#'   \item `nodes`   — named list of data frames (node‑level metrics per layer);
#'   \item `paths`   — list of absolute paths: `run_dir`, `summary_csv`, and either
#'         `nodes_xlsx` (when **writexl** is available) or `per_layer_dir` (CSV folder).
#' }
#'
#' @examples
#' \dontrun{
#' # Optionally set your project results folder once:
#' # options(mlnet.results_dir = "/path/to/results")
#'
#' out <- layer_metrics(net)  # writes to omicsDNA_results/layer_metrics_<stamp>/
#' out$summary[, c("layer","n_nodes","n_edges","density")]
#' head(out$nodes[[ out$summary$layer[1] ]])
#' out$paths
#' }
#'
#' @seealso
#'   \code{\link{detectCom}} for community detection prior to metric summaries;
#'   \code{\link{build_multiNet}} for assembling multilayer networks.
#'
#' @importFrom multinet edges_ml
#' @importFrom igraph graph_from_data_frame degree betweenness eigen_centrality closeness page_rank transitivity coreness components average.path.length diameter vcount ecount V edge_density
#' @importFrom utils write.csv
#' @importFrom stats sd median
#' @export
layer_metrics <- function(net,
                          layers      = NULL,
                          directed    = FALSE,
                          file_prefix = "layer_metrics",
                          verbose     = TRUE) {
  # --------- resolve run folder under the global results dir ----------
  base_dir <- getOption("mlnet.results_dir", "omicsDNA_results")
  if (!dir.exists(base_dir)) dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)
  stamp   <- format(Sys.time(), "%Y-%m-%d_%H%M%S")
  run_dir <- file.path(base_dir, paste0(file_prefix, "_", stamp))
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)

  # --------- capability check ----------
  ok_layers <- !inherits(try(multinet::layers_ml(net), silent = TRUE), "try-error")
  ok_edges  <- !inherits(try(multinet::edges_ml(net),  silent = TRUE), "try-error")
  if (!(ok_layers || ok_edges)) {
    stop("`net` does not behave like a multinet object: both layers_ml() and edges_ml() failed.")
  }

  # --------- fetch edges once ----------
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

  # --------- identify endpoint columns (robustly) ----------
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

  # --------- normalize `layer` & keep intra-layer ----------
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
  if (!nrow(Eall)) stop("No intra-layer edges found in the network.")

  # --------- choose layers ----------
  all_layers <- sort(unique(Eall$layer))
  if (is.null(layers)) layers <- all_layers else layers <- intersect(as.character(layers), all_layers)
  if (!length(layers)) stop("No valid layers to analyze after intersection.")

  # --------- helpers ----------
  .skewness <- function(x) {
    x <- x[is.finite(x)]
    n <- length(x); if (n < 3) return(NA_real_)
    m <- mean(x); s <- stats::sd(x); if (isTRUE(all.equal(s, 0))) return(NA_real_)
    sum((x - m)^3) / ((n - 1) * s^3)
  }
  .summ <- function(x) {
    x <- x[is.finite(x)]
    c(mean = mean(x), sd = stats::sd(x), median = stats::median(x),
      min = min(x), max = max(x))
  }
  .sanitize <- function(s) gsub("[^[:alnum:]_.-]+", "_", s)

  per_layer_nodes <- list()
  summary_rows    <- list()

  # --------- per-layer graphs & metrics ----------
  for (ly in layers) {
    if (verbose) message("• Processing layer: ", ly)
    ed <- Eall[Eall$layer == ly, , drop = FALSE]
    vs <- unique(c(as.character(ed[[a_col]]), as.character(ed[[b_col]])))

    g <- igraph::graph_from_data_frame(
      d = data.frame(from = ed[[a_col]], to = ed[[b_col]], stringsAsFactors = FALSE),
      directed = directed,
      vertices = if (length(vs)) data.frame(name = vs, stringsAsFactors = FALSE) else NULL
    )

    # -- node metrics
    if (igraph::vcount(g) > 0) {
      deg  <- igraph::degree(g, mode = "all")
      bet  <- suppressWarnings(igraph::betweenness(g, directed = directed, normalized = TRUE))
      eigv <- tryCatch(igraph::eigen_centrality(g, directed = directed, scale = TRUE)$vector,
                       error = function(e) rep(NA_real_, igraph::vcount(g)))
      clo  <- suppressWarnings(igraph::closeness(g, mode = "all"))
      clo[!is.finite(clo)] <- 0
      pr   <- igraph::page_rank(g, directed = directed)$vector
      clus <- igraph::transitivity(g, type = "local", isolates = "zero")
      core <- igraph::coreness(g, mode = "all")

      nodes_tbl <- data.frame(
        actor       = igraph::V(g)$name,
        degree      = as.numeric(deg),
        betweenness = as.numeric(bet),
        eigenvector = as.numeric(eigv),
        closeness   = as.numeric(clo),
        pagerank    = as.numeric(pr),
        clustering  = as.numeric(clus),
        coreness    = as.numeric(core),
        stringsAsFactors = FALSE
      )
    } else {
      nodes_tbl <- data.frame(
        actor       = character(),
        degree      = numeric(),
        betweenness = numeric(),
        eigenvector = numeric(),
        closeness   = numeric(),
        pagerank    = numeric(),
        clustering  = numeric(),
        coreness    = numeric(),
        stringsAsFactors = FALSE
      )
    }
    per_layer_nodes[[ly]] <- nodes_tbl

    # -- layer metrics
    n_nodes <- igraph::vcount(g); n_edges <- igraph::ecount(g)
    dens    <- if (n_nodes > 1) igraph::edge_density(g, loops = FALSE) else NA_real_
    comps   <- if (n_nodes > 0) igraph::components(g, mode = "weak")$no else 0
    apl     <- tryCatch(igraph::average.path.length(g, directed = directed, unconnected = TRUE),
                        error = function(e) NA_real_)
    diam    <- tryCatch(igraph::diameter(g, directed = directed, unconnected = TRUE),
                        error = function(e) NA_real_)

    if (nrow(nodes_tbl)) {
      deg_s  <- .summ(nodes_tbl$degree)
      deg_cv <- if (isTRUE(all.equal(deg_s["mean"], 0))) NA_real_ else deg_s["sd"] / deg_s["mean"]
      deg_sk <- .skewness(nodes_tbl$degree)
      dmin <- deg_s["min"]; dmax <- deg_s["max"]
    } else {
      deg_s  <- c(mean = NA_real_, sd = NA_real_, median = NA_real_, min = NA_real_, max = NA_real_)
      deg_cv <- NA_real_; deg_sk <- NA_real_; dmin <- NA_real_; dmax <- NA_real_
    }

    sum_row <- data.frame(
      layer         = ly,
      n_nodes       = n_nodes,
      n_edges       = n_edges,
      density       = dens,
      n_components  = comps,
      avg_path_len  = apl,
      diameter      = diam,
      degree_mean   = deg_s["mean"], degree_sd = deg_s["sd"],
      degree_median = deg_s["median"], degree_min = dmin, degree_max = dmax,
      degree_cv     = deg_cv, degree_skew = deg_sk,
      bet_mean      = if (nrow(nodes_tbl)) .summ(nodes_tbl$betweenness)["mean"] else NA_real_,
      eig_mean      = if (nrow(nodes_tbl)) .summ(nodes_tbl$eigenvector)["mean"] else NA_real_,
      close_mean    = if (nrow(nodes_tbl)) .summ(nodes_tbl$closeness)["mean"] else NA_real_,
      pr_mean       = if (nrow(nodes_tbl)) .summ(nodes_tbl$pagerank)["mean"] else NA_real_,
      clust_mean    = if (nrow(nodes_tbl)) .summ(nodes_tbl$clustering)["mean"] else NA_real_,
      core_mean     = if (nrow(nodes_tbl)) .summ(nodes_tbl$coreness)["mean"] else NA_real_,
      stringsAsFactors = FALSE
    )
    summary_rows[[ly]] <- sum_row

    # (optional) write per-layer node CSV immediately (only if you prefer streaming)
    # We keep all writing to the end to produce either XLSX or per-layer CSV set.
  }

  summary_df <- do.call(rbind, summary_rows)

  # --------- save outputs inside the run folder ----------
  summary_csv <- file.path(run_dir, "layer_metrics_summary.csv")
  utils::write.csv(summary_df, summary_csv, row.names = FALSE)

  nodes_xlsx    <- NULL
  per_layer_dir <- NULL
  if (requireNamespace("writexl", quietly = TRUE)) {
    nodes_xlsx <- file.path(run_dir, "layer_metrics_nodes.xlsx")
    writexl::write_xlsx(per_layer_nodes, nodes_xlsx)
  } else {
    per_layer_dir <- file.path(run_dir, "nodes_per_layer_csv")
    dir.create(per_layer_dir, recursive = TRUE, showWarnings = FALSE)
    for (ly in names(per_layer_nodes)) {
      f <- file.path(per_layer_dir, paste0(.sanitize(ly), "_nodes.csv"))
      utils::write.csv(per_layer_nodes[[ly]], f, row.names = FALSE)
    }
  }

  # --------- messages with absolute paths ----------
  run_dir_abs    <- normalizePath(run_dir,    winslash = "/", mustWork = FALSE)
  summary_csv_abs<- normalizePath(summary_csv,winslash = "/", mustWork = FALSE)
  if (!is.null(nodes_xlsx))    nodes_xlsx    <- normalizePath(nodes_xlsx,    winslash = "/", mustWork = FALSE)
  if (!is.null(per_layer_dir)) per_layer_dir <- normalizePath(per_layer_dir, winslash = "/", mustWork = FALSE)

  if (verbose) {
    message("📂 Run folder: ", run_dir_abs)
    message("✅ Summary CSV: ", summary_csv_abs)
    if (!is.null(nodes_xlsx)) {
      message("✅ Node metrics workbook: ", nodes_xlsx)
    } else {
      message("ℹ️ Per-layer node CSVs: ", per_layer_dir)
    }
  }

  invisible(list(
    summary = summary_df,
    nodes   = per_layer_nodes,
    paths   = list(
      run_dir      = run_dir_abs,
      summary_csv  = summary_csv_abs,
      nodes_xlsx   = nodes_xlsx,
      per_layer_dir= per_layer_dir
    )
  ))
}


# ---- helper (same as your other utilities) ----------------------------------
.normalize_keys <- function(x, steps = c("strip_version","trim","tolower")) {
  x <- as.character(x)
  for (s in steps) {
    x <- switch(s,
                "strip_version" = sub("\\.\\d+$", "", x),
                "trim"          = trimws(x),
                "tolower"       = tolower(x),
                "toupper"       = toupper(x),
                "rm_dash"       = gsub("[-_]", "", x),
                "rm_punct"      = gsub("[[:punct:]]+", "", x),
                "alnum"         = gsub("[^[:alnum:]]+", "", x),
                x
    )
  }
  x
}
# ---------------------------------------------------------------------------
# 14) Degrees for any feature list across layers
# ---------------------------------------------------------------------------


#' Degrees for a selected feature set across layers (robust, \pkg{igraph}-based)
#'
#' @description
#' Compute per‑layer **degrees** for an arbitrary set of actors (e.g., genes of
#' interest). The function extracts all **intra‑layer** edges from a
#' \pkg{multinet} object, builds one \pkg{igraph} graph per layer, and reports
#' each requested actor’s degree in that layer. Actors not present in a given
#' layer are added as **isolates** so their degree is reported as `0` rather than
#' `NA`. Results can be returned in **wide** (default) or **long** shape and
#' optionally written to a timestamped CSV in your project results folder.
#'
#' @details
#' **Input normalisation and scope**
#' - Edges are retrieved once via `multinet::edges_ml(net)` and coerced to a data
#'   frame. If `from_layer` and `to_layer` exist, only **intra‑layer** edges are
#'   kept (`from_layer == to_layer`) and a single `layer` column is created from
#'   `from_layer`. Otherwise an existing `layer`/`Layer` column is used; if none
#'   exists, a single layer `"L1"` is assumed.
#' - Endpoint columns are identified robustly by name (tries common pairs such as
#'   `from`/`to`, `source`/`target`, `actor1`/`actor2`). If these are absent, the
#'   **first two character columns** (excluding layer‑like columns) are used.
#'
#' **Actor matching and normalisation**
#' - Actor names in `featureList` are matched to vertex names per layer using the
#'   same normalisation pipeline you use elsewhere (default:
#'   `c("strip_version","trim","tolower")`), applied via the helper
#'   `.normalize_keys()`. This improves robustness to minor ID formatting
#'   differences.
#'
#' **Degree definition**
#' - Degrees are computed with `igraph::degree()` and `loops = FALSE`.
#' - When `directed = FALSE` (default), mode is effectively `"all"`. When
#'   `directed = TRUE`, `mode` can be `"all"`, `"in"`, or `"out"`.
#'
#' **Output shape and persistence**
#' - By default (`return_long = FALSE`), a **wide** table is returned with one
#'   row per `Layer` and one column per `Feature`. If \pkg{tidyr} is available,
#'   reshaping is done with `pivot_wider`; otherwise a base‑R fallback is used.
#' - With `return_long = TRUE`, a **long** table is returned with columns
#'   `Layer`, `Feature`, `Degree`.
#' - If `write_csv = TRUE`, results are written under
#'   `getOption("mlnet.results_dir","omicsDNA_results")` when `output_file` is a
#'   relative path or left as the default placeholder (a timestamped name is
#'   then generated automatically). Absolute paths are respected. The absolute
#'   file path is attached as `attr(x, "file")`.
#'
#' **Caveats**
#' - The heuristic for endpoint detection (first two character columns) can be
#'   fooled if your edge table contains other character columns before the
#'   endpoints; reorder columns upstream when in doubt.
#'
#' @param net A \pkg{multinet} object (must work with `multinet::edges_ml()`).
#' @param featureList Character vector of actor names to measure (duplicates are
#'   ignored; input order is preserved in the output).
#' @param layers Optional character vector of layer names to include. Default:
#'   all layers found in `edges_ml(net)` after intra‑layer filtering.
#' @param directed Logical; build directed per‑layer graphs? Default `FALSE`.
#' @param mode Degree mode used when `directed = TRUE`: one of `"all"` (default),
#'   `"in"`, or `"out"`. Ignored when `directed = FALSE`.
#' @param actor_normalize Character vector of normalisation steps applied to both
#'   sides before matching (`"strip_version"`, `"trim"`, `"tolower"`, `"toupper"`,
#'   `"rm_dash"`, `"rm_punct"`, `"alnum"`). Default
#'   `c("strip_version","trim","tolower")`.
#' @param write_csv Logical; if `TRUE`, write a CSV (wide by default; long if
#'   `return_long = TRUE`). Default `FALSE`.
#' @param output_file File name for the CSV. If relative, it is saved under
#'   `getOption("mlnet.results_dir","omicsDNA_results")`. If `NULL` or left as
#'   the default placeholder (`"features-degree_byLayer.csv"`), an informative,
#'   timestamped name is generated automatically in that folder.
#' @param results_dir Directory used when `output_file` is relative or
#'   auto‑generated. Default `getOption("mlnet.results_dir","omicsDNA_results")`.
#' @param return_long Logical; return long format (`Layer`, `Feature`, `Degree`)?
#'   Default `FALSE` (wide).
#' @param verbose Logical; print the saved file path when writing. Default `TRUE`.
#'
#' @return A data frame of degrees by layer × feature:
#' - **Wide** (default): columns `Layer`, `<Feature1>`, `<Feature2>`, …
#' - **Long** (`return_long = TRUE`): columns `Layer`, `Feature`, `Degree`.
#' If `write_csv = TRUE`, the absolute file path is attached as `attr(x, "file")`.
#'
#' @examples
#' \dontrun{
#' # 1) lncRNAs across all layers, undirected (wide) + save under results dir
#' lnc_list  <- unique(subset(comm_annot, GeneType == "lncRNA")$actor)
#' all_layers <- multinet::layers_ml(net)
#' deg_lnc <- get_FeatureDeg(
#'   net, featureList = lnc_list, layers = all_layers,
#'   directed   = FALSE, mode = "all",
#'   write_csv  = TRUE,
#'   output_file = "lncRNA_degree_byLayer.csv"
#' )
#'
#' # 2) Protein‑coding, long shape (no file)
#' pc_list <- unique(subset(comm_annot, GeneType == "protein_coding")$actor)
#' deg_pc  <- get_FeatureDeg(
#'   net, featureList = pc_list, layers = all_layers,
#'   return_long = TRUE
#' )
#'
#' # 3) Directed in‑degree for a custom list on selected layers
#' sel_layers <- c("Young","Old")
#' deg_in <- get_FeatureDeg(
#'   net, featureList = c("TP53","MYC","EGFR"),
#'   layers = sel_layers, directed = TRUE, mode = "in"
#' )
#' }
#'
#' @seealso
#'   \code{\link{layer_metrics}} for broader per‑layer summaries;
#'   \code{\link{sumComFeat}} and \code{\link{annotateCom}} for feature‑aware
#'   community summaries and annotations.
#'
#' @importFrom multinet edges_ml
#' @importFrom igraph graph_from_data_frame degree
#' @importFrom utils write.csv
#' @export
get_FeatureDeg <- function(
    net,
    featureList,
    layers          = NULL,
    directed        = FALSE,
    mode            = "all",
    actor_normalize = c("strip_version","trim","tolower"),
    write_csv       = FALSE,
    output_file     = "features-degree_byLayer.csv",
    results_dir     = getOption("mlnet.results_dir","omicsDNA_results"),
    return_long     = FALSE,
    verbose         = TRUE
) {
  stopifnot(length(featureList) >= 1)
  featureList <- as.character(featureList)
  featureList <- unique(featureList)  # avoid dup columns
  feat_key    <- .normalize_keys(featureList, steps = actor_normalize)

  if (directed && !(mode %in% c("all","in","out"))) {
    stop('When directed=TRUE, `mode` must be one of "all", "in", or "out".')
  }

  # --- read all edges once (robust to multinet versions)
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
  # --- identify endpoints
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
    if (length(char_cols) < 2) stop("Could not identify two endpoint columns in the edges table.")
    a_col <- nm[char_cols[1]]; b_col <- nm[char_cols[2]]
  }

  # --- normalize single 'layer' column and drop inter-layer edges
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

  # --- choose layers
  all_layers <- unique(Eall$layer)
  if (is.null(layers)) {
    layers <- all_layers
  } else {
    layers <- intersect(as.character(layers), all_layers)
  }
  if (!length(layers)) stop("No valid layers after intersection.")
  # preserve requested order
  layers <- as.character(layers)

  # --- split edges by layer
  edges_by_layer <- lapply(layers, function(ly) {
    ed <- Eall[Eall$layer == ly, , drop = FALSE]
    if (!nrow(ed)) {
      data.frame(from = character(), to = character(), stringsAsFactors = FALSE)
    } else {
      data.frame(from = as.character(ed[[a_col]]),
                 to   = as.character(ed[[b_col]]),
                 stringsAsFactors = FALSE)
    }
  })
  names(edges_by_layer) <- layers

  # --- compute degrees per layer (include features as isolates if needed)
  long_rows <- vector("list", length(layers))
  for (i in seq_along(layers)) {
    ly <- layers[i]
    ed <- edges_by_layer[[ly]]

    verts   <- unique(c(ed$from, ed$to, featureList))
    verts_df <- if (length(verts)) data.frame(name = verts, stringsAsFactors = FALSE) else NULL

    g <- igraph::graph_from_data_frame(ed, directed = directed, vertices = verts_df)

    v_names <- igraph::V(g)$name
    v_key   <- .normalize_keys(v_names, steps = actor_normalize)
    hit     <- match(feat_key, v_key)

    vals <- rep(NA_real_, length(featureList))
    present <- which(!is.na(hit))
    if (length(present)) {
      ids <- hit[present]
      vals[present] <- igraph::degree(g, v = ids, mode = if (directed) mode else "all", loops = FALSE)
    }
    # Any not matched (should be none, because we inserted isolates) remain NA -> set to 0
    vals[is.na(vals)] <- 0

    long_rows[[i]] <- data.frame(
      Layer   = ly,
      Feature = featureList,
      Degree  = as.numeric(vals),
      stringsAsFactors = FALSE
    )
  }
  long_df <- do.call(rbind, long_rows)

  # --- shape output (wide or long), preserve orders
  long_df$Layer   <- factor(long_df$Layer,   levels = layers)
  long_df$Feature <- factor(long_df$Feature, levels = featureList)

  out_df <- NULL
  if (!return_long) {
    if (requireNamespace("tidyr", quietly = TRUE)) {
      out_df <- tidyr::pivot_wider(long_df, names_from = Feature, values_from = Degree)
      out_df <- as.data.frame(out_df, stringsAsFactors = FALSE)
    } else {
      # base-R fallback
      w <- reshape(long_df, idvar = "Layer", timevar = "Feature", direction = "wide")
      # columns look like Degree.<feature>; clean names
      names(w) <- sub("^Degree\\.", "", names(w))
      # order by Layer factor
      w <- w[order(as.integer(w$Layer)), , drop = FALSE]
      rownames(w) <- NULL
      out_df <- w
    }
  } else {
    out_df <- long_df[order(as.integer(long_df$Layer)), , drop = FALSE]
    rownames(out_df) <- NULL
  }

  # --- save (to same results dir used elsewhere) -----------------------------
  if (isTRUE(write_csv)) {
    .is_abs <- function(p) grepl("^(/|[A-Za-z]:[\\/])", p)

    # Decide filename: if user kept the default placeholder or set NULL -> auto timestamp
    auto_name <- is.null(output_file) || identical(output_file, "features-degree_byLayer.csv")
    if (auto_name) {
      stamp <- format(Sys.time(), "%Y-%m-%d_%H%M%S")
      shape <- if (return_long) "long" else "wide"
      nL    <- length(layers)
      nF    <- length(featureList)
      output_file <- sprintf("features-degree_byLayer_%s_%dlayers_%dfeat_%s.csv", shape, nL, nF, stamp)
    }

    # If relative -> place under results_dir; if absolute -> use as is
    if (!.is_abs(output_file)) {
      if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
      out_path <- file.path(results_dir, output_file)
    } else {
      out_path <- output_file
      # ensure parent exists
      out_parent <- dirname(out_path)
      if (!dir.exists(out_parent)) dir.create(out_parent, recursive = TRUE, showWarnings = FALSE)
    }

    utils::write.csv(out_df, out_path, row.names = FALSE)
    out_abs <- normalizePath(out_path, winslash = "/", mustWork = FALSE)
    if (isTRUE(verbose)) message("Saved CSV: ", out_abs)
    attr(out_df, "file") <- out_abs
  }

  out_df
}

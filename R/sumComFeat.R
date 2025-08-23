# ---- small helper for consistent key matching (unchanged) -------------------
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
# 13) Summarize community membership for any feature type (short name)
# ---------------------------------------------------------------------------

#' Summarise community membership for selected feature types (pure data‑frame workflow)
#'
#' @description
#' Filter a community assignment table (e.g., from `detectCom()` plus your
#' node‑metadata merge) to one or more **feature types** (e.g., `"lncRNA"`,
#' `"protein_coding"`, `"TF"`), then produce three concise summaries:
#' (1) **per‑actor** membership across communities and layers;
#' (2) **per‑community** composition (actors and layers); and
#' (3) **counts** of unique actors per `CID × layer`.
#' This function does **not** call \pkg{multinet}: it operates purely on data
#' frames, making it fast, reproducible, and script‑friendly.
#'
#' @details
#' **Input expectations and harmonisation**
#' - `communities` must contain:
#'   - an actor/ID column: preferably `actor`; if missing, `GeneName` is copied to `actor`;
#'   - a `layer` column; and
#'   - either `cid` **or** `com`. If only `com` is present, a numeric `cid` is
#'     generated via `as.integer(factor(com))`. If only `cid` is present, a label
#'     `com = paste0("C", cid)` is generated.
#' - A feature attribute column named in `feature_col` (e.g., `"GeneType"`) must
#'   be present; it is used for the type filter.
#'
#' **Filtering and normalisation**
#' - Rows are kept only if their `feature_col` value is in `feature_type`
#'   (optionally case‑insensitive via `ignore_case`). You may further restrict to
#'   a subset of `layers`.
#' - An auxiliary `actor_clean` is created via `.normalize_keys()` using
#'   `actor_normalize` (default `c("strip_version","trim","tolower")`) **only to
#'   provide a deterministic display order**; it does not affect grouping or counts.
#'
#' **Outputs produced**
#' 1) **Per‑actor summary** (`by_actor`): columns
#'    `actor, CIDs, Layers, n_cids, n_layers`, ordered by `actor_clean`.
#' 2) **Per‑community summary** (`by_cid`): columns
#'    `cid, com, Features, Layers, n_features, n_layers`, where `Features`
#'    is a comma‑separated list of unique actor IDs belonging to the community
#'    (after filtering), and `Layers` lists the layers in which those actors occur.
#' 3) **Counts by community × layer** (`counts_cid_layer`): columns
#'    `cid, com, layer, n`, where `n` is the number of **unique actors** in that
#'    community‑layer slice.
#'
#' **Saving to disk**
#' - When `write_csv = TRUE`, three timestamped CSVs are written under
#'   `results_dir` (default:
#'   `getOption("mlnet.results_dir", "omicsDNA_results")`), with file stems
#'   derived from `prefix` (or, if missing, from `feature_col` and `feature_type`):
#'   `"<prefix>_by-actor_<timestamp>.csv"`,
#'   `"<prefix>_by-cid_<timestamp>.csv"`,
#'   `"<prefix>_counts_cid-layer_<timestamp>.csv"`.
#' - The returned list then carries an attribute `"files"` with the absolute
#'   paths of the three CSVs.
#'
#' **Behaviour on empty selections**
#' - If no rows remain after the type/layer filters or required columns are
#'   missing, the function raises an informative error.
#'
#' @param communities Data frame with at least: `actor` (or `GeneName`),
#'   `layer`, one of `cid` or `com`, and the feature attribute column named by
#'   `feature_col`.
#' @param feature_col Character scalar; name of the feature attribute column in
#'   `communities` used for filtering (e.g., `"GeneType"`). Default `"GeneType"`.
#' @param feature_type Character vector of feature values to keep
#'   (e.g., `c("lncRNA", "TF")`). Required.
#' @param layers Optional character vector of layer names to retain. Default `NULL`
#'   (use all layers present).
#' @param ignore_case Logical; if `TRUE`, compare `feature_type` to `feature_col`
#'   values case‑insensitively. Default `FALSE`.
#' @param actor_normalize Character vector of normalisation steps used to build
#'   `actor_clean` for display ordering (does **not** affect grouping or counts).
#'   Supported steps: `"strip_version"`, `"trim"`, `"tolower"`, `"toupper"`,
#'   `"rm_dash"`, `"rm_punct"`, `"alnum"`. Default `c("strip_version","trim","tolower")`.
#' @param write_csv Logical; if `TRUE`, write the three CSVs described above.
#'   Default `FALSE`.
#' @param prefix Optional filename stem for the CSVs. If `NULL`, a stem is
#'   auto‑generated from `feature_col` and `feature_type` (sanitised).
#' @param results_dir Output directory for CSVs when `write_csv = TRUE`. Default
#'   `getOption("mlnet.results_dir", "omicsDNA_results")`.
#' @param verbose Logical; print short status messages (e.g., saved file paths).
#'   Default `TRUE`.
#'
#' @return (Invisibly) a list with three data frames:
#' \itemize{
#'   \item `by_actor` — per‑actor membership summary,
#'   \item `by_cid` — per‑community composition summary,
#'   \item `counts_cid_layer` — unique‑actor counts by community × layer.
#' }
#' If `write_csv = TRUE`, the list has an attribute `"files"` containing absolute
#' paths to the three CSVs.
#'
#' @examples
#' \dontrun{
#' # Keep TFs and lncRNAs across all layers, then write summaries to disk
#' out <- sumComFeat(
#'   communities   = comm_annot,        # e.g., detectCom() + annotateCom()
#'   feature_col   = "GeneType",
#'   feature_type  = c("TF", "lncRNA"),
#'   write_csv     = TRUE
#' )
#' out$by_actor[1:5, ]
#' out$by_cid[1:5, ]
#' out$counts_cid_layer[1:5, ]
#' attr(out, "files")
#' }
#'
#' @seealso
#'   \code{\link{annotateCom}} to add feature annotations before summarising;
#'   \code{\link{detectCom}} to generate community assignments.
#'
#' @importFrom utils write.csv
#' @importFrom stats aggregate
#' @export
sumComFeat <- function(
    communities,
    feature_col     = "GeneType",
    feature_type,
    layers          = NULL,
    ignore_case     = FALSE,
    actor_normalize = c("strip_version","trim","tolower"),
    write_csv       = FALSE,
    prefix          = NULL,
    results_dir     = getOption("mlnet.results_dir","omicsDNA_results"),
    verbose         = TRUE
) {
  stopifnot(is.data.frame(communities))

  # Ensure required columns are present / harmonized
  if (!("actor" %in% names(communities))) {
    if ("GeneName" %in% names(communities)) {
      communities$actor <- as.character(communities$GeneName)
    } else {
      stop("`communities` must contain 'actor' (or 'GeneName').")
    }
  }
  if (!("layer" %in% names(communities))) {
    stop("`communities` must contain a 'layer' column.")
  }
  if (!(feature_col %in% names(communities))) {
    stop("`communities` must contain feature column '", feature_col, "'.")
  }
  if (!("cid" %in% names(communities)) && ("com" %in% names(communities))) {
    communities$cid <- as.integer(factor(communities$com))
  }
  if (!("com" %in% names(communities)) && ("cid" %in% names(communities))) {
    communities$com <- paste0("C", communities$cid)
  }
  if (!all(c("cid","com") %in% names(communities))) {
    stop("`communities` must contain either 'cid' or 'com' (or both).")
  }

  dat <- communities

  # Optional layer restriction
  if (!is.null(layers)) {
    layers <- as.character(layers)
    dat <- dat[dat$layer %in% layers, , drop = FALSE]
  }

  # Filter by feature type
  fvals <- as.character(dat[[feature_col]])
  keep <- if (isTRUE(ignore_case)) tolower(fvals) %in% tolower(feature_type) else fvals %in% feature_type
  dat <- dat[keep, , drop = FALSE]
  if (!nrow(dat)) stop("No rows remain after filtering by ", feature_col, " in feature_type.")

  # Deduplicate rows on (actor, layer, cid) to avoid double counting
  dat <- unique(dat[, c("actor","layer","cid","com", feature_col), drop = FALSE])

  # Cleaned label for display order only
  dat$actor_clean <- .normalize_keys(dat$actor, steps = actor_normalize)

  ## 1) per-actor summary
  split_actor <- split(dat, dat$actor, drop = TRUE)
  by_actor <- data.frame(
    actor     = names(split_actor),
    CIDs      = vapply(split_actor, function(df) paste(sort(unique(df$cid)), collapse = ", "), character(1)),
    Layers    = vapply(split_actor, function(df) paste(sort(unique(df$layer)), collapse = ", "), character(1)),
    n_cids    = vapply(split_actor, function(df) length(unique(df$cid)),   integer(1)),
    n_layers  = vapply(split_actor, function(df) length(unique(df$layer)), integer(1)),
    stringsAsFactors = FALSE
  )
  # order by normalized actor label
  ord <- order(.normalize_keys(by_actor$actor, steps = actor_normalize))
  by_actor <- by_actor[ord, , drop = FALSE]

  ## 2) per-community summary
  key_cid <- paste(dat$cid, dat$com, sep = "\r")
  split_cid <- split(dat, key_cid, drop = TRUE)
  # extract cid/com back from the key
  cid_vals <- vapply(split_cid, function(df) unique(df$cid)[1], integer(1))
  com_vals <- vapply(split_cid, function(df) unique(df$com)[1], character(1))
  by_cid <- data.frame(
    cid        = as.integer(cid_vals),
    com        = as.character(com_vals),
    Features   = vapply(split_cid, function(df) paste(sort(unique(df$actor)), collapse = ", "), character(1)),
    Layers     = vapply(split_cid, function(df) paste(sort(unique(df$layer)), collapse = ", "), character(1)),
    n_features = vapply(split_cid, function(df) length(unique(df$actor)),   integer(1)),
    n_layers   = vapply(split_cid, function(df) length(unique(df$layer)),   integer(1)),
    stringsAsFactors = FALSE
  )
  by_cid <- by_cid[order(by_cid$cid), , drop = FALSE]

  ## 3) counts by CID × layer
  agg <- aggregate(actor ~ cid + com + layer, data = dat, FUN = function(v) length(unique(v)))
  names(agg)[names(agg) == "actor"] <- "n"
  counts_cid_layer <- agg[order(agg$cid, agg$layer), , drop = FALSE]

  # Optional CSV export (to the same place as previous files)
  files <- NULL
  if (isTRUE(write_csv)) {
    if (is.null(prefix)) {
      ft_tag <- paste0(feature_col, "-", paste(feature_type, collapse = "_"))
      ft_tag <- gsub("[^A-Za-z0-9._-]+", "_", ft_tag)
      prefix <- paste0("summ_", ft_tag)
    }
    if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
    stamp <- format(Sys.time(), "%Y-%m-%d_%H%M%S")

    f_actor  <- file.path(results_dir, paste0(prefix, "_by-actor_",          stamp, ".csv"))
    f_cid    <- file.path(results_dir, paste0(prefix, "_by-cid_",            stamp, ".csv"))
    f_counts <- file.path(results_dir, paste0(prefix, "_counts_cid-layer_",  stamp, ".csv"))

    utils::write.csv(by_actor,         f_actor,  row.names = FALSE)
    utils::write.csv(by_cid,           f_cid,    row.names = FALSE)
    utils::write.csv(counts_cid_layer, f_counts, row.names = FALSE)

    f_actor  <- normalizePath(f_actor,  winslash = "/", mustWork = FALSE)
    f_cid    <- normalizePath(f_cid,    winslash = "/", mustWork = FALSE)
    f_counts <- normalizePath(f_counts, winslash = "/", mustWork = FALSE)

    if (isTRUE(verbose)) {
      message("Saved CSVs:")
      message("  • by_actor         → ", f_actor)
      message("  • by_cid           → ", f_cid)
      message("  • counts_cid_layer → ", f_counts)
    }

    files <- list(by_actor = f_actor, by_cid = f_cid, counts_cid_layer = f_counts)
  }

  out <- list(by_actor = by_actor, by_cid = by_cid, counts_cid_layer = counts_cid_layer)
  if (!is.null(files)) attr(out, "files") <- files
  invisible(out)
}

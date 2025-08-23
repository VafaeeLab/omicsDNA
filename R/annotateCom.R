

# --- helpers (same as before, safe) ------------------------------------------
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
# 12 - annotateCom: shorter name, comprehensive docs, robust saving
# ---------------------------------------------------------------------------

#' Annotate multilayer community rows with a node attribute (pure data-frame join)
#'
#' @description
#' Add a node attribute (e.g., a gene annotation such as `"GeneType"`) to a
#' **community assignment table** by matching actor IDs in `communities` to an
#' identifier column in `nodesMetadata`. The operation is **metadata‑only**:
#' it never calls \pkg{multinet}; it normalises IDs, performs a vectorised join,
#' reports coverage, can fill unmatched rows with a sentinel value, and can write
#' the annotated table to a timestamped CSV under your project’s results folder.
#'
#' @details
#' **Expected inputs**
#' - `communities` must have an actor/ID column:
#'   - Preferably named `actor`; if absent, the function looks for `GeneName` and
#'     copies it into `actor`. If neither is present, it errors.
#'   - Any other columns (e.g., `cid`, `com`, `layer`) are preserved verbatim.
#' - `nodesMetadata` must provide:
#'   - a key column `featureID_col` (e.g., `"GeneName"`) whose values match
#'     `communities$actor` *after normalisation*, and
#'   - the attribute column `attribute` to be joined (e.g., `"GeneType"`).
#'
#' **How matching works**
#' 1. Both sides are normalised with `.normalize_keys()` using the sequence in
#'    `actor_normalize` (default: `c("strip_version","trim","tolower")`).
#' 2. `nodesMetadata` is deduplicated by the normalised key (first occurrence wins).
#' 3. A vectorised `match()` maps community keys to metadata keys and collects
#'    the corresponding attribute values.
#' 4. If `fill_missing` is not `NULL`, unmatched values are replaced by this
#'    sentinel (default `"Unknown"`).
#'
#' **Saving behaviour**
#' - If `write_csv = TRUE`, the result is written to
#'   `file.path(getOption("mlnet.results_dir", "omicsDNA_results"), <basename>_<timestamp>.csv)`.
#' - The file stem defaults to `paste0("communities_with_", attribute)` and can be
#'   overridden via `output_basename` (or its alias `output_prefix` for
#'   backward compatibility). A timestamp in the form `YYYY-mm-dd_HHMMSS` is appended.
#'
#' **Type note**
#' - The appended attribute column is **character** (internal coercion via
#'   `as.character()`), even if the source column in `nodesMetadata` was numeric.
#'
#' @param communities Data frame of community assignments. Must contain `actor`
#'   (or `GeneName`, which will be copied into `actor`). Other columns are
#'   preserved.
#' @param nodesMetadata Data frame containing the join key and the attribute to
#'   merge (e.g., a gene annotation table).
#' @param featureID_col Character scalar; name of the identifier column in
#'   `nodesMetadata` that corresponds to `communities$actor`. Default `"GeneName"`.
#' @param attribute Character scalar; name of the attribute column in
#'   `nodesMetadata` to attach (e.g., `"GeneType"`).
#' @param actor_normalize Character vector of normalisation steps applied to both
#'   sides before matching. Supported steps: `"strip_version"`, `"trim"`,
#'   `"tolower"`, `"toupper"`, `"rm_dash"`, `"rm_punct"`, `"alnum"`. Default
#'   `c("strip_version","trim","tolower")`.
#' @param fill_missing Value used for unmatched rows; set `NULL` to leave `NA`.
#'   Default `"Unknown"`.
#' @param write_csv Logical; if `TRUE`, write a timestamped CSV of the annotated
#'   table to disk. Default `FALSE`.
#' @param results_dir Directory where the CSV is written when `write_csv = TRUE`.
#'   Default `getOption("mlnet.results_dir", "omicsDNA_results")`.
#' @param output_basename Basename (without timestamp or extension) for the CSV.
#'   If `NULL`, uses `paste0("communities_with_", attribute)`.
#' @param output_prefix Deprecated alias of `output_basename` (kept for backward
#'   compatibility).
#' @param verbose Logical; print a coverage summary and, if applicable, the saved
#'   file path. Default `TRUE`.
#'
#' @return The input `communities` with one new column named exactly as
#'   `attribute` (e.g., `"GeneType"`). When `write_csv = TRUE`, the absolute file
#'   path is also attached as `attr(x, "file")`. Coverage (fraction of
#'   non‑missing, non‑empty annotations before filling) is reported to the console
#'   when `verbose = TRUE`.
#'
#' @examples
#' \dontrun{
#' # Return-only (no file written)
#' comm_annot <- annotateCom(
#'   communities   = comm,
#'   nodesMetadata = genes_info,
#'   featureID_col = "GeneName",
#'   attribute     = "GeneType"
#' )
#'
#' # Write a CSV into your project results folder
#' # options(mlnet.results_dir = "/path/to/results")
#' comm_annot <- annotateCom(
#'   communities     = comm,
#'   nodesMetadata   = genes_info,
#'   featureID_col   = "GeneName",
#'   attribute       = "GeneType",
#'   write_csv       = TRUE,
#'   output_basename = "communities_with_GeneType"
#' )
#'
#' # Backward-compatible alias (same as output_basename)
#' comm_annot <- annotateCom(
#'   communities   = comm,
#'   nodesMetadata = genes_info,
#'   featureID_col = "GeneName",
#'   attribute     = "GeneType",
#'   write_csv     = TRUE,
#'   output_prefix = "communities_with_GeneType"
#' )
#' }
#'
#' @seealso
#'   \code{\link{add_network_attributes}} for attaching attributes inside a
#'   \pkg{multinet} object; \code{\link{build_multiNet}} for assembling layers.
#'
#' @importFrom utils write.csv
#' @export
annotateCom <- function(
    communities,
    nodesMetadata,
    featureID_col    = "GeneName",
    attribute        = "GeneType",
    actor_normalize  = c("strip_version","trim","tolower"),
    fill_missing     = "Unknown",
    write_csv        = FALSE,
    results_dir      = getOption("mlnet.results_dir","omicsDNA_results"),
    output_basename  = NULL,
    output_prefix    = NULL,  # alias for backward compatibility
    verbose          = TRUE
) {
  stopifnot(is.data.frame(communities), is.data.frame(nodesMetadata))

  # Ensure 'actor' column is present
  if (!("actor" %in% names(communities))) {
    if ("GeneName" %in% names(communities)) {
      communities$actor <- as.character(communities$GeneName)
    } else {
      stop("`communities` must contain 'actor' (or 'GeneName').")
    }
  }

  # Basic checks on metadata
  if (!(featureID_col %in% names(nodesMetadata)))
    stop("`nodesMetadata` must contain key column: ", featureID_col)
  if (!(attribute %in% names(nodesMetadata)))
    stop("`nodesMetadata` must contain attribute column: ", attribute)

  # Normalize keys on both sides
  comm_key <- .normalize_keys(communities$actor, steps = actor_normalize)
  meta_key <- .normalize_keys(as.character(nodesMetadata[[featureID_col]]), steps = actor_normalize)
  meta_val <- as.character(nodesMetadata[[attribute]])

  # Deduplicate metadata by normalized key (first occurrence wins)
  keep <- !duplicated(meta_key)
  meta_key <- meta_key[keep]
  meta_val <- meta_val[keep]

  # Vectorized match and value extraction
  hit <- match(comm_key, meta_key)
  out_vals <- rep(NA_character_, length(hit))
  ok <- which(!is.na(hit))
  if (length(ok)) out_vals[ok] <- meta_val[hit[ok]]

  # Coverage before filling
  coverage_raw <- mean(!is.na(out_vals) & out_vals != "")

  # Optional fill for missing annotations
  if (!is.null(fill_missing)) {
    miss <- is.na(out_vals) | out_vals == ""
    if (any(miss)) out_vals[miss] <- as.character(fill_missing)
  }

  # Append the attribute column
  out <- communities
  out[[attribute]] <- out_vals

  # Reporting
  if (isTRUE(verbose)) {
    message(sprintf("Coverage from metadata (pre-fill): %.1f%% (%d / %d)",
                    100 * coverage_raw,
                    sum(!is.na(out_vals) & out_vals != "" & out_vals != fill_missing),
                    length(out_vals)))
  }

  # Optional CSV export
  file_path <- NULL
  if (isTRUE(write_csv)) {
    if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
    if (is.null(output_basename) && !is.null(output_prefix)) {
      # backward-compatible alias
      output_basename <- output_prefix
    }
    if (is.null(output_basename)) {
      output_basename <- paste0("communities_with_", attribute)
    }
    stamp <- format(Sys.time(), "%Y-%m-%d_%H%M%S")
    file_path <- file.path(results_dir, paste0(output_basename, "_", stamp, ".csv"))
    utils::write.csv(out, file_path, row.names = FALSE)
    if (isTRUE(verbose)) message("Saved CSV: ", normalizePath(file_path, winslash = "/", mustWork = FALSE))
    attr(out, "file") <- file_path
  }

  out
}

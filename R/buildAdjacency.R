#' Build correlation‑filtered adjacency matrices by group:
#'
#' Construct one or more feature‑by‑feature networks for each group of samples
#' by correlating feature profiles across samples in that group and keeping only
#' edges that satisfy **both** an absolute correlation threshold and a p‑value
#' cutoff. Optionally, repeat the procedure with balanced draws per group to
#' produce multiple resampled networks (useful for consensus or stability work).
#'
#' @description
#' The function is domain‑agnostic:
#' - *Omics:* features = genes/transcripts; samples = subjects/cells; groups = condition/age.
#' - *Social/behavioral:* features = variables; samples = respondents; groups = cohorts.
#' - *Sensors/engineering:* features = channels; samples = runs/windows; groups = modes.
#'
#' @details
#' **Input layout**
#' - `dataMatrix` must be numeric with **features in rows** and **samples in columns**.
#' - `sample_metadata` must contain a grouping column (`group_col`) and a sample ID column
#'   (`sample_col`) whose values match `colnames(dataMatrix)`. Samples not found are ignored.
#' - The analysis is restricted to `feature_ids` (order is preserved). Requested features
#'   not present in `dataMatrix` are zero‑padded in the final adjacency so all outputs
#'   share identical dimensions and ordering.
#'
#' **Computation per group**
#' 1. Assemble the group’s samples (optionally sub‑sampled if resampling is enabled).
#' 2. (If available) run `WGCNA::goodSamplesGenes()` to drop problematic samples/features.
#' 3. Compute correlations and p‑values with `WGCNA::corAndPvalue()` when **WGCNA** is
#'    installed; otherwise fall back to `Hmisc::rcorr()`. Diagonals are set to 0 (corr) and 1 (p).
#' 4. Form the adjacency by zeroing entries that fail either filter:
#'    `abs(corr) < corr_threshold` **or** `p > pval_cutoff`.
#' 5. Expand to the full `feature_ids × feature_ids` grid; missing features remain zero.
#'
#' **Resampling**
#' - When `resample = TRUE`, each repeat draws up to `samples_per_group` samples
#'   **without replacement** within each group (`min(samples_per_group, n_in_group)`).
#' - Groups with too few samples still contribute (using all available samples).
#'
#' **Size requirements & NA handling**
#' - Requires **≥ 3 samples** per group and **≥ 2 features** after QC; otherwise a zero
#'   matrix is returned for that group.
#' - Non‑numeric columns are coerced to numeric; coercion may introduce `NA`s.
#'   If **WGCNA** is present, QC removes egregious rows/columns; otherwise the fallback
#'   correlation uses pairwise information via `Hmisc::rcorr()`.
#'
#' **Output shape**
#' - If `resample = FALSE`: a **named list** of adjacency matrices (one per group).
#' - If `resample = TRUE`: a **list of length `n_repeats`**, where each element is a
#'   named sub‑list of groups (`out[[repeat]][[group]]`).
#'
#' **Persistence**
#' - When `save_rds = TRUE`, the object is written to
#'   `<out_dir>/<file_prefix>_(single|boot)_YYYYMMDD-HHMMSS.rds` and also returned.
#'   The absolute file path is attached as `attr(result, "rds_file")`.
#'
#' @param dataMatrix Numeric matrix/data.frame; **features in rows**, **samples in columns**.
#' @param sample_metadata Data frame with at least `group_col` and `sample_col`.
#' @param feature_ids Character vector of feature IDs to include (should match `rownames(dataMatrix)`);
#'   missing features are zero‑padded in outputs so dimensions are consistent. IDs should be unique.
#' @param feature_metadata Optional feature annotations (currently not used; kept for API symmetry).
#' @param group_col Name of the grouping column in `sample_metadata`. Default `"group"`.
#' @param sample_col Name of the sample ID column in `sample_metadata`. Default `"sample"`.
#' @param feature_id_col Column name of IDs in `feature_metadata` (reserved; currently unused).
#' @param cor_method Correlation type: `"spearman"` or `"pearson"`. Default `"spearman"`.
#' @param pval_cutoff P‑value cutoff for keeping an edge. Default `0.05`.
#' @param corr_threshold Absolute correlation threshold for keeping an edge. Default `0.7`.
#' @param resample Logical; enable balanced draws per group. Default `TRUE`.
#' @param samples_per_group Target samples drawn per group when resampling. Default `10`.
#' @param n_repeats Number of resampling repeats. Default `50`.
#' @param verbose Logical; print progress/messages. Default `TRUE`.
#' @param save_rds Logical; write result to an `.rds` file. Default `TRUE`.
#' @param out_dir Directory for saved files. Default `file.path(getwd(), "omicsDNA_results")`.
#' @param file_prefix Basename for saved files (no extension). Default `"adjacency"`.
#' @param compress Compression passed to `saveRDS()`: `TRUE/FALSE` or `"xz"`, `"gzip"`, `"bzip2"`.
#'   Default `"xz"`.
#'
#' @return A list of adjacency matrices as described under **Output shape**. The return
#' object carries `attr(x, "rds_file")` with the saved file path (or `NULL` if `save_rds = FALSE`).
#'
#' @section Tips
#' - Call `set.seed()` for reproducible resampling.
#' - For very large `feature_ids`, pre‑filter (e.g., top DE features) to keep the
#'   number of pairwise correlations manageable.
#'
#' @examples
#' # Minimal toy example (no resampling)
#' set.seed(1)
#' X <- matrix(rnorm(5 * 6), nrow = 5,
#'             dimnames = list(paste0("g", 1:5), paste0("s", 1:6)))
#' meta <- data.frame(sample = colnames(X),
#'                    group  = rep(c("A","B"), each = 3))
#' adj <- buildAdjacency(
#'   dataMatrix        = X,
#'   sample_metadata   = meta,
#'   feature_ids       = rownames(X),
#'   group_col         = "group",
#'   sample_col        = "sample",
#'   resample          = FALSE,
#'   save_rds          = FALSE,
#'   verbose           = FALSE
#' )
#' lapply(adj, dim)  # one matrix per group
#'
#' @seealso \code{\link{edgesFromAdjacency}} to convert adjacency matrices to edge lists.
#' @export
buildAdjacency <- function(
    dataMatrix,
    sample_metadata,
    feature_ids,
    feature_metadata   = NULL,
    group_col          = "group",
    sample_col         = "sample",
    feature_id_col     = "feature_id",
    cor_method         = c("spearman", "pearson"),
    pval_cutoff        = 0.05,
    corr_threshold     = 0.7,
    resample           = TRUE,
    samples_per_group  = 10,
    n_repeats          = 50,
    verbose            = TRUE,
    save_rds           = TRUE,
    out_dir            = file.path(getwd(), "omicsDNA_results"),
    file_prefix        = "adjacency",
    compress           = "xz"
) {
  cor_method <- match.arg(cor_method)

  # ---- coerce to numeric matrix ----
  M <- as.matrix(dataMatrix)
  if (is.null(rownames(M))) stop("`dataMatrix` must have row names (feature IDs).")
  if (is.null(colnames(M))) stop("`dataMatrix` must have column names (sample IDs).")
  if (!is.numeric(M)) {
    df <- as.data.frame(M, check.names = FALSE, stringsAsFactors = FALSE)
    for (j in seq_len(ncol(df))) df[[j]] <- suppressWarnings(as.numeric(df[[j]]))
    M <- as.matrix(df); rownames(M) <- rownames(dataMatrix)
  }

  # ---- restrict to requested features (preserve order) ----
  found   <- intersect(feature_ids, rownames(M))
  missing <- setdiff(feature_ids, found)
  if (!length(found)) stop("None of `feature_ids` found in `dataMatrix` rownames.")
  if (length(missing) && verbose) {
    message("buildAdjacency: ", length(missing), " feature(s) missing; they will be zero-padded.")
  }
  M <- M[found, , drop = FALSE]

  # ---- validate metadata columns ----
  if (!all(c(group_col, sample_col) %in% names(sample_metadata))) {
    stop("`sample_metadata` must contain columns: ", group_col, " and ", sample_col, ".")
  }
  groups <- unique(as.character(sample_metadata[[group_col]]))
  groups <- groups[!is.na(groups)]
  if (!length(groups)) stop("No groups found in `sample_metadata` column ", group_col, ".")

  # ---- correlation + p-values helper (samples x features) ----
  .cor_with_p <- function(X) {
    Cp <- try({
      if (requireNamespace("WGCNA", quietly = TRUE)) {
        method_map <- c(pearson = "p", spearman = "s")
        WGCNA::corAndPvalue(X, method = method_map[[cor_method]], alternative = "two.sided")
      } else stop("noWGCNA")
    }, silent = TRUE)
    if (!inherits(Cp, "try-error")) {
      list(cor = Cp$cor, p = Cp$p)
    } else {
      if (!requireNamespace("Hmisc", quietly = TRUE)) {
        stop("Need either WGCNA or Hmisc installed to compute correlation p-values.")
      }
      rc <- Hmisc::rcorr(as.matrix(X), type = cor_method)
      list(cor = rc$r, p = rc$P)
    }
  }

  # ---- per-group adjacency (optionally with resampling) ----
  .adj_for_group <- function(group_id, draw_n = NULL, rep_index = NA_integer_) {
    sub_meta   <- sample_metadata[sample_metadata[[group_col]] == group_id, , drop = FALSE]
    sample_ids <- intersect(colnames(M), as.character(sub_meta[[sample_col]]))

    if (length(sample_ids) < 3L) {
      if (verbose) message("Group '", group_id, "' has <3 samples; returning zeros.")
      return(matrix(0, nrow = length(feature_ids), ncol = length(feature_ids),
                    dimnames = list(feature_ids, feature_ids)))
    }

    if (!is.null(draw_n)) {
      draw_n     <- min(draw_n, length(sample_ids))
      sample_ids <- sample(sample_ids, draw_n)
    }

    sub_M <- M[, sample_ids, drop = FALSE]  # features x samples
    X     <- t(sub_M)                       # samples x features

    # ---- Data quality check (prints removed features/samples) ----
    if (requireNamespace("WGCNA", quietly = TRUE)) {
      gsg <- WGCNA::goodSamplesGenes(X, verbose = 0)
      if (!gsg$allOK) {
        if (verbose) {
          if (sum(!gsg$goodGenes) > 0) {
            removed_features <- colnames(X)[!gsg$goodGenes]
            message("Removing features: ", paste(removed_features, collapse = ", "),
                    " | group = ", group_id,
                    if (!is.na(rep_index)) paste0(" | repeat = ", rep_index) else "")
          }
          if (sum(!gsg$goodSamples) > 0) {
            removed_samples <- rownames(X)[!gsg$goodSamples]
            message("Removing samples: ", paste(removed_samples, collapse = ", "),
                    " | group = ", group_id,
                    if (!is.na(rep_index)) paste0(" | repeat = ", rep_index) else "")
          }
        }
        X <- X[gsg$goodSamples, gsg$goodGenes, drop = FALSE]
      }
    }

    if (nrow(X) < 3L || ncol(X) < 2L) {
      if (verbose) message("After QC, group '", group_id, "' too small; returning zeros.")
      return(matrix(0, nrow = length(feature_ids), ncol = length(feature_ids),
                    dimnames = list(feature_ids, feature_ids)))
    }

    cp <- .cor_with_p(X)
    cor_mat <- cp$cor; p_mat <- cp$p
    diag(cor_mat) <- 0; diag(p_mat) <- 1

    keep <- (abs(cor_mat) >= corr_threshold) & (p_mat <= pval_cutoff)
    adj  <- cor_mat; adj[!keep] <- 0

    # expand to full feature_ids order
    present  <- colnames(X)
    full_adj <- matrix(0, nrow = length(feature_ids), ncol = length(feature_ids),
                       dimnames = list(feature_ids, feature_ids))
    common <- intersect(feature_ids, present)
    if (length(common)) full_adj[common, common] <- adj[common, common, drop = FALSE]
    full_adj
  }

  result <- NULL
  if (resample) {
    if (verbose) message("Resampling: ", n_repeats, " repeats; ", samples_per_group, " samples/group.")
    result <- vector("list", n_repeats)
    for (i in seq_len(n_repeats)) {
      mats <- lapply(groups, function(g) .adj_for_group(g, draw_n = samples_per_group, rep_index = i))
      names(mats) <- groups
      result[[i]] <- mats
      if (verbose) message("sampling ", i, " completed")
    }
  } else {
    if (verbose) message("Computing adjacency once per group (no resampling).")
    mats <- lapply(groups, function(g) .adj_for_group(g, draw_n = NULL, rep_index = NA_integer_))
    names(mats) <- groups
    result <- mats
  }

  # ---- save to RDS (optional) ----
  rds_path <- NULL
  if (isTRUE(save_rds)) {
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    stamp <- format(Sys.time(), "%Y%m%d-%H%M%S")
    kind  <- if (resample) "boot" else "single"
    prefix <- if (is.null(file_prefix) || !nzchar(file_prefix)) "adjacency" else file_prefix
    rds_path <- file.path(out_dir, sprintf("%s_%s_%s.rds", prefix, kind, stamp))
    comp <- if (isTRUE(compress) || isFALSE(compress)) compress else as.character(compress)
    saveRDS(result, rds_path, compress = comp)
    if (verbose) message("💾 Saved adjacency to: ", rds_path)
  }

  attr(result, "rds_file") <- rds_path
  return(result)
}

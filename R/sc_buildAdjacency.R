
# ---------------------------------------------------------------------
#  1.1 - sc_buildAdjacency(): Seurat → per–cell type gene–gene adjacency
# ---------------------------------------------------------------------

#' Build single-cell (per cell-type) correlation-filtered adjacency matrices.
#'
#' This helper constructs **one correlation network per cell type** from a Seurat
#' object by thresholding pairwise gene–gene correlations using both an absolute
#' correlation cutoff and a (optionally adjusted) p-value cutoff. It can perform
#' **balanced resampling** (bootstrap-like) so each repeat draws the same number
#' of cells per cell type, stabilising downstream edges/metrics. The output shape
#' mirrors `omicsDNA::buildAdjacency()` (per-layer matrices, optionally repeated),
#' so you can pipe results directly into `omicsDNA` edge extraction / consensus /
#' multilayer assembly workflows. See the `omicsDNA` overview for the layered
#' network workflow and typical end-to-end usage.  # (overview mirrors buildAdjacency)
#'
#' @examples
#' \donttest{
#'   # Minimal single-cell example using Seurat's toy dataset
#'   library(Seurat)
#'   data("pbmc_small")
#'
#'   # Create a simple cell-type label for demonstration (use your real labels)
#'   pbmc_small$celltype <- as.character(Idents(pbmc_small))
#'
#'   # Build one adjacency matrix per celltype (no resampling for speed)
#'   adj_list <- sc_buildAdjacency(
#'     seurat_obj          = pbmc_small,
#'     assay               = DefaultAssay(pbmc_small),
#'     slot                = "data",
#'     cell_type_col       = "celltype",
#'     cor_method          = "spearman",
#'     corr_threshold      = 0.6,
#'     pval_adjust         = "fdr",
#'     pval_cutoff         = 0.1,
#'     resample            = FALSE,
#'     save_rds            = FALSE,
#'     verbose             = TRUE
#'   )
#'
#'   # Each element is a (genes x genes) adjacency for a cell type
#'   str(adj_list, max.level = 1)
#' }
#'
#' @details
#' **What it does (per layer/cell type):**
#' 1. Selects a fixed set of features (genes) shared across layers/repeats:
#'    uses `feature_ids` if supplied; otherwise uses Seurat's `VariableFeatures()`
#'    if present; otherwise picks the top-variance genes (`n_var_features`).
#' 2. Optionally draws a **balanced** sample of cells per cell type
#'    (`samples_per_group`) for each repeat (`n_repeats`).
#' 3. Optionally filters bad samples/genes via `WGCNA::goodSamplesGenes()`
#'    (if WGCNA is installed); logs removals when `verbose = TRUE`.
#' 4. Computes pairwise gene–gene correlations and associated p-values using
#'    **WGCNA** (`corAndPvalue`) when available, otherwise **Hmisc** (`rcorr`).
#' 5. Adjusts the **upper-triangular** p-values only (to avoid double counting),
#'    mirrors them to the lower triangle, and sets diagonals appropriately.
#' 6. Builds an adjacency by keeping edges where `|correlation| ≥ corr_threshold`
#'    **and** `p ≤ pval_cutoff` (after adjustment if requested); others are zeroed.
#' 7. **Pads** back to the requested `feature_ids` order so every returned
#'    matrix has identical dimensions, even if some features are absent in a layer.
#'
#' **Notes & tips**
#' - **Reproducibility:** This function does not set a seed internally. For
#'   reproducible resampling, call `set.seed(<integer>)` before running.
#' - **Performance:** Correlation scales roughly with O(G^2) in the number of
#'   genes `G`. Supplying a targeted `feature_ids` (or relying on variable
#'   features) is strongly recommended for speed and memory.
#' - **QC fallback:** If after QC there are too few cells or genes in a layer,
#'   the function returns a zero matrix of the right shape and sets the
#'   `"n_cells"` attribute accordingly.
#'
#' Mirrors `omicsDNA::buildAdjacency()` semantics and return shape (per-layer matrices,
#' with optional repeated balanced draws). See the omicsDNA overview for the
#' layered network workflow.
#'
#' @param seurat_obj A **Seurat/SeuratObject** containing your single-cell data.
#'   Must have (i) an assay with a gene (rows) × cell (columns) matrix and
#'   (ii) a `@meta.data` column with cell-type labels. The function supports both
#'   `SeuratObject` and `Seurat` namespaces for `GetAssayData`, `DefaultAssay`,
#'   and `VariableFeatures`.
#' @param assay Character; which assay to read from. Defaults to
#'   `DefaultAssay(seurat_obj)`. Ensure it matches the modality you want to
#'   correlate (e.g., "RNA" for gene expression).
#' @param slot Character; which slot to pull: `"data"`, `"counts"`, or
#'   `"scale.data"`. Default `"data"`. Use `"data"` after normalization/log-
#'   transform, `"counts"` for raw counts (less typical for correlations), or
#'   `"scale.data"` if you pre-scaled/centered features.
#' @param cell_type_col Character; name of the `@meta.data` column that holds
#'   **cell-type (layer) labels**. Required unless `group_col` (alias) is used.
#' @param group_col Character; **alias** for `cell_type_col`. If provided and
#'   `cell_type_col` is `NULL`, the function uses this and prints a note.
#' @param feature_ids Character vector of **gene IDs** (row names of the assay
#'   matrix) to include. If `NULL`, the function tries `VariableFeatures()`;
#'   if none are set, it picks the top-variance genes (see `n_var_features`).
#'   Returned matrices are always ordered by this vector, padding missing genes
#'   with zero rows/columns to ensure identical shapes across layers/repeats.
#' @param n_var_features Integer; how many top-variance features to pick when
#'   `VariableFeatures()` is absent. Default `3000`.
#' @param cor_method Correlation method: `"spearman"` (rank-based; default) or
#'   `"pearson"` (linear). Choose `"spearman"` for robustness to outliers and
#'   nonlinearity; `"pearson"` for linear associations on scaled data.
#' @param pval_adjust P-value adjustment method. Accepts `"none"`, `"fdr"`,
#'   `"BH"`, `"bonferroni"`, `"bf"` (synonym for `"bonferroni"`). The method is
#'   applied to the **upper triangle only** and mirrored. Default `"none"`.
#' @param pval_cutoff Numeric; p-value cutoff **after** adjustment (if any).
#'   Edges with p-values above this threshold are zeroed. Default `0.05`.
#' @param corr_threshold Numeric; absolute correlation threshold. Only edges
#'   with `|r| ≥ corr_threshold` survive. Default `0.7`.
#' @param resample Logical; if `TRUE` (default), perform **balanced** sampling
#'   of `samples_per_group` cells from each cell type, repeated `n_repeats` times.
#'   If `FALSE`, use all available cells per cell type exactly once.
#' @param samples_per_group Integer; **target** number of cells per cell type
#'   in each repeat when `resample = TRUE`. If a cell type has fewer cells, the
#'   function draws as many as available (no upsampling).
#' @param n_repeats Integer; number of repeated draws when `resample = TRUE`.
#'   Default `5`.
#' @param min_cells_per_layer Integer; minimum number of cells required in a
#'   cell type for attempting an adjacency. Layers with fewer cells return a
#'   zero matrix of the correct shape. Default `20`.
#' @param group_order Character vector; **alias** for `layer_order`. If used,
#'   the function prints a note and applies it as `layer_order`.
#' @param layer_order Character vector; desired **ordering of layers** (cell
#'   types) in the output. Any unspecified layers are appended after the given
#'   order.
#' @param save_rds Logical; whether to save the return object to **RDS**
#'   on disk. Default `TRUE`. The path is stored in the top-level attribute
#'   `"rds_file"`.
#' @param out_dir Directory where the RDS file is written when `save_rds = TRUE`.
#'   Defaults to `file.path(getwd(), "omicsDNA_sc_results")`. Created if missing.
#' @param file_prefix Basename (no extension) for the saved file. Default
#'   `"sc_adjacency"`. A timestamp and a tag indicating resampling/no-resampling
#'   are appended automatically.
#' @param compress Compression for `saveRDS()`. Accepts logical (`TRUE`/`FALSE`)
#'   or a method string (`"gzip"`, `"bzip2"`, `"xz"`). Default `"xz"` (smallest
#'   size, slower write/read).
#' @param verbose Logical; print progress/QC notes. Default `TRUE`.
#'
#' @return
#' If `resample = TRUE`: a **list of length `n_repeats`**, each element being a
#' **named sub-list** of per-layer adjacency matrices (one matrix per cell type,
#' ordered by `layer_order` if supplied).
#' If `resample = FALSE`: a **single named list** of per-layer adjacency matrices.
#'
#' For every adjacency matrix:
#' - dimension is `length(feature_ids) × length(feature_ids)` and row/column
#'   names are exactly `feature_ids` (padded with zeros for missing genes).
#' - diagonal is zero; nonzero entries are the (signed) correlation estimates.
#' - attribute `"n_cells"` records how many cells contributed for that layer/repeat.
#'
#' The **top-level** return object has attribute `"rds_file"` with the saved path
#' (or `NULL` if `save_rds = FALSE`).
#'
#' @seealso
#' `omicsDNA::buildAdjacency()`, `omicsDNA::edgesFromAdjacency()`,
#' `omicsDNA::consensusEdges()`, `omicsDNA::build_multiNet()`
#'
#' @references
#' Langfelder P, Horvath S (2008) WGCNA: an R package for weighted correlation
#' network analysis. *BMC Bioinformatics* 9, 559.
#' Harrell FE Jr et al. (2024) Hmisc: Harrell Miscellaneous. R package version.
#'
#' @note
#' - For reproducible draws, set a seed (e.g., `set.seed(1)`) **before** calling.
#' - The function prefers `WGCNA` for correlation + p-values and falls back to
#'   `Hmisc` when `WGCNA` is not installed.
#'
#' @export
sc_buildAdjacency <- function(
    seurat_obj,
    assay                 = NULL,
    slot                  = c("data","counts","scale.data"),
    cell_type_col         = NULL,
    group_col             = NULL,      # alias
    feature_ids           = NULL,
    n_var_features        = 3000,
    cor_method            = c("spearman","pearson"),
    pval_adjust           = c("none","fdr","BH","bonferroni","bf"),
    pval_cutoff           = 0.05,
    corr_threshold        = 0.7,
    resample              = TRUE,
    samples_per_group     = 100,
    n_repeats             = 5,
    min_cells_per_layer   = 20,
    group_order           = NULL,      # alias
    layer_order           = NULL,
    save_rds              = TRUE,
    out_dir               = file.path(getwd(), "omicsDNA_results"),
    file_prefix           = "sc_adjacency",
    compress              = "xz",
    verbose               = TRUE
) {

  slot <- match.arg(slot)
  cor_method <- match.arg(cor_method)

  # -- Resolve p-value adjustment (accept synonyms) --
  pval_adjust <- if (length(pval_adjust) > 1L) pval_adjust[1L] else pval_adjust
  method_key  <- tolower(pval_adjust)
  allowed     <- c("none","fdr","bh","bonferroni","bf")
  if (!method_key %in% allowed) {
    stop("Invalid `pval_adjust`. Choose one of: ", paste(allowed, collapse=", "))
  }
  adj_method <- switch(method_key,
                       "none"="none","fdr"="fdr","bh"="BH",
                       "bonferroni"="bonferroni","bf"="bonferroni")
  if (verbose && adj_method != "none") {
    message("Using adjusted p-values: method = ", adj_method)
  }

  # -- Seurat helpers (support SeuratObject or Seurat) --
  .haspkg <- function(pkg) requireNamespace(pkg, quietly = TRUE)
  if (.haspkg("SeuratObject")) {
    GetAssayData   <- SeuratObject::GetAssayData
    DefaultAssay   <- SeuratObject::DefaultAssay
    VariableFeat   <- SeuratObject::VariableFeatures
  } else if (.haspkg("Seurat")) {
    GetAssayData   <- Seurat::GetAssayData
    DefaultAssay   <- Seurat::DefaultAssay
    VariableFeat   <- Seurat::VariableFeatures
  } else {
    stop("Please install 'SeuratObject' or 'Seurat'.")
  }

  # -- Assay and cell-type column resolution (alias messages) --
  if (is.null(assay)) assay <- DefaultAssay(seurat_obj)
  if (!is.null(group_col) && is.null(cell_type_col)) {
    cell_type_col <- group_col
    if (verbose) message("Using `group_col` as alias of `cell_type_col`: ", group_col)
  }
  if (is.null(cell_type_col)) {
    stop("Provide `cell_type_col` (or `group_col` alias) pointing to cell-type labels in @meta.data.")
  }

  if (!is.null(group_order) && is.null(layer_order)) {
    layer_order <- group_order
    if (verbose) message("Using `group_order` as alias of `layer_order`.")
  }

  # -- Pull the matrix (features x cells) --
  Mat <- GetAssayData(seurat_obj, assay = assay, slot = slot)
  if (is.null(rownames(Mat)) || is.null(colnames(Mat)))
    stop("Assay matrix must have rownames (features) and colnames (cells).")
  if (is.data.frame(Mat)) Mat <- as.matrix(Mat)  # <<< FIX: guard against data.frame

  # -- Determine feature set (consistent across layers/repeats) --
  if (is.null(feature_ids)) {
    vf <- tryCatch(VariableFeat(seurat_obj, assay = assay), error = function(e) character(0))
    vf <- vf[!is.na(vf)]
    if (length(vf) == 0L) {
      if (verbose) message("No VariableFeatures; selecting top-", n_var_features, " by variance.")
      .row_vars <- function(m) {
        n <- ncol(m)
        if (inherits(m, "dgCMatrix")) {
          m2 <- m; m2@x <- m2@x^2
          s1 <- Matrix::rowSums(m)
          s2 <- Matrix::rowSums(m2)
        } else {
          s1 <- rowSums(m)
          s2 <- rowSums(m * m)
        }
        (s2 - (s1 * s1) / n) / max(1, n - 1)
      }
      rv <- .row_vars(Mat)
      o  <- order(rv, decreasing = TRUE)
      keep <- head(rownames(Mat)[o], n = min(n_var_features, length(o)))
      feature_ids <- unique(keep)
    } else {
      feature_ids <- unique(as.character(vf))
    }
  } else {
    feature_ids <- unique(as.character(feature_ids))
  }

  # Keep only features that actually exist (pad later)
  present_feats <- intersect(feature_ids, rownames(Mat))
  missing_feats <- setdiff(feature_ids, present_feats)
  if (!length(present_feats)) stop("None of `feature_ids` found in assay rows.")
  Mat <- Mat[present_feats, , drop = FALSE]

  # -- Layers (cell types) --
  md <- seurat_obj@meta.data
  if (!cell_type_col %in% names(md)) {
    stop("`cell_type_col` (", cell_type_col, ") not found in seurat_obj@meta.data.")
  }
  layers <- as.character(md[[cell_type_col]])
  names(layers) <- rownames(md)  # cell barcodes
  uniq_layers <- unique(layers[!is.na(layers)])

  # Apply user-specified order, append any remainder
  if (!is.null(layer_order)) {
    in_both   <- intersect(layer_order, uniq_layers)
    remainder <- setdiff(uniq_layers, in_both)
    uniq_layers <- c(in_both, remainder)
  }

  # -- Correlation + P helper (X = cells x genes) --
  .cor_with_p <- function(X) {
    Cp <- try({
      if (requireNamespace("WGCNA", quietly = TRUE)) {
        method_map <- c(pearson="p", spearman="s")
        WGCNA::corAndPvalue(as.matrix(X), method = method_map[[cor_method]], alternative = "two.sided")
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

  # -- Build one layer adjacency (optional draw) --
  .adj_for_layer <- function(layer_id, draw_n = NULL, rep_index = NA_integer_) {
    cells_all <- names(layers)[layers == layer_id]
    n_all     <- length(cells_all)

    if (n_all < min_cells_per_layer) {
      if (verbose) message("Skipping '", layer_id, "': only ", n_all, " cells (< ", min_cells_per_layer, ").")
      Z <- matrix(0, nrow = length(feature_ids), ncol = length(feature_ids),
                  dimnames = list(feature_ids, feature_ids))
      attr(Z, "n_cells") <- n_all
      return(Z)
    }

    if (!is.null(draw_n)) {
      draw_n <- min(draw_n, n_all)
      set.seed(NULL)
      cells <- sample(cells_all, draw_n)
    } else {
      cells <- cells_all
      draw_n <- length(cells)
    }

    # <<< FIX: ensure the sampled cells actually exist in the assay matrix
    cells <- intersect(cells, colnames(Mat))
    if (length(cells) < 2L) {
      if (verbose) message("Layer '", layer_id, "' has < 2 matching cells in assay after intersect; returning zeros.")
      Z <- matrix(0, nrow = length(feature_ids), ncol = length(feature_ids),
                  dimnames = list(feature_ids, feature_ids))
      attr(Z, "n_cells") <- length(cells)
      return(Z)
    }

    subM <- Mat[, cells, drop = FALSE]   # genes x cells

    # <<< FIX: robust transpose for sparse/base matrices; avoid t.default() on non-matrix
    if (inherits(subM, "Matrix")) {
      X <- Matrix::t(subM)               # cells x genes (sparse-friendly)
    } else {
      X <- t(as.matrix(subM))            # force base matrix then transpose
    }

    # QC using WGCNA if available
    if (requireNamespace("WGCNA", quietly = TRUE)) {
      gsg <- WGCNA::goodSamplesGenes(as.matrix(X), verbose = 0)
      if (!gsg$allOK) {
        if (verbose) {
          if (sum(!gsg$goodGenes) > 0) {
            removed_features <- colnames(X)[!gsg$goodGenes]
            message(
              paste0(
                paste(removed_features, collapse = ", "),
                " removed in iteration ", ifelse(is.na(rep_index), "NA", rep_index),
                " from layer ", layer_id
              )
            )
          }
          if (sum(!gsg$goodSamples) > 0) {
            removed_samples <- rownames(X)[!gsg$goodSamples]
            message("Removing cells: ", paste(removed_samples, collapse = ", "),
                    " | layer = ", layer_id,
                    if (!is.na(rep_index)) paste0(" | repeat = ", rep_index))
          }
        }
        X <- as.matrix(X)[gsg$goodSamples, gsg$goodGenes, drop = FALSE]  # ensure base matrix after QC
      } else {
        X <- as.matrix(X)  # ensure base matrix for downstream stats either way
      }
    } else {
      X <- as.matrix(X)    # ensure base matrix when WGCNA not present
    }

    if (nrow(X) < 3L || ncol(X) < 2L) {
      if (verbose) message("After QC, layer '", layer_id, "' too small; returning zeros.")
      Z <- matrix(0, nrow = length(feature_ids), ncol = length(feature_ids),
                  dimnames = list(feature_ids, feature_ids))
      attr(Z, "n_cells") <- nrow(X)
      return(Z)
    }

    cp <- .cor_with_p(X)
    cor_mat <- cp$cor; p_mat <- cp$p
    diag(cor_mat) <- 0; diag(p_mat) <- 1

    # Adjust p-values on unique tests (upper triangle), mirror back
    p_use <- p_mat
    if (adj_method != "none") {
      ut <- upper.tri(p_use, diag = FALSE)
      p_use[ut] <- p.adjust(p_use[ut], method = adj_method)
      p_use[lower.tri(p_use, diag = FALSE)] <- t(p_use)[lower.tri(p_use, diag = FALSE)]
      diag(p_use) <- 1
    }

    keep <- (abs(cor_mat) >= corr_threshold) & (p_use <= pval_cutoff)
    keep[is.na(keep)] <- FALSE

    adj <- cor_mat
    adj[!keep] <- 0

    # Expand back to the full feature_ids order (pad zeros for missing)
    present <- colnames(X)
    full_adj <- matrix(0, nrow = length(feature_ids), ncol = length(feature_ids),
                       dimnames = list(feature_ids, feature_ids))
    common <- intersect(feature_ids, present)
    if (length(common)) {
      full_adj[common, common] <- adj[common, common, drop = FALSE]
    }

    if (verbose) message("Built adjacency for '", layer_id, "' (genes=", length(feature_ids),
                         ", cells=", draw_n, ").")
    attr(full_adj, "n_cells") <- draw_n
    full_adj
  }

  # -- Driver: resampling or single pass --
  result <- NULL
  if (resample) {
    if (verbose) message("Resampling: ", n_repeats, " repeats; ", samples_per_group, " cells/type.")
    result <- vector("list", n_repeats)
    for (i in seq_len(n_repeats)) {
      mats <- lapply(uniq_layers, function(L) .adj_for_layer(L, draw_n = samples_per_group, rep_index = i))
      names(mats) <- uniq_layers
      result[[i]] <- mats
      if (verbose) message("sampling ", i, " completed")
    }
  } else {
    if (verbose) message("Computing adjacency once per layer (no resampling).")
    mats <- lapply(uniq_layers, function(L) .adj_for_layer(L, draw_n = NULL, rep_index = NA_integer_))
    names(mats) <- uniq_layers
    result <- mats
  }

  # -- Save to RDS (optional) --
  rds_path <- NULL
  if (isTRUE(save_rds)) {
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    stamp <- format(Sys.time(), "%Y%m%d-%H%M%S")
    tag   <- if (resample) "bootstrap" else "no_bootstrap"
    prefix <- if (is.null(file_prefix) || !nzchar(file_prefix)) "sc_adjacency" else file_prefix
    rds_path <- file.path(out_dir, sprintf("%s_%s_%s.rds", prefix, tag, stamp))
    comp <- if (isTRUE(compress) || isFALSE(compress)) compress else as.character(compress)
    saveRDS(result, rds_path, compress = comp)
    if (verbose) message("💾 Saved adjacency to: ", rds_path)
  }

  attr(result, "rds_file") <- rds_path
  return(result)
}

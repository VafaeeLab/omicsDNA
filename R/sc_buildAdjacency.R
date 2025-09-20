
# --------------------------------------------------------------------------------
#  1b - sc_buildAdjacency(): Seurat → per–cell type gene–gene adjacency
# --------------------------------------------------------------------------------

#' Build per–cell type gene–gene adjacency matrices from a Seurat object (no aggregation)
#'
#' @title Single‑cell coexpression adjacency per cell type (layers) from Seurat
#'
#' @description
#' Given a \pkg{Seurat} object with one or more assays (e.g., \code{"RNA"}, \code{"SCT"},
#' \code{"integrated"}), this function extracts an assay/slot matrix (default: \code{slot="data"}),
#' pulls the cell metadata, and for each **cell type** (chosen via a metadata column) computes
#' a gene–gene correlation matrix **across individual cells** of that type (no aggregation).
#' Edges are kept when they satisfy both an absolute correlation threshold and a p‑value cutoff.
#' The result is a list of adjacency matrices—one per cell type—sharing an identical gene order,
#' ready to be used as layers in a multilayer network.
#'
#' @details
#' \strong{Assay and slot.}
#' You can point to any assay present in the object (e.g., \code{"RNA"}, \code{"SCT"},
#' \code{"integrated"}) and choose the slot to read (\code{"data"} is recommended).
#' For typical Seurat workflows:
#' \itemize{
#'   \item \code{assay="RNA", slot="data"} → log‑normalized expression (sparse).
#'   \item \code{assay="SCT", slot="data"} → SCT residuals (already normalized).
#'   \item \code{assay="integrated", slot="data"} → corrected expression used for integration.
#' }
#'
#' \strong{Cell types (layers).}
#' The function reads \code{object@meta.data}. Set \code{cell_type_col} to the column
#' that contains the cell‑type labels (these labels become the layer names). If omitted,
#' the function tries common candidates (e.g., \code{"cell_type"}, \code{"celltype"},
#' \code{"CellType"}, \code{"seurat_clusters"}), and errors if none is found.
#'
#' \strong{Gene universe.}
#' Use \code{feature_ids} to fix the genes used across all layers. If \code{feature_ids} is
#' \code{NULL}, the function (by default) uses \code{\link[Seurat]{VariableFeatures}} from
#' the selected assay if available; otherwise it picks the top \code{n_features} genes by
#' row variance (computed efficiently on sparse matrices). All adjacency matrices are then
#' padded to this common gene universe and share the same row/column order.
#'
#' \strong{Correlation & statistics.}
#' Correlations are computed with \code{cor_method = "spearman"} (default), or
#' \code{"pearson"}, or robust \code{"bicor"} (requires \pkg{WGCNA}). P‑values come from
#' \pkg{WGCNA} (\code{corAndPvalue}) when available; otherwise the function falls back to
#' \pkg{Hmisc} (\code{rcorr}) for Pearson/Spearman. Adjusted p‑values are supported via
#' \code{pval_adjust = "fdr"|"BH"|"bonferroni"|"none"}.
#'
#' \strong{Scale considerations.}
#' Computing correlations across cells is feasible provided you keep the gene set modest
#' (e.g., 1–3k genes). For very large datasets, consider limiting cells per type using
#' \code{max_cells_per_type}, or pass a curated \code{feature_ids}. Note that with large
#' numbers of cells, **p‑values can be extremely small**; if you want to keep edges based
#' on correlation only, set \code{pval_cutoff = 1}.
#'
#' @param seurat_obj A \pkg{Seurat} object.
#' @param assay Character. Assay to use. If \code{NULL}, the function tries, in order,
#'   \code{"integrated"}, \code{"SCT"}, \code{"RNA"}, then falls back to the first assay present.
#' @param slot Character; which slot to read from the assay. One of \code{"data"},
#'   \code{"scale.data"}, \code{"counts"}. Default \code{"data"}.
#' @param cell_type_col Character; name of the metadata column with cell‑type labels.
#'   If \code{NULL}, the function searches common candidates and errors if none are present.
#' @param feature_ids Optional character vector of genes to include. If \code{NULL},
#'   the function uses \code{VariableFeatures(seurat_obj[[assay]])} when available,
#'   otherwise selects the top \code{n_features} by row variance.
#' @param use_variable_features Logical; if \code{TRUE} and \code{feature_ids} is \code{NULL},
#'   attempt to use \code{VariableFeatures}. Default \code{TRUE}.
#' @param n_features Integer; when \code{feature_ids} is \code{NULL} and no variable features
#'   are stored, pick this many top‑variance genes. Default \code{2000}.
#' @param min_cells_per_type Integer; minimum number of cells required to build an adjacency
#'   for a cell type. Types with fewer cells return a zero matrix. Default \code{20}.
#' @param max_cells_per_type Optional integer; if provided and a cell type has more cells than this,
#'   randomly subsample to this many cells to control time/memory. Default \code{NULL} (use all).
#' @param cor_method Correlation method: \code{"spearman"} (default), \code{"pearson"},
#'   or \code{"bicor"} (requires \pkg{WGCNA}).
#' @param pval_adjust P‑value adjustment: \code{"fdr"}, \code{"BH"}, \code{"bonferroni"},
#'   or \code{"none"}. Default \code{"fdr"}.
#' @param pval_cutoff Numeric; (adjusted) p‑value cutoff for edges. Default \code{0.05}.
#' @param corr_threshold Numeric; absolute correlation threshold for edges. Default \code{0.25}.
#' @param layer_order Optional character vector to control the order of cell types in the output.
#' @param verbose Logical; print progress messages. Default \code{TRUE}.
#'
#' @return
#' A named list of adjacency matrices \code{adj[[cell_type]] = matrix} (genes × genes), all with
#' identical row/column order. Attributes:
#' \itemize{
#'   \item \code{attr(x, "gene_order")} — the gene order used across matrices.
#'   \item \code{attr(x, "cell_types")} — the cell types included (order reflects output).
#'   \item \code{attr(x, "assay")} and \code{attr(x, "slot")} — provenance.
#' }
#'
#' @section Tips:
#' \itemize{
#'   \item If you see “too few cells” messages, lower \code{min_cells_per_type} or merge rare types.
#'   \item To keep only correlation‑based edges (ignore p‑values), set \code{pval_cutoff = 1}.
#'   \item For very large cell counts, set \code{max_cells_per_type} (e.g., 5000–10000).
#' }
#'
#' @seealso \code{\link[Seurat]{GetAssayData}}, \code{\link[Seurat]{VariableFeatures}}
#'
#' @importFrom Seurat GetAssayData VariableFeatures
#' @importFrom stats p.adjust
#' @examples
#' \dontrun{
#' library(Seurat)
#' seu <- readRDS("my_seurat.rds")
#'
#' adj <- sc_buildAdjacency(
#'   seurat_obj         = seu,
#'   assay              = "SCT",        # or "RNA", "integrated", etc.
#'   slot               = "data",
#'   cell_type_col      = "celltype",   # choose your metadata column
#'   feature_ids        = NULL,         # use VariableFeatures if present; else top-variance
#'   use_variable_features = TRUE,
#'   n_features         = 3000,
#'   min_cells_per_type = 30,
#'   max_cells_per_type = 8000,         # subsample if needed; NULL = use all
#'   cor_method         = "spearman",
#'   pval_adjust        = "fdr",
#'   pval_cutoff        = 0.05,
#'   corr_threshold     = 0.25,
#'   layer_order        = c("T cell","B cell","NK","Myeloid","Endothelial"),
#'   verbose            = TRUE
#' )
#'
#' names(adj)                 # cell types (layers)
#' dim(adj[[1]])              # genes × genes (same order across layers)
#' attr(adj, "gene_order")[1:10]
#' }
#' @export
sc_buildAdjacency <- function(
    seurat_obj,
    assay               = NULL,
    slot                = c("data","scale.data","counts"),
    cell_type_col       = NULL,
    feature_ids         = NULL,
    use_variable_features = TRUE,
    n_features          = 2000,
    min_cells_per_type  = 20,
    max_cells_per_type  = NULL,
    cor_method          = c("spearman","pearson","bicor"),
    pval_adjust         = c("fdr","BH","bonferroni","none"),
    pval_cutoff         = 0.05,
    corr_threshold      = 0.25,
    layer_order         = NULL,
    verbose             = TRUE
) {
  # ---- deps ----
  if (!requireNamespace("Seurat", quietly = TRUE))
    stop("Please install Seurat.")
  if (!requireNamespace("Matrix", quietly = TRUE))
    stop("Please install Matrix.")

  slot       <- match.arg(slot)
  cor_method <- match.arg(cor_method)
  pval_adjust<- match.arg(pval_adjust)

  # ---- pick assay (priority: integrated → SCT → RNA → first) ----
  assays_avail <- names(seurat_obj@assays)
  if (is.null(assay)) {
    assay <- if ("integrated" %in% assays_avail) "integrated" else
      if ("SCT" %in% assays_avail) "SCT" else
        if ("RNA" %in% assays_avail) "RNA" else assays_avail[1]
    if (verbose) message("Assay auto-selected: ", assay)
  } else {
    if (!assay %in% assays_avail)
      stop("Assay '", assay, "' not found. Available: ", paste(assays_avail, collapse = ", "))
  }

  # ---- pull matrix and metadata ----
  mat  <- Seurat::GetAssayData(seurat_obj, assay = assay, slot = slot)
  meta <- seurat_obj@meta.data
  if (is.null(colnames(mat)) || is.null(rownames(meta)))
    stop("Assay matrix needs colnames (cells) and meta.data needs rownames (cells).")

  # align metadata to matrix columns
  common_cells <- intersect(colnames(mat), rownames(meta))
  if (length(common_cells) < 3L)
    stop("Fewer than 3 shared cells between assay and metadata.")
  mat  <- mat[, common_cells, drop = FALSE]
  meta <- meta[common_cells, , drop = FALSE]

  # ---- choose cell-type column ----
  if (is.null(cell_type_col)) {
    candidates <- c("cell_type","celltype","CellType","celltype_l1","celltype.l1",
                    "major_cell_type","annotation","annot","seurat_clusters","cellType")
    cell_type_col <- candidates[candidates %in% colnames(meta)][1]
    if (is.na(cell_type_col))
      stop("Please provide `cell_type_col` (no common candidates found in meta.data).")
    if (verbose) message("cell_type_col auto-selected: ", cell_type_col)
  } else if (!cell_type_col %in% colnames(meta)) {
    stop("`cell_type_col` not found in meta.data: ", cell_type_col)
  }

  cell_types <- unique(as.character(meta[[cell_type_col]]))
  cell_types <- cell_types[!is.na(cell_types)]
  if (length(cell_types) == 0L) stop("No cell types found in column: ", cell_type_col)

  # optional ordering
  if (!is.null(layer_order)) {
    keep <- intersect(layer_order, cell_types)
    cell_types <- c(keep, setdiff(cell_types, keep))
  }

  # ---- decide gene universe ----
  all_genes <- rownames(mat)
  if (!is.null(feature_ids)) {
    genes <- intersect(all_genes, unique(as.character(feature_ids)))
    if (length(genes) < 2L)
      stop("After intersecting with the assay, <2 genes remain in `feature_ids`.")
    if (verbose) message("Using ", length(genes), " genes from `feature_ids`.")
  } else if (isTRUE(use_variable_features)) {
    # try stored variable features for this assay
    vf <- try(Seurat::VariableFeatures(seurat_obj[[assay]]), silent = TRUE)
    vf <- if (inherits(vf, "try-error")) character(0) else vf
    vf <- intersect(all_genes, vf)
    if (length(vf) > 0L) {
      genes <- if (!is.null(n_features) && is.finite(n_features)) head(vf, n_features) else vf
      if (verbose) message("Using ", length(genes), " VariableFeatures from assay '", assay, "'.")
    } else {
      # fallback to top-variance genes across all cells (sparse-friendly)
      if (verbose) message("No VariableFeatures stored; selecting top-", n_features, " by variance.")
      genes <- .top_variance_genes(mat, n = n_features)
    }
  } else {
    if (verbose) message("Using top-", n_features, " genes by variance (no VariableFeatures requested).")
    genes <- .top_variance_genes(mat, n = n_features)
  }

  # ---- helpers ----
  .rc_fast_rowvars <- function(smat) {
    # row variance for genes × cells sparse or dense matrix
    n <- ncol(smat)
    # sums and sums of squares (sparse-friendly)
    rs  <- Matrix::rowSums(smat)
    rss <- Matrix::rowSums(smat^2)
    mu  <- rs / pmax(n, 1L)
    v   <- pmax(rss / pmax(n, 1L) - mu^2, 0)
    as.numeric(v)
  }

  .cor_with_p <- function(X, method) {
    # X: cells × genes (observations × variables)
    # prefer WGCNA for p-values (supports bicor)
    Cp <- try({
      if (requireNamespace("WGCNA", quietly = TRUE)) {
        method_map <- c(pearson = "p", spearman = "s", bicor = "bicor")
        WGCNA::corAndPvalue(X, method = method_map[[method]], alternative = "two.sided")
      } else stop("noWGCNA")
    }, silent = TRUE)
    if (!inherits(Cp, "try-error")) return(list(cor = Cp$cor, p = Cp$p))

    # fallback to Hmisc::rcorr for pearson/spearman
    if (!requireNamespace("Hmisc", quietly = TRUE))
      stop("Need WGCNA or Hmisc to compute correlation p-values.")
    if (method == "bicor") stop("bicor requires WGCNA; install WGCNA or use 'spearman'/'pearson'.")
    rc <- Hmisc::rcorr(as.matrix(X), type = method)
    list(cor = rc$r, p = rc$P)
  }

  .pad_adjust <- function(pmat, method) {
    if (method == "none") return(pmat)
    meth <- if (tolower(method) == "fdr") "BH" else method
    ut <- upper.tri(pmat, diag = FALSE)
    pvec <- pmat[ut]
    pmat[ut] <- p.adjust(pvec, method = meth)
    pmat[lower.tri(pmat, diag = FALSE)] <- t(pmat)[lower.tri(pmat, diag = FALSE)]
    diag(pmat) <- 1
    pmat
  }

  # expose top-variance selector
  .top_variance_genes <- function(smat, n = 2000) {
    v <- .rc_fast_rowvars(smat)
    ord <- order(v, decreasing = TRUE)
    head(rownames(smat)[ord], n = min(n, nrow(smat)))
  }

  # ---- per-type adjacency ----
  gene_universe <- genes
  adj_list <- setNames(vector("list", length(cell_types)), cell_types)

  for (ct in cell_types) {
    cells_ct <- rownames(meta)[meta[[cell_type_col]] == ct]
    cells_ct <- intersect(cells_ct, colnames(mat))
    n_ct <- length(cells_ct)

    if (n_ct < max(3L, min_cells_per_type)) {
      if (verbose) message("Skipping '", ct, "': only ", n_ct, " cells (< ", max(3L, min_cells_per_type), ").")
      adj_list[[ct]] <- matrix(0, nrow = length(gene_universe), ncol = length(gene_universe),
                               dimnames = list(gene_universe, gene_universe))
      next
    }

    if (!is.null(max_cells_per_type) && is.finite(max_cells_per_type) && n_ct > max_cells_per_type) {
      set.seed(1)
      cells_ct <- sample(cells_ct, max_cells_per_type)
      n_ct <- length(cells_ct)
      if (verbose) message("Subsampled '", ct, "' to ", n_ct, " cells (max_cells_per_type).")
    }

    # subset to gene_universe × cells_ct
    M <- mat[gene_universe, cells_ct, drop = FALSE]   # genes × cells
    # drop zero-variance genes within this cell type (avoid NA correlations)
    v_ct <- .rc_fast_rowvars(M)
    gkeep <- gene_universe[v_ct > 0]
    if (length(gkeep) < 2L) {
      if (verbose) message("Cell type '", ct, "': <2 informative genes; returning zeros.")
      adj_list[[ct]] <- matrix(0, nrow = length(gene_universe), ncol = length(gene_universe),
                               dimnames = list(gene_universe, gene_universe))
      next
    }
    M2 <- M[gkeep, , drop = FALSE]

    # compute correlation & p across cells
    X <- t(as.matrix(M2))  # cells × genes
    cp <- .cor_with_p(X, method = cor_method)
    R <- cp$cor; P <- cp$p
    # ensure row/colnames
    rownames(R) <- colnames(R) <- gkeep
    rownames(P) <- colnames(P) <- gkeep
    diag(R) <- 0; diag(P) <- 1

    # adjust p-values if requested
    P2 <- .pad_adjust(P, pval_adjust)

    # threshold to adjacency
    keep <- (abs(R) >= corr_threshold) & (P2 <= pval_cutoff)
    keep[is.na(keep)] <- FALSE
    A <- R; A[!keep] <- 0

    # pad to full gene_universe order
    fullA <- matrix(0, nrow = length(gene_universe), ncol = length(gene_universe),
                    dimnames = list(gene_universe, gene_universe))
    fullA[gkeep, gkeep] <- A[gkeep, gkeep, drop = FALSE]
    adj_list[[ct]] <- fullA
    if (verbose) message("Built adjacency for '", ct, "' (genes=", nrow(fullA), ", cells=", n_ct, ").")
  }

  attr(adj_list, "gene_order") <- gene_universe
  attr(adj_list, "cell_types") <- cell_types
  attr(adj_list, "assay")      <- assay
  attr(adj_list, "slot")       <- slot

  # ---- NEW: save the output list to a timestamped .rds (no behavior change) ----
  # Uses getOption("mlnet.results_dir","omicsDNA_results") as the base folder.
  {
    .ensure_dir <- function(d) if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
    .san        <- function(x) gsub("[^A-Za-z0-9._-]+", "_", as.character(x))
    results_dir <- getOption("mlnet.results_dir", "omicsDNA_results")
    .ensure_dir(results_dir)
    stamp <- format(Sys.time(), "%Y-%m-%d_%H%M%S")
    aslot <- paste0("_", .san(assay), "_", .san(slot))
    nl    <- paste0("_Nlayers-", length(adj_list))
    fname <- paste0("sc_adjacency", aslot, nl, "_", stamp, ".rds")
    path  <- file.path(results_dir, fname)
    saveRDS(adj_list, file = path, compress = "xz")
    attr(adj_list, "rds_file") <- normalizePath(path, winslash = "/", mustWork = FALSE)
    if (isTRUE(verbose)) message("💾 Saved adjacency list: ", attr(adj_list, "rds_file"))
  }
  # ------------------------------------------------------------------------------

  adj_list
}

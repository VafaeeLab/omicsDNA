
# -------------------------------------------------------------------
# 2 - Convert adjacency matrix (or nested list of matrices) to edge list(s)
# -------------------------------------------------------------------

#' Turn adjacency matrices (single or nested) into tidy edge lists
#'
#' Convert one adjacency matrix, a list of matrices, or a nested list of matrices
#' into data frame(s) of edges with columns `from`, `to`, and `weight`. Works with
#' dense base matrices and sparse matrices from the **Matrix** package. Supports
#' directed/undirected output, filtering, and optional flattening of nested inputs.
#'
#' @details
#' **Directed vs. undirected**
#' - `directed = FALSE` (default): return one edge per unordered pair (upper triangle),
#'   excluding self‑loops.
#' - `directed = TRUE`: return all off‑diagonal pairs.
#'
#' **Filtering**
#' - Use `drop_zeros`, `min_abs`, and `na_rm` to remove edges. The sign of `weight`
#'   is preserved.
#'
#' **Sparse input**
#' - Sparse matrices (`dgCMatrix`, `dsCMatrix`, …) are handled efficiently by
#'   enumerating only stored (non‑zero) entries. Note that explicit zeros are not
#'   stored in sparse formats; therefore `drop_zeros = FALSE` will not add zeros.
#' - Matrices are expected to be square (enforced for dense matrices and assumed
#'   for sparse matrices).
#'
#' **Lists and nesting**
#' - If `x` is a list (possibly nested), the output can either preserve the list
#'   shape (default) or be flattened into one tidy data frame. For flattened output,
#'   supply `id_cols` to name the list levels (e.g., `c("repeat","group")`). If
#'   `id_cols` is omitted, generic names `level1..levelK` are used. Ragged nesting
#'   is supported; missing levels are filled with `NA`.
#'
#' **Persistence**
#' - When `save_to_rds = TRUE`, the result is saved under `results_dir`. If `x` is
#'   a matrix the file name is `edges_df_<timestamp>.rds`; if a flattened list,
#'   `edges_flat_<timestamp>.rds`; otherwise `edges_list_<timestamp>.rds`. If you
#'   provide `rds_file`, a relative path is placed under `results_dir`, while an
#'   absolute path is respected.
#'
#' @param x Numeric square matrix (adjacency), a sparse Matrix, or a (possibly nested)
#'   list of such matrices.
#' @param directed Logical; if `FALSE` use upper triangle (no self‑loops); if `TRUE`
#'   include all off‑diagonal pairs. Default `FALSE`.
#' @param drop_zeros Logical; drop edges with weight exactly 0 (dense matrices).
#'   Default `TRUE`. For sparse matrices, zeros are not stored and thus never appear.
#' @param min_abs Optional numeric; keep edges with `abs(weight) >= min_abs`. Default `NULL`.
#' @param na_rm Logical; drop edges with `NA` weights. Default `TRUE`.
#' @param flatten Logical; when `x` is a list, return a single data frame with
#'   identifier columns instead of preserving the list structure. Default `FALSE`.
#' @param id_cols Character vector naming the list levels when `flatten = TRUE`.
#'   Length must equal the maximum nesting depth of `x`. If `NULL`, defaults to
#'   `level1..levelK`. Default `NULL`.
#' @param save_to_rds Logical; save the returned edge list(s) as an RDS file. Default `TRUE`.
#' @param rds_file File name to use when saving. If `NULL`, an informative name is
#'   created automatically. Relative names are placed under `results_dir`; absolute
#'   paths are respected. Default `NULL`.
#' @param results_dir Output directory for saved files. Default
#'   `getOption("mlnet.results_dir", "omicsDNA_results")`.
#' @param verbose Logical; print progress/messages. Default `TRUE`.
#'
#' @return
#' - If `x` is a matrix: a data frame with columns `from`, `to`, `weight`.
#' - If `x` is a list: a list (same shape) of data frames, or a single data frame
#'   when `flatten = TRUE` (with `id_cols` prepended).
#'
#' @examples
#' # Dense matrix
#' A <- matrix(c(0, 0.8, 0.8, 0), 2, 2, dimnames = list(c("a","b"), c("a","b")))
#' edgesFromAdjacency(A)
#'
#' # Per‑group list, flattened with a label
#' L <- list(G1 = A, G2 = A * 0.5)
#' edgesFromAdjacency(L, flatten = TRUE, id_cols = "group")
#'
#' # Nested list: repeat -> group (ragged OK)
#' NL <- list(`1` = list(E1 = A, E2 = A), `2` = list(E1 = A))
#' edgesFromAdjacency(NL, flatten = TRUE, id_cols = c("repeat","group"))
#'
#' @seealso \code{\link{buildAdjacency}} to build group‑wise adjacency matrices.
#' @export
edgesFromAdjacency <- function(
    x,
    directed    = FALSE,
    drop_zeros  = TRUE,
    min_abs     = NULL,
    na_rm       = TRUE,
    flatten     = FALSE,
    id_cols     = NULL,
    save_to_rds = TRUE,                                      # <- default save on
    rds_file    = NULL,
    results_dir = getOption("mlnet.results_dir","omicsDNA_results"),  # <- new default
    verbose     = TRUE
) {
  .ensure_dir <- function(d) if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  .is_abs     <- function(p) grepl("^(/|[A-Za-z]:[\\/])", p)

  # ---- dispatch on input type ----
  if (is.matrix(x) || inherits(x, "Matrix")) {
    res <- .adj_matrix_to_edges(
      x,
      directed   = directed,
      drop_zeros = drop_zeros,
      min_abs    = min_abs,
      na_rm      = na_rm
    )
    if (isTRUE(save_to_rds)) {
      .ensure_dir(results_dir)
      if (is.null(rds_file)) {
        stamp   <- format(Sys.time(), "%Y-%m-%d_%H%M%S")
        rds_file <- file.path(results_dir, sprintf("edges_df_%s.rds", stamp))
      } else if (!.is_abs(rds_file)) {
        rds_file <- file.path(results_dir, rds_file)
      }
      saveRDS(res, rds_file)
      if (verbose) message("Saved edge list to: ", normalizePath(rds_file, FALSE))
    }
    return(res)
  }

  if (is.list(x)) {
    depth <- .list_depth(x)
    if (isTRUE(flatten)) {
      # default id_cols: level1..levelK
      if (is.null(id_cols)) id_cols <- paste0("level", seq_len(depth))
      if (length(id_cols) != depth) {
        stop("`id_cols` length (", length(id_cols), ") must match list nesting depth (", depth, ").")
      }
      rows <- .walk_list_collect(
        x,
        id_cols = id_cols,
        path    = character(),
        f_leaf  = function(mat, path) {
          df <- .adj_matrix_to_edges(
            mat,
            directed   = directed,
            drop_zeros = drop_zeros,
            min_abs    = min_abs,
            na_rm      = na_rm
          )
          # Fill missing levels with NA to support ragged lists
          add <- as.list(rep(NA_character_, length(id_cols)))
          add[seq_along(path)] <- as.character(path)
          names(add) <- id_cols
          cbind(as.data.frame(add, optional = TRUE, stringsAsFactors = FALSE), df)
        }
      )
      res <- if (!length(rows)) {
        data.frame(from = character(), to = character(), weight = numeric(), stringsAsFactors = FALSE)
      } else {
        do.call(rbind, rows)
      }
      if (isTRUE(save_to_rds)) {
        .ensure_dir(results_dir)
        if (is.null(rds_file)) {
          stamp   <- format(Sys.time(), "%Y-%m-%d_%H%M%S")
          rds_file <- file.path(results_dir, sprintf("edges_flat_%s.rds", stamp))
        } else if (!.is_abs(rds_file)) {
          rds_file <- file.path(results_dir, rds_file)
        }
        saveRDS(res, rds_file)
        if (verbose) message("Saved flattened edge list to: ", normalizePath(rds_file, FALSE))
      }
      return(res)
    } else {
      # preserve original list structure
      res <- .walk_list_shape(x, f_leaf = function(mat, path) {
        .adj_matrix_to_edges(
          mat,
          directed   = directed,
          drop_zeros = drop_zeros,
          min_abs    = min_abs,
          na_rm      = na_rm
        )
      })
      if (isTRUE(save_to_rds)) {
        .ensure_dir(results_dir)
        if (is.null(rds_file)) {
          stamp   <- format(Sys.time(), "%Y-%m-%d_%H%M%S")
          rds_file <- file.path(results_dir, sprintf("edges_list_%s.rds", stamp))
        } else if (!.is_abs(rds_file)) {
          rds_file <- file.path(results_dir, rds_file)
        }
        saveRDS(res, rds_file)
        if (verbose) message("Saved edge-list (list) to: ", normalizePath(rds_file, FALSE))
      }
      return(res)
    }
  }

  stop("`x` must be a matrix/sparse Matrix or a (nested) list of such matrices.")
}

# ---- internal helpers (unchanged except for sparse support) -------------------

.empty_edges_df <- function() {
  data.frame(from = character(), to = character(), weight = numeric(), stringsAsFactors = FALSE)
}

.adj_matrix_to_edges <- function(adjMatrix,
                                 directed   = FALSE,
                                 drop_zeros = TRUE,
                                 min_abs    = NULL,
                                 na_rm      = TRUE) {
  # Sparse path
  if (inherits(adjMatrix, "Matrix")) {
    if (!methods::is(adjMatrix, "dsCMatrix") && !methods::is(adjMatrix, "dgCMatrix")) {
      adjMatrix <- as(adjMatrix, "dgCMatrix")
    }
    s <- Matrix::summary(adjMatrix) # i,j,x (1-based)
    if (!nrow(s)) return(.empty_edges_df())

    sel <- s$i != s$j
    if (!directed) sel <- sel & (s$i < s$j)
    if (na_rm) sel <- sel & !is.na(s$x)
    if (isTRUE(drop_zeros)) sel <- sel & (s$x != 0)
    if (!is.null(min_abs))  sel <- sel & (abs(s$x) >= min_abs)
    s <- s[sel, , drop = FALSE]
    if (!nrow(s)) return(.empty_edges_df())

    rn <- rownames(adjMatrix); if (is.null(rn)) rn <- paste0("V", seq_len(nrow(adjMatrix)))
    cn <- colnames(adjMatrix); if (is.null(cn)) cn <- paste0("V", seq_len(ncol(adjMatrix)))
    return(data.frame(from = rn[s$i], to = cn[s$j], weight = as.numeric(s$x), stringsAsFactors = FALSE))
  }

  # Dense path
  adj <- as.matrix(adjMatrix)
  if (!is.numeric(adj)) {
    df <- as.data.frame(adj, check.names = FALSE, stringsAsFactors = FALSE)
    for (j in seq_len(ncol(df))) df[[j]] <- suppressWarnings(as.numeric(df[[j]]))
    adj <- as.matrix(df)
  }
  if (nrow(adj) != ncol(adj)) stop("Adjacency matrix must be square.")

  if (is.null(rownames(adj))) rownames(adj) <- paste0("V", seq_len(nrow(adj)))
  if (is.null(colnames(adj))) colnames(adj) <- paste0("V", seq_len(ncol(adj)))

  sel <- if (!directed) upper.tri(adj, diag = FALSE) else row(adj) != col(adj)
  if (!any(sel)) return(.empty_edges_df())

  rr <- row(adj)[sel]; cc <- col(adj)[sel]; w <- adj[sel]

  keep <- rep(TRUE, length(w))
  if (na_rm)              keep <- keep & !is.na(w)
  if (isTRUE(drop_zeros)) keep <- keep & (w != 0)
  if (!is.null(min_abs))  keep <- keep & (abs(w) >= min_abs)
  if (!any(keep)) return(.empty_edges_df())

  data.frame(
    from   = rownames(adj)[rr][keep],
    to     = colnames(adj)[cc][keep],
    weight = as.numeric(w[keep]),
    stringsAsFactors = FALSE
  )
}

.list_depth <- function(x) {
  if (!is.list(x) || !length(x)) return(0L)
  1L + max(vapply(x, .list_depth, integer(1L)))
}

.walk_list_shape <- function(x, f_leaf, path = character()) {
  if (is.null(x)) return(NULL)
  if (is.matrix(x) || inherits(x, "Matrix") || is.data.frame(x)) return(f_leaf(x, path))
  stopifnot(is.list(x))
  out <- vector("list", length(x)); names(out) <- names(x)
  for (i in seq_along(x)) {
    nm <- names(x)[i]
    out[[i]] <- .walk_list_shape(
      x[[i]],
      f_leaf = f_leaf,
      path   = c(path, if (is.null(nm) || nm == "") as.character(i) else nm)
    )
  }
  out
}

.walk_list_collect <- function(x, id_cols, path, f_leaf) {
  res <- list()
  .recurse <- function(node, path) {
    if (is.null(node)) return(invisible(NULL))
    if (is.matrix(node) || inherits(node, "Matrix") || is.data.frame(node)) {
      res[[length(res) + 1L]] <<- f_leaf(node, path)
    } else if (is.list(node) && length(node)) {
      for (i in seq_along(node)) {
        nm <- names(node)[i]
        .recurse(node[[i]], c(path, if (is.null(nm) || nm == "") as.character(i) else nm))
      }
    }
  }
  .recurse(x, path)
  res
}

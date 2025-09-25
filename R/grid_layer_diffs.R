
# -------------------------------------------------------
# 24 — Grid of consecutive layer differences
# -------------------------------------------------------

#' Grid of consecutive layer differences (E1→E2, E2→E3, …) with timestamped outputs
#'
#' @description
#' Builds a grid of **consecutive** layer difference plots from a multilayer network
#' (e.g., E1→E2, E2→M1, …), saves **one per‑pair figure** in a user‑selected format
#' (\code{"png"}, \code{"pdf"}, \code{"jpg"}) **with a timestamped name**, and writes a
#' **single combined CSV** stacking edges for all pairs (including per‑layer weights).
#' A **grid figure** summarizing all pairs is also saved with a timestamped filename.
#'
#' @details
#' The grid is assembled from PNG panels; if \code{pair_format != "png"}, an additional
#' temporary PNG per pair is rendered for layout, while the per‑pair image is saved in
#' the requested format.
#'
#' @param net A multilayer network (usable with \code{multinet::layers_ml()}).
#' @param layer_order Optional explicit layer order (e.g., \code{names(cons_list)}).
#'   Defaults to \code{multinet::layers_ml(net)}.
#' @param ncol Integer; number of columns in the grid layout.
#' @param results_dir Output directory for all artifacts. Default:
#'   \code{getOption("mlnet.results_dir","omicsDNA_results")}.
#' @param grid_file Optional base filename for the grid figure (e.g., \code{"diffs_all_pairs.png"}).
#'   If relative, it is created under \code{results_dir}. A timestamp is appended
#'   before the extension. If \code{NULL}, a timestamped name is auto‑generated.
#' @param grid_format Format for the grid figure: \code{"png"}, \code{"pdf"}, or \code{"jpg"}.
#'   Default \code{"png"}.
#' @param width,height,dpi Device size and resolution for the grid figure (\code{dpi} ignored for PDF).
#' @param pair_format Format for **per‑pair** figures: \code{"png"}, \code{"pdf"}, or \code{"jpg"}.
#'   Default \code{"png"}.
#' @param pair_file_prefix Basename prefix for per‑pair files (e.g., \code{"diff_pair"}).
#'   Files are named \code{"<prefix>_<L1>_vs_<L2>_<timestamp>.<ext>"} under \code{results_dir}.
#' @param panel_width,panel_height,dpi_panel Geometry used when rendering panel PNGs for the grid
#'   (and per‑pair PNGs). \code{dpi_panel} is ignored for PDF.
#' @param combined_csv_file Optional base filename for the combined CSV (e.g.,
#'   \code{"diff_edges_all_pairs.csv"}). If relative, it is created under \code{results_dir}.
#'   A timestamp is appended before \code{.csv}. If \code{NULL}, a default name is used.
#' @param weight_attr Candidate edge attribute names for weights (forwarded to \code{plot_layer_diff()}).
#' @param weight_aggregate How to combine duplicates per edge when weights exist (forwarded to \code{plot_layer_diff()}).
#'   One of \code{"sum","mean","median","first"}. Default \code{"sum"}.
#' @param bg Background color for grid and JPEG conversions. Default \code{"white"}.
#' @param ... Additional aesthetics forwarded to \code{plot_layer_diff()} (layout, colors, sizes, labels, seed, etc.).
#'
#' @return (Invisibly) a list with:
#' \itemize{
#'   \item \code{pairs} — \code{data.frame} of consecutive pairs (\code{L1}, \code{L2});
#'   \item \code{pair_files} — named character vector of per‑pair file paths (names \code{"L1_vs_L2"});
#'   \item \code{grid_file} — absolute path to the timestamped grid figure;
#'   \item \code{combined_csv} — absolute path to the timestamped combined CSV;
#'   \item \code{edges} — combined tidy \code{data.frame} across all pairs;
#'   \item \code{counts} — per‑pair summary (\code{only_<L1>}, \code{only_<L2>}, \code{common}, \code{total}).
#' }
#'
#' @importFrom png readPNG
#' @export
grid_layer_diffs <- function(
    net,
    layer_order        = NULL,
    ncol               = 3,
    results_dir        = getOption("mlnet.results_dir","omicsDNA_results"),
    grid_file          = NULL,
    grid_format        = c("png","pdf","jpg"),
    width              = 12,
    height             = 10,
    dpi                = 300,
    pair_format        = c("png","pdf","jpg"),
    pair_file_prefix   = "diff_pair",
    panel_width        = 6,
    panel_height       = 5,
    dpi_panel          = dpi,
    combined_csv_file  = NULL,
    weight_attr        = c("weight","w","value","prop_present"),
    weight_aggregate   = c("sum","mean","median","first"),
    bg                 = "white",
    ...
) {
  # ---- deps ----
  if (!requireNamespace("png", quietly = TRUE)) {
    stop("Package 'png' is required. Install with install.packages('png').")
  }

  grid_format      <- match.arg(grid_format)
  pair_format      <- match.arg(pair_format)
  weight_aggregate <- match.arg(weight_aggregate)

  # ---- helpers ----
  .is_abs <- function(p) grepl("^(/|[A-Za-z]:[\\/])", p)
  .ensure_dir <- function(d) if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  .stamp <- function() format(Sys.time(), "%Y-%m-%d_%H%M%S")
  .append_stamp <- function(path, stamp, default_ext = NULL) {
    if (is.null(path) || !nzchar(path)) return(NULL)
    if (!grepl("\\.[A-Za-z0-9]+$", path) && !is.null(default_ext)) path <- paste0(path, ".", default_ext)
    sub("(\\.[A-Za-z0-9]+)$", paste0("_", stamp, "\\1"), path)
  }
  `%||%` <- function(x, y) if (!is.null(x)) x else y

  # ---- layers & pairs ----
  Ls <- try(multinet::layers_ml(net), silent = TRUE)
  if (inherits(Ls, "try-error") || !length(Ls))
    stop("Could not retrieve layers via multinet::layers_ml(net).")
  Ls <- as.character(Ls)

  if (is.null(layer_order)) {
    layers <- Ls
  } else {
    layers <- as.character(layer_order)
    layers <- layers[layers %in% Ls]
    if (!length(layers)) stop("No valid layers after applying `layer_order`.")
  }
  if (length(layers) < 2L) stop("Need at least two layers to build consecutive pairs.")

  pairs_df <- data.frame(
    L1 = head(layers, -1L),
    L2 = tail(layers, -1L),
    stringsAsFactors = FALSE
  )
  n_pairs <- nrow(pairs_df)

  # ---- I/O prep ----
  .ensure_dir(results_dir)
  # *** keep results_dir ABSOLUTE so downstream file paths are absolute ***
  results_dir <- normalizePath(results_dir, winslash = "/", mustWork = FALSE)

  stamp <- .stamp()

  # NEW: timestamped run_dir that will hold all outputs for this call -----------
  run_dir <- file.path(results_dir, paste0("layer_diffs_", stamp))   # NEW: run_dir
  .ensure_dir(run_dir)                                               # NEW: create it
  # ---------------------------------------------------------------------------

  # grid output path (always timestamped)
  if (is.null(grid_file) || !nzchar(grid_file)) {
    grid_file <- file.path(run_dir, paste0("diffs_all_pairs_", stamp, ".", grid_format))   # (run_dir)
  } else {
    if (!.is_abs(grid_file)) grid_file <- file.path(run_dir, grid_file)                    # (run_dir)
    grid_file <- .append_stamp(grid_file, stamp, default_ext = grid_format)
  }

  # combined CSV output path (always timestamped)
  if (is.null(combined_csv_file) || !nzchar(combined_csv_file)) {
    combined_csv_file <- file.path(run_dir, paste0("diff_edges_all_pairs_", stamp, ".csv"))  # (run_dir)
  } else {
    if (!.is_abs(combined_csv_file)) combined_csv_file <- file.path(run_dir, combined_csv_file)  # (run_dir)
    combined_csv_file <- .append_stamp(combined_csv_file, stamp, default_ext = "csv")
  }

  # ---- per-pair processing ----
  panel_pngs <- character(n_pairs)
  pair_files <- character(n_pairs)
  names(pair_files) <- paste0(pairs_df$L1, "_vs_", pairs_df$L2)

  edges_list  <- vector("list", n_pairs)
  counts_list <- vector("list", n_pairs)

  for (i in seq_len(n_pairs)) {
    L1i <- pairs_df$L1[i]
    L2i <- pairs_df$L2[i]
    pair_base <- file.path(run_dir, sprintf("%s_%s_vs_%s", pair_file_prefix, L1i, L2i))  # (run_dir)
    out_edges <- NULL

    if (identical(pair_format, "png")) {
      # per-pair PNG (absolute path) used both as deliverable and as grid panel
      pair_png <- paste0(pair_base, "_", stamp, ".png")
      out <- do.call(
        plot_layer_diff,
        c(
          list(
            net              = net,
            L1               = L1i,
            L2               = L2i,
            file             = pair_png,
            format           = "png",
            results_dir      = run_dir,   # (run_dir) safe; file is absolute anyway
            save_edge_csv    = FALSE,
            weight_attr      = weight_attr,
            weight_aggregate = weight_aggregate
          ),
          list(...)
        )
      )
      panel_pngs[i] <- pair_png
      pair_files[i] <- normalizePath(pair_png, winslash = "/", mustWork = FALSE)
      out_edges <- out$edges

    } else if (identical(pair_format, "pdf")) {
      # deliverable PDF + extra PNG (absolute temp path) for the grid
      pair_pdf <- paste0(pair_base, "_", stamp, ".pdf")
      do.call(
        plot_layer_diff,
        c(
          list(
            net              = net, L1 = L1i, L2 = L2i,
            file             = pair_pdf, format = "pdf",
            results_dir      = run_dir,   # (run_dir)
            save_edge_csv    = FALSE,
            weight_attr      = weight_attr, weight_aggregate = weight_aggregate
          ),
          list(...)
        )
      )
      pair_files[i] <- normalizePath(pair_pdf, winslash = "/", mustWork = FALSE)

      panel_png <- tempfile(pattern = sprintf("panel_%s_vs_%s_", L1i, L2i), fileext = ".png")
      out <- do.call(
        plot_layer_diff,
        c(
          list(
            net              = net, L1 = L1i, L2 = L2i,
            file             = panel_png, format = "png",
            results_dir      = tempdir(),      # temp is absolute
            save_edge_csv    = FALSE,
            weight_attr      = weight_attr, weight_aggregate = weight_aggregate
          ),
          list(...)
        )
      )
      panel_pngs[i] <- panel_png
      out_edges <- out$edges

    } else { # pair_format == "jpg"
      # render PNG panel in temp (absolute), then convert to JPG deliverable
      panel_png <- tempfile(pattern = sprintf("panel_%s_vs_%s_", L1i, L2i), fileext = ".png")
      out <- do.call(
        plot_layer_diff,
        c(
          list(
            net              = net, L1 = L1i, L2 = L2i,
            file             = panel_png, format = "png",
            results_dir      = tempdir(),      # temp is absolute
            save_edge_csv    = FALSE,
            weight_attr      = weight_attr, weight_aggregate = weight_aggregate
          ),
          list(...)
        )
      )
      panel_pngs[i] <- panel_png
      out_edges <- out$edges

      pair_jpg <- paste0(pair_base, "_", stamp, ".jpg")
      grDevices::jpeg(pair_jpg, width = panel_width, height = panel_height, units = "in",
                      res = dpi_panel, quality = 95, bg = bg)
      op <- par(mar = c(0,0,0,0), xaxs = "i", yaxs = "i"); on.exit(par(op), add = TRUE)
      plot.new()
      img <- png::readPNG(panel_png)
      rasterImage(img, 0, 0, 1, 1)
      grDevices::dev.off()
      pair_files[i] <- normalizePath(pair_jpg, winslash = "/", mustWork = FALSE)
    }

    # collect edges + counts
    if (!is.null(out_edges) && nrow(out_edges)) {
      e <- out_edges
      e$pair_L1 <- L1i; e$pair_L2 <- L2i
      e <- e[, c("pair_L1","pair_L2","from","to","category","weight_L1","weight_L2"), drop = FALSE]
      edges_list[[i]] <- e

      ct <- table(e$category)
      counts_list[[i]] <- data.frame(
        pair_L1  = L1i,
        pair_L2  = L2i,
        only_L1  = as.integer(ct[paste0("only_", L1i)] %||% 0L),
        only_L2  = as.integer(ct[paste0("only_", L2i)] %||% 0L),
        common   = as.integer(ct["common"]             %||% 0L),
        total    = nrow(e),
        stringsAsFactors = FALSE
      )
    } else {
      counts_list[[i]] <- data.frame(
        pair_L1  = L1i, pair_L2 = L2i,
        only_L1  = 0L,  only_L2 = 0L, common = 0L, total = 0L,
        stringsAsFactors = FALSE
      )
    }
  }

  # ---- one combined CSV only ----
  combined_edges <- if (any(vapply(edges_list, function(x) !is.null(x) && nrow(x) > 0, logical(1)))) {
    do.call(rbind, edges_list)
  } else {
    data.frame(pair_L1 = character(), pair_L2 = character(), from = character(),
               to = character(), category = character(),
               weight_L1 = numeric(), weight_L2 = numeric(),
               stringsAsFactors = FALSE)
  }
  utils::write.csv(combined_edges, combined_csv_file, row.names = FALSE)
  message("Saved combined edges CSV: ", normalizePath(combined_csv_file, winslash = "/", mustWork = FALSE))

  counts_df <- do.call(rbind, counts_list)

  # ---- draw the grid from PNG panels ----
  nr <- ceiling(n_pairs / ncol)
  draw_grid <- function() {
    op <- par(mfrow = c(nr, ncol), mar = c(0, 0, 2, 0), bg = bg); on.exit(par(op), add = TRUE)
    for (i in seq_len(n_pairs)) {
      if (!nzchar(panel_pngs[i]) || !file.exists(panel_pngs[i])) {
        plot.new(); title(sprintf("%s → %s (no panel)", pairs_df$L1[i], pairs_df$L2[i]), cex.main = 0.9)
        next
      }
      img <- png::readPNG(panel_pngs[i])
      plot.new()
      rasterImage(img, 0, 0, 1, 1)
      title(sprintf("%s → %s", pairs_df$L1[i], pairs_df$L2[i]), cex.main = 0.9)
    }
  }

  if (identical(grid_format, "png")) {
    grDevices::png(grid_file, width = width, height = height, units = "in", res = dpi, bg = bg)
    on.exit(grDevices::dev.off(), add = TRUE)
    draw_grid()
  } else if (identical(grid_format, "jpg")) {
    grDevices::jpeg(grid_file, width = width, height = height, units = "in", res = dpi, bg = bg, quality = 95)
    on.exit(grDevices::dev.off(), add = TRUE)
    draw_grid()
  } else {
    grDevices::pdf(grid_file, width = width, height = height, onefile = FALSE, paper = "special", bg = bg)
    on.exit(grDevices::dev.off(), add = TRUE)
    draw_grid()
  }
  message("Saved grid of diffs: ", normalizePath(grid_file, winslash = "/", mustWork = FALSE))

  # clean up temporary PNG panels (keep per-pair outputs)
  tmp_idx <- !startsWith(panel_pngs, results_dir)
  try(unlink(panel_pngs[tmp_idx & file.exists(panel_pngs)]), silent = TRUE)

  invisible(list(
    pairs        = pairs_df,
    pair_files   = pair_files,
    grid_file    = grid_file,
    combined_csv = combined_csv_file,
    edges        = combined_edges,
    counts       = counts_df
  ))
}

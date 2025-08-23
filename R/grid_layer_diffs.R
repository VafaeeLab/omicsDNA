#' Grid of consecutive layer differences (E1→E2, E2→E3, …)
#'
#' @param net multinet object.
#' @param layer_order Optional explicit order; default = multinet::layers_ml(net).
#' @param ncol Number of columns in the grid.
#' @param results_dir, file Optional PNG output. If file = NULL, draws to the current device.
#' @param width,height,dpi Image size when saving.
grid_layer_diffs <- function(net, layer_order = NULL, ncol = 3,
                             results_dir = getOption("mlnet.results_dir","omicsDNA_results"),
                             file = NULL, width = 12, height = 10, dpi = 300) {
  layers <- try(multinet::layers_ml(net), silent = TRUE)
  if (inherits(layers, "try-error") || !length(layers)) stop("Could not get layers.")
  layers <- as.character(layers)
  if (!is.null(layer_order)) layers <- intersect(as.character(layer_order), layers)
  if (length(layers) < 2) stop("Need at least two layers.")
  pairs <- Map(c, head(layers, -1L), tail(layers, -1L))

  # precompute plots into temporary pngs; then arrange them
  png_paths <- character(length(pairs))
  for (i in seq_along(pairs)) {
    f <- tempfile(fileext = ".png")
    plot_layer_diff(net, L1 = pairs[[i]][1], L2 = pairs[[i]][2],
                    file = f, width = 6, height = 5, dpi = dpi)
    png_paths[i] <- f
  }

  # draw them in a grid to the current device or a final PNG
  draw_grid <- function() {
    op <- par(mfrow = c(ceiling(length(pairs)/ncol), ncol), mar = c(0,0,2,0))
    on.exit(par(op), add = TRUE)
    for (i in seq_along(png_paths)) {
      img <- png::readPNG(png_paths[i])
      plot.new(); rasterImage(img, 0, 0, 1, 1)
      title(sprintf("%s → %s", pairs[[i]][1], pairs[[i]][2]), cex.main = 0.9)
    }
  }

  if (!is.null(file) && nzchar(file)) {
    if (!grepl("\\.png$", file, ignore.case = TRUE)) file <- paste0(file, ".png")
    if (!grepl("^(/|[A-Za-z]:[\\/])", file)) file <- file.path(results_dir, file)
    grDevices::png(file, width = width, height = height, units = "in", res = dpi, bg = "white")
    on.exit(grDevices::dev.off(), add = TRUE)
    draw_grid()
    message("Saved grid of diffs: ", normalizePath(file, winslash = "/", mustWork = FALSE))
  } else {
    draw_grid()  # show in RStudio
  }

  invisible(list(pairs = pairs, files = png_paths))
}

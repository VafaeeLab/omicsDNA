#' Plot multilayer communities with an interactive preview and file exports:
#'
#' @description
#' This routine visualises community assignments on a multilayer network and
#' simultaneously supports reproducible export of the figure and the underlying
#' long table of assignments. By default it:
#' (i) draws to the current graphics device (e.g., the RStudio *Plots* pane) so
#' you can inspect or manually *Export*; (ii) writes a high‑resolution image to
#' disk; and (iii) saves a tidy `(actor, layer, cid)` table for downstream use.
#'
#' @details
#' **What it does**
#' 1. Validates and normalises the `communities` table to ensure one community ID
#'    (`cid`) per `(actor, layer)` pair. If multiple rows exist for the same pair,
#'    the first is kept (a message is printed when `verbose = TRUE`).
#' 2. Intersects the requested layers with those actually present in the network
#'    (queried via `multinet::edges_ml()`), and drops community rows whose actors
#'    do not appear in the corresponding layer. This prevents plotting nodes that
#'    are absent from the chosen layers.
#' 3. Computes a multilayer layout using either
#'    `multinet::layout_multiforce_ml()` (force‑directed; parameterised by
#'    `gravity`) or `multinet::layout_circular_ml()` (arranged on circles).
#' 4. Calls the multilayer network `plot()` method once, passing the community
#'    mapping as a three‑column data frame (`actor`, `layer`, `cid`) and laying
#'    out the selected layers on a user‑controlled grid of panels.
#' 5. Optionally duplicates the on‑screen plot to a file with
#'    `grDevices::dev.copy()` (when `show_in_rstudio = TRUE`) or renders directly
#'    to an off‑screen device (when `show_in_rstudio = FALSE`). It also saves the
#'    long community table in `.rds` and/or `.csv` format without printing it.
#'
#' **How community labels are handled**
#' - The function accepts one of `cid`, `com`, or `community` in `communities`.
#'   If the column is numeric, it is used as‑is; otherwise labels are converted
#'   to integer IDs via factor coding for plotting purposes. This does not alter
#'   your original labels upstream (e.g., `"C1"`, `"C2"`)—only the numeric
#'   encoding passed to the plotting backend.
#'
#' **Layer selection and grid arrangement**
#' - If `layersToPlot = NULL`, the function visualises the intersection of the
#'   layers present in `communities` and in the network. Otherwise, the supplied
#'   subset is intersected with available layers.
#' - When `grid = NULL`, a near‑square arrangement is chosen automatically
#'   (number of rows ≈ `sqrt(#layers)`); you may override with `grid = c(nr, nc)`.
#'
#' **Persistence and file names**
#' - Figures are written under `results_dir` using an informative stem that
#'   includes the layout and number of layers, e.g.,
#'   `plotCom_multiforce_4layers_<timestamp>.png`. Provide `file` to override.
#' - The long table is saved as RDS and/or CSV. When `df_file` is provided, it
#'   is honoured **for the RDS** (absolute paths are respected; relative paths
#'   are created under `results_dir`). For CSV, the base name of `df_file`
#'   (without extension) is used under `results_dir`.
#'
#' @param net A `multinet::ml.network` to be plotted (must be compatible with
#'   `multinet::edges_ml()` and the multilayer `plot()` method).
#' @param communities A data frame containing at least `actor` and `layer`, plus
#'   one of `cid`, `com`, or `community` indicating community membership.
#' @param layout Character; multilayer layout to use: `"multiforce"` (force‑
#'   directed; supports `gravity`) or `"circular"`. Default `c("multiforce","circular")`
#'   (matched to `"multiforce"`).
#' @param layersToPlot Optional character vector selecting which layers to draw.
#'   By default, plots the intersection of layers present in `communities` and
#'   in the network.
#' @param grid `NULL` (automatic near‑square arrangement) or an integer vector
#'   `c(nrow, ncol)` specifying the panel grid. Default `NULL`.
#' @param vertex.size Numeric; node size passed to the plotting backend.
#'   Default `5`.
#' @param vertex.cex Numeric; node size multiplier passed to the plotting
#'   backend. Default `1.2`.
#' @param show.labels Logical; show vertex labels (`TRUE`) or suppress them
#'   (`FALSE`, default). Internally this toggles `vertex.labels`.
#' @param gravity Numeric in roughly `[0, 1]`; attraction strength for the
#'   multiforce layout. Ignored when `layout = "circular"`. Default `0.3`.
#' @param seed Optional integer seed for reproducible layout initialisation.
#'   Default `NULL`.
#' @param results_dir Directory where outputs are written. Default
#'   `getOption("mlnet.results_dir", "omicsDNA_results")`.
#' @param show_in_rstudio Logical; if `TRUE` (default), draw to the current
#'   device (e.g., RStudio *Plots* pane) for interactive inspection.
#' @param save_plot Logical; if `TRUE`, save the figure to a file. Default `TRUE`.
#' @param file Optional output path for the figure. If `NULL`, an informative
#'   name is constructed under `results_dir`. Relative paths are resolved under
#'   `results_dir`; absolute paths are respected.
#' @param format Image format for the saved figure: `"png"` or `"pdf"`. Default
#'   `c("png","pdf")` (matched to `"png"`).
#' @param width,height Numeric dimensions of the saved figure. Interpreted in
#'   `units` for PNG and in inches for PDF. Defaults `10` × `8`.
#' @param units Character; units for PNG dimensions (`"in"`, `"cm"`, or `"mm"`).
#'   Default `"in"`. Ignored for PDF.
#' @param dpi Numeric; raster resolution (PNG only). Default `300`.
#' @param save_df Logical; if `TRUE`, save the long `(actor, layer, cid)` table.
#'   Default `TRUE`.
#' @param df_format Either `"rds"`, `"csv"`, or both (e.g., `c("rds","csv")`).
#'   Default `"rds"`.
#' @param df_file Optional file name for the RDS output. If relative, it is
#'   created under `results_dir`; absolute paths are respected. (For CSV, the
#'   base name of `df_file` is used under `results_dir`.)
#' @param verbose Logical; print informative messages (dropped rows, file paths).
#'   Default `TRUE`.
#' @param ... Additional arguments forwarded to the multilayer network `plot()`
#'   method (e.g., colour palettes or edge styling parameters supported by your
#'   plotting backend).
#'
#' @return (Invisibly) the normalised community table used for plotting with
#'   columns `actor`, `layer`, `cid`. The returned object carries attributes:
#'   - `file`: path to the saved figure (if `save_plot = TRUE`);
#'   - `df_file`: character vector of paths to the saved RDS/CSV (if any);
#'   - `grid`: the panel grid used;
#'   - `layers`: the layers actually plotted.
#'
#' @section Practical notes
#' - When `show_in_rstudio = TRUE` and `save_plot = TRUE`, the function clones
#'   the current plot to a file using `grDevices::dev.copy()`. To render only to
#'   an off‑screen file (without drawing in the Plots pane), set
#'   `show_in_rstudio = FALSE`.
#' - The function only plots actors that appear in the chosen layers of the
#'   network; assignments for absent actors/layers are dropped (reported when
#'   `verbose = TRUE`).
#'
#' @examples
#' \dontrun{
#' # Minimal usage: preview in RStudio, save PNG and the long table as RDS
#' plotCom(
#'   net, comm,
#'   layout          = "multiforce",
#'   gravity         = 0.3,
#'   show_in_rstudio = TRUE,
#'   save_plot       = TRUE,
#'   format          = "png",
#'   save_df         = TRUE,
#'   df_format       = "rds",
#'   seed            = 1
#' )
#'
#' # Custom grid and both table formats
#' plotCom(
#'   net, comm,
#'   layersToPlot    = c("Young","Old"),
#'   grid            = c(1, 2),
#'   save_df         = TRUE,
#'   df_format       = c("rds","csv")
#' )
#' }
#'
#' @seealso
#'   \code{\link{detectCom}} for obtaining community assignments;
#'   \code{\link{build_multiNet}} for assembling multilayer networks.
#'
#' @importFrom multinet layout_multiforce_ml layout_circular_ml edges_ml
#' @importFrom grDevices png pdf dev.off dev.copy
#' @importFrom utils write.csv
#' @export
plotCom <- function(net,
                    communities,
                    layout        = c("multiforce","circular"),
                    layersToPlot  = NULL,
                    grid          = NULL,
                    vertex.size   = 5,
                    vertex.cex    = 1.2,
                    show.labels   = FALSE,
                    gravity       = 0.3,
                    seed          = NULL,
                    results_dir   = getOption("mlnet.results_dir","omicsDNA_results"),
                    # NEW:
                    show_in_rstudio = TRUE,
                    save_plot     = TRUE,
                    file          = NULL,
                    format        = c("png","pdf"),
                    width         = 10,
                    height        = 8,
                    units         = "in",
                    dpi           = 300,
                    # NEW (save the long df, don't print):
                    save_df       = TRUE,
                    df_format     = "rds",  # one or both of c("rds","csv")
                    df_file       = NULL,
                    verbose       = TRUE,
                    ...) {
  layout <- match.arg(layout)
  format <- match.arg(format)

  .ensure_dir <- function(d) if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  .is_abs     <- function(p) grepl("^(/|[A-Za-z]:[\\/])", p)
  .pick       <- function(cands, nms) { z <- cands[cands %in% nms]; if (length(z)) z[1] else NA_character_ }

  if (!is.data.frame(communities))
    stop("`communities` must be a data.frame.")
  if (!all(c("actor","layer") %in% names(communities)))
    stop("`communities` must contain columns `actor` and `layer`.")

  # Accept 'cid' or 'com' or 'community' and coerce to numeric cid starting at 1
  comm_col <- intersect(c("cid","com","community"), names(communities))
  if (!length(comm_col))
    stop("`communities` must include either `cid`, `com`, or `community`.")
  comm_col <- comm_col[1L]
  cid_raw  <- communities[[comm_col]]
  cid_num  <- if (is.numeric(cid_raw)) as.integer(cid_raw) else as.integer(factor(cid_raw))

  df <- data.frame(
    actor = as.character(communities$actor),
    layer = as.character(communities$layer),
    cid   = as.integer(cid_num),
    stringsAsFactors = FALSE
  )
  # De‑duplicate: one cid per (actor, layer)
  if (any(duplicated(df[, c("actor","layer")]))) {
    if (verbose) message("Duplicate (actor,layer) rows in `communities`: keeping the first for each pair.")
    df <- df[!duplicated(df[, c("actor","layer")]), , drop = FALSE]
  }

  # --------- Intersect with what actually exists in the network ----------
  Eraw <- try(multinet::edges_ml(net), silent = TRUE)
  if (inherits(Eraw, "try-error") || is.null(Eraw))
    stop("Could not retrieve edges from `net` via multinet::edges_ml().")
  Eall <- if (is.data.frame(Eraw)) Eraw else {
    tmp <- try(as.data.frame(Eraw, stringsAsFactors = FALSE), silent = TRUE)
    if (inherits(tmp, "try-error") || is.null(tmp) || !nrow(tmp)) {
      stop("`edges_ml(net)` did not return a coercible table of edges.")
    }
    tmp
  }
  nm <- names(Eall)

  # Robust endpoint/layer detection
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
      stop("Could not identify two endpoint columns in the network edge table.")
    a_col <- nm[char_cols[1]]; b_col <- nm[char_cols[2]]
  }

  if ("from_layer" %in% nm && "to_layer" %in% nm) {
    Eall <- Eall[Eall$from_layer == Eall$to_layer, , drop = FALSE]
    Eall$layer <- as.character(Eall$from_layer)
  } else if ("layer" %in% nm) {
    Eall$layer <- as.character(Eall$layer)
  } else {
    Eall$layer <- "L1"
  }
  if (!nrow(Eall)) stop("No intra-layer edges found in the network.")

  # Layers to plot (intersection of requested, communities, and network)
  layers_in_comm <- unique(df$layer)
  layers_in_net  <- sort(unique(Eall$layer))
  if (is.null(layersToPlot)) {
    layersToPlot <- intersect(layers_in_comm, layers_in_net)
  } else {
    layersToPlot <- intersect(as.character(layersToPlot), layers_in_net)
  }
  if (!length(layersToPlot)) stop("No valid layers to plot (after intersection with the network).")

  # Filter community rows to actors present in each selected layer of the network
  actors_by_layer <- lapply(layersToPlot, function(ly) {
    ed <- Eall[Eall$layer == ly, , drop = FALSE]
    unique(c(as.character(ed[[a_col]]), as.character(ed[[b_col]])))
  })
  names(actors_by_layer) <- layersToPlot

  keep_idx <- df$layer %in% names(actors_by_layer) &
    mapply(function(act, ly) act %in% actors_by_layer[[ly]], df$actor, df$layer)
  if (verbose) {
    dropped <- sum(!keep_idx)
    if (dropped > 0) message("Dropping ", dropped, " community rows (actors not present in selected layers).")
  }
  df <- df[keep_idx, , drop = FALSE]
  if (!nrow(df)) stop("No community assignments left after filtering to actors present in selected layers.")

  # --------- Layout ----------
  if (!is.null(seed)) set.seed(as.integer(seed))
  lay_coords <- if (layout == "multiforce") {
    res <- try(multinet::layout_multiforce_ml(net, gravity = gravity), silent = TRUE)
    if (inherits(res, "try-error")) {
      if (verbose) message("`layout_multiforce_ml` failed with gravity=", gravity, "; trying default.")
      multinet::layout_multiforce_ml(net)
    } else res
  } else {
    multinet::layout_circular_ml(net)
  }

  # --------- Grid ----------
  if (is.null(grid)) {
    k <- length(layersToPlot)
    nr <- max(1, floor(sqrt(k)))
    nc <- ceiling(k / nr)
    grid <- c(nr, nc)
  } else {
    if (!(is.numeric(grid) && length(grid) == 2L && all(is.finite(grid)) && all(grid >= 1)))
      stop("`grid` must be NULL or a numeric length-2 vector (nrow, ncol) >= 1.")
    grid <- as.integer(grid)
  }

  # --------- Labels on/off ----------
  vlabels <- if (show.labels) NULL else NA

  # --------- Ensure results dir ----------
  .ensure_dir(results_dir)
  stamp <- format(Sys.time(), "%Y-%m-%d_%H%M%S")

  # --------- 1) Draw to RStudio (current device) so user can Export ----------
  do_plot <- function() {
    plot(net,
         layout        = lay_coords,
         grid          = grid,
         layers        = layersToPlot,
         com           = df,                # actor, layer, cid
         vertex.size   = vertex.size,
         vertex.cex    = vertex.cex,
         vertex.labels = vlabels,
         ...)
  }
  if (isTRUE(show_in_rstudio)) {
    do_plot()  # shows in the Plots pane; user can Export from RStudio
  }

  # --------- 2) Save the plot under results/ ----------
  if (isTRUE(save_plot)) {
    if (is.null(file)) {
      base <- sprintf("plotCom_%s_%dlayers_%s", tolower(layout), length(layersToPlot), stamp)
      file <- paste0(base, ".", format)
    }
    if (!.is_abs(file)) file <- file.path(results_dir, file)

    if (isTRUE(show_in_rstudio)) {
      # copy what you see to the chosen device (keeps the RStudio plot visible)
      if (format == "png") {
        dev_id <- grDevices::dev.copy(grDevices::png, filename = file,
                                      width = width, height = height, units = units, res = dpi)
        grDevices::dev.off(dev_id)
      } else {
        dev_id <- grDevices::dev.copy(grDevices::pdf, file = file,
                                      width = width, height = height)  # inches
        grDevices::dev.off(dev_id)
      }
    } else {
      # open off‑screen device, plot, close
      if (format == "png") {
        grDevices::png(filename = file, width = width, height = height, units = units, res = dpi)
        on.exit(grDevices::dev.off(), add = TRUE)
        do_plot()
      } else {
        grDevices::pdf(file, width = width, height = height)  # inches
        on.exit(grDevices::dev.off(), add = TRUE)
        do_plot()
      }
    }
    if (verbose) message("Saved ", toupper(format), ": ", normalizePath(file, FALSE))
  }

  # --------- 3) Save the (actor, layer, cid) table under results/ ----------
  df_paths <- character(0)
  if (isTRUE(save_df)) {
    df_formats <- unique(match.arg(df_format, several.ok = TRUE))
    if (is.null(df_file)) {
      base <- sprintf("plotCom_communities_%s_%dlayers_%s", tolower(layout), length(layersToPlot), stamp)
    } else {
      base <- tools::file_path_sans_ext(basename(df_file))
    }
    if ("rds" %in% df_formats) {
      rds_path <- if (is.null(df_file)) file.path(results_dir, paste0(base, ".rds")) else
        if (.is_abs(df_file)) df_file else file.path(results_dir, df_file)
      saveRDS(df, rds_path)
      df_paths <- c(df_paths, rds_path)
      if (verbose) message("Saved community table (RDS): ", normalizePath(rds_path, FALSE))
    }
    if ("csv" %in% df_formats) {
      csv_path <- file.path(results_dir, paste0(base, ".csv"))
      utils::write.csv(df, csv_path, row.names = FALSE)
      df_paths <- c(df_paths, csv_path)
      if (verbose) message("Saved community table (CSV): ", normalizePath(csv_path, FALSE))
    }
  }

  # Return invisibly (so nothing prints to the Console)
  res <- structure(df,
                   file    = if (isTRUE(save_plot)) file else NULL,
                   df_file = if (length(df_paths)) df_paths else NULL,
                   grid    = grid,
                   layers  = layersToPlot)
  return(invisible(res))
}

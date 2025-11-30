
# ---------------------------------------
# 12 - plotCom
# ---------------------------------------

#' Plot multilayer communities with an interactive preview and file exports
#'
#' @description
#' Visualise community assignments on a multilayer network while **also**
#' exporting (i) a high‑resolution image and (ii) a tidy long table
#' `(actor, layer, cid)` for downstream analyses. By default the function
#' (1) draws to the current graphics device (e.g., the RStudio *Plots* pane),
#' (2) writes an image to disk, and (3) saves the normalised assignment table.
#'
#' @details
#' **What this function guarantees**
#'
#' 1. **Normalisation** — The `communities` input is validated and reduced to a
#'    single `cid` per `(actor, layer)` pair. If duplicates exist, the first is
#'    kept (message shown when `verbose = TRUE`). The community column can be
#'    any of `cid`, `com`, or `community`. Non‑numeric labels are factor‑coded to
#'    integers for plotting; your original labels are unaffected upstream.
#' 2. **Layer & actor intersection** — Only layers present in the network are
#'    plotted, and only actors that actually appear in those layers are kept.
#'    Intra‑layer edges are detected via `multinet::edges_ml()`; cross‑layer
#'    edges are discarded for the purpose of per‑layer drawing.
#' 3. **Layout** — Choose either a force‑directed layout
#'    (`"multiforce"`, accepts `gravity`) or a `"circular"` arrangement via
#'    `multinet::layout_multiforce_ml()` / `layout_circular_ml()`.
#' 4. **Grid arrangement** — Layers are laid out on a near‑square panel grid
#'    unless you supply `grid = c(nrow, ncol)`.
#' 5. **Reproducible export** — If `show_in_rstudio = TRUE`, you see the plot on
#'    screen and the same content is copied to file via `grDevices::dev.copy()`.
#'    If `show_in_rstudio = FALSE`, the figure is rendered directly to an
#'    off‑screen device. The long `(actor, layer, cid)` table is saved as RDS
#'    and/or CSV without printing to the console.
#'
#' **File naming & persistence**
#'
#' - Output files live under `results_dir` unless you pass absolute paths.
#' - If `file` is `NULL`, a descriptive stem is generated, e.g.:
#'   `plotCom_multiforce_4layers_<timestamp>.png`.
#' - For the table, if `df_file` is provided, it is used for the **RDS** path
#'   (relative paths resolved under `results_dir`, absolute respected). For CSV,
#'   the **base name** of `df_file` (sans extension) is used under `results_dir`.
#'
#' @param net A `multinet::ml.network`. Must be compatible with
#'   `multinet::edges_ml()` and the multilayer `plot()` method.
#' @param communities `data.frame` with at least `actor` and `layer`, plus one of
#'   `cid`, `com`, or `community` describing community membership.
#' @param layout Character; `"multiforce"` (force‑directed; supports `gravity`)
#'   or `"circular"`. Default `c("multiforce","circular")` (matched to `"multiforce"`).
#' @param layersToPlot Optional character vector of layer names to draw. If
#'   `NULL`, plots the intersection of layers present in `communities` and in
#'   the network.
#' @param grid `NULL` (auto near‑square) or integer vector `c(nrow, ncol)`.
#' @param vertex.size Numeric; node size passed to the plotting backend. Default `5`.
#' @param vertex.cex Numeric; node size multiplier. Default `1.2`.
#' @param show.labels Logical; whether to show vertex labels. Default `FALSE`.
#' @param gravity Numeric (roughly `[0, 1]`); attraction for multiforce layout.
#'   Ignored when `layout = "circular"`. Default `0.3`.
#' @param seed Optional integer seed for reproducible layout initialisation. Default `NULL`.
#' @param results_dir Directory where outputs are written. Default
#'   `getOption("mlnet.results_dir", "omicsDNA_results")`.
#' @param show_in_rstudio Logical; if `TRUE` (default) draw to current device
#'   (interactive preview).
#' @param save_plot Logical; if `TRUE`, save the figure to a file. Default `TRUE`.
#' @param file Optional output path for the figure. If `NULL`, constructed under
#'   `results_dir`. Relative paths are resolved under `results_dir`; absolute paths are respected.
#' @param format Image format for the saved figure: `"png"` or `"pdf"`. Default
#'   `c("png","pdf")` (matched to `"png"`).
#' @param width,height Numeric dimensions of the saved figure. Interpreted in
#'   `units` for PNG; inches for PDF. Defaults `10` × `8`.
#' @param units Character units for PNG dimensions (`"in"`, `"cm"`, or `"mm"`).
#'   Ignored for PDF. Default `"in"`.
#' @param dpi Numeric; raster resolution (PNG only). Default `300`.
#' @param save_df Logical; if `TRUE`, save the long `(actor, layer, cid)` table.
#'   Default `TRUE`.
#' @param df_format `"rds"`, `"csv"`, or both (e.g., `c("rds","csv")`). Default `"rds"`.
#' @param df_file Optional file name for the **RDS** output. Relative paths are
#'   created under `results_dir`; absolute paths are respected. For CSV, the base
#'   name of `df_file` is used under `results_dir`.
#' @param verbose Logical; print informative messages (dropped rows, file paths).
#'   Default `TRUE`.
#' @param ... Additional arguments forwarded to the multilayer network `plot()`
#'   method (e.g., colour palette or edge styling supported by your plotting backend).
#'
#' @return (Invisibly) the normalised `data.frame` used for plotting with
#'   columns `actor`, `layer`, `cid`. The object carries attributes:
#'   - `file`: path to the saved figure (if `save_plot = TRUE`)
#'   - `df_file`: character vector of saved RDS/CSV paths (if any)
#'   - `grid`: the panel grid used
#'   - `layers`: the layers actually plotted
#'
#' @section Practical notes:
#' - To draw **only** to file (without preview), set `show_in_rstudio = FALSE`.
#' - Assignments for actors absent from the chosen layers are dropped
#'   (reported when `verbose = TRUE`).
#' - `vertex.labels` is suppressed by passing `NA` unless `show.labels = TRUE`.
#'
#' @examples
#' \dontrun{
#' # Minimal usage: preview in RStudio, save PNG and the long table as RDS
#' out <- plotCom(
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
#' attr(out, "file")     # path to the image
#' attr(out, "df_file")  # path(s) to saved table(s)
#'
#' # Custom layer subset, custom grid, and both table formats
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
#'   \code{\link{detectCom}} for community detection;
#'   \code{\link{build_multiNet}} for assembling multilayer networks.
#'
#' @importFrom multinet layout_multiforce_ml layout_circular_ml edges_ml
#' @importFrom grDevices png pdf dev.off dev.copy
#' @importFrom utils write.csv
#' @export
plotCom <- function(net,
                    communities,
                    layout          = c("multiforce","circular"),
                    layersToPlot    = NULL,
                    grid            = NULL,
                    vertex.size     = 5,
                    vertex.cex      = 1.2,
                    show.labels     = FALSE,
                    gravity         = 0.3,
                    seed            = NULL,
                    results_dir     = getOption("mlnet.results_dir","omicsDNA_results"),
                    show_in_rstudio = TRUE,
                    save_plot       = TRUE,
                    file            = NULL,
                    format          = c("png","pdf"),
                    width           = 10,
                    height          = 8,
                    units           = "in",
                    dpi             = 300,
                    save_df         = TRUE,
                    df_format       = "rds",
                    df_file         = NULL,
                    verbose         = TRUE,
                    ...) {

  layout <- match.arg(layout)
  format <- match.arg(format)

  # ---- small internal helpers ------------------------------------------------
  .dir_create <- function(d) if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  .is_abs     <- function(p) grepl("^(/|[A-Za-z]:[\\/])", p)  # unix or windows
  .msg        <- function(...) if (isTRUE(verbose)) message(...)

  # sanity checks on numeric inputs
  stopifnot(is.numeric(vertex.size), is.numeric(vertex.cex), vertex.size > 0, vertex.cex > 0)
  stopifnot(is.numeric(width), is.numeric(height), width > 0, height > 0)
  stopifnot(is.numeric(dpi), dpi > 0)
  stopifnot(is.logical(show_in_rstudio), is.logical(save_plot), is.logical(save_df), is.logical(show.labels))

  # validate communities
  if (!is.data.frame(communities))
    stop("`communities` must be a data.frame.")
  if (!all(c("actor","layer") %in% names(communities)))
    stop("`communities` must contain columns `actor` and `layer`.")

  # discover community column and normalise to integer `cid`
  comm_col <- intersect(c("cid","com","community"), names(communities))
  if (!length(comm_col))
    stop("`communities` must include one of: `cid`, `com`, or `community`.")
  comm_col <- comm_col[1L]

  cid_raw <- communities[[comm_col]]
  cid_num <- if (is.numeric(cid_raw)) as.integer(cid_raw) else as.integer(factor(cid_raw))

  df <- data.frame(
    actor = as.character(communities$actor),
    layer = as.character(communities$layer),
    cid   = as.integer(cid_num),
    stringsAsFactors = FALSE
  )

  # de-duplicate to one row per (actor, layer)
  if (any(duplicated(df[, c("actor","layer")]))) {
    .msg("Duplicate (actor, layer) rows in `communities`: keeping the first for each pair.")
    df <- df[!duplicated(df[, c("actor","layer")]), , drop = FALSE]
  }

  # ---- network interrogation & layer/actor intersection ----------------------
  Eraw <- try(multinet::edges_ml(net), silent = TRUE)
  if (inherits(Eraw, "try-error") || is.null(Eraw))
    stop("Could not retrieve edges from `net` via multinet::edges_ml()`.")

  Eall <- try(as.data.frame(Eraw, stringsAsFactors = FALSE), silent = TRUE)
  if (inherits(Eall, "try-error") || is.null(Eall) || !nrow(Eall))
    stop("`edges_ml(net)` did not yield a coercible, non-empty edge table.")

  nm <- names(Eall)

  # endpoint detection: common source/target name pairs
  endpoint_pairs <- list(
    c("from_actor","to_actor"), c("from","to"), c("source","target"),
    c("actor1","actor2"), c("i","j"), c("v1","v2")
  )
  a_col <- b_col <- NA_character_
  for (p in endpoint_pairs) if (all(p %in% nm)) { a_col <- p[1]; b_col <- p[2]; break }
  if (is.na(a_col)) {
    # fallback: choose two character/factor columns that are not layer columns
    layer_like <- c("layer","Layer","from_layer","to_layer","l1","l2")
    char_cols  <- which(vapply(Eall, function(x) is.character(x) || is.factor(x), logical(1)))
    char_cols  <- setdiff(char_cols, match(layer_like, nm, nomatch = 0))
    if (length(char_cols) < 2)
      stop("Could not identify two endpoint columns in the edge table.")
    a_col <- nm[char_cols[1]]; b_col <- nm[char_cols[2]]
  }

  # construct a single 'layer' column for intra-layer edges only
  if (all(c("from_layer","to_layer") %in% nm)) {
    Eall <- Eall[Eall$from_layer == Eall$to_layer, , drop = FALSE]
    Eall$layer <- as.character(Eall$from_layer)
  } else if ("layer" %in% nm) {
    Eall$layer <- as.character(Eall$layer)
  } else {
    # if the network has no explicit layer info, treat it as a single layer
    Eall$layer <- "L1"
  }
  if (!nrow(Eall))
    stop("No intra-layer edges found in the network (after filtering).")

  layers_in_net  <- sort(unique(Eall$layer))
  layers_in_comm <- sort(unique(df$layer))

  # choose layers to plot
  if (is.null(layersToPlot)) {
    layersToPlot <- intersect(layers_in_comm, layers_in_net)
  } else {
    layersToPlot <- intersect(as.character(layersToPlot), layers_in_net)
  }
  if (!length(layersToPlot))
    stop("No valid layers to plot after intersecting with the network.")

  # actors present per selected layer (based on endpoints)
  actors_by_layer <- lapply(layersToPlot, function(ly) {
    ed <- Eall[Eall$layer == ly, , drop = FALSE]
    unique(c(as.character(ed[[a_col]]), as.character(ed[[b_col]])))
  })
  names(actors_by_layer) <- layersToPlot

  # keep only (actor, layer) pairs that exist in the network layers selected
  keep_idx <- mapply(function(act, ly) {
    if (!ly %in% names(actors_by_layer)) return(FALSE)
    act %in% actors_by_layer[[ly]]
  }, df$actor, df$layer)

  dropped <- sum(!keep_idx)
  if (dropped > 0) .msg("Dropping ", dropped, " community rows (actors not present in selected layers).")
  df <- df[keep_idx, , drop = FALSE]
  if (!nrow(df))
    stop("No community assignments remain after filtering to actors present in selected layers.")

  # ---- layout computation ----------------------------------------------------
  if (!is.null(seed)) set.seed(as.integer(seed))
  lay_coords <- if (layout == "multiforce") {
    out <- try(multinet::layout_multiforce_ml(net, gravity = gravity), silent = TRUE)
    if (inherits(out, "try-error")) {
      .msg("`layout_multiforce_ml()` failed with gravity = ", gravity, "; falling back to default parameters.")
      multinet::layout_multiforce_ml(net)
    } else out
  } else {
    multinet::layout_circular_ml(net)
  }

  # ---- grid arrangement ------------------------------------------------------
  if (is.null(grid)) {
    k  <- length(layersToPlot)
    nr <- max(1L, round(sqrt(k)))  # near-square
    nc <- ceiling(k / nr)
    grid <- c(as.integer(nr), as.integer(nc))
  } else {
    if (!(is.numeric(grid) && length(grid) == 2L && all(is.finite(grid)) && all(grid >= 1)))
      stop("`grid` must be NULL or a numeric length-2 vector (nrow, ncol) >= 1.")
    grid <- as.integer(grid)
  }

  # ---- labels on/off ---------------------------------------------------------
  vlabels <- if (isTRUE(show.labels)) NULL else NA

  # ---- ensure output directory & stem ---------------------------------------
  .dir_create(results_dir)
  stamp <- format(Sys.time(), "%Y-%m-%d_%H%M%S")

  # function that performs a single draw (used for both preview and export)
  .do_plot <- function() {
    plot(
      net,
      layout        = lay_coords,
      grid          = grid,
      layers        = layersToPlot,
      com           = df,                # three columns: actor, layer, cid
      vertex.size   = vertex.size,
      vertex.cex    = vertex.cex,
      vertex.labels = vlabels,
      ...
    )
  }

  # ---- 1) Interactive preview (if requested) --------------------------------
  if (isTRUE(show_in_rstudio)) {
    .do_plot()
  }

  # ---- 2) Save plot to file (if requested) ----------------------------------
  saved_plot_path <- NULL
  if (isTRUE(save_plot)) {
    if (is.null(file)) {
      stem <- sprintf("plotCom_%s_%dlayers_%s", tolower(layout), length(layersToPlot), stamp)
      file <- paste0(stem, ".", format)
    }
    if (!.is_abs(file)) file <- file.path(results_dir, file)

    if (isTRUE(show_in_rstudio)) {
      # copy the current device content to disk
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
      # draw directly to an off-screen device
      if (format == "png") {
        grDevices::png(filename = file, width = width, height = height, units = units, res = dpi)
        on.exit(grDevices::dev.off(), add = TRUE)
        .do_plot()
      } else {
        grDevices::pdf(file = file, width = width, height = height)  # inches
        on.exit(grDevices::dev.off(), add = TRUE)
        .do_plot()
      }
    }
    saved_plot_path <- normalizePath(file, mustWork = FALSE)
    .msg("Saved ", toupper(format), ": ", saved_plot_path)
  }

  # ---- 3) Save the long (actor, layer, cid) table (if requested) ------------
  df_paths <- character(0)
  if (isTRUE(save_df)) {
    df_formats <- unique(match.arg(df_format, c("rds","csv"), several.ok = TRUE))

    base_name <- if (is.null(df_file)) {
      sprintf("plotCom_communities_%s_%dlayers_%s", tolower(layout), length(layersToPlot), stamp)
    } else {
      tools::file_path_sans_ext(basename(df_file))
    }

    # RDS path: honour df_file exactly (absolute respected; relative into results_dir)
    if ("rds" %in% df_formats) {
      rds_path <- if (is.null(df_file)) {
        file.path(results_dir, paste0(base_name, ".rds"))
      } else if (.is_abs(df_file)) {
        df_file
      } else {
        file.path(results_dir, df_file)
      }
      saveRDS(df, rds_path)
      df_paths <- c(df_paths, normalizePath(rds_path, mustWork = FALSE))
      .msg("Saved community table (RDS): ", df_paths[length(df_paths)])
    }

    # CSV path: always under results_dir, using base name
    if ("csv" %in% df_formats) {
      csv_path <- file.path(results_dir, paste0(base_name, ".csv"))
      utils::write.csv(df, csv_path, row.names = FALSE)
      df_paths <- c(df_paths, normalizePath(csv_path, mustWork = FALSE))
      .msg("Saved community table (CSV): ", df_paths[length(df_paths)])
    }
  }

  # ---- return invisibly with useful attributes -------------------------------
  res <- structure(
    df,
    file    = saved_plot_path,
    df_file = if (length(df_paths)) df_paths else NULL,
    grid    = grid,
    layers  = layersToPlot
  )
  invisible(res)
}


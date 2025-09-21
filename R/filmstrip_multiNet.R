
# ----------------------------------------------------------------------
# 21 — Filmstrip (grid of time‑ordered snapshots) for multi-layer networks
# ----------------------------------------------------------------------

#' Filmstrip of multilayer snapshots with a single, stable layout
#'
#' @description
#' This function renders a *filmstrip* of a multilayer network—one panel per
#' layer—using a **single set of vertex coordinates** computed on the union of
#' all layers. The result is a small multiples display in which spatial position
#' is comparable across layers (conditions, ages, time points), making visual
#' differences easier to interpret.
#'
#' If a community assignment is provided (a data frame with columns `actor`,
#' `layer`, and either `com` or `cid`), vertices in each panel are coloured by
#' their community membership for that layer. Colours are drawn from
#' **RColorBrewer** when available, otherwise a perceptible HCL palette is used.
#'
#' The function can draw into the current device (e.g., the RStudio Plots pane)
#' *and* save the same page to file. Layouts are computed with **igraph**; panels
#' are drawn with **network** so that label and point aesthetics are familiar to
#' many R network users.
#'
#' @details
#' **Stable coordinates.** We build an *undirected* union graph from all edges in
#' the selected layers, compute a single layout once (Kamada–Kawai by default),
#' and reuse those coordinates for every panel. This prevents the common “nodes
#' jumping around” problem across slices.
#'
#' **Label size.** The visible label size is controlled by `label.cex` (passed to
#' `plot.network()` as `label.cex`). A soft‑deprecated alias `cex` is accepted and
#' mapped to `label.cex` to keep older scripts working.
#'
#' **Robust I/O.** The function is defensive about column names in
#' `multinet::edges_ml()` and will detect reasonable synonyms for endpoints and
#' layer columns. If both `from_layer`/`to_layer` are present, we show only
#' **intra‑layer** edges (`from_layer == to_layer`) and store that layer name.
#'
#' @param net A multilayer network usable with `multinet::layers_ml()` and
#'   `multinet::edges_ml()`.
#' @param communities Optional data frame with vertex community assignment per
#'   layer. Must contain columns `actor`, `layer`, and either `com` **or** `cid`
#'   (the latter will be converted to a string community label).
#' @param layer_order Optional character vector specifying the subset and order
#'   of layers to display. Default `NULL` = all layers in the order returned by
#'   `layers_ml(net)`.
#' @param layout One of `c("kamadakawai","mds","force","fr","circle")`. `"force"`
#'   and `"fr"` are convenient aliases to `"kamadakawai"`. Default `"kamadakawai"`.
#' @param actor_normalize Character vector of simple normalisers applied when
#'   matching community `actor` IDs to network vertices. Any subset of
#'   `c("strip_version","trim","tolower","toupper","rm_dash","rm_punct","alnum")`.
#'   Default `c("strip_version","trim","tolower")`.
#' @param ncol,nrow Grid dimensions. If both are `NULL`, a near‑square grid is
#'   chosen automatically.
#' @param width,height Size of the saved figure (in inches). Default `12 × 8`.
#' @param dpi Resolution for PNG output. Ignored for PDF. Default `300`.
#' @param format File format for saving: `"png"` or `"pdf"`. Default `"png"`.
#' @param vertex.cex Numeric multiplier for vertex symbol size (passed to
#'   `plot.network()` as `vertex.cex`). Default `0.9`.
#' @param label.cex Numeric multiplier for **label** size (passed to
#'   `plot.network()` as `label.cex`). Default `0.6`.
#' @param cex Soft‑deprecated alias for `label.cex`. If supplied, it overrides
#'   `label.cex` and a message is emitted. Default `NULL`.
#' @param label.col Single colour for labels. Default `"black"`.
#' @param displaylabels Logical; whether to draw vertex labels. Default `TRUE`.
#' @param edge.col Edge colour (single hex/name) for every panel, with optional
#'   transparency (e.g. `"#55555555"`). Default `"#55555555"`.
#' @param bg Background colour for the page. Default `"white"`.
#' @param slice.par Kept for API symmetry (currently unused). Default is a list
#'   with harmless placeholders.
#' @param seed Optional random seed for reproducible layout coordinates.
#' @param results_dir Output directory for the saved file. Default
#'   `getOption("mlnet.results_dir","omicsDNA_results")`.
#' @param file Optional output filename; if relative, it is placed inside
#'   `results_dir`. If `NULL`, a timestamped name is generated.
#' @param show_in_rstudio Logical; if `TRUE`, draw into the current device (e.g.
#'   RStudio Plots) and then copy the exact page to file. If `FALSE`, draw
#'   off‑screen directly to the file device. Default `TRUE`.
#' @param verbose Logical; print progress messages. Default `TRUE`.
#'
#' @return (Invisibly) a list with:
#' \itemize{
#'   \item `layers` – character vector of layers drawn (in panel order);
#'   \item `coords` – numeric matrix of vertex coordinates (rows named by actor);
#'   \item `file` – absolute path of the saved image (PNG/PDF).
#' }
#'
#' @section Practical notes:
#' - Panels are undirected for positioning (union layout), but your network can
#'   contain directed edges—the drawing simply suppresses arrows for clarity.
#' - If your network has very long labels, consider setting a smaller
#'   `label.cex`, hiding labels (`displaylabels = FALSE`), or exporting a larger
#'   image (`width`/`height`).
#'
#' @examples
#' \dontrun{
#' fs <- filmstrip_multiNet(
#'   net,
#'   communities     = comm,                 # data.frame: actor, layer, com/cid
#'   layer_order     = c("E1","E2","M1","M2"),
#'   layout          = "kamadakawai",
#'   ncol            = 4,
#'   vertex.cex      = 0.8,                  # node size
#'   label.cex       = 0.12,                 # label size (use this, not `cex`)
#'   displaylabels   = TRUE,
#'   edge.col        = "#55555555",
#'   format          = "png",
#'   width           = 12, height = 8, dpi = 300,
#'   seed            = 1,
#'   show_in_rstudio = TRUE,
#'   verbose         = TRUE
#' )
#' fs$file  # path to the saved image
#' }
#'
#' @seealso
#'   \code{\link[multinet]{layers_ml}}, \code{\link[multinet]{edges_ml}},
#'   \code{\link[igraph]{layout_with_kk}}, \code{\link[network]{plot.network}}
#'
#' @importFrom multinet layers_ml edges_ml actors_ml
#' @importFrom igraph graph_from_data_frame layout_with_kk layout_with_mds layout_in_circle
#' @importFrom network network.initialize set.vertex.attribute add.edges network.edgecount
#' @importFrom grDevices png pdf dev.copy dev.off
#' @importFrom graphics par
#' @export
filmstrip_multiNet <- function(
    net,
    communities     = NULL,
    layer_order     = NULL,
    layout          = c("kamadakawai","mds","force","fr","circle"),
    actor_normalize = c("strip_version","trim","tolower"),
    ncol            = NULL,
    nrow            = NULL,
    width           = 12,
    height          = 8,
    dpi             = 300,
    format          = c("png","pdf"),
    vertex.cex      = 0.9,
    label.cex       = 0.6,
    cex             = NULL,   # deprecated alias for label.cex
    label.col       = "black",
    displaylabels   = TRUE,
    edge.col        = "#55555555",
    bg              = "white",
    slice.par       = list(start = 0, end = NULL, interval = 1, aggregate.dur = 1, rule = "any"),
    seed            = NULL,
    results_dir     = getOption("mlnet.results_dir","omicsDNA_results"),
    file            = NULL,
    show_in_rstudio = TRUE,
    verbose         = TRUE
) {
  # ---- small helpers ---------------------------------------------------------
  .need <- function(pkgs) {
    miss <- pkgs[!vapply(pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
    if (length(miss)) stop("Missing packages: ", paste(miss, collapse = ", "),
                           ". install.packages(c(", paste(sprintf('\"%s\"', miss), collapse = ", "), ")).")
  }
  .is_abs     <- function(p) grepl("^(/|[A-Za-z]:[\\/])", p)
  .ensure_dir <- function(d) if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  .norm <- function(x, steps) {
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
                  x)
    }
    x
  }
  .fix_layout <- function(x) {
    v <- c("kamadakawai","mds","force","fr","circle")
    x0 <- tolower(x[1]); if (x0 %in% c("force","fr")) return("kamadakawai")
    if (x0 %in% v) x0 else "kamadakawai"
  }

  # ---- dependencies ----------------------------------------------------------
  .need(c("multinet","network","igraph"))
  has_brewer <- requireNamespace("RColorBrewer", quietly = TRUE)

  # ---- handle deprecated alias ----------------------------------------------
  if (!is.null(cex)) {
    label.cex <- cex
    if (isTRUE(verbose)) message("`cex` is deprecated; use `label.cex` instead.")
  }

  # ---- choose layout / format ------------------------------------------------
  layout  <- .fix_layout(match.arg(layout))
  format  <- match.arg(format)
  if (!is.null(seed)) set.seed(as.integer(seed))

  # ---- layers to show --------------------------------------------------------
  Ls <- try(multinet::layers_ml(net), silent = TRUE)
  if (inherits(Ls, "try-error") || is.null(Ls) || !length(Ls))
    stop("Could not retrieve layers via multinet::layers_ml(net).")
  Ls <- as.character(Ls)
  layers <- if (is.null(layer_order)) Ls else {
    keep <- intersect(as.character(layer_order), Ls)
    if (!length(keep)) stop("No layers to plot after applying `layer_order`.")
    keep
  }
  K <- length(layers)

  # ---- edges (robust extraction) --------------------------------------------
  Eraw <- try(multinet::edges_ml(net), silent = TRUE)
  if (inherits(Eraw, "try-error") || is.null(Eraw))
    stop("Could not retrieve edges via multinet::edges_ml(net).")
  E <- if (is.data.frame(Eraw)) Eraw else {
    tmp <- try(as.data.frame(Eraw, stringsAsFactors = FALSE), silent = TRUE)
    if (inherits(tmp, "try-error")) stop("`edges_ml(net)` is not coercible to data.frame in this multinet version.")
    tmp
  }
  nm <- names(E)
  # endpoint columns (several synonyms)
  pairs <- list(c("from_actor","to_actor"), c("from","to"), c("source","target"),
                c("actor1","actor2"), c("i","j"), c("v1","v2"))
  a_col <- b_col <- NA_character_
  for (p in pairs) if (all(p %in% nm)) { a_col <- p[1]; b_col <- p[2]; break }
  if (is.na(a_col)) {
    # fallback: pick two character columns that are not layer columns
    layer_like <- c("layer","Layer","from_layer","to_layer","l1","l2","layer1","layer2")
    chr <- which(vapply(E, function(x) is.character(x) || is.factor(x), logical(1)))
    chr <- setdiff(chr, match(layer_like, nm, nomatch = 0))
    if (length(chr) < 2) stop("Could not identify two endpoint columns in edges table.")
    a_col <- nm[chr[1]]; b_col <- nm[chr[2]]
  }
  # layer column
  if ("from_layer" %in% nm && "to_layer" %in% nm) {
    E <- E[E$from_layer == E$to_layer, , drop = FALSE]
    E$layer <- as.character(E$from_layer)
  } else if ("layer" %in% nm) {
    E$layer <- as.character(E$layer)
  } else if ("Layer" %in% nm) {
    E$layer <- as.character(E$Layer)
  } else {
    # single-layer fallback if no layer information is present
    E$layer <- Ls[1]
  }
  E <- E[E$layer %in% layers, , drop = FALSE]

  # ---- vertex universe -------------------------------------------------------
  verts <- sort(unique(c(as.character(E[[a_col]]), as.character(E[[b_col]]))))
  if (!length(verts)) {
    A <- try(multinet::actors_ml(net), silent = TRUE)
    if (!inherits(A, "try-error") && !is.null(A)) {
      if (is.data.frame(A)) {
        cand <- intersect(c("actor","name","vertex","node","id"), names(A))
        if (length(cand)) verts <- sort(unique(as.character(A[[cand[1]]])))
      } else if (is.character(A)) verts <- sort(unique(A))
    }
  }
  if (!length(verts)) stop("No vertices found in the selected layers.")

  # ---- build one `network` per layer (consistent vertex order) ---------------
  nets <- vector("list", K); names(nets) <- layers
  for (i in seq_len(K)) {
    ly <- layers[i]
    ed <- E[E$layer == ly, , drop = FALSE]
    gR <- network::network.initialize(n = length(verts), directed = FALSE, loops = FALSE)
    network::set.vertex.attribute(gR, "vertex.names", verts)
    if (nrow(ed)) {
      ti <- match(as.character(ed[[a_col]]), verts)
      tj <- match(as.character(ed[[b_col]]), verts)
      keep <- which(!is.na(ti) & !is.na(tj) & ti != tj)
      if (length(keep)) network::add.edges(gR, tail = ti[keep], head = tj[keep])
    }
    nets[[i]] <- gR
  }

  # ---- compute stable coordinates on the union graph (via igraph) ------------
  el_all <- E[, c(a_col, b_col), drop = FALSE]
  ig <- igraph::graph_from_data_frame(el_all, directed = FALSE, vertices = data.frame(name = verts))
  coords <- switch(layout,
                   "kamadakawai" = igraph::layout_with_kk(ig),
                   "mds"         = igraph::layout_with_mds(ig),
                   "circle"      = igraph::layout_in_circle(ig),
                   igraph::layout_with_kk(ig))
  coords <- as.matrix(coords)
  rownames(coords) <- verts

  # ---- choose grid shape ------------------------------------------------------
  .ensure_dir(results_dir)
  stamp <- format(Sys.time(), "%Y-%m-%d_%H%M%S")
  if (is.null(file) || !nzchar(file)) file <- sprintf("multiNet_filmstrip_%dlayers_%s.%s", K, stamp, format)
  if (!.is_abs(file)) file <- file.path(results_dir, file)

  if (is.null(ncol) && is.null(nrow)) {
    nr <- max(1L, floor(sqrt(K))); nc <- ceiling(K / nr)
  } else if (is.null(ncol)) {
    nr <- as.integer(nrow); nc <- ceiling(K / nr)
  } else if (is.null(nrow)) {
    nc <- as.integer(ncol); nr <- ceiling(K / nc)
  } else { nr <- as.integer(nrow); nc <- as.integer(ncol) }

  # ---- panel drawing function -------------------------------------------------
  draw_panels <- function() {
    op <- par(no.readonly = TRUE); on.exit(par(op), add = TRUE)
    par(mfrow = c(nr, nc), mar = c(0.6, 0.6, 1.6, 0.6), xaxs = "i", yaxs = "i", bg = bg)

    # palette for communities (if provided)
    col_map <- NULL
    if (!is.null(communities)) {
      dfc <- communities
      if (!("com" %in% names(dfc)) && ("cid" %in% names(dfc))) dfc$com <- paste0("C", dfc$cid)
      if (!all(c("actor","layer","com") %in% names(dfc)))
        stop("`communities` must contain columns: actor, layer, and com (or cid).")
      uniq_com <- sort(unique(as.character(dfc$com)))
      col_map <- if (has_brewer) {
        pal <- RColorBrewer::brewer.pal(max(3, min(12, length(uniq_com))), "Set3")
        stats::setNames(rep(pal, length.out = length(uniq_com)), uniq_com)
      } else {
        hues <- grDevices::hcl(h = seq(0, 360, length.out = length(uniq_com) + 1L)[-1L], c = 60, l = 65)
        stats::setNames(hues, uniq_com)
      }
    }

    for (i in seq_len(K)) {
      gR <- nets[[i]]
      vnames <- network::get.vertex.attribute(gR, "vertex.names")

      # vertex colours: default grey; overwrite by community if provided
      vcols <- rep("grey85", length(vnames))
      if (!is.null(communities)) {
        dfly <- communities[communities$layer == layers[i], , drop = FALSE]
        if (nrow(dfly)) {
          idx <- match(.norm(as.character(dfly$actor), actor_normalize),
                       .norm(vnames, actor_normalize))
          ok  <- which(!is.na(idx))
          if (length(ok)) vcols[idx[ok]] <- col_map[as.character(dfly$com[ok])]
        }
      }

      ne <- network::network.edgecount(gR)
      ecols <- if (ne > 0L) rep(as.character(edge.col), ne) else NULL

      # draw one panel
      plot(gR,
           coord          = coords[vnames, , drop = FALSE],
           displaylabels  = displaylabels,
           vertex.cex     = vertex.cex,
           label.cex      = label.cex,
           label.col      = label.col,
           vertex.col     = vcols,
           edge.col       = ecols,
           usearrows      = FALSE,
           main           = layers[i])
    }
  }

  # ---- draw and save ----------------------------------------------------------
  if (isTRUE(show_in_rstudio)) {
    # 1) draw to current device (e.g., RStudio Plots)
    draw_panels()
    # 2) copy the exact page to file
    if (identical(format, "png")) {
      dev_id <- grDevices::dev.copy(grDevices::png, filename = file,
                                    width = width, height = height, units = "in", res = dpi, bg = bg)
      grDevices::dev.off(dev_id)
    } else {
      dev_id <- grDevices::dev.copy(grDevices::pdf, file = file,
                                    width = width, height = height, onefile = TRUE, paper = "special")
      grDevices::dev.off(dev_id)
    }
  } else {
    # off‑screen directly to file
    if (identical(format, "png")) {
      grDevices::png(file, width = width, height = height, units = "in", res = dpi, bg = bg)
      on.exit(grDevices::dev.off(), add = TRUE)
    } else {
      grDevices::pdf(file, width = width, height = height, onefile = TRUE, paper = "special", bg = bg)
      on.exit(grDevices::dev.off(), add = TRUE)
    }
    draw_panels()
  }

  if (isTRUE(verbose)) {
    message("Saved filmstrip: ", normalizePath(file, winslash = "/", mustWork = FALSE))
    message("Frames: ", K, " | Grid: ", nr, "×", nc, " | Layout: ", layout)
  }

  invisible(list(layers = layers, coords = coords, file = file))
}


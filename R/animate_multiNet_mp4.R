
# -------------------------------------------------
# 25 - Animate multi-layer networks to GIF/MP4
# -------------------------------------------------

#' Animate multilayer networks to GIF/MP4 (one layer = one frame)
#'
#' @description
#' Builds a consistent 2D layout from the union of all layers (so node positions
#' are stable), then renders a frame for each layer and animates transitions.
#' If a communities table is provided (`actor`, `layer`, and `com` or `cid`),
#' nodes are coloured by community per layer; otherwise they are drawn uniformly.
#'
#' @param net A multilayer object usable with multinet::layers_ml() and edges_ml().
#' @param communities Optional data.frame with columns `actor`, `layer`, and either
#'   `com` or `cid` (if only `cid`, labels `com = paste0("C", cid)` are generated).
#' @param layer_order Optional character vector to specify the order of layers in the animation.
#'   Default: order returned by layers_ml(net).
#' @param layout One of c("kamadakawai","mds","circle"). Default "kamadakawai".
#' @param actor_normalize Character vector of normalization steps used to match
#'   community actor IDs to network actors. Default c("strip_version","trim","tolower").
#' @param format Output format: "gif" (default) or "mp4".
#' @param fps Frames per second (default 12).
#' @param frames_per_layer How many frames to allocate to each layer state (default 12).
#'   Higher values slow the animation and give longer dwell time.
#' @param width,height Output size in *pixels* (default 1000 x 800).
#' @param edge.col Edge colour (hex or name). Default "#555555".
#' @param edge.alpha Edge transparency in [0,1] (default 0.5).
#' @param point.size Node point size (default 2.5).
#' @param show.labels Logical; draw node labels. Default FALSE.
#' @param results_dir Output directory (default getOption("mlnet.results_dir","omicsDNA_results")).
#' @param file Optional filename. If NULL, an auto-stamped name is used under results_dir.
#'             If you pass a relative path, it will be placed under results_dir.
#' @param seed Optional integer seed for reproducible layouts.
#' @param verbose Logical; print progress.
#'
#' @return Invisibly returns a list with:
#'   - `file`: absolute path to the rendered GIF/MP4,
#'   - `layers`: the layer order used,
#'   - `anim`: the gganimate object (so you can re-render with different parameters).
#'
#' @examples
#' \dontrun{
#' gif <- animate_multiNet_mp4(net, communities = comm, format = "gif")
#' mp4 <- animate_multiNet_mp4(net, communities = comm, format = "mp4", fps = 15)
#' }
animate_multiNet_mp4 <- function(
    net,
    communities     = NULL,
    layer_order     = NULL,
    layout          = c("kamadakawai","mds","circle"),
    actor_normalize = c("strip_version","trim","tolower"),
    format          = c("gif","mp4"),
    fps             = 12,
    frames_per_layer= 12,
    width           = 1000,
    height          = 800,
    edge.col        = "#555555",
    edge.alpha      = 0.5,
    point.size      = 2.5,
    show.labels     = FALSE,
    results_dir     = getOption("mlnet.results_dir","omicsDNA_results"),
    file            = NULL,
    seed            = NULL,
    verbose         = TRUE
) {
  ## ---- utilities
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
                  x
      )
    }
    x
  }
  .fix_layout <- function(x) {
    x0 <- tolower(x[1])
    if (x0 %in% c("kamadakawai","mds","circle")) x0 else "kamadakawai"
  }

  ## ---- deps
  .need(c("multinet","igraph","ggplot2","gganimate"))
  has_brewer <- requireNamespace("RColorBrewer", quietly = TRUE)
  has_repel  <- requireNamespace("ggrepel", quietly = TRUE)
  format <- match.arg(format)
  layout_arg <- .fix_layout(layout)

  ## ---- layers
  Ls <- try(multinet::layers_ml(net), silent = TRUE)
  if (inherits(Ls, "try-error") || is.null(Ls) || !length(Ls))
    stop("Could not retrieve layers via multinet::layers_ml(net).")
  Ls <- as.character(Ls)
  layers <- if (is.null(layer_order)) Ls else {
    keep <- intersect(as.character(layer_order), Ls)
    if (!length(keep)) stop("No layers to animate after applying `layer_order`.")
    keep
  }
  K <- length(layers)

  ## ---- edges (robust)
  Eraw <- try(multinet::edges_ml(net), silent = TRUE)
  if (inherits(Eraw, "try-error") || is.null(Eraw))
    stop("Could not retrieve edges via multinet::edges_ml(net).")
  E <- if (is.data.frame(Eraw)) Eraw else {
    tmp <- try(as.data.frame(Eraw, stringsAsFactors = FALSE), silent = TRUE)
    if (inherits(tmp, "try-error") || !nrow(tmp)) stop("edges_ml(net) not coercible to data.frame.")
    tmp
  }
  nm <- names(E)
  pairs <- list(c("from_actor","to_actor"), c("from","to"), c("source","target"),
                c("actor1","actor2"), c("i","j"), c("v1","v2"))
  a_col <- b_col <- NA_character_
  for (p in pairs) if (all(p %in% nm)) { a_col <- p[1]; b_col <- p[2]; break }
  if (is.na(a_col)) {
    layer_like <- c("layer","Layer","from_layer","to_layer","l1","l2")
    chr <- which(vapply(E, function(x) is.character(x) || is.factor(x), logical(1)))
    chr <- setdiff(chr, match(layer_like, nm, nomatch = 0))
    if (length(chr) < 2) stop("Could not identify two endpoint columns in edges table.")
    a_col <- nm[chr[1]]; b_col <- nm[chr[2]]
  }
  if ("from_layer" %in% nm && "to_layer" %in% nm) {
    E <- E[E$from_layer == E$to_layer, , drop = FALSE]; E$layer <- as.character(E$from_layer)
  } else if ("layer" %in% nm) {
    E$layer <- as.character(E$layer)
  } else if ("Layer" %in% nm) {
    E$layer <- as.character(E$Layer)
  } else {
    E$layer <- "L1"
  }
  E <- E[E$layer %in% layers, , drop = FALSE]

  ## ---- vertex universe & union layout
  verts <- sort(unique(c(as.character(E[[a_col]]), as.character(E[[b_col]]))))
  if (!length(verts)) {
    A <- try(multinet::actors_ml(net), silent = TRUE)
    if (!inherits(A, "try-error") && !is.null(A)) {
      if (is.data.frame(A)) {
        cand <- intersect(c("actor","name","vertex","node","id"), names(A))
        if (length(cand)) verts <- sort(unique(as.character(A[[cand[1]]])))
      } else if (is.character(A)) {
        verts <- sort(unique(A))
      }
    }
  }
  if (!length(verts)) stop("No vertices found in the selected layers.")

  if (!is.null(seed)) set.seed(as.integer(seed))
  ig <- igraph::graph_from_data_frame(E[, c(a_col, b_col)], directed = FALSE,
                                      vertices = data.frame(name = verts))
  coords <- switch(layout_arg,
                   "kamadakawai" = igraph::layout_with_kk(ig),
                   "mds"         = igraph::layout_with_mds(ig),
                   "circle"      = igraph::layout_in_circle(ig)
  )
  coords <- as.data.frame(coords); names(coords) <- c("x","y")
  coords$actor <- verts

  ## ---- per-layer nodes/edges data frames for ggplot/gganimate
  # nodes present per layer
  actors_by_layer <- lapply(layers, function(ly) {
    ed <- E[E$layer == ly, , drop = FALSE]
    unique(c(as.character(ed[[a_col]]), as.character(ed[[b_col]])))
  })
  names(actors_by_layer) <- layers

  nodes <- do.call(rbind, lapply(layers, function(ly) {
    a <- actors_by_layer[[ly]]
    if (length(a)) merge(data.frame(actor = a, layer = ly, stringsAsFactors = FALSE),
                         coords, by = "actor", all.x = TRUE)
  }))
  if (is.null(nodes) || !nrow(nodes)) stop("No nodes found to plot.")

  # colour by community per layer (if provided)
  if (!is.null(communities)) {
    if (!("com" %in% names(communities)) && ("cid" %in% names(communities)))
      communities$com <- paste0("C", communities$cid)
    if (!all(c("actor","layer","com") %in% names(communities)))
      stop("`communities` must contain columns: actor, layer, and com (or cid).")

    uniq_com <- sort(unique(as.character(communities$com)))
    col_map <- if (has_brewer) {
      pal <- RColorBrewer::brewer.pal(max(3, min(12, length(uniq_com))), "Set3")
      stats::setNames(rep(pal, length.out = length(uniq_com)), uniq_com)
    } else {
      hues <- grDevices::hcl(h = seq(0, 360, length.out = length(uniq_com) + 1L)[-1L], c = 60, l = 65)
      stats::setNames(hues, uniq_com)
    }

    # join per layer
    key <- paste(.norm(communities$actor, actor_normalize), communities$layer, sep = "\r")
    com_vec <- setNames(as.character(communities$com), key)
    nodes$key <- paste(.norm(nodes$actor, actor_normalize), nodes$layer, sep = "\r")
    nodes$col <- col_map[com_vec[nodes$key]]
    nodes$col[is.na(nodes$col)] <- "grey70"
  } else {
    nodes$col <- "steelblue"
  }

  # edges per layer as segments (two points per edge)
  edges <- do.call(rbind, lapply(layers, function(ly) {
    ed <- E[E$layer == ly, c(a_col,b_col), drop = FALSE]
    if (!nrow(ed)) return(NULL)
    ed$layer <- ly
    ed
  }))
  if (!is.null(edges) && nrow(edges)) {
    edges <- merge(edges, coords, by.x = a_col, by.y = "actor", all.x = TRUE)
    names(edges)[names(edges) %in% c("x","y")] <- c("x0","y0")
    edges <- merge(edges, coords, by.x = b_col, by.y = "actor", all.x = TRUE)
    names(edges)[names(edges) %in% c("x","y")] <- c("x1","y1")
    edges <- edges[is.finite(edges$x0) & is.finite(edges$y0) &
                     is.finite(edges$x1) & is.finite(edges$y1), , drop = FALSE]
  }

  ## ---- ggplot + gganimate
  suppressPackageStartupMessages(requireNamespace("ggplot2"))
  p <- ggplot2::ggplot() +
    {
      if (!is.null(edges) && nrow(edges)) {
        ggplot2::geom_segment(
          data = edges,
          ggplot2::aes(x = x0, y = y0, xend = x1, yend = y1),
          colour = edge.col, alpha = edge.alpha, linewidth = 0.3
        )
      }
    } +
    ggplot2::geom_point(
      data = nodes,
      ggplot2::aes(x = x, y = y, fill = col),
      shape = 21, colour = NA, size = point.size
    ) +
    {
      if (isTRUE(show.labels)) {
        if (has_repel) {
          ggrepel::geom_text_repel(
            data = nodes, ggplot2::aes(x = x, y = y, label = actor),
            size = 2.6, max.overlaps = 50, min.segment.length = 0
          )
        } else {
          ggplot2::geom_text(
            data = nodes, ggplot2::aes(x = x, y = y, label = actor),
            size = 2.6, alpha = 0.8
          )
        }
      }
    } +
    ggplot2::scale_fill_identity() +
    ggplot2::coord_equal() +
    ggplot2::theme_void(base_size = 12) +
    ggplot2::labs(title = "Layer: {closest_state}") +
    gganimate::transition_states(
      states = nodes$layer,
      state_length = 1, transition_length = 1, wrap = FALSE
    ) +
    gganimate::ease_aes("cubic-in-out")

  nframes <- max(1L, K * as.integer(frames_per_layer))

  .ensure_dir(results_dir)
  stamp <- format(Sys.time(), "%Y-%m-%d_%H%M%S")
  if (is.null(file) || !nzchar(file)) {
    file <- sprintf("multiNet_anim_%dlayers_%s.%s", K, stamp, format)
  }
  if (!.is_abs(file)) file <- file.path(results_dir, file)

  renderer <- if (identical(format, "gif")) {
    if (!requireNamespace("gifski", quietly = TRUE))
      stop("`format='gif'` requires the {gifski} package. install.packages('gifski')")
    gganimate::gifski_renderer(file)
  } else {
    if (!requireNamespace("av", quietly = TRUE))
      stop("`format='mp4'` requires the {av} package. install.packages('av')")
    gganimate::av_renderer(file)
  }

  if (isTRUE(verbose)) message("Rendering ", toupper(format), " → ", normalizePath(file, winslash = "/", mustWork = FALSE))
  anim <- gganimate::animate(
    p, nframes = nframes, fps = fps, width = width, height = height,
    renderer = renderer, bg = "white"
  )

  invisible(list(file = normalizePath(file, winslash = "/", mustWork = FALSE),
                 layers = layers, anim = anim))
}


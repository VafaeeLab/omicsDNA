#' Filmstrip of multilayer snapshots (shows in RStudio + saves to file):
#'
#' Builds one `network` per layer from your multinet object, computes a single
#' stable layout from the union graph (so positions are consistent across panels),
#' and draws a grid of panels (one per layer). If `communities` is provided
#' (`actor`, `layer`, and `com` or `cid`), vertices are coloured by community
#' per slice. When `show_in_rstudio = TRUE`, the filmstrip is drawn to the
#' current device (RStudio Plots pane), and then copied to PNG/PDF.
#'
#' @param net A multilayer network usable with multinet::layers_ml() and edges_ml().
#' @param communities Optional data.frame with `actor`, `layer`, and `com` or `cid`.
#' @param layer_order Optional explicit order of layers; default = order from layers_ml(net).
#' @param layout One of c("kamadakawai","mds","force","fr","circle"); "force"/"fr" map to "kamadakawai".
#' @param actor_normalize Normalization steps for matching community actors to network actors.
#' @param ncol,nrow Grid dimensions; if both NULL a near-square grid is chosen.
#' @param width,height Figure size (inches) for the saved file. @param dpi PNG resolution.
#' @param format "png" or "pdf" (for saving).
#' @param vertex.cex Node size multiplier. @param displaylabels Show vertex labels?
#' @param edge.col Edge colour (single hex/name) for every panel.
#' @param bg Background colour. @param slice.par Kept for API symmetry (unused here).
#' @param seed Optional RNG seed for the layout. @param results_dir Output directory for the file.
#' @param file Optional filename; if relative it is placed in `results_dir`. Auto‑named if NULL.
#' @param show_in_rstudio Logical; draw to the Plots pane before saving. Default TRUE.
#' @param verbose Print progress.
#' @return (Invisibly) list with `layers`, `coords` (matrix), `file` (abs path).
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
  ## ---- helpers
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

  ## ---- deps
  .need(c("multinet","network","igraph"))
  has_brewer <- requireNamespace("RColorBrewer", quietly = TRUE)

  ## ---- layers
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

  ## ---- vertex universe
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

  ## ---- one network per layer (consistent vertex order)
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

  ## ---- stable layout from the union graph (via igraph)
  if (!is.null(seed)) set.seed(as.integer(seed))
  el_all <- E[, c(a_col, b_col)]
  ig <- igraph::graph_from_data_frame(el_all, directed = FALSE, vertices = data.frame(name = verts))
  layout_arg <- .fix_layout(layout)
  coords <- switch(layout_arg,
                   "kamadakawai" = igraph::layout_with_kk(ig),
                   "mds"         = igraph::layout_with_mds(ig),
                   "circle"      = igraph::layout_in_circle(ig),
                   igraph::layout_with_kk(ig)
  )
  coords <- as.matrix(coords)
  rownames(coords) <- verts

  ## ---- grid shape
  .ensure_dir(results_dir)
  format  <- match.arg(format)
  stamp   <- format(Sys.time(), "%Y-%m-%d_%H%M%S")
  if (is.null(file) || !nzchar(file)) file <- sprintf("multiNet_filmstrip_%dlayers_%s.%s", K, stamp, format)
  if (!.is_abs(file)) file <- file.path(results_dir, file)

  if (is.null(ncol) && is.null(nrow)) {
    nr <- max(1L, floor(sqrt(K))); nc <- ceiling(K / nr)
  } else if (is.null(ncol)) {
    nr <- as.integer(nrow); nc <- ceiling(K / nr)
  } else if (is.null(nrow)) {
    nc <- as.integer(ncol); nr <- ceiling(K / nc)
  } else { nr <- as.integer(nrow); nc <- as.integer(ncol) }

  draw_panels <- function() {
    op <- par(no.readonly = TRUE); on.exit(par(op), add = TRUE)
    par(mfrow = c(nr, nc), mar = c(0.6, 0.6, 1.5, 0.6), xaxs = "i", yaxs = "i", bg = bg)
    ## community palette (if any)
    col_map <- NULL
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
    }
    for (i in seq_len(K)) {
      gR <- nets[[i]]
      vnames <- network::get.vertex.attribute(gR, "vertex.names")
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
      plot(gR,
           coord          = coords[vnames, , drop = FALSE],
           displaylabels  = displaylabels,
           vertex.cex     = vertex.cex,
           vertex.col     = vcols,
           edge.col       = ecols,
           usearrows      = FALSE,
           main           = layers[i])
    }
  }

  if (isTRUE(show_in_rstudio)) {
    # 1) draw to the current device (RStudio Plots)
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
    # Off-screen draw directly to file device (previous behaviour)
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
    message("Frames: ", K, " | Grid: ", nr, "×", nc, " | Layout: ", layout_arg)
  }

  invisible(list(layers = layers, coords = coords, file = file))
}

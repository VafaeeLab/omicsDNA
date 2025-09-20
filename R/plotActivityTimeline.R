
# ----------------------------------------------------------------
# 22 - Edge/Vertex Activity Timelines (Gantt-style) for multilayer networks
# ----------------------------------------------------------------

#' Edge/Vertex Activity Timelines (Gantt-style) for multilayer networks
#'
#' @description
#' Build a temporal network from a multilayer \pkg{multinet} object and plot
#' activity “spells” as a timeline using \pkg{ndtv}. Each horizontal bar
#' represents either a single **edge** (co-expression relationship) or a single
#' **vertex** (gene), with bar length covering the age-group layers in which the
#' element is present. Contiguous presence across adjacent layers is merged into
#' a single spell so that “stable” elements appear as long bars.
#'
#' @details
#' The function:
#' 1) Reads all intra-layer edges from `net` and normalizes to a single `layer`.
#' 2) Maps layers to discrete time points (t = 0,1,2,...) using `slice.par`.
#' 3) Computes **contiguous spells** of presence for either:
#'    - `type = "edge"`: undirected pairs (u,v) present in one or more layers.
#'    - `type = "vertex"`: genes that appear in at least one edge in a layer.
#' 4) Constructs a base `network` with the full vertex set, wraps it as a
#'    `networkDynamic`, and **activates** the relevant edges/vertices over
#'    the computed time spells.
#' 5) Calls `ndtv::timeline()` to draw a Gantt‑style plot.
#'
#' **Biological interpretation**
#' - **Stable edges** (long bars) can suggest core, age‑independent regulatory
#'   relationships.
#' - **Transient edges** (short bars) may reflect age‑specific responses.
#' - **Late‑appearing vertices/edges** indicate genes/relationships activated
#'   only at older age groups.
#'
#' @param net A multilayer network compatible with `multinet::edges_ml()` and
#'   `multinet::layers_ml()`.
#' @param type `"edge"` (default) or `"vertex"`. What to plot as timeline rows.
#' @param layer_order Optional character vector to fix the chronological order
#'   of layers (e.g., `c("E1","E2","M1",...)`). Default: all layers in their
#'   natural order from `layers_ml(net)`.
#' @param slice.par List controlling the layer→time mapping:
#'   `list(start=0, interval=1, aggregate.dur=1)`. Each layer *i* occupies
#'   `[start + (i-1)*interval, start + (i-1)*interval + aggregate.dur)`.
#' @param top_n Optional integer; plot only the `top_n` longest‑duration edges
#'   (for `type="edge"`) or vertices (`type="vertex"`). Default `NULL` = all.
#' @param min_duration Minimum total duration (in time units) to keep an item.
#'   Default `1` (keeps elements present in at least one layer).
#' @param edge.col Colour for edge bars (e.g., `"#2c7fb8"`). Used when
#'   `type="edge"`.
#' @param vertex.col Colour for vertex bars (e.g., `"#d95f02"`). Used when
#'   `type="vertex"`.
#' @param width,height,dpi Figure size if saving to PNG (inches; DPI). For PDF,
#'   `dpi` is ignored.
#' @param format `"png"` or `"pdf"` for saving. Default `"png"`.
#' @param results_dir Output directory (default:
#'   `getOption("mlnet.results_dir","omicsDNA_results")`).
#' @param file Optional output filename. If `NULL`, an informative name is made.
#'   If relative, it is written under `results_dir`.
#' @param show_in_rstudio Draw to the current device (RStudio Plots pane).
#'   Default `TRUE`.
#' @param save_plot Also save to file. Default `TRUE`.
#' @param seed Optional integer seed used only for reproducible ordering ties.
#' @param verbose Print progress messages. Default `TRUE`.
#'
#' @return Invisibly returns a list:
#'   - `nd`: the `networkDynamic` object (with activated spells),
#'   - `file`: absolute path to the saved plot (if written),
#'   - `layers`: layer order used,
#'   - `spells`: a data.frame of computed spells (element, onset, terminus).
#'
#' @examples
#' \dontrun{
#' # Edges timeline (default)
#' plotActivityTimeline(net, type="edge", top_n=300)
#'
#' # Vertex timeline, only longest 100 genes
#' plotActivityTimeline(net, type="vertex", top_n=100, vertex.col="tomato")
#' }
plotActivityTimeline <- function(
    net,
    type            = c("edge","vertex"),
    layer_order     = NULL,
    slice.par       = list(start = 0, interval = 1, aggregate.dur = 1),
    top_n           = NULL,
    min_duration    = 1,
    edge.col        = "#2c7fb8",
    vertex.col      = "#d95f02",
    width           = 12, height = 8, dpi = 300,
    format          = c("png","pdf"),
    results_dir     = getOption("mlnet.results_dir","omicsDNA_results"),
    file            = NULL,
    show_in_rstudio = TRUE,
    save_plot       = TRUE,
    seed            = NULL,
    verbose         = TRUE
) {
  ## --- helpers
  .need <- function(pkgs) {
    miss <- pkgs[!vapply(pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
    if (length(miss)) stop("Missing packages: ", paste(miss, collapse = ", "),
                           ". install.packages(c(", paste(sprintf('\"%s\"', miss), collapse=", "), ")).")
  }
  .is_abs     <- function(p) grepl("^(/|[A-Za-z]:[\\/])", p)
  .ensure_dir <- function(d) if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  .canon_pair <- function(u, v) {
    u <- as.character(u); v <- as.character(v)
    paste(pmin(u, v), pmax(u, v), sep = "\t")
  }
  .runs_to_spells <- function(idx, onsets, dur) {
    # idx: integer layer indices (sorted unique)
    if (!length(idx)) return(list())
    idx <- sort(unique(idx))
    brk <- c(1L, which(diff(idx) != 1L) + 1L)
    res <- vector("list", length(brk))
    for (k in seq_along(brk)) {
      start_pos <- brk[k]
      end_pos   <- if (k < length(brk)) brk[k+1] - 1L else length(idx)
      i1 <- idx[start_pos]
      i2 <- idx[end_pos]
      onset   <- onsets[i1]
      terminus<- onsets[i2] + dur
      res[[k]] <- c(onset = onset, terminus = terminus)
    }
    res
  }

  ## --- deps
  .need(c("multinet","network","networkDynamic","ndtv"))
  type   <- match.arg(type)
  format <- match.arg(format)
  sp     <- modifyList(list(start = 0, interval = 1, aggregate.dur = 1), slice.par)

  if (!is.null(seed)) set.seed(as.integer(seed))

  ## --- layers
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
  onsets  <- sp$start + (seq_len(K) - 1L) * sp$interval
  dur     <- sp$aggregate.dur

  ## --- global edge table (robust)
  Eraw <- try(multinet::edges_ml(net), silent = TRUE)
  if (inherits(Eraw, "try-error") || is.null(Eraw))
    stop("Could not retrieve edges via multinet::edges_ml(net).")
  E <- if (is.data.frame(Eraw)) Eraw else {
    tmp <- try(as.data.frame(Eraw, stringsAsFactors = FALSE), silent = TRUE)
    if (inherits(tmp, "try-error") || !nrow(tmp))
      stop("`edges_ml(net)` not coercible to data.frame.")
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
    E <- E[E$from_layer == E$to_layer, , drop = FALSE]
    E$layer <- as.character(E$from_layer)
  } else if ("layer" %in% nm) {
    E$layer <- as.character(E$layer)
  } else if ("Layer" %in% nm) {
    E$layer <- as.character(E$Layer)
  } else {
    E$layer <- "L1"
  }
  # keep only selected layers
  E <- E[E$layer %in% layers, , drop = FALSE]
  if (!nrow(E)) stop("No intra-layer edges found in the selected layers.")

  ## --- vertex universe
  verts <- sort(unique(c(as.character(E[[a_col]]), as.character(E[[b_col]]))))
  if (!length(verts)) stop("No vertices found in the selected layers.")
  nv <- length(verts)
  idmap <- setNames(seq_len(nv), verts)

  ## --- base network + wrap as dynamic
  g0 <- network::network.initialize(n = nv, directed = FALSE, loops = FALSE)
  network::set.vertex.attribute(g0, "vertex.names", verts)
  nd <- networkDynamic::networkDynamic(base.net = g0)

  ## --- compute spells & activate
  spells_df <- NULL

  if (type == "edge") {
    # Build index: for each canonical pair, which layer indices contain it?
    idx_map <- new.env(parent = emptyenv())
    for (i in seq_len(K)) {
      ly <- layers[i]
      ed <- E[E$layer == ly, , drop = FALSE]
      if (!nrow(ed)) next
      u <- as.character(ed[[a_col]]); v <- as.character(ed[[b_col]])
      key <- .canon_pair(u, v)
      if (length(key)) {
        tab <- table(key) # unique within layer
        for (k in names(tab)) {
          vv <- if (exists(k, envir = idx_map, inherits = FALSE)) get(k, envir = idx_map) else integer(0)
          assign(k, c(vv, i), envir = idx_map)
        }
      }
    }
    all_pairs <- ls(idx_map)
    if (!length(all_pairs)) stop("No edges to plot.")
    # Preorder pairs by total duration (desc)
    total_dur <- vapply(all_pairs, function(k) length(get(k, envir = idx_map)) * dur, numeric(1))
    ord <- order(total_dur, decreasing = TRUE)
    all_pairs <- all_pairs[ord]
    total_dur <- total_dur[ord]

    # Optional filters
    keep <- rep(TRUE, length(all_pairs))
    if (!is.null(top_n)) keep <- keep & seq_along(all_pairs) <= as.integer(top_n)
    if (!is.null(min_duration)) keep <- keep & total_dur >= as.numeric(min_duration)
    all_pairs <- all_pairs[keep]; total_dur <- total_dur[keep]
    if (!length(all_pairs)) stop("No edges remain after filtering by `top_n`/`min_duration`.")

    # Add one edge per unique pair, then activate multiple spells on that edge id
    spells_list <- vector("list", length(all_pairs))
    for (j in seq_along(all_pairs)) {
      k <- all_pairs[j]
      parts <- strsplit(k, "\t", fixed = TRUE)[[1]]
      tail_id <- idmap[parts[1]]; head_id <- idmap[parts[2]]

      # create the edge, remember its id
      prev <- network::network.edgecount(nd)
      network::add.edges(nd, tail = tail_id, head = head_id)
      e_id <- prev + 1L

      idx <- get(k, envir = idx_map)
      rs <- .runs_to_spells(idx, onsets, dur)
      if (length(rs)) {
        for (r in rs) {
          networkDynamic::activate.edges(nd, e = e_id, onset = r["onset"], terminus = r["terminus"])
        }
        # collect spells for return
        spells_list[[j]] <- data.frame(
          element  = k,  # "u\tv"
          onset    = vapply(rs, `[[`, numeric(1), "onset"),
          terminus = vapply(rs, `[[`, numeric(1), "terminus"),
          stringsAsFactors = FALSE
        )
      }
    }
    spells_df <- do.call(rbind, Filter(Negate(is.null), spells_list))

  } else {  # type == "vertex"
    # Build index: for each actor, which layer indices contain it?
    idx_map <- new.env(parent = emptyenv())
    for (i in seq_len(K)) {
      ly <- layers[i]
      ed <- E[E$layer == ly, , drop = FALSE]
      if (!nrow(ed)) next
      acts <- unique(c(as.character(ed[[a_col]]), as.character(ed[[b_col]])))
      for (a in acts) {
        vv <- if (exists(a, envir = idx_map, inherits = FALSE)) get(a, envir = idx_map) else integer(0)
        assign(a, c(vv, i), envir = idx_map)
      }
    }
    actors <- ls(idx_map)
    if (!length(actors)) stop("No vertices to plot (none appear in any selected layer).")

    total_dur <- vapply(actors, function(a) length(get(a, envir = idx_map)) * dur, numeric(1))
    ord <- order(total_dur, decreasing = TRUE)
    actors <- actors[ord]; total_dur <- total_dur[ord]
    keep <- rep(TRUE, length(actors))
    if (!is.null(top_n)) keep <- keep & seq_along(actors) <= as.integer(top_n)
    if (!is.null(min_duration)) keep <- keep & total_dur >= as.numeric(min_duration)
    actors <- actors[keep]; total_dur <- total_dur[keep]
    if (!length(actors)) stop("No vertices remain after filtering by `top_n`/`min_duration`.")

    # Activate vertex spells
    spells_list <- vector("list", length(actors))
    for (j in seq_along(actors)) {
      a <- actors[j]
      v_id <- idmap[a]
      idx <- get(a, envir = idx_map)
      rs <- .runs_to_spells(idx, onsets, dur)
      if (length(rs)) {
        # Note: we *add* spells; vertices are otherwise inactive
        for (r in rs) {
          networkDynamic::activate.vertices(nd, v = v_id, onset = r["onset"], terminus = r["terminus"])
        }
        spells_list[[j]] <- data.frame(
          element  = a,
          onset    = vapply(rs, `[[`, numeric(1), "onset"),
          terminus = vapply(rs, `[[`, numeric(1), "terminus"),
          stringsAsFactors = FALSE
        )
      }
    }
    spells_df <- do.call(rbind, Filter(Negate(is.null), spells_list))
  }

  ## --- draw (to Plots pane) and/or save
  .ensure_dir(results_dir)
  stamp <- format(Sys.time(), "%Y-%m-%d_%H%M%S")
  if (is.null(file) || !nzchar(file)) {
    base <- sprintf("timeline_%s_%dlayers_%s", type, K, stamp)
    file <- paste0(base, ".", format)
  }
  if (!.is_abs(file)) file <- file.path(results_dir, file)

  # 1) Show in RStudio device (Plots pane)
  if (isTRUE(show_in_rstudio)) {
    if (type == "edge") {
      ndtv::timeline(nd, plot.edge.spells = TRUE,  plot.vertex.spells = FALSE,
                     edge.col = edge.col, main = sprintf("Edge activity timeline (%d layers)", K))
    } else {
      ndtv::timeline(nd, plot.edge.spells = FALSE, plot.vertex.spells = TRUE,
                     vertex.col = vertex.col, main = sprintf("Vertex activity timeline (%d layers)", K))
    }
  }

  # 2) Save to file
  out_path <- NULL
  if (isTRUE(save_plot)) {
    if (identical(format, "png")) {
      grDevices::dev.copy(grDevices::png, filename = file, width = width, height = height,
                          units = "in", res = dpi, bg = "white")
      grDevices::dev.off()
    } else {
      grDevices::dev.copy(grDevices::pdf, file = file, width = width, height = height, onefile = TRUE)
      grDevices::dev.off()
    }
    out_path <- normalizePath(file, winslash = "/", mustWork = FALSE)
    if (isTRUE(verbose)) message("Saved timeline: ", out_path)
  }

  invisible(list(nd = nd, file = out_path, layers = layers, spells = spells_df))
}

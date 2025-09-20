
# -------------------------------------------------------------
# 23 — Plot differences between two layers on a shared layout
#        (customisable figure + CSV with per-layer weights)
# -------------------------------------------------------------

#' Plot differences between two multilayer slices with a shared layout
#' and export a tidy CSV that includes the correct per-layer edge weights.
#'
#' @description
#' Compares two layers \code{L1} and \code{L2} from a \pkg{multinet} object by
#' computing a layout on the \emph{union} graph and colouring edges by category:
#' edges present only in \code{L1}, only in \code{L2}, or in both (common).
#' A tidy CSV is written with columns \code{from}, \code{to}, \code{category},
#' \code{weight_L1}, \code{weight_L2}. The weights are pulled from the per-layer
#' igraphs stored inside the multilayer object (via \code{as.list(net)}), so
#' they match exactly what you built with \code{build_multiNet()}.
#'
#' @details
#' - If a weight attribute exists on layer edges (default name \code{"weight"}),
#'   those values are used. You can specify alternatives via \code{weight_attr}.
#' - Duplicate edges in a layer (if any) are aggregated with \code{weight_aggregate}.
#' - For undirected comparisons (\code{directed=FALSE}), edge keys use
#'   \code{min(from,to)}–\code{max(from,to)} canonicalization; for directed,
#'   keys are \code{"from->to"}.
#' - The CSV \code{category} values are \code{only_<L1>}, \code{common},
#'   \code{only_<L2>} to match your requested naming.
#'
#' @param net multinet object (usable with \code{as.list(net)} and
#'   \code{multinet::layers_ml()}).
#' @param L1,L2 Character layer identifiers to compare.
#' @param layout igraph layout for the union graph: one of
#'   \code{"fr"}, \code{"kk"}, \code{"mds"}, \code{"lgl"}, \code{"graphopt"}.
#' @param directed Logical; if \code{TRUE}, treat edges as directed. Default \code{FALSE}.
#' @param col_only_L1,col_only_L2,col_common Edge colours for categories.
#' @param edge_alpha  Edge transparency [0,1].
#' @param edge_width_only,edge_width_common Numeric edge widths for only/common.
#' @param vertex_size,vertex_col,vertex_frame_col Vertex style.
#' @param vertex_label FALSE/TRUE/character vector for labels; label_* control fonts.
#' @param show_legend,legend_pos,legend_cex Legend style.
#' @param title Optional plot title; if NULL a summary is used. @param bg background.
#' @param seed Optional integer seed for layout.
#' @param results_dir,file,format,width,height,dpi Output controls for figure.
#' @param save_edge_csv,csv_file Whether/where to write the edge CSV.
#' @param weight_attr Candidate edge attribute names to use as weights.
#' @param weight_aggregate How to combine duplicates (\code{"sum"}, \code{"mean"},
#'   \code{"median"}, \code{"first"}).
#'
#' @return Invisibly: list(g_union, g1, g2, counts, edges, file, csv_file)
#' @export
plot_layer_diff <- function(
    net, L1, L2,
    layout = c("fr","kk","mds","lgl","graphopt"),
    directed = FALSE,
    col_only_L1 = "#e74c3c",
    col_only_L2 = "#1f9d55",
    col_common  = "#95a5a6",
    edge_alpha  = 0.7,
    edge_width_only   = 1.8,
    edge_width_common = 1.2,
    vertex_size      = 4,
    vertex_col       = "grey85",
    vertex_frame_col = NA,
    vertex_label     = FALSE,
    label_cex        = 0.6,
    label_col        = "black",
    label_family     = NULL,
    label_font       = 1,
    show_legend      = TRUE,
    legend_pos       = "topleft",
    legend_cex       = 0.9,
    title            = NULL,
    bg               = "white",
    seed             = NULL,
    results_dir      = getOption("mlnet.results_dir","omicsDNA_results"),
    file             = NULL,
    format           = c("png","pdf"),
    width            = 9, height = 7, dpi = 300,
    save_edge_csv    = TRUE,
    csv_file         = NULL,
    weight_attr      = c("weight","w","value"),
    weight_aggregate = c("sum","mean","median","first")
) {
  layout <- match.arg(layout)
  format <- match.arg(format)
  weight_aggregate <- match.arg(weight_aggregate)

  .is_abs <- function(p) grepl("^(/|[A-Za-z]:[\\/])", p)
  .ensure_dir <- function(d) if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

  # ---- per-layer igraphs ----
  glist <- try(as.list(net), silent = TRUE)
  if (inherits(glist, "try-error") || !is.list(glist))
    stop("Could not coerce `net` to list of igraphs via as.list(net).")
  nm_g <- names(glist)
  if (!is.null(nm_g) && any(nm_g == "_flat_")) glist <- glist[nm_g != "_flat_"]
  if (!all(c(L1, L2) %in% names(glist)))
    stop("Requested layers not found in as.list(net): ",
         paste(setdiff(c(L1,L2), names(glist)), collapse = ", "))

  g1 <- glist[[L1]]
  g2 <- glist[[L2]]
  stopifnot(inherits(g1,"igraph"), inherits(g2,"igraph"))

  # ---- key + weight map helpers ----
  key_fun <- if (!directed) {
    function(el) paste(pmin(el[,1], el[,2]), pmax(el[,1], el[,2]), sep = "\t")
  } else {
    function(el) paste(el[,1], el[,2], sep = "->")
  }
  pick_weight_attr <- function(g, cands) {
    ea <- igraph::edge_attr_names(g)
    for (cand in cands) {
      if (cand %in% ea) return(cand)
      if (tolower(cand) %in% ea) return(tolower(cand))
      if (toupper(cand) %in% ea) return(toupper(cand))
    }
    NULL
  }
  agg_fun <- switch(weight_aggregate,
                    sum    = function(z) sum(z, na.rm = TRUE),
                    mean   = function(z) mean(z, na.rm = TRUE),
                    median = function(z) stats::median(z, na.rm = TRUE),
                    first  = function(z) z[1])

  layer_weight_map <- function(g) {
    if (igraph::ecount(g) == 0) return(setNames(numeric(0), character(0)))
    el <- igraph::as_edgelist(g, names = TRUE)
    kk <- key_fun(el)
    nm <- pick_weight_attr(g, weight_attr)
    ww <- if (!is.null(nm)) igraph::edge_attr(g, nm) else rep(1, nrow(el))
    ww <- suppressWarnings(as.numeric(ww))
    stats::aggregate(ww, by = list(key = kk), FUN = agg_fun) |>
      (\(df) setNames(df$x, df$key))()
  }

  w1_map <- layer_weight_map(g1)
  w2_map <- layer_weight_map(g2)

  # ---- union graph + shared layout ----
  gU <- igraph::simplify(igraph::union(g1, g2))
  if (igraph::vcount(gU) == 0) stop("Both layers have no edges.")
  if (!is.null(seed)) set.seed(as.integer(seed))
  Lcoords <- switch(layout,
                    fr      = igraph::layout_with_fr(gU),
                    kk      = igraph::layout_with_kk(gU),
                    mds     = igraph::layout_with_mds(gU),
                    lgl     = igraph::layout_with_lgl(gU),
                    graphopt= igraph::layout_with_graphopt(gU)
  )
  rownames(Lcoords) <- igraph::V(gU)$name

  # keys & categories
  elU <- igraph::as_edgelist(gU, names = TRUE)
  kU  <- key_fun(elU)
  k1  <- names(w1_map)
  k2  <- names(w2_map)
  only1_keys  <- setdiff(k1, k2)
  only2_keys  <- setdiff(k2, k1)
  common_keys <- intersect(k1, k2)

  # edge aesthetics
  edge_col <- rep(col_common, length(kU))
  edge_col[kU %in% only1_keys] <- col_only_L1
  edge_col[kU %in% only2_keys] <- col_only_L2
  edge_col <- grDevices::adjustcolor(edge_col, alpha.f = edge_alpha)

  edge_w <- rep(edge_width_common, length(kU))
  edge_w[(kU %in% only1_keys) | (kU %in% only2_keys)] <- edge_width_only

  # vertex labels
  vnames <- igraph::V(gU)$name
  vlab <- NA
  if (isTRUE(vertex_label)) {
    vlab <- vnames
  } else if (is.character(vertex_label)) {
    if (!is.null(names(vertex_label))) {
      idx <- match(vnames, names(vertex_label))
      vlab <- vertex_label[idx]
    } else if (length(vertex_label) == length(vnames)) {
      vlab <- vertex_label
    } else {
      warning("`vertex_label` length does not match vertices; labels suppressed.")
      vlab <- NA
    }
  }

  # ---- draw ----
  draw <- function() {
    op <- par(mar=c(1.2,1.2,3.2,1.2), xaxs="i", yaxs="i", bg=bg); on.exit(par(op), add=TRUE)
    main_title <- if (!is.null(title) && nzchar(title)) title else
      sprintf("Edges: %s vs %s  |  only %s: %d   only %s: %d   common: %d",
              L1, L2, L1, length(only1_keys), L2, length(only2_keys), length(common_keys))
    plot(gU,
         layout             = Lcoords,
         edge.color         = edge_col,
         edge.width         = edge_w,
         vertex.size        = vertex_size,
         vertex.color       = vertex_col,
         vertex.frame.color = vertex_frame_col,
         vertex.label       = vlab,
         vertex.label.cex   = label_cex,
         vertex.label.color = label_col,
         vertex.label.family= label_family,
         vertex.label.font  = label_font,
         main               = main_title
    )
    if (isTRUE(show_legend)) {
      legend(legend_pos,
             legend = c(sprintf("only %s", L1), "common", sprintf("only %s", L2)),
             col    = c(col_only_L1, col_common, col_only_L2),
             lwd    = 2, bty = "n", cex = legend_cex)
    }
  }

  fig_file <- NULL
  if (!is.null(file) && nzchar(file)) {
    if (!.is_abs(file)) { .ensure_dir(results_dir); file <- file.path(results_dir, file) }
    if (identical(format, "png")) {
      grDevices::png(file, width = width, height = height, units = "in", res = dpi, bg = bg)
      on.exit(grDevices::dev.off(), add = TRUE); draw()
    } else {
      grDevices::pdf(file, width = width, height = height, onefile = FALSE, paper = "special", bg = bg)
      on.exit(grDevices::dev.off(), add = TRUE); draw()
    }
    fig_file <- normalizePath(file, winslash = "/", mustWork = FALSE)
    message("Saved diff figure: ", fig_file)
  } else {
    draw()
  }

  # ---- CSV with per-layer weights (correct) ----
  cat_vec <- ifelse(kU %in% common_keys, "common",
                    ifelse(kU %in% only1_keys, paste0("only_", L1), paste0("only_", L2)))
  edges_out <- data.frame(
    from      = if (!directed) pmin(elU[,1], elU[,2]) else elU[,1],
    to        = if (!directed) pmax(elU[,1], elU[,2]) else elU[,2],
    category  = cat_vec,
    weight_L1 = as.numeric(unname(w1_map[kU])),
    weight_L2 = as.numeric(unname(w2_map[kU])),
    stringsAsFactors = FALSE
  )

  csv_out <- NULL
  if (isTRUE(save_edge_csv)) {
    .ensure_dir(results_dir)
    stamp <- format(Sys.time(), "%Y-%m-%d_%H%M%S")
    if (is.null(csv_file) || !nzchar(csv_file)) {
      csv_file <- file.path(results_dir, sprintf("edges_diff_%s_vs_%s_%s.csv", L1, L2, stamp))
    } else if (!.is_abs(csv_file)) {
      csv_file <- file.path(results_dir, csv_file)
    }
    utils::write.csv(edges_out, csv_file, row.names = FALSE)
    csv_out <- normalizePath(csv_file, winslash = "/", mustWork = FALSE)
    message("Saved edge list CSV: ", csv_out)
  }

  counts_vec <- setNames(
    c(length(only1_keys), length(only2_keys), length(common_keys)),
    c(paste0("only_", L1), paste0("only_", L2), "common")
  )

  invisible(list(
    g_union  = gU, g1 = g1, g2 = g2,
    counts   = counts_vec,
    edges    = edges_out,
    file     = fig_file,
    csv_file = csv_out
  ))
}


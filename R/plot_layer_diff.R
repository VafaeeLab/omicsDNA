#' Plot differences between two layers on a shared layout:
#'
#' Shows edges present only in layer A (red), only in layer B (green),
#' and edges common to both (grey). Positions are computed on the union graph.
#'
#' @param net   multinet object (usable with multinet::edges_ml()).
#' @param L1,L2 Character layer ids to compare.
#' @param layout Name of igraph layout to use on the union ("fr","kk","mds","lgl","graphopt").
#' @param directed Logical; build directed graphs? Default FALSE.
#' @param edge_alpha Numeric in [0,1]; edge transparency.
#' @param vertex_cex Vertex point size.
#' @param results_dir Where to save when file is relative.
#' @param file Optional filename to save (PNG). If NULL, just draws to the current device.
#' @param width,height,dpi PNG device size (inches) and resolution when saving.
#' @return Invisibly returns a list with igraphs and counts.
plot_layer_diff <- function(net, L1, L2,
                            layout = c("fr","kk","mds","lgl","graphopt"),
                            directed = FALSE,
                            edge_alpha = 0.7,
                            vertex_cex = 0.6,
                            results_dir = getOption("mlnet.results_dir","omicsDNA_results"),
                            file = NULL,
                            width = 9, height = 7, dpi = 300) {
  layout <- match.arg(layout)
  .is_abs <- function(p) grepl("^(/|[A-Za-z]:[\\/])", p)
  .ensure_dir <- function(d) if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

  # ---- Pull edges once and normalise layer column (same heuristics you used) ----
  Eraw <- try(multinet::edges_ml(net), silent = TRUE)
  if (inherits(Eraw, "try-error") || is.null(Eraw)) stop("Could not read edges via multinet::edges_ml().")
  E <- if (is.data.frame(Eraw)) Eraw else {
    tmp <- try(as.data.frame(Eraw, stringsAsFactors = FALSE), silent = TRUE)
    if (inherits(tmp, "try-error") || !nrow(tmp)) stop("edges table not coercible.")
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
    if (length(chr) < 2) stop("Could not find two endpoint columns in edges.")
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

  L1 <- as.character(L1); L2 <- as.character(L2)
  if (!all(c(L1,L2) %in% unique(E$layer))) stop("Layers not found in the data.")

  # ---- Build undirected graphs for each layer ----
  mk_g <- function(ed) {
    if (!nrow(ed)) igraph::make_empty_graph(directed = directed)
    else igraph::graph_from_data_frame(
      d = data.frame(from = as.character(ed[[a_col]]),
                     to   = as.character(ed[[b_col]]),
                     stringsAsFactors = FALSE),
      directed = directed
    )
  }
  g1 <- mk_g(E[E$layer == L1, , drop = FALSE])
  g2 <- mk_g(E[E$layer == L2, , drop = FALSE])

  # ---- Keys for undirected edges ----
  key <- function(g) if (igraph::ecount(g)==0) character(0) else {
    a <- igraph::as_edgelist(g, names = TRUE); paste(pmin(a[,1],a[,2]), pmax(a[,1],a[,2]), sep="\t")
  }
  k1 <- key(g1); k2 <- key(g2)
  only1   <- setdiff(k1, k2)
  only2   <- setdiff(k2, k1)
  common  <- intersect(k1, k2)

  # ---- Union graph + consistent layout ----
  gU <- igraph::simplify(igraph::union(g1, g2))
  if (igraph::vcount(gU) == 0) stop("Both layers have no edges.")
  L <- switch(layout,
              fr = igraph::layout_with_fr(gU),
              kk = igraph::layout_with_kk(gU),
              mds = igraph::layout_with_mds(gU),
              lgl = igraph::layout_with_lgl(gU),
              graphopt = igraph::layout_with_graphopt(gU))
  rownames(L) <- igraph::V(gU)$name

  # Map categories onto union edges
  kU <- key(gU)
  colU <- rep("#95a5a6", length(kU))              # common (grey)
  colU[kU %in% only1] <- "#e74c3c"                # only in L1 (red)
  colU[kU %in% only2] <- "#2ecc71"                # only in L2 (green)
  colU <- grDevices::adjustcolor(colU, alpha.f = edge_alpha)

  # ---- Draw (to device; optionally also save) ----
  draw <- function() {
    op <- par(mar=c(1.2,1.2,3.2,1.2))
    on.exit(par(op), add=TRUE)
    plot(gU, layout = L,
         edge.color = colU,
         edge.width = 1.2,
         vertex.size = 3.5 * vertex_cex,
         vertex.label = NA,
         main = sprintf("Edges: %s vs %s  |  only %s: %d   only %s: %d   common: %d",
                        L1, L2, L1, length(only1), L2, length(only2), length(common)))
    legend("topleft",
           legend = c(sprintf("only %s", L1), "common", sprintf("only %s", L2)),
           col    = c("#e74c3c", "#95a5a6", "#2ecc71"),
           lwd = 2, bty = "n", cex = 0.9)
  }

  if (!is.null(file) && nzchar(file)) {
    if (!.is_abs(file)) { .ensure_dir(results_dir); file <- file.path(results_dir, file) }
    grDevices::png(file, width = width, height = height, units = "in", res = dpi, bg = "white")
    on.exit(grDevices::dev.off(), add = TRUE)
    draw()
    message("Saved diff figure: ", normalizePath(file, winslash = "/", mustWork = FALSE))
  } else {
    draw()  # show in RStudio Plots pane
  }

  invisible(list(g_union = gU, g1 = g1, g2 = g2,
                 counts = c(only_L1 = length(only1), only_L2 = length(only2), common = length(common))))
}

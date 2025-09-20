

# --------------------------------------------------------
# 20 - Animate multi-layer networks
# --------------------------------------------------------

#' Animate a multilayer network across layers using \pkg{ndtv}
#'
#' @description
#' \code{animate_multiNet()} converts a multilayer network (built with
#' \pkg{multinet}) into a time‑aware \code{networkDynamic} object and renders an
#' interactive HTML animation with \pkg{ndtv}. Each **layer** is treated as a
#' discrete **time slice**; nodes persist across slices (isolates are kept) and
#' edges are taken from the corresponding layer only. Optionally, the function
#' colours vertices by community membership per slice (e.g., the output of
#' \code{detectCom()}), so you can watch communities evolve across age groups
#' or conditions.
#'
#' The conversion is robust to minor schema differences in \code{multinet} edge
#' tables (it auto‑detects endpoint and layer columns), and it writes a
#' timestamped HTML file under your project results folder for easy sharing.
#'
#' @details
#' **Workflow / What the function does**
#'
#' 1. **Layer order.** Retrieves the available layers via
#'    \code{multinet::layers_ml(net)} and optionally reorders/subsets them using
#'    \code{layer_order}. If not supplied, all layers are used in their original
#'    order.
#' 2. **Edge normalization.** Pulls a single global edge table via
#'    \code{multinet::edges_ml(net)} and harmonizes columns to one
#'    \code{layer} field while keeping **intra‑layer** edges only
#'    (if \code{from_layer}/\code{to_layer} exist, they are filtered to
#'    \code{from_layer == to_layer}).
#' 3. **Vertex universe.** Creates the union of all actors observed in the
#'    selected layers (falling back to \code{multinet::actors_ml(net)} if
#'    necessary) and constructs one static \pkg{network} object per layer with
#'    that common vertex set (so isolates are preserved).
#' 4. **Time mapping.** Wraps the list of per‑layer networks into a single
#'    \code{networkDynamic} object, assigning onsets/termini for each slice
#'    (layer) based on \code{slice.par}. With \code{K} layers, the default
#'    slices are \eqn{[0,1)}, \eqn{[1,2)}, …, \eqn{[K-1,K)}.
#' 5. **Optional colours by community.** If a \code{communities} data frame is
#'    provided (columns \code{actor}, \code{layer}, and \code{com} or \code{cid}),
#'    the function assigns a per‑slice dynamic vertex attribute
#'    \code{"vertex.col"} so that \code{render.d3movie()} automatically colours
#'    nodes by community on each layer.
#' 6. **Layout and rendering.** Computes an animation layout with
#'    \code{ndtv::compute.animation()} (default Kamada–Kawai) and writes an
#'    interactive HTML animation via \code{ndtv::render.d3movie()}.
#'
#' **Why this design?**
#'
#' - Treating layers as time slices is a natural way to visualize temporal or
#'   age‑stratified multilayer graphs without requiring continuous timestamps.
#' - Keeping isolates ensures that vertices can reappear after disappearing in
#'   intermediate slices, making trajectories easy to follow.
#' - Storing colours in the dynamic attribute \code{"vertex.col"} exploits
#'   \pkg{ndtv}'s native mechanism for per‑slice styling.
#'
#' **Common interpretation**
#'
#' - Clusters of nodes that remain close across slices indicate communities that
#'   persist across conditions/ages. Sudden movements or fragmentations suggest
#'   community reorganization. Colour stability (when \code{communities} are
#'   supplied) helps track membership changes.
#'
#' **Dependencies**
#'
#' Requires \pkg{ndtv}, \pkg{networkDynamic}, and \pkg{network}; \pkg{RColorBrewer}
#' is optional (used for palettes; otherwise a perceptually uniform HCL palette
#' is generated).
#'
#' @param net A multilayer network object compatible with
#'   \code{multinet::layers_ml()}, \code{multinet::edges_ml()}, and (optionally)
#'   \code{multinet::actors_ml()}.
#' @param communities Optional \code{data.frame} with per‑slice community labels
#'   used for vertex colours. Must contain columns:
#'   \itemize{
#'     \item \code{actor}: vertex identifiers;
#'     \item \code{layer}: layer name matching those in \code{net};
#'     \item \code{com} or \code{cid}: community label. If only \code{cid} is
#'           present, it is converted to \code{com = paste0("C", cid)}.
#'   }
#'   Actor IDs are matched robustly using \code{actor_normalize}.
#' @param layer_order Optional character vector of layer names to use and their
#'   order in the animation. Defaults to all layers in the order returned by
#'   \code{layers_ml(net)}.
#' @param layout Character; animation layout mode passed to
#'   \code{ndtv::compute.animation()}. One of \code{"kamadakawai"}, \code{"mds"},
#'   \code{"force"}, \code{"fr"}, or \code{"circle"}. The synonyms \code{"force"}
#'   and \code{"fr"} are mapped internally to \code{"kamadakawai"}. Minor typos
#'   are tolerated and corrected when possible.
#' @param slice.par A named list controlling time slicing (mirrors \pkg{ndtv}):
#'   \describe{
#'     \item{\code{start}}{Numeric start time for the first slice (default 0).}
#'     \item{\code{end}}{Optional numeric end time. If \code{NULL} (default), it
#'       is set to \code{start + (K - 1) * interval}, where \code{K} is the
#'       number of layers.}
#'     \item{\code{interval}}{Width between successive slice onsets (default 1).}
#'     \item{\code{aggregate.dur}}{Duration of each slice (default 1).}
#'     \item{\code{rule}}{Aggregation rule for edges within a slice (default
#'       \code{"any"}).}
#'   }
#' @param directed Logical; whether to build directed \pkg{network} objects for
#'   each slice (default \code{FALSE}). Note that most visual narratives benefit
#'   from undirected layouts even if the underlying graph is directed.
#' @param actor_normalize Character vector of normalization steps used to match
#'   \code{communities$actor} to network vertex names. Defaults to
#'   \code{c("strip_version","trim","tolower")}. Supported steps are:
#'   \code{"strip_version"}, \code{"trim"}, \code{"tolower"}, \code{"toupper"},
#'   \code{"rm_dash"}, \code{"rm_punct"}, \code{"alnum"}.
#' @param results_dir Output directory for the HTML file. Defaults to
#'   \code{getOption("mlnet.results_dir","omicsDNA_results")}.
#' @param html_file Optional filename for the HTML animation. If relative, it is
#'   created under \code{results_dir}. If \code{NULL}, a timestamped name is
#'   generated (e.g., \code{"multiNet_animation_YYYY-mm-dd_HHMMSS.html"}).
#' @param vertex.cex Numeric vertex size passed to \code{render.d3movie()}.
#'   Default \code{0.9}.
#' @param displaylabels Logical; show vertex labels in the animation
#'   (\code{TRUE} by default). For large graphs, consider \code{FALSE}.
#' @param edge.col Edge colour (can include alpha), e.g. \code{"#55555555"}.
#' @param bg Background colour for the movie device (default \code{"white"}).
#' @param launchBrowser Logical; open the resulting HTML file in the browser
#'   after writing (default \code{TRUE}).
#' @param seed Optional integer seed for reproducible layouts.
#' @param verbose Logical; print progress messages (default \code{TRUE}).
#'
#' @return
#' Invisibly returns a list with:
#' \itemize{
#'   \item \code{nd}: the \code{networkDynamic} object containing all slices;
#'   \item \code{layers}: the character vector of layers used (in order);
#'   \item \code{html_file}: the absolute path to the saved animation;
#'   \item \code{anim}: the result of \code{ndtv::compute.animation()}.
#' }
#' The function also writes an interactive HTML file to \code{results_dir}.
#'
#' @section Interpretation guide:
#' \itemize{
#'   \item \emph{Stable positions and colours} across adjacent slices suggest
#'         persistent communities or roles.
#'   \item \emph{Nodes changing colour} (with \code{communities} supplied)
#'         indicate community reassignments; cohesive colour blocks splitting or
#'         merging reflect structural transitions.
#'   \item \emph{Isolates appearing/disappearing} signal entry/exit of actors in
#'         specific layers rather than total removal from the study.
#' }
#'
#' @section Troubleshooting:
#' \itemize{
#'   \item \strong{“invalid color name 'ndtv_col'”} — occurs when colours are not
#'         stored under the dynamic attribute \code{"vertex.col"} that
#'         \pkg{ndtv} expects. This function writes to \code{"vertex.col"}
#'         internally to avoid that error.
#'   \item \strong{No layers found} — ensure \code{net} responds to
#'         \code{multinet::layers_ml()} and that \code{layer_order}, if supplied,
#'         matches existing layer names exactly (after normalization).
#' }
#'
#' @examples
#' \dontrun{
#' ## Minimal animation (no colouring), all layers in default order
#' anim <- animate_multiNet(net)
#'
#' ## With per-slice colours from community assignments (detectCom output)
#' anim <- animate_multiNet(
#'   net,
#'   communities = comm,                   # has actor, layer, com (or cid)
#'   layout      = "kamadakawai",
#'   slice.par   = list(start = 0, interval = 1, aggregate.dur = 1),
#'   seed        = 1
#' )
#'
#' ## Restrict to a subset / custom order of layers and save with a custom name
#' anim <- animate_multiNet(
#'   net,
#'   layer_order = c("E1","E2","M1","M2"),
#'   html_file   = "my_movie.html",
#'   displaylabels = FALSE
#' )
#' }
#'
#' @seealso
#' \code{\link[ndtv]{compute.animation}}, \code{\link[ndtv]{render.d3movie}},
#' \code{\link[networkDynamic]{networkDynamic}}, and the pipeline functions
#' \code{detectCom()} (community detection) and \code{plotCom()} (static plots).
#'
#' @importFrom ndtv compute.animation render.d3movie
#' @importFrom network network.initialize set.vertex.attribute add.edges
#' @importFrom networkDynamic networkDynamic
#' @export
animate_multiNet <- function(
    net,
    communities     = NULL,
    layer_order     = NULL,
    layout          = c("kamadakawai","mds","force","fr","circle"),
    slice.par       = list(start = 0, end = NULL, interval = 1, aggregate.dur = 1, rule = "any"),
    directed        = FALSE,
    actor_normalize = c("strip_version","trim","tolower"),
    results_dir     = getOption("mlnet.results_dir","omicsDNA_results"),
    html_file       = NULL,
    vertex.cex      = 0.9,
    displaylabels   = TRUE,
    edge.col        = "#55555555",
    bg              = "white",
    launchBrowser   = TRUE,
    seed            = NULL,
    verbose         = TRUE
) {
  ## --- tiny utilities ---
  .need <- function(pkgs) {
    miss <- pkgs[!vapply(pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
    if (length(miss)) stop("Missing packages: ", paste(miss, collapse = ", "),
                           ". Install with install.packages(c(", paste(sprintf('\"%s\"', miss), collapse = ", "), ")).")
  }
  .is_abs <- function(p) grepl("^(/|[A-Za-z]:[\\/])", p)
  .ensure_dir <- function(d) if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  .normalize_keys <- function(x, steps) {
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
  .fix_layout <- function(x, verbose=FALSE) {
    valid <- c("kamadakawai","mds","force","fr","circle")
    x0 <- tolower(x[1])
    if (x0 %in% valid) return(if (x0 %in% c("force","fr")) "kamadakawai" else x0)
    cand <- agrep(x0, valid, max.distance = 1, value = TRUE)
    if (length(cand)) {
      if (verbose) message("Interpreting layout='", x[1], "' as '", cand[1], "'.")
      return(if (cand[1] %in% c("force","fr")) "kamadakawai" else cand[1])
    }
    if (verbose) message("Unknown layout '", x[1], "'. Using 'kamadakawai'.")
    "kamadakawai"
  }

  ## --- deps ---
  .need(c("ndtv","networkDynamic","network"))
  has_brewer <- requireNamespace("RColorBrewer", quietly = TRUE)

  ## --- layers in order ---
  Ls <- try(multinet::layers_ml(net), silent = TRUE)
  if (inherits(Ls, "try-error") || is.null(Ls) || !length(Ls))
    stop("Could not retrieve layers via multinet::layers_ml(net).")
  Ls <- as.character(Ls)
  layers <- if (is.null(layer_order)) Ls else intersect(as.character(layer_order), Ls)
  if (!length(layers)) stop("No layers to animate after applying `layer_order`.")
  K <- length(layers)

  ## --- robust edge table ---
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
    if (length(chr) < 2) stop("Could not identify two endpoint columns.")
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
  E <- E[E$layer %in% layers, , drop = FALSE]

  ## --- vertex universe ---
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

  ## --- build one network per slice (keep isolates) ---
  nets <- vector("list", K); names(nets) <- layers
  for (i in seq_len(K)) {
    ly <- layers[i]
    ed <- E[E$layer == ly, , drop = FALSE]
    gR <- network::network.initialize(n = length(verts), directed = directed, loops = FALSE)
    network::set.vertex.attribute(gR, "vertex.names", verts)
    if (nrow(ed)) {
      ti <- match(as.character(ed[[a_col]]), verts)
      tj <- match(as.character(ed[[b_col]]), verts)
      keep <- which(!is.na(ti) & !is.na(tj) & ti != tj)
      if (length(keep)) network::add.edges(gR, tail = ti[keep], head = tj[keep])
    }
    nets[[i]] <- gR
  }

  ## --- time mapping -> networkDynamic (exactly K slices) ---
  sp <- utils::modifyList(list(start=0, end=NULL, interval=1, aggregate.dur=1, rule="any"), slice.par)
  onsets  <- sp$start + (seq_len(K) - 1L) * sp$interval
  termini <- onsets + sp$aggregate.dur
  end_for_slices <- sp$start + (K - 1L) * sp$interval

  nd <- networkDynamic::networkDynamic(
    network.list = nets,
    base.net     = nets[[1]],
    onsets       = onsets,
    termini      = termini,
    vertex.pid   = "vertex.names"
  )

  ## --- optional: per-slice colours from communities ---
  if (!is.null(communities)) {
    if (!is.data.frame(communities) || !all(c("actor","layer") %in% names(communities)))
      stop("`communities` must have columns `actor` and `layer`, plus `com` or `cid`.")
    if (!("com" %in% names(communities)) && ("cid" %in% names(communities)))
      communities$com <- paste0("C", communities$cid)
    if (!("com" %in% names(communities)))
      stop("`communities` must contain either `com` or `cid`.")

    uniq_com <- sort(unique(as.character(communities$com)))
    if (has_brewer) {
      pal <- RColorBrewer::brewer.pal(max(3, min(12, length(uniq_com))), "Set3")
      col_map <- setNames(rep(pal, length.out = length(uniq_com)), uniq_com)
    } else {
      hues <- grDevices::hcl(h = seq(0, 360, length.out = length(uniq_com) + 1L)[-1L], c = 60, l = 65)
      col_map <- setNames(hues, uniq_com)
    }

    v_keys <- .normalize_keys(verts, actor_normalize)
    use_activate <- exists("activate.vertex.attribute", where = asNamespace("networkDynamic"), inherits = FALSE)
    setter <- if (use_activate)
      get("activate.vertex.attribute", asNamespace("networkDynamic"))
    else
      get("set.vertex.attribute.active", asNamespace("networkDynamic"))

    for (i in seq_len(K)) {
      ly   <- layers[i]
      dfly <- communities[communities$layer == ly, , drop = FALSE]
      vcols <- rep("grey85", length(verts))
      if (nrow(dfly)) {
        a_keys <- .normalize_keys(as.character(dfly$actor), actor_normalize)
        idx    <- match(a_keys, v_keys)
        ok     <- which(!is.na(idx))
        if (length(ok)) vcols[idx[ok]] <- col_map[as.character(dfly$com[ok])]
      }
      ## Store colours in dynamic attribute **named "vertex.col"** (what ndtv expects)
      do.call(setter, list(x = nd, prefix = "vertex.col", value = unname(vcols),
                           onset = onsets[i], terminus = termini[i], v = seq_along(verts)))
    }
  }

  ## --- layout & render ---
  if (!is.null(seed)) set.seed(as.integer(seed))
  layout_arg <- .fix_layout(layout, verbose = verbose)

  anim <- ndtv::compute.animation(
    nd,
    animation.mode = layout_arg,
    slice.par      = list(start = sp$start, end = end_for_slices, interval = sp$interval,
                          aggregate.dur = sp$aggregate.dur, rule = sp$rule),
    verbose        = verbose
  )

  .ensure_dir(results_dir)
  stamp <- format(Sys.time(), "%Y-%m-%d_%H%M%S")
  if (is.null(html_file) || !nzchar(html_file)) {
    html_file <- file.path(results_dir, paste0("multiNet_animation_", stamp, ".html"))
  } else if (!.is_abs(html_file)) {
    html_file <- file.path(results_dir, html_file)
  }

  ## Do NOT pass vertex.col=; ndtv will read dynamic vertex attribute "vertex.col"
  ndtv::render.d3movie(
    nd,
    filename      = html_file,
    usearrows     = FALSE,
    displaylabels = displaylabels,
    vertex.cex    = vertex.cex,
    edge.col      = edge.col,
    launchBrowser = launchBrowser,
    plot.par      = list(bg = bg),
    slice.par     = list(start = sp$start, end = end_for_slices, interval = sp$interval,
                         aggregate.dur = sp$aggregate.dur, rule = sp$rule)
  )
  if (isTRUE(verbose)) message("Saved animation: ", normalizePath(html_file, winslash = "/", mustWork = FALSE))

  invisible(list(nd = nd, layers = layers, html_file = html_file, anim = anim))
}

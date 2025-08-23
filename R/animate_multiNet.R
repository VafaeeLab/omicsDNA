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

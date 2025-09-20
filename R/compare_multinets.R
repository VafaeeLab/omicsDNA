
# ------------------------------------------------------------------------------
# 8 - Robust comparison of two ml.network objects via per-layer igraphs
# ------------------------------------------------------------------------------

#' Compare two `ml.network` objects and **print the differences**
#'
#' Uses multinet's `as.list()` to obtain one `igraph` per layer and then checks
#' that both multilayer networks have the same layers, the same **actor sets**,
#' and the same **edge sets/weights** (to a tolerance). When a mismatch is
#' found, the function **prints a human-readable diff** (which layers differ,
#' how they differ, and sample items).
#'
#' Differences reported per layer:
#' - **Directedness** mismatch (one directed, the other undirected)
#' - **Actors only in A** / **only in B** (vertex `name`)
#' - **Edges only in A** / **only in B** (canonical undirected keys `u||v` or
#'   directed keys `u->v`)
#' - **Weight mismatches** on common edges beyond `tol` (shows `wA`, `wB`, `Δ`)
#'
#' @param net_a,net_b Two `multinet::ml.network` objects to compare.
#' @param tol Numeric tolerance for weight comparison. Default `1e-9`.
#' @param show_max Maximum number of mismatched items to print per category
#'   (actors / edges / weights) for each layer. Default `10`.
#' @param verbose Logical; if `TRUE`, prints a summary and per-layer differences.
#'   Default `TRUE`.
#'
#' @return (Invisibly) a list with:
#'   - `ok`: logical, `TRUE` if fully equivalent, else `FALSE`
#'   - `layers_all`: character vector of all layer names considered
#'   - `missing_in_a`, `missing_in_b`: layers missing from A or B
#'   - `per_layer`: named list; for each layer:
#'       - `ok_layer`: logical; `TRUE` if that layer matches
#'       - `directed_mismatch`: logical
#'       - `actors_only_a`, `actors_only_b`: character vectors
#'       - `edges_only_a`, `edges_only_b`: character vectors of edge keys
#'       - `weight_mismatches`: data.frame with `key`, `wA`, `wB`, `delta`
#'
#' @examples
#' \dontrun{
#' # Assuming you've already built both multilayer networks:
#' # net_edges <- build_multiNet(...); net_graphs <- build_multinet_from_graphs(...)
#'
#' diff <- compare_multinets(
#'   net_a   = net_edges,
#'   net_b   = net_graphs,
#'   tol     = 1e-9,
#'   show_max= 10,
#'   verbose = TRUE
#' )
#'
#' if (!diff$ok) {
#'   # Programmatic access to details:
#'   names(diff$per_layer)                # layers compared
#'   diff$per_layer[["E1"]]$edges_only_a  # edges only in A for layer E1
#' }
#' }
#' @export
compare_multinets <- function(net_a, net_b, tol = 1e-9, show_max = 10, verbose = TRUE) {
  # ---- helpers ---------------------------------------------------------------
  .print_head <- function(x, n = show_max) {
    x <- unique(x)
    if (!length(x)) return(invisible(NULL))
    cat(paste(utils::head(x, n), collapse = ", "))
    if (length(x) > n) cat(" ...")
    cat("\n")
  }
  .agg_edges <- function(g) {
    # aggregate parallel edges by endpoint pair; sum weights
    el <- igraph::as_edgelist(g, names = TRUE)
    if (nrow(el) == 0L) {
      return(data.frame(key = character(), weight = numeric(), stringsAsFactors = FALSE))
    }
    w <- igraph::edge_attr(g, "weight")
    if (is.null(w)) w <- rep(1, nrow(el)) else {
      w <- suppressWarnings(as.numeric(w))
      w[is.na(w)] <- 1
    }
    df <- data.frame(from = el[, 1], to = el[, 2], weight = w, stringsAsFactors = FALSE)
    if (!igraph::is_directed(g)) {
      key <- ifelse(df$from <= df$to,
                    paste(df$from, df$to, sep = "||"),
                    paste(df$to, df$from, sep = "||"))
    } else {
      key <- paste(df$from, "->", df$to, sep = "")
    }
    stats::aggregate(weight ~ key, data = transform(df, key = key),
                     FUN = function(z) sum(z, na.rm = TRUE))
  }
  .to_layer_graphs <- function(net) {
    lst <- try(as.list(net), silent = TRUE)
    if (inherits(lst, "try-error") || !is.list(lst))
      stop("Could not coerce ml.network to a list of igraphs via as.list().")
    # keep only igraph elements and drop '_flat_' if present
    keep <- vapply(lst, inherits, logical(1), what = "igraph")
    lst <- lst[keep]
    nm <- names(lst)
    if (!is.null(nm) && "_flat_" %in% nm) lst <- lst[nm != "_flat_"]
    lst
  }

  # ---- layer set reconciliation ---------------------------------------------
  la <- multinet::layers_ml(net_a)
  lb <- multinet::layers_ml(net_b)
  layers_all <- sort(unique(c(la, lb)))
  missing_in_a <- setdiff(layers_all, la)
  missing_in_b <- setdiff(layers_all, lb)

  if (verbose) {
    cat("Comparing multilayer networks...\n")
    cat("  Layers in A:", length(la), " | Layers in B:", length(lb), "\n")
    if (length(missing_in_a)) cat("  Missing in A:", paste(missing_in_a, collapse = ", "), "\n")
    if (length(missing_in_b)) cat("  Missing in B:", paste(missing_in_b, collapse = ", "), "\n")
  }

  # if layer sets differ, we still compute diffs for the intersection and report
  layers <- sort(intersect(la, lb))
  if (!length(layers)) {
    if (verbose) cat("No common layers to compare.\n")
    return(invisible(list(
      ok = FALSE,
      layers_all = layers_all,
      missing_in_a = missing_in_a,
      missing_in_b = missing_in_b,
      per_layer = list()
    )))
  }

  # ---- convert to per-layer igraphs -----------------------------------------
  gla <- .to_layer_graphs(net_a)
  glb <- .to_layer_graphs(net_b)

  # ---- per-layer comparisons -------------------------------------------------
  details <- setNames(vector("list", length(layers)), layers)
  matched_layers <- 0L

  for (ln in layers) {
    g1 <- gla[[ln]]; g2 <- glb[[ln]]
    if (is.null(g1) || is.null(g2)) next

    # directedness
    dir_mis <- igraph::is_directed(g1) != igraph::is_directed(g2)

    # actors
    a1 <- sort(as.character(igraph::V(g1)$name))
    a2 <- sort(as.character(igraph::V(g2)$name))
    act_only_a <- setdiff(a1, a2)
    act_only_b <- setdiff(a2, a1)

    # edges & weights
    e1 <- .agg_edges(g1); e2 <- .agg_edges(g2)
    keys_a <- e1$key; keys_b <- e2$key
    edges_only_a <- setdiff(keys_a, keys_b)
    edges_only_b <- setdiff(keys_b, keys_a)

    # common edges weight comparison
    common <- intersect(keys_a, keys_b)
    wm <- data.frame()
    if (length(common)) {
      e1c <- e1[match(common, e1$key), , drop = FALSE]
      e2c <- e2[match(common, e2$key), , drop = FALSE]
      delta <- e1c$weight - e2c$weight
      bad <- which(abs(delta) > tol)
      if (length(bad)) {
        wm <- data.frame(
          key   = common[bad],
          wA    = e1c$weight[bad],
          wB    = e2c$weight[bad],
          delta = delta[bad],
          stringsAsFactors = FALSE
        )
        # sort by largest absolute delta
        wm <- wm[order(-abs(wm$delta)), , drop = FALSE]
      }
    }

    ok_layer <- !dir_mis && !length(act_only_a) && !length(act_only_b) &&
      !length(edges_only_a) && !length(edges_only_b) && !nrow(wm)

    details[[ln]] <- list(
      ok_layer          = ok_layer,
      directed_mismatch = dir_mis,
      actors_only_a     = act_only_a,
      actors_only_b     = act_only_b,
      edges_only_a      = edges_only_a,
      edges_only_b      = edges_only_b,
      weight_mismatches = wm
    )

    if (ok_layer) {
      matched_layers <- matched_layers + 1L
      if (verbose) cat(sprintf("✓ Layer %-12s — OK\n", ln))
    } else if (verbose) {
      cat(sprintf("✗ Layer %-12s — differences detected:\n", ln))
      if (dir_mis) {
        cat("  • Directedness mismatch: ",
            if (igraph::is_directed(g1)) "A=directed, " else "A=undirected, ",
            if (igraph::is_directed(g2)) "B=directed\n" else "B=undirected\n", sep = "")
      }
      if (length(act_only_a)) {
        cat("  • Actors only in A (", length(act_only_a), "): ", sep = "")
        .print_head(act_only_a, show_max)
      }
      if (length(act_only_b)) {
        cat("  • Actors only in B (", length(act_only_b), "): ", sep = "")
        .print_head(act_only_b, show_max)
      }
      if (length(edges_only_a)) {
        cat("  • Edges only in A (", length(edges_only_a), "): ", sep = "")
        .print_head(edges_only_a, show_max)
      }
      if (length(edges_only_b)) {
        cat("  • Edges only in B (", length(edges_only_b), "): ", sep = "")
        .print_head(edges_only_b, show_max)
      }
      if (nrow(wm)) {
        cat("  • Weight mismatches beyond tol (", nrow(wm), "), showing up to ", show_max, ":\n", sep = "")
        to_show <- utils::head(wm, show_max)
        print(to_show, row.names = FALSE)
      }
    }
  }

  all_ok <- matched_layers == length(layers) && !length(missing_in_a) && !length(missing_in_b)

  if (verbose) {
    cat(sprintf("\nSummary: %d/%d common layers matched exactly (tol=%.2g).\n",
                matched_layers, length(layers), tol))
    if (!all_ok && (length(missing_in_a) || length(missing_in_b))) {
      cat("Note: layer sets differ between A and B; see 'Missing in A/B' above.\n")
    }
  }

  invisible(list(
    ok            = all_ok,
    layers_all    = layers_all,
    missing_in_a  = missing_in_a,
    missing_in_b  = missing_in_b,
    per_layer     = details
  ))
}

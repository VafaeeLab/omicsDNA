
# ---------------------------------------------------------------------------
# detectCom.R
# Community detection wrapper for multilayer networks (multinet)
# ---------------------------------------------------------------------------

#' Community detection for multilayer networks
#'
#' @md
#' @description
#' `detectCom()` identifies communities (modules) in a multilayer network of class
#' \code{multinet::ml.network}. The function supports two analysis modes:
#'
#' * **Multilayer community detection (default; \code{supra_graph = FALSE})**:
#'   community detection is performed directly on the multilayer structure using
#'   the corresponding \pkg{multinet} routine (e.g., generalized Louvain, Infomap,
#'   clique percolation, ABACUS). In this mode the function returns the native
#'   \pkg{multinet} community object.
#'
#' * **Supra-graph community detection (\code{supra_graph = TRUE})**:
#'   selected layers are first aggregated into a single undirected actor–actor graph
#'   (a “supra-graph”). A standard \pkg{igraph} community algorithm is then applied,
#'   and community memberships are mapped back to each layer where the actor appears.
#'
#' In supra-graph mode, communities can be filtered by minimum size (\code{min.actors})
#' and minimum layer coverage (\code{min.layers}), and optionally relabelled by decreasing
#' community size (\code{relabel_by_size}).
#'
#' @details
#' ## Input to the function
#' The primary input is a multilayer network object (\code{net}) where actors correspond
#' to node identifiers (e.g., gene symbols) and edges represent within-layer interactions.
#' In supra-graph mode, edges are read from \code{multinet::edges_ml(net)} and (when available)
#' restricted to intra-layer edges before aggregation.
#'
#' ## Edge aggregation in supra-graph mode
#' When \code{supra_graph = TRUE}, edge aggregation across layers is controlled by \code{edgeWeight}:
#' * \code{"count"}: supra-edge weight equals the number of layers containing that edge;
#' * \code{"sum"}: supra-edge weight equals the sum of edge weights across layers (when available).
#'
#' ## Output files written to disk
#' Output files are optional and are written under \code{results_dir} (created if missing).
#'
#' * *CSV outputs (two files; controlled by \code{write_csv})*
#'   When \code{write_csv = TRUE} (default), the function writes two CSV files:
#'   \enumerate{
#'     \item {A CSV representation of the returned object:}
#'           In multilayer mode, this file is produced via \code{as.data.frame()} on the
#'           \pkg{multinet} community object. If coercion fails, a small placeholder CSV is written.
#'
#'     \item {The main membership table (CSV):}
#'           containing (at minimum) \code{actor}, \code{layer}, \code{com}, and \code{method}.
#'           In multilayer mode, memberships are extracted from the \pkg{multinet} result when possible;
#'           otherwise an empty membership table is written.
#'   }
#'
#'
#' @param net A multilayer network of class \code{multinet::ml.network}.
#'
#' @param method Community detection method. One of \code{"glouvain"} (or \code{"louvain"} as an alias),
#' \code{"infomap"}, \code{"clique"}, \code{"abacus"}.
#'
#' @param supra_graph Logical. If \code{FALSE} (default), communities are detected using \pkg{multinet}.
#' If \code{TRUE}, communities are detected on a layer-aggregated supra-graph using \pkg{igraph}.
#'
#' @param layers Optional character vector of layer names to include when \code{supra_graph = TRUE}.
#' Ignored when \code{supra_graph = FALSE}.
#'
#' @param gamma Numeric (>0). Resolution parameter for \code{multinet::glouvain_ml()} (multilayer mode).
#' @param omega Numeric (>=0). Inter-layer coupling for \code{multinet::glouvain_ml()} (multilayer mode).
#'
#' @param overlapping,directed,self.links Passed to \code{multinet::infomap_ml()} (multilayer mode).
#'
#' @param k Integer (>=3). Clique size threshold for clique percolation (multilayer mode) and
#' minimum clique size in supra-graph mode.
#'
#' @param m Integer (>=1). Clique percolation parameter (multilayer mode only).
#'
#' @param edgeWeight Aggregation rule for supra-graph edges: \code{"count"} or \code{"sum"}.
#'
#' @param min.actors Integer (>=1). In supra-graph mode, communities with fewer than this many
#' unique actors are removed. In multilayer mode, this is passed to \code{multinet::abacus_ml()}.
#'
#' @param min.layers Integer (>=1). In supra-graph mode, communities represented in fewer than
#' this many layers are removed. In multilayer mode, this is passed to \code{multinet::abacus_ml()}.
#'
#' @param relabel_by_size Logical. In supra-graph mode, optionally relabel retained communities
#' to \code{C1, C2, ...} by decreasing size. In multilayer mode, this affects only the *saved*
#' membership CSV (the returned \pkg{multinet} object is not modified).
#'
#' @param seed Optional integer seed for reproducibility.
#'
#' @param results_dir Output directory for optional files. Created if missing.
#'
#' @param save_to_rds Logical. If \code{TRUE}, saves the returned object (multilayer mode) or the
#' membership table (supra-graph mode) as an \code{.rds}.
#'
#' @param rds_file Optional filename for the RDS output. If \code{NULL}, a timestamped name is used.
#'
#' @param write_csv Logical. If \code{TRUE} (default), writes \strong{two} CSV files under \code{results_dir}:
#' (i) a CSV representation of the returned object (\code{<base>.csv}) and (ii) the main community
#' membership table (\code{<base>_membership.csv}).
#'
#' @param csv_prefix Character prefix used to build output filenames.
#'
#' @param verbose Logical. If \code{TRUE}, prints progress messages and output paths.
#'
#' @param ... Additional arguments passed to the underlying \code{multinet::*_ml()} function in
#' multilayer mode. Also accepts legacy alias names for compatibility.
#'
#' @return
#' * If \code{supra_graph = FALSE}: returns the object returned by the selected \code{multinet::*_ml()} method.
#' * If \code{supra_graph = TRUE}: returns a data frame with columns \code{actor}, \code{com}, \code{layer}, \code{method},
#'   with attributes \code{community_sizes}, \code{layer_span}, and \code{files}.
#'
#' @importFrom multinet ml_empty add_igraph_layer_ml edges_ml
#' @importFrom multinet glouvain_ml infomap_ml clique_percolation_ml abacus_ml
#' @importFrom igraph graph_from_data_frame simplify cliques components cluster_louvain cluster_infomap membership V E vcount ecount
#' @importFrom stats aggregate
#' @importFrom utils write.csv
#' @export
detectCom <- function(
    net,
    method = c("glouvain", "louvain", "infomap", "clique", "abacus"),
    supra_graph = FALSE,
    layers = NULL,

    # --- legacy / package-style names (kept for compatibility) ---
    glouvain_gamma = 1,
    glouvain_omega = 1,
    clique.k = 3,
    clique.m = 1,
    infomap_overlapping = FALSE,
    infomap_directed    = FALSE,
    infomap_self_links  = TRUE,

    # --- user-facing aliases ---
    gamma       = NULL,
    omega       = NULL,
    k           = NULL,
    m           = NULL,
    overlapping = NULL,
    directed    = NULL,
    self.links  = NULL,

    # --- supra-graph options ---
    edgeWeight = c("count", "sum"),

    # --- post-processing (applied only when supra_graph = TRUE) ---
    min.actors      = 15,
    min.layers      = 2,
    relabel_by_size = TRUE,

    # --- reproducibility ---
    seed = NULL,

    # --- outputs ---
    results_dir       = getOption("mlnet.results_dir", "omicsDNA_results"),
    save_to_rds       = FALSE,
    rds_file          = NULL,
    write_csv         = TRUE,
    csv_prefix        = "communities",
    verbose           = TRUE,

    ...
) {
  # ---------------------------
  # normalise / validate
  # ---------------------------
  method <- match.arg(method)
  if (identical(method, "louvain")) {
    warning("`method = \"louvain\"` is deprecated; using `method = \"glouvain\"`.", call. = FALSE)
    method <- "glouvain"
  }
  edgeWeight <- match.arg(edgeWeight)

  stopifnot(is.logical(supra_graph), length(supra_graph) == 1L)
  stopifnot(is.logical(verbose),     length(verbose)     == 1L)
  stopifnot(is.numeric(min.actors),  length(min.actors)  == 1L, min.actors >= 1)
  stopifnot(is.numeric(min.layers),  length(min.layers)  == 1L, min.layers >= 1)
  stopifnot(is.logical(relabel_by_size), length(relabel_by_size) == 1L)

  # resolve aliases (user-facing names override legacy)
  gamma_use <- if (!is.null(gamma)) gamma else glouvain_gamma
  omega_use <- if (!is.null(omega)) omega else glouvain_omega

  k_use <- if (!is.null(k)) k else clique.k
  m_use <- if (!is.null(m)) m else clique.m

  overlapping_use <- if (!is.null(overlapping)) overlapping else infomap_overlapping
  directed_use    <- if (!is.null(directed))    directed    else infomap_directed
  self_links_use  <- if (!is.null(self.links))  self.links  else infomap_self_links

  stopifnot(is.numeric(gamma_use), length(gamma_use) == 1L, gamma_use > 0)
  stopifnot(is.numeric(omega_use), length(omega_use) == 1L, omega_use >= 0)
  stopifnot(is.numeric(k_use),     length(k_use)     == 1L, k_use >= 3)
  stopifnot(is.numeric(m_use),     length(m_use)     == 1L, m_use >= 1)

  stopifnot(is.logical(overlapping_use), length(overlapping_use) == 1L)
  stopifnot(is.logical(directed_use),    length(directed_use)    == 1L)
  stopifnot(is.logical(self_links_use),  length(self_links_use)  == 1L)

  # ABACUS is only implemented here via multinet
  use_multilayer_method <- (!isTRUE(supra_graph)) || method == "abacus"

  # ---------------------------
  # helpers
  # ---------------------------
  .ensure_dir <- function(d) {
    if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
    invisible(d)
  }
  .is_abs <- function(p) {
    is.character(p) && length(p) == 1L && grepl("^(/|[A-Za-z]:[\\/])", p)
  }
  .stamp <- function() format(Sys.time(), "%Y-%m-%d_%H%M%S")

  .pick_ci <- function(cands, nms) {
    nms0 <- as.character(nms)
    nms1 <- tolower(trimws(nms0))
    for (c in cands) {
      w <- which(nms1 == tolower(trimws(c)))
      if (length(w)) return(nms0[w[1]])
    }
    NA_character_
  }

  # ---- robust coercion to membership table: actor / layer / com ------------
  .coerce_membership_table <- function(obj, net = NULL) {
    df <- NULL
    if (is.data.frame(obj)) {
      df <- obj
    } else {
      tmp <- try(as.data.frame(obj, stringsAsFactors = FALSE), silent = TRUE)
      if (!inherits(tmp, "try-error")) df <- tmp
    }

    # If we have a table, standardise column selection
    if (is.data.frame(df) && nrow(df)) {
      nm <- names(df)

      actor_col <- .pick_ci(c("actor","node","vertex","name"), nm)
      layer_col <- .pick_ci(c("layer"), nm)
      com_col   <- .pick_ci(c("com","cid","community","cluster","community_id","module","group"), nm)

      # Fallback: if column names are not informative but the object is clearly 3-column,
      # interpret as actor / layer / community-id (multinet’s conventional layout).
      if (any(is.na(c(actor_col, layer_col, com_col))) && ncol(df) >= 3L) {
        actor_col <- nm[1]
        layer_col <- nm[2]
        com_col   <- nm[3]
      }

      if (!any(is.na(c(actor_col, layer_col, com_col)))) {
        out <- unique(data.frame(
          actor = as.character(df[[actor_col]]),
          layer = as.character(df[[layer_col]]),
          com   = as.character(df[[com_col]]),
          stringsAsFactors = FALSE
        ))
        out$actor <- trimws(out$actor)
        out$layer <- trimws(out$layer)
        out$com   <- trimws(out$com)
        out <- out[nzchar(out$actor) & nzchar(out$layer) & nzchar(out$com), , drop = FALSE]
        return(out)
      }
    }

    # Fallback: if obj is not directly tabular, attempt multinet-native conversion
    # using get_community_list_ml() + vertices_ml().
    if (!is.null(net)) {
      cl <- try(multinet::get_community_list_ml(obj, net), silent = TRUE)
      vt <- try(multinet::vertices_ml(net), silent = TRUE)

      if (!inherits(cl, "try-error") && is.list(cl) && length(cl) &&
          !inherits(vt, "try-error") && is.data.frame(vt) && nrow(vt) >= 1L) {

        # vertices_ml() is documented as a 2-column table: actor + layer
        vnm <- names(vt)
        v_actor <- .pick_ci(c("actor"), vnm)
        v_layer <- .pick_ci(c("layer"), vnm)
        if (is.na(v_actor) || is.na(v_layer)) {
          # last-resort: take first two columns
          if (ncol(vt) >= 2L) {
            v_actor <- vnm[1]
            v_layer <- vnm[2]
          }
        }

        rows <- list()
        for (i in seq_along(cl)) {
          el <- cl[[i]]

          cid_i <- NA_character_
          verts_i <- NULL

          if (is.list(el)) {
            # community id
            if ("cid" %in% names(el))       cid_i <- as.character(el$cid)
            if ("com" %in% names(el))       cid_i <- as.character(el$com)
            if ("community" %in% names(el)) cid_i <- as.character(el$community)

            # vertices
            for (kk in c("vertices","vertex","v","nodes","ids")) {
              if (kk %in% names(el)) { verts_i <- el[[kk]]; break }
            }
            if (is.null(verts_i)) {
              # try first numeric-like component
              for (jj in seq_along(el)) {
                if (is.numeric(el[[jj]]) || is.integer(el[[jj]])) { verts_i <- el[[jj]]; break }
              }
            }
          } else if (is.numeric(el) || is.integer(el)) {
            verts_i <- el
            # sometimes cid may be encoded in the list element name
            nm_i <- names(cl)[i]
            if (!is.null(nm_i) && nzchar(nm_i)) {
              # extract digits after "cid" when possible
              m <- regexpr("cid\\s*=?\\s*([0-9]+)", nm_i, perl = TRUE)
              if (m[1] != -1) {
                cid_i <- sub(".*cid\\s*=?\\s*([0-9]+).*", "\\1", nm_i, perl = TRUE)
              }
            }
          }

          if (is.null(verts_i) || !length(verts_i)) next
          verts_i <- as.integer(verts_i)
          verts_i <- verts_i[is.finite(verts_i) & verts_i >= 1L & verts_i <= nrow(vt)]
          if (!length(verts_i)) next

          subvt <- vt[verts_i, , drop = FALSE]
          if (is.na(cid_i) || !nzchar(cid_i)) next

          dd <- data.frame(
            actor = as.character(subvt[[v_actor]]),
            layer = as.character(subvt[[v_layer]]),
            com   = as.character(cid_i),
            stringsAsFactors = FALSE
          )
          dd$actor <- trimws(dd$actor)
          dd$layer <- trimws(dd$layer)
          dd$com   <- trimws(dd$com)
          dd <- dd[nzchar(dd$actor) & nzchar(dd$layer) & nzchar(dd$com), , drop = FALSE]

          if (nrow(dd)) rows[[length(rows) + 1L]] <- dd
        }

        if (length(rows)) {
          out <- unique(do.call(rbind, rows))
          if (nrow(out)) return(out)
        }
      }
    }

    NULL
  }

  # ---- WRITE MAIN MEMBERSHIP CSV (legacy function name retained) ------------
  # NOTE: This no longer writes a com/size/span summary. It writes the *membership table*.
  .write_summary_csv <- function(membership_df, base) {
    .ensure_dir(results_dir)
    mem_file <- file.path(results_dir, paste0(base, "_membership.csv"))

    # standardise empty output
    if (is.null(membership_df) || !is.data.frame(membership_df) || !nrow(membership_df)) {
      out <- data.frame(
        actor  = character(),
        com    = character(),
        layer  = character(),
        method = character(),
        stringsAsFactors = FALSE
      )
      utils::write.csv(out, mem_file, row.names = FALSE)
      if (isTRUE(verbose)) message("Saved membership CSV (empty): ",
                                   normalizePath(mem_file, winslash = "/", mustWork = FALSE))
      return(mem_file)
    }

    out <- membership_df
    # ensure required columns exist
    need <- c("actor","layer","com")
    miss <- setdiff(need, names(out))
    if (length(miss)) {
      stop("Membership table is missing required columns: ", paste(miss, collapse = ", "))
    }

    out$actor <- trimws(as.character(out$actor))
    out$layer <- trimws(as.character(out$layer))
    out$com   <- trimws(as.character(out$com))
    out <- out[nzchar(out$actor) & nzchar(out$layer) & nzchar(out$com), , drop = FALSE]
    out <- unique(out[, intersect(c("actor","com","layer","method"), names(out)), drop = FALSE])

    if (!"method" %in% names(out)) out$method <- NA_character_

    # stable column order
    out <- out[, c("actor","com","layer","method"), drop = FALSE]

    utils::write.csv(out, mem_file, row.names = FALSE)
    if (isTRUE(verbose)) message("Saved membership CSV: ",
                                 normalizePath(mem_file, winslash = "/", mustWork = FALSE))
    mem_file
  }

  .save_outputs_raw <- function(obj, method_used, mode_tag, weight_tag = NULL) {
    if (!isTRUE(save_to_rds) && !isTRUE(write_csv)) return(invisible(NULL))

    .ensure_dir(results_dir)
    base <- if (is.null(weight_tag)) {
      sprintf("%s_%s_%s_%s", csv_prefix, tolower(method_used), mode_tag, .stamp())
    } else {
      sprintf("%s_%s_%s_%s_%s", csv_prefix, tolower(method_used), mode_tag, tolower(weight_tag), .stamp())
    }

    # RDS (save raw object unchanged)
    if (isTRUE(save_to_rds)) {
      rds_use <- if (is.null(rds_file)) paste0(base, ".rds") else rds_file
      if (!.is_abs(rds_use)) rds_use <- file.path(results_dir, rds_use)
      saveRDS(obj, rds_use)
      if (isTRUE(verbose)) message("Saved RDS: ",
                                   normalizePath(rds_use, winslash = "/", mustWork = FALSE))
    }

    # CSV outputs (two files)
    if (isTRUE(write_csv)) {

      # (1) CSV representation of returned object (coerce copy only; placeholder if coercion fails)
      csv_file <- file.path(results_dir, paste0(base, ".csv"))
      obj_df <- if (is.data.frame(obj)) obj else try(as.data.frame(obj, stringsAsFactors = FALSE), silent = TRUE)
      if (inherits(obj_df, "try-error") || is.null(obj_df)) {
        obj_df <- data.frame(
          note  = "Output object could not be coerced to a data.frame; writing placeholder CSV.",
          class = paste(class(obj), collapse = ";"),
          stringsAsFactors = FALSE
        )
        warning("write_csv=TRUE but output could not be coerced to a data.frame; wrote placeholder CSV.", call. = FALSE)
      }
      utils::write.csv(obj_df, csv_file, row.names = FALSE)
      if (isTRUE(verbose)) message("Saved CSV: ",
                                   normalizePath(csv_file, winslash = "/", mustWork = FALSE))

      # (2) MAIN MEMBERSHIP CSV
      mem <- .coerce_membership_table(obj, net = net)
      if (!is.null(mem) && is.data.frame(mem) && nrow(mem)) {
        # optional relabel-by-size for file output only (multilayer mode)
        if (isTRUE(relabel_by_size)) {
          size_by <- tapply(mem$actor, mem$com, function(v) length(unique(v)))
          ord <- order(as.integer(size_by), decreasing = TRUE, na.last = TRUE)
          old_ids <- names(size_by)[ord]
          old_ids <- old_ids[!is.na(old_ids) & nzchar(old_ids)]
          if (length(old_ids)) {
            map <- stats::setNames(paste0("C", seq_along(old_ids)), old_ids)
            mem$com <- unname(map[mem$com])
          }
        }
        mem$method <- as.character(method_used)
      }
      .write_summary_csv(mem, base)
    }

    invisible(NULL)
  }

  .finalize_supra <- function(res_df, method_used, weight_tag = NULL) {
    .ensure_dir(results_dir)
    files <- list()

    base <- if (is.null(weight_tag)) {
      sprintf("%s_%s_supra_%s", csv_prefix, tolower(method_used), .stamp())
    } else {
      sprintf("%s_%s_supra_%s_%s", csv_prefix, tolower(method_used), tolower(weight_tag), .stamp())
    }

    if (isTRUE(save_to_rds)) {
      rds_use <- if (is.null(rds_file)) paste0(base, ".rds") else rds_file
      if (!.is_abs(rds_use)) rds_use <- file.path(results_dir, rds_use)
      saveRDS(res_df, rds_use)
      files$rds <- rds_use
      if (isTRUE(verbose)) message("Saved RDS: ",
                                   normalizePath(rds_use, winslash = "/", mustWork = FALSE))
    }

    # CSV outputs (two files)
    if (isTRUE(write_csv)) {
      csv_file <- file.path(results_dir, paste0(base, ".csv"))
      utils::write.csv(res_df, csv_file, row.names = FALSE)
      files$csv <- csv_file
      if (isTRUE(verbose)) message("Saved CSV: ",
                                   normalizePath(csv_file, winslash = "/", mustWork = FALSE))

      mem_file <- .write_summary_csv(res_df, base)
      files$membership_csv <- mem_file
    }

    # attributes (always)
    if (nrow(res_df)) {
      size_by <- tapply(res_df$actor, res_df$com, function(v) length(unique(v)))
      span_by <- tapply(res_df$layer, res_df$com, function(v) length(unique(v)))
      ord <- order(as.integer(size_by), decreasing = TRUE)
      size_by <- as.integer(size_by[ord])
      span_by <- as.integer(span_by[names(size_by)])
    } else {
      size_by <- integer(0)
      span_by <- integer(0)
    }

    attr(res_df, "community_sizes") <- size_by
    attr(res_df, "layer_span")      <- span_by
    attr(res_df, "files")           <- files

    res_df
  }

  # ---------------------------
  # RNG seed
  # ---------------------------
  if (!is.null(seed)) set.seed(as.integer(seed))

  # =======================================================================
  # 1) Multilayer community detection via multinet (supra_graph = FALSE OR abacus)
  # =======================================================================
  if (use_multilayer_method) {
    if (isTRUE(verbose)) {
      msg <- paste0("Running multilayer community detection via multinet, method = \"", method, "\".")
      if (!is.null(layers) && length(layers) && isFALSE(supra_graph)) {
        msg <- paste0(msg, " (`layers` is only used when supra_graph = TRUE.)")
      }
      message(msg)
    }

    cm_raw <- switch(
      method,
      glouvain = multinet::glouvain_ml(net, gamma = gamma_use, omega = omega_use, ...),
      infomap  = multinet::infomap_ml(net,
                                      overlapping = overlapping_use,
                                      directed    = directed_use,
                                      self.links  = self_links_use, ...),
      clique   = multinet::clique_percolation_ml(net, k = as.integer(k_use), m = as.integer(m_use), ...),
      abacus   = multinet::abacus_ml(net, min.actors = min.actors, min.layers = min.layers, ...),
      stop("Unsupported method: ", method, call. = FALSE)
    )

    # side-effect outputs only; return object unchanged
    .save_outputs_raw(cm_raw, method_used = method, mode_tag = "ml", weight_tag = NULL)

    return(cm_raw)
  }

  # =======================================================================
  # 2) Supra-graph mode (supra_graph = TRUE; glouvain/infomap/clique)
  # =======================================================================
  if (!method %in% c("glouvain", "infomap", "clique")) {
    stop("supra_graph=TRUE is supported only for method = 'glouvain', 'infomap', or 'clique'.", call. = FALSE)
  }

  if (isTRUE(verbose)) {
    message("Running supra-graph community detection (layer-aggregated), method = \"", method, "\".")
  }

  # ---- pull multilayer edges ----
  Eraw <- try(multinet::edges_ml(net), silent = TRUE)
  if (inherits(Eraw, "try-error") || is.null(Eraw)) {
    stop("Could not retrieve edges from `net` via multinet::edges_ml().", call. = FALSE)
  }
  Eall <- if (is.data.frame(Eraw)) Eraw else {
    tmp <- try(as.data.frame(Eraw, stringsAsFactors = FALSE), silent = TRUE)
    if (inherits(tmp, "try-error") || is.null(tmp)) stop("edges_ml(net) did not return a coercible table.", call. = FALSE)
    tmp
  }
  if (!nrow(Eall)) stop("No edges found in `net`.", call. = FALSE)

  nm <- names(Eall)

  # endpoints
  pairs <- list(
    c("from_actor","to_actor"), c("from","to"), c("source","target"),
    c("actor1","actor2"), c("i","j"), c("v1","v2")
  )
  a_col <- b_col <- NA_character_
  for (p in pairs) {
    if (all(p %in% nm)) { a_col <- p[1]; b_col <- p[2]; break }
  }
  if (is.na(a_col)) {
    layer_like <- c("layer","Layer","from_layer","to_layer","l1","l2")
    char_cols <- which(vapply(Eall, function(x) is.character(x) || is.factor(x), logical(1)))
    char_cols <- setdiff(char_cols, match(layer_like, nm, nomatch = 0))
    if (length(char_cols) < 2) stop("Could not identify two endpoint columns in edges table.", call. = FALSE)
    a_col <- nm[char_cols[1]]
    b_col <- nm[char_cols[2]]
  }

  # layer + intra-layer filtering
  if ("from_layer" %in% nm && "to_layer" %in% nm) {
    Eall <- Eall[Eall$from_layer == Eall$to_layer, , drop = FALSE]
    Eall$layer <- as.character(Eall$from_layer)
  } else if ("layer" %in% nm) {
    Eall$layer <- as.character(Eall$layer)
  } else if ("Layer" %in% nm) {
    Eall$layer <- as.character(Eall$Layer)
  } else {
    Eall$layer <- "L1"
  }
  if (!nrow(Eall)) stop("No intra-layer edges found.", call. = FALSE)

  # choose layers
  if (is.null(layers)) {
    layers_use <- sort(unique(Eall$layer))
  } else {
    layers_use <- intersect(as.character(layers), unique(Eall$layer))
  }
  if (!length(layers_use)) stop("No layers available after applying `layers` selection.", call. = FALSE)
  Eall <- Eall[Eall$layer %in% layers_use, , drop = FALSE]
  if (!nrow(Eall)) stop("No edges left after applying `layers` selection.", call. = FALSE)

  # weight col (optional)
  w_col <- .pick_ci(c("weight","w","w_","value","score"), nm)

  edges_all <- data.frame(
    layer  = as.character(Eall$layer),
    from   = as.character(Eall[[a_col]]),
    to     = as.character(Eall[[b_col]]),
    weight = if (!is.na(w_col) && w_col %in% names(Eall)) suppressWarnings(as.numeric(Eall[[w_col]])) else 1,
    stringsAsFactors = FALSE
  )

  # canonicalize undirected
  aa <- pmin(edges_all$from, edges_all$to)
  bb <- pmax(edges_all$from, edges_all$to)
  edges_all$from <- aa
  edges_all$to   <- bb

  # aggregate supra-edges
  if (edgeWeight == "count") {
    per_layer_unique <- unique(edges_all[, c("layer","from","to")])
    agg <- stats::aggregate(rep(1L, nrow(per_layer_unique)),
                            by = per_layer_unique[c("from","to")], FUN = sum)
    names(agg) <- c("from","to","weight")
  } else {
    agg <- stats::aggregate(weight ~ from + to, data = edges_all, FUN = sum, na.rm = TRUE)
  }

  g <- igraph::graph_from_data_frame(agg, directed = FALSE)
  if (igraph::ecount(g) == 0) stop("Supra-graph has no edges after aggregation.", call. = FALSE)
  g <- igraph::simplify(g, edge.attr.comb = list(weight = "sum", "ignore"))

  if (isTRUE(verbose)) {
    message("Built supra-graph: |V|=", igraph::vcount(g), " |E|=", igraph::ecount(g),
            " (edgeWeight=", edgeWeight, ")")
  }

  # detect communities on supra graph
  memb <- NULL
  if (method == "glouvain") {
    co <- igraph::cluster_louvain(g, weights = igraph::E(g)$weight)
    memb <- igraph::membership(co)
  } else if (method == "infomap") {
    co <- igraph::cluster_infomap(g, e.weights = igraph::E(g)$weight)
    memb <- igraph::membership(co)
  } else if (method == "clique") {
    cls <- igraph::cliques(g, min = as.integer(k_use))
    if (!length(cls)) stop("No cliques of size >= ", as.integer(k_use), " found in supra-graph.", call. = FALSE)

    pairs_df <- do.call(rbind, lapply(cls, function(cl) {
      vs <- igraph::V(g)$name[cl]
      if (length(vs) < 2) return(NULL)
      mm <- utils::combn(vs, 2)
      data.frame(from = mm[1, ], to = mm[2, ], stringsAsFactors = FALSE)
    }))
    pairs_df <- unique(pairs_df)

    h <- igraph::graph_from_data_frame(pairs_df, directed = FALSE, vertices = igraph::V(g)$name)
    comp <- igraph::components(h)
    memb <- comp$membership
    names(memb) <- igraph::V(h)$name
  }

  if (is.null(memb) || !length(memb)) stop("Supra-graph community detection failed.", call. = FALSE)

  memb_vec <- as.integer(memb)
  names(memb_vec) <- names(memb)
  memb_df <- data.frame(actor = names(memb_vec), com = as.character(memb_vec), stringsAsFactors = FALSE)

  # map back to layers (actors that appear in each layer)
  actors_by_layer <- lapply(layers_use, function(ly) {
    ed <- edges_all[edges_all$layer == ly, , drop = FALSE]
    unique(c(ed$from, ed$to))
  })
  names(actors_by_layer) <- layers_use

  res <- do.call(rbind, lapply(layers_use, function(ly) {
    a <- actors_by_layer[[ly]]
    df <- memb_df[memb_df$actor %in% a, , drop = FALSE]
    if (!nrow(df)) return(NULL)
    df$layer  <- ly
    df$method <- method
    df
  }))

  if (is.null(res) || !nrow(res)) {
    empty <- data.frame(actor = character(), com = character(), layer = character(), method = character(),
                        stringsAsFactors = FALSE)
    return(.finalize_supra(empty, method_used = method, weight_tag = edgeWeight))
  }

  # post-filter (supra only)
  size_by <- tapply(res$actor, res$com, function(v) length(unique(v)))
  span_by <- tapply(res$layer, res$com, function(v) length(unique(v)))
  keep <- names(size_by)[as.integer(size_by) >= min.actors & as.integer(span_by) >= min.layers]
  res <- res[res$com %in% keep, , drop = FALSE]

  if (!length(keep) || !nrow(res)) {
    empty <- data.frame(actor = character(), com = character(), layer = character(), method = character(),
                        stringsAsFactors = FALSE)
    return(.finalize_supra(empty, method_used = method, weight_tag = edgeWeight))
  }

  # relabel (supra only)
  if (isTRUE(relabel_by_size)) {
    size2 <- tapply(res$actor, res$com, function(v) length(unique(v)))
    ord <- order(as.integer(size2), decreasing = TRUE)
    old_ids <- names(size2)[ord]
    new_ids <- paste0("C", seq_along(old_ids))
    map <- stats::setNames(new_ids, old_ids)
    res$com <- unname(map[res$com])
  }

  # standardize column order/types
  res <- res[, c("actor","com","layer","method"), drop = FALSE]
  res$actor  <- as.character(res$actor)
  res$com    <- as.character(res$com)
  res$layer  <- as.character(res$layer)
  res$method <- as.character(res$method)

  .finalize_supra(res, method_used = method, weight_tag = edgeWeight)
}

#------------------------------------------------------
# 9 - Add network-level attributes
#------------------------------------------------------

#' Attach node (actor) and edge attributes to a multilayer network, with match reporting
#'
#' @description
#' Adds external metadata to an existing multilayer network (\code{multinet::ml.network}).
#' The function can attach:
#' \itemize{
#'   \item \strong{actor (node) attributes} from \code{nodesMetadata}, and/or
#'   \item \strong{edge attributes} from \code{edgesMetadata}.
#' }
#' It is designed to be robust to common identifier formatting issues (e.g., case differences,
#' version suffixes, punctuation) by supporting simple ID normalisation steps and an optional
#' auto-detection procedure that selects the normalisation strategy with the highest match coverage.
#' A compact report of match rates and attached attributes can be saved to disk.
#'
#' @details
#' \section{What is the input to this function?}{
#' The main input is an existing \code{multinet::ml.network} object (\code{net}).
#'
#' Actor attributes are supplied via \code{nodesMetadata}, a data frame in which
#' \code{nodesMetadata[[featureID_col]]} contains the IDs that should match actors in the network.
#'
#' Edge attributes are supplied via \code{edgesMetadata}, which can be either:
#' \itemize{
#'   \item a \strong{data frame} containing \code{from}, \code{to}, a layer column, and attribute columns; or
#'   \item a \strong{named list} of per-layer data frames (one data frame per layer).
#' }
#' If \code{map_edges_by_actor_key = TRUE}, edge endpoints in \code{edgesMetadata} are first mapped to
#' network actors using the same ID key logic used for node matching.
#' }
#'
#' \section{How matching and attachment are performed (simple overview)}{
#' \enumerate{
#'   \item \strong{Extract network structure:} edges and actors are read from \code{net}. If the network
#'         provides two layer columns (\code{layer1}, \code{layer2}), only intra-layer edges are used
#'         for attachment (\code{layer1 == layer2}).
#'   \item \strong{Match actors:} when \code{nodesMetadata} is provided, actor IDs are normalised
#'         (via \code{actor_key_fun} or \code{actor_key_normalize}) and matched to metadata IDs. When
#'         \code{auto_detect_keys = TRUE}, multiple normalisation recipes are tried and the best
#'         (highest coverage) is selected.
#'   \item \strong{Attach actor attributes:} selected columns (\code{nodeAttrCols}) are stored as actor
#'         attributes in the multilayer network; missing values can be filled when
#'         \code{fillMissingNodes = TRUE}.
#'   \item \strong{Attach edge attributes:} when \code{edgesMetadata} is provided, endpoints are mapped
#'         to network actors, unmappable rows are dropped, and duplicate rows per edge are aggregated
#'         (controlled by \code{edgeAggregate}). For undirected networks, endpoints are canonicalised so
#'         \code{A–B} and \code{B–A} map to the same edge.
#'   \item \strong{Report:} a concise summary is attached to the returned object and can be saved to disk.
#' }
#' }
#'
#' \section{What files does this function write?}{
#' File outputs are optional and controlled by \code{save_report} and \code{export_edge_join_csv}.
#' When enabled, files are written under \code{results_dir}:
#' \itemize{
#'   \item \strong{Report (machine-readable):} \code{addattrs_report_<timestamp>.rds}
#'   \item \strong{Report (human-readable):} \code{addattrs_report_<timestamp>.txt}
#'   \item \strong{Optional edge-join table (debugging):} if \code{export_edge_join_csv = TRUE},
#'         the final edge join table is saved as \code{<edge_join_prefix>_<timestamp>.csv}.
#' }
#' If \code{save_report = FALSE} and \code{export_edge_join_csv = FALSE}, the function writes no files.
#' }
#'
#' @param net A \code{multinet::ml.network} object to which attributes will be attached.
#' @param nodesMetadata Optional data frame of actor (node) metadata.
#' @param featureID_col Column name in \code{nodesMetadata} containing the actor identifier used for
#'   matching (default \code{"feature_id"}).
#' @param nodeAttrCols Columns from \code{nodesMetadata} to attach as actor attributes. If \code{NULL},
#'   all columns except \code{featureID_col} are attached.
#' @param fillMissingNodes Logical; if \code{TRUE}, missing node attributes are filled using
#'   \code{nodeFillValue}.
#' @param nodeFillValue Fill value used when \code{fillMissingNodes = TRUE}.
#' @param actor_key_fun Optional function that transforms IDs before matching (overrides
#'   \code{actor_key_normalize}).
#' @param actor_key_normalize Optional character vector of built-in normalisation steps applied to IDs
#'   before matching (e.g., \code{c("strip_version","trim","tolower")}).
#' @param auto_detect_keys Logical; if \code{TRUE}, the function tests multiple normalisation recipes and
#'   uses the one that yields the highest match coverage.
#' @param min_match_warn If match coverage is below this fraction, the function prints examples of
#'   unmatched actors.
#' @param report_n_examples Maximum number of unmatched actor examples to show in messages/reports.
#' @param map_edges_by_actor_key Logical; if \code{TRUE}, edge endpoints in \code{edgesMetadata} are
#'   normalised/mapped using the actor key logic before joining to network edges.
#' @param edgesMetadata Optional edge metadata, either a data frame (with endpoints + layer column) or a
#'   named list of per-layer data frames.
#' @param edge_from,edge_to Column names for endpoints in \code{edgesMetadata}.
#' @param edge_layer_col Layer column name for data-frame \code{edgesMetadata}. If \code{NULL}, common
#'   names are auto-detected.
#' @param edgeAttrCols Edge metadata columns to attach. If \code{NULL}, all columns except endpoints and
#'   (when present) the layer column are attached.
#' @param edgeAggregate How to combine duplicate metadata rows for the same edge. One of
#'   \code{"first"}, \code{"mean"}, \code{"median"}, \code{"sum"}.
#' @param directed_network Logical; controls whether edges are treated as directed or undirected when
#'   matching edge metadata.
#' @param fillMissingEdges Logical; if \code{TRUE}, missing edge attribute values are filled before
#'   attachment.
#' @param edgeFillValue Fill value used when \code{fillMissingEdges = TRUE}.
#' @param overwrite_existing Logical; if \code{FALSE}, existing attributes are not re-declared, though
#'   values for specified actors/edges may still be set.
#' @param verbose Logical; print progress messages and file locations.
#' @param results_dir Output directory for reports and optional debug CSV.
#' @param save_report Logical; if \code{TRUE}, write the RDS and TXT report files to \code{results_dir}.
#' @param export_edge_join_csv Logical; if \code{TRUE}, export the final edge join table as CSV for debugging.
#' @param edge_join_prefix Prefix used for the edge join CSV filename.
#'
#' @return
#' The modified \code{multinet::ml.network} (returned invisibly). The object also receives two
#' attached R attributes:
#' \describe{
#'   \item{\code{attr(net, "actor_match_report")}}{Summary of actor matching coverage and unmatched examples.}
#'   \item{\code{attr(net, "edge_attach_report")}}{Summary of edge-attribute attachment and dropped rows.}
#' }
#'
#' @export
add_network_attributes <- function(
    net,
    nodesMetadata       = NULL,
    featureID_col       = "feature_id",
    nodeAttrCols        = NULL,
    fillMissingNodes    = TRUE,
    nodeFillValue       = "Unknown",
    actor_key_fun       = NULL,
    actor_key_normalize = NULL,
    auto_detect_keys    = TRUE,
    min_match_warn      = 0.7,
    report_n_examples   = 10,
    map_edges_by_actor_key = TRUE,
    edgesMetadata       = NULL,
    edge_from           = "from",
    edge_to             = "to",
    edge_layer_col      = NULL,
    edgeAttrCols        = NULL,
    edgeAggregate       = c("first","mean","median","sum"),
    directed_network    = FALSE,
    fillMissingEdges    = FALSE,
    edgeFillValue       = NA,
    overwrite_existing  = TRUE,
    verbose             = TRUE,
    results_dir         = getOption("mlnet.results_dir","omicsDNA_results"),
    save_report         = TRUE,
    export_edge_join_csv = FALSE,
    edge_join_prefix     = "edge_attributes_join"
) {
  edgeAggregate <- match.arg(edgeAggregate)

  # ---------- helpers ----------
  .collapse_fun <- function(which) switch(which,
                                          mean   = function(z) mean(z, na.rm = TRUE),
                                          median = function(z) stats::median(z, na.rm = TRUE),
                                          sum    = function(z) sum(z, na.rm = TRUE),
                                          first  = function(z) z[1]
  )
  .infer_type <- function(v) if (is.numeric(v)) "numeric" else "string"
  .ensure_dir <- function(d) if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  .is_abs     <- function(p) grepl("^(/|[A-Za-z]:[\\/])", p)
  .pick <- function(cands, nms) { z <- cands[cands %in% nms]; if (length(z)) z[1] else NA_character_ }

  .norm_catalog <- list(
    identity      = function(x) x,
    trim          = function(x) trimws(x),
    strip_version = function(x) sub("\\.\\d+$", "", x),
    tolower       = function(x) tolower(x),
    toupper       = function(x) toupper(x),
    rm_dash       = function(x) gsub("[-_]", "", x),
    rm_punct      = function(x) gsub("[[:punct:]]+", "", x),
    alnum         = function(x) gsub("[^[:alnum:]]+", "", x)
  )
  .compose <- function(steps) {
    if (is.null(steps) || !length(steps)) return(.norm_catalog$identity)
    funs <- lapply(steps, function(s) {
      if (!s %in% names(.norm_catalog)) stop("Unknown normalizer step: ", s)
      .norm_catalog[[s]]
    })
    function(x) { x <- as.character(x); for (f in funs) x <- f(x); x }
  }
  .frac_match <- function(a, b) { if (!length(a)) return(0); mean(as.character(a) %in% as.character(b)) }

  # ---------- pull edges and actors from net ----------
  E_raw <- try(multinet::edges_ml(net), silent = TRUE)
  if (inherits(E_raw, "try-error") || is.null(E_raw)) stop("Could not retrieve edges via multinet::edges_ml(net).")
  Edf <- if (is.data.frame(E_raw)) E_raw else as.data.frame(E_raw, stringsAsFactors = FALSE)

  a1 <- .pick(c("actor1","from_actor","from","source","node1","i","v1"), names(Edf))
  a2 <- .pick(c("actor2","to_actor","to","target","node2","j","v2"),     names(Edf))
  l1 <- .pick(c("layer1","from_layer","layer","Layer","l1"),             names(Edf))
  l2 <- .pick(c("layer2","to_layer","l2","layer","Layer"),               names(Edf))

  if (any(is.na(c(a1, a2, l1)))) {
    stop("Could not recognize edge columns in `edges_ml(net)`. Found: ", paste(names(Edf), collapse = ", "))
  }

  # Keep intra-layer only if a l2 exists
  if (!is.na(l2) && (l2 %in% names(Edf))) {
    Edf <- Edf[Edf[[l1]] == Edf[[l2]], , drop = FALSE]
  }
  Edf$layer <- as.character(Edf[[l1]])
  e_from <- as.character(Edf[[a1]])
  e_to   <- as.character(Edf[[a2]])
  Edf$u <- if (directed_network) e_from else pmin(e_from, e_to)
  Edf$v <- if (directed_network) e_to   else pmax(e_from, e_to)

  # Actors present in the network
  actors_tbl <- try(multinet::actors_ml(net), silent = TRUE)
  actors <- if (!inherits(actors_tbl, "try-error") && !is.null(actors_tbl)) {
    # try to pick a reasonable column
    nms <- names(actors_tbl)
    cand <- .pick(c("actor","actors","name","vertex","node","id"), nms)
    if (is.na(cand)) {
      w <- which(vapply(actors_tbl, function(col) is.character(col) || is.factor(col), logical(1)))
      if (length(w)) unique(as.character(actors_tbl[[w[1]]])) else character()
    } else unique(as.character(actors_tbl[[cand]]))
  } else unique(c(Edf$u, Edf$v))

  # ---------- attach ACTOR attributes ----------
  actor_report <- NULL
  key_fun_used <- NULL

  if (!is.null(nodesMetadata)) {
    if (!is.data.frame(nodesMetadata)) stop("`nodesMetadata` must be a data.frame.")
    if (!(featureID_col %in% names(nodesMetadata))) stop("`nodesMetadata` must contain key column: ", featureID_col)

    key_raw <- as.character(nodesMetadata[[featureID_col]])
    if (anyDuplicated(key_raw)) {
      nodesMetadata <- nodesMetadata[!duplicated(key_raw), , drop = FALSE]
      key_raw <- as.character(nodesMetadata[[featureID_col]])
    }

    # Build key function
    if (is.function(actor_key_fun)) {
      key_fun <- function(x) actor_key_fun(as.character(x)); key_fun_used <- "user-supplied function"
    } else {
      if (is.null(actor_key_normalize)) { key_fun <- .compose("identity"); key_fun_used <- "identity" }
      else { key_fun <- .compose(actor_key_normalize); key_fun_used <- paste(actor_key_normalize, collapse = " -> ") }
    }

    # Try auto-detect recipes if coverage looks poor
    coverage0 <- .frac_match(key_fun(actors), key_fun(key_raw))
    if (isTRUE(auto_detect_keys) && (!is.finite(coverage0) || coverage0 < min_match_warn)) {
      candidates <- list(
        c("identity"),
        c("trim"),
        c("strip_version"),
        c("strip_version","trim"),
        c("strip_version","tolower"),
        c("trim","strip_version","tolower"),
        c("tolower"),
        c("toupper")
      )
      best_cov <- coverage0; best_fun <- key_fun; best_desc <- key_fun_used
      for (steps in candidates) {
        f <- .compose(steps)
        cov <- .frac_match(f(actors), f(key_raw))
        if (is.finite(cov) && cov > best_cov) { best_cov <- cov; best_fun <- f; best_desc <- paste(steps, collapse = " -> ") }
      }
      key_fun <- best_fun; key_fun_used <- best_desc
    }

    # Build keyed tables
    A <- data.frame(actor = actors, key = key_fun(actors), stringsAsFactors = FALSE)
    A <- A[!duplicated(A$key), , drop = FALSE]

    if (is.null(nodeAttrCols)) nodeAttrCols <- setdiff(names(nodesMetadata), featureID_col)
    if (length(nodeAttrCols)) {
      meta_keyed <- cbind(
        data.frame(key = key_fun(key_raw), stringsAsFactors = FALSE),
        nodesMetadata[, nodeAttrCols, drop = FALSE]
      )
      merged <- merge(A, meta_keyed, by = "key", all.x = TRUE, sort = FALSE)
      rownames(merged) <- NULL

      # coverage
      row_has <- if (length(nodeAttrCols) == 1L) !is.na(merged[[nodeAttrCols[1]]])
      else rowSums(!is.na(merged[, nodeAttrCols, drop = FALSE])) > 0
      coverage <- if (length(row_has)) mean(row_has) else 0
      matched  <- sum(row_has); total <- nrow(merged)

      if (fillMissingNodes) {
        for (a in nodeAttrCols) {
          nas <- is.na(merged[[a]]) | (is.character(merged[[a]]) & merged[[a]] == "")
          if (!any(nas)) next
          if (is.numeric(merged[[a]])) {
            if (is.numeric(nodeFillValue)) merged[[a]][nas] <- as.numeric(nodeFillValue)
          } else merged[[a]][nas] <- as.character(nodeFillValue)
        }
      }

      # declare + set attributes
      existing <- try(multinet::attributes_ml(net), silent = TRUE)
      existing_names <- if (!inherits(existing, "try-error") && !is.null(existing)) as.character(existing$name) else character()

      for (a in nodeAttrCols) {
        if (overwrite_existing || !(a %in% existing_names)) {
          multinet::add_attributes_ml(net, attributes = a, target = "actor", type = .infer_type(merged[[a]]))
        }
        multinet::set_values_ml(net, attribute = a, actors = merged$actor, values = merged[[a]])
      }

      if (verbose) {
        msg <- sprintf("Actor attributes: matched %.1f%% (%d / %d) via: %s",
                       100 * (if (is.finite(coverage)) coverage else 0), matched, total, key_fun_used)
        message(msg)
        if (!is.finite(coverage) || coverage < min_match_warn) {
          bad <- merged$actor[!row_has]
          nshow <- min(report_n_examples, length(bad))
          if (nshow > 0) message("Unmatched actors (", nshow, " examples): ", paste(utils::head(bad, nshow), collapse = ", "))
        }
      }

      actor_report <- list(
        coverage = coverage, matched = matched, total = total,
        normalizer = key_fun_used,
        unmatched_examples = utils::head(merged$actor[!row_has], report_n_examples)
      )
    } else if (verbose) {
      message("No actor attributes to attach (no columns besides `", featureID_col, "`).")
    }
  }

  # ---------- attach EDGE attributes ----------
  edge_report <- NULL
  if (!is.null(edgesMetadata)) {
    # Build a normalized metadata table A_meta: (layer, u, v, <attrs>)
    # Helper to make (layer,u,v) from a df (either per-layer df or a flat df with layer column)
    # mapping endpoints through the same key function if requested.
    if (exists("key_fun")) {
      key_fun_final <- key_fun
    } else {
      # Default identity mapping when no actor key work was done
      key_fun_final <- .compose("identity")
      key_fun_used  <- "identity"
    }

    actor_map <- data.frame(actor = actors, key = key_fun_final(actors), stringsAsFactors = FALSE)
    actor_map <- actor_map[!duplicated(actor_map$key), , drop = FALSE]

    .build_attr_df <- function(df, layer_name = NULL, layer_col = edge_layer_col) {
      if (!all(c(edge_from, edge_to) %in% names(df))) {
        stop("Edge metadata is missing endpoint columns `", edge_from, "` and/or `", edge_to, "`.")
      }
      attrs <- edgeAttrCols
      if (is.null(attrs)) attrs <- setdiff(names(df), c(edge_from, edge_to, layer_col))
      if (!length(attrs)) stop("No `edgeAttrCols` to attach in edge metadata.")

      from_ext <- as.character(df[[edge_from]])
      to_ext   <- as.character(df[[edge_to]])

      if (isTRUE(map_edges_by_actor_key)) {
        # map external endpoints -> actor names in the network
        k_from <- key_fun_final(from_ext); k_to <- key_fun_final(to_ext)
        mfrom <- merge(data.frame(key = k_from, idx = seq_along(k_from), stringsAsFactors = FALSE),
                       actor_map, by = "key", all.x = TRUE, sort = FALSE)
        mto   <- merge(data.frame(key = k_to,   idx = seq_along(k_to),   stringsAsFactors = FALSE),
                       actor_map, by = "key", all.x = TRUE, sort = FALSE)
        from <- mfrom$actor[order(mfrom$idx)]
        to   <- mto$actor[order(mto$idx)]
      } else {
        from <- from_ext; to <- to_ext
      }

      uu <- if (directed_network) from else pmin(from, to)
      vv <- if (directed_network) to   else pmax(from, to)

      out <- data.frame(
        layer = if (!is.null(layer_name)) layer_name else as.character(df[[layer_col]]),
        u = uu, v = vv,
        stringsAsFactors = FALSE
      )
      cbind(out, df[, attrs, drop = FALSE])
    }

    if (is.list(edgesMetadata)) {
      nm <- names(edgesMetadata); if (is.null(nm)) nm <- as.character(seq_along(edgesMetadata))
      pieces <- vector("list", length(edgesMetadata))
      for (i in seq_along(edgesMetadata)) {
        pieces[[i]] <- .build_attr_df(edgesMetadata[[i]], layer_name = nm[i], layer_col = NULL)
      }
      A_meta <- do.call(rbind, pieces)
    } else if (is.data.frame(edgesMetadata)) {
      if (is.null(edge_layer_col)) {
        edge_layer_col <- .pick(c("layer","Layer","group","Group"), names(edgesMetadata))
        if (is.na(edge_layer_col)) stop("Could not detect layer column in `edgesMetadata`; specify `edge_layer_col`.")
      }
      A_meta <- .build_attr_df(edgesMetadata, layer_col = edge_layer_col)
    } else {
      stop("`edgesMetadata` must be a named list of data frames or a data frame.")
    }

    # Drop rows that failed to map (NA endpoints)
    bad_endpoints <- is.na(A_meta$u) | is.na(A_meta$v)
    dropped_meta  <- sum(bad_endpoints)
    if (dropped_meta && verbose) message("Edge metadata: dropped ", dropped_meta, " row(s) with unmapped endpoints.")
    A_meta <- A_meta[!bad_endpoints, , drop = FALSE]

    # Aggregate duplicates within metadata (layer,u,v)
    if (nrow(A_meta)) {
      if (edgeAggregate == "first") {
        A_meta <- A_meta[!duplicated(A_meta[, c("layer","u","v")]), , drop = FALSE]
      } else {
        FUN  <- .collapse_fun(edgeAggregate)
        keep <- setdiff(names(A_meta), c("layer","u","v"))
        A_meta <- stats::aggregate(. ~ layer + u + v,
                                   data = A_meta[, c("layer","u","v", keep), drop = FALSE],
                                   FUN  = function(col) if (is.numeric(col)) FUN(col) else col[1]
        )
      }
    }

    # Deterministic left join: ensure the joined table preserves Edf order
    join_cols <- c("layer","u","v")
    A_meta <- A_meta[order(A_meta$layer, A_meta$u, A_meta$v), , drop = FALSE]
    Edf_ord <- Edf[, join_cols, drop = FALSE]
    Edf_ord$.row_id <- seq_len(nrow(Edf_ord))
    join <- merge(Edf_ord, A_meta, by = join_cols, all.x = TRUE, sort = FALSE)
    join <- join[order(join$.row_id), , drop = FALSE]
    rownames(join) <- NULL
    join$.row_id <- NULL

    # Which attributes do we actually set?
    attr_cols <- setdiff(names(join), join_cols)
    if (length(attr_cols) == 0L) {
      if (verbose) message("No edge attributes to attach.")
    } else {
      # Fill missing edge attrs (optional)
      if (isTRUE(fillMissingEdges)) {
        for (a in attr_cols) {
          nas <- is.na(join[[a]])
          if (!any(nas)) next
          if (is.numeric(join[[a]])) join[[a]][nas] <- if (is.numeric(edgeFillValue)) as.numeric(edgeFillValue) else NA_real_
          else                       join[[a]][nas] <- as.character(edgeFillValue)
        }
      }

      # Declare attributes (edge-level)
      existing <- try(multinet::attributes_ml(net), silent = TRUE)
      existing_names <- if (!inherits(existing, "try-error") && !is.null(existing)) as.character(existing$name) else character()
      add_formals <- names(formals(multinet::add_attributes_ml))

      .declare_edge_attr <- function(attr, type, layers_vec) {
        if ("layers" %in% add_formals) {
          ok <- try(multinet::add_attributes_ml(net, attributes = attr, target = "edge", type = type, layers = layers_vec),
                    silent = TRUE)
          if (inherits(ok, "try-error")) {
            for (L in layers_vec) multinet::add_attributes_ml(net, attributes = attr, target = "edge", type = type, layers = L)
          }
        } else if ("layer" %in% add_formals) {
          for (L in layers_vec) multinet::add_attributes_ml(net, attributes = attr, target = "edge", type = type, layer = L)
        } else {
          multinet::add_attributes_ml(net, attributes = attr, target = "edge", type = type)
        }
      }

      layers_vec <- unique(Edf$layer)
      for (a in attr_cols) {
        if (overwrite_existing || !(a %in% existing_names)) .declare_edge_attr(a, .infer_type(join[[a]]), layers_vec)
      }

      # Build edge spec in the exact row order of Edf (to match `join`)
      if (!is.na(l2) && (l2 %in% names(Edf))) {
        edge_spec <- data.frame(actor1 = Edf[[a1]], layer1 = Edf[[l1]],
                                actor2 = Edf[[a2]], layer2 = Edf[[l2]],
                                stringsAsFactors = FALSE)
      } else {
        edge_spec <- data.frame(actor1 = Edf[[a1]], layer1 = Edf[[l1]],
                                actor2 = Edf[[a2]], layer2 = Edf[[l1]],
                                stringsAsFactors = FALSE)
      }

      # Set values per attribute (skip NAs)
      for (a in attr_cols) {
        vals <- join[[a]]
        keep <- !is.na(vals)
        if (any(keep)) multinet::set_values_ml(net, attribute = a, edges = edge_spec[keep, , drop = FALSE], values = vals[keep])
      }

      if (verbose) cat("Attached", length(attr_cols), "edge attribute(s) to", nrow(join), "edges.\n")
    }

    edge_report <- list(
      n_edges_in_net   = nrow(Edf),
      n_meta_rows_in   = if (!exists("A_meta")) 0L else nrow(A_meta),
      n_meta_dropped   = dropped_meta,
      attrs_attached   = if (exists("attr_cols")) attr_cols else character()
    )

    # Optional CSV of the join view (for debugging)
    if (isTRUE(export_edge_join_csv) && length(attr_cols)) {
      .ensure_dir(results_dir)
      stamp <- format(Sys.time(), "%Y-%m-%d_%H%M%S")
      csv_file <- file.path(results_dir, sprintf("%s_%s.csv", edge_join_prefix, stamp))
      utils::write.csv(join, csv_file, row.names = FALSE)
      if (verbose) message("Saved edge-join CSV: ", normalizePath(csv_file, FALSE))
      edge_report$edge_join_csv <- csv_file
    }
  }

  # ---------- final report & save ----------
  if (verbose) {
    cat("Attributes present in network:\n")
    print(try(multinet::attributes_ml(net), silent = TRUE))
  }

  report <- list(
    actor_match_report = actor_report,
    edge_attach_report = edge_report,
    timestamp          = Sys.time()
  )
  attr(net, "actor_match_report") <- actor_report
  attr(net, "edge_attach_report") <- edge_report

  if (isTRUE(save_report)) {
    .ensure_dir(results_dir)
    stamp <- format(Sys.time(), "%Y-%m-%d_%H%M%S")
    rds_file <- file.path(results_dir, sprintf("addattrs_report_%s.rds", stamp))
    txt_file <- file.path(results_dir, sprintf("addattrs_report_%s.txt", stamp))
    saveRDS(report, rds_file)
    # small human-readable summary
    cat(
      sprintf("Actor report: %s\n", if (is.null(actor_report)) "<none>" else
        sprintf("coverage=%.1f%% matched=%s/%s normalizer=%s",
                100*(actor_report$coverage %||% 0),
                actor_report$matched %||% 0, actor_report$total %||% 0,
                actor_report$normalizer %||% "NA")),
      sprintf("Edge report:  %s\n", if (is.null(edge_report)) "<none>" else
        sprintf("n_edges=%s n_meta_in=%s dropped_meta=%s attrs=[%s]",
                edge_report$n_edges_in_net %||% 0,
                edge_report$n_meta_rows_in %||% 0,
                edge_report$n_meta_dropped %||% 0,
                paste(edge_report$attrs_attached %||% character(), collapse = ", "))),
      file = txt_file, sep = ""
    )
    if (verbose) {
      message("Saved attribute attach report RDS: ", normalizePath(rds_file, FALSE))
      message("Saved attribute attach report TXT: ", normalizePath(txt_file, FALSE))
    }
  }

  invisible(net)
}


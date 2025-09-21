
#------------------------------------------------------
# 9 - Add network-level attributes
#------------------------------------------------------

#' Attach actor- and edge-level attributes to a multilayer network, with coverage reporting
#'
#' @description
#' This procedure enriches an existing `multinet::ml.network` with **actor** (node)
#' and **edge** attributes drawn from external metadata tables. It is designed to
#' be robust to common ID-format differences by offering simple, composable
#' normalisers (e.g., `tolower`, `strip_version`, `trim`) and an optional
#' auto-detection step that chooses the recipe yielding the highest coverage.
#' A compact report (RDS + TXT) is generated to document match rates and the
#' attributes attached.
#'
#' @details
#' **What the function does**
#' 1. Extracts the edges and actors from `net`. If the network stores two layer
#'    columns (`layer1`, `layer2`), edges are restricted to **intra-layer**
#'    interactions (i.e., `layer1 == layer2`).
#' 2. If `nodesMetadata` is provided, actor IDs in `net` are reconciled against
#'    `nodesMetadata[[featureID_col]]` via either a user-supplied `actor_key_fun`
#'    or a pipeline of built-in normalisers listed in `actor_key_normalize`.
#'    When `auto_detect_keys = TRUE`, several normaliser pipelines are tried and
#'    the one with the highest actor-coverage is used. Selected columns
#'    (`nodeAttrCols`, or all except `featureID_col` when `NULL`) are then
#'    attached as actor attributes. Missing values can be filled with
#'    `nodeFillValue` when `fillMissingNodes = TRUE`.
#' 3. If `edgesMetadata` is supplied, the external endpoints are mapped to actors
#'    in `net` using the same key function (when `map_edges_by_actor_key = TRUE`).
#'    Rows that cannot be mapped to a pair of actors are **dropped** (with a
#'    message). When multiple metadata rows correspond to the same `(layer, u, v)`
#'    edge, they are aggregated according to `edgeAggregate`. Numeric columns are
#'    collapsed using the chosen summary (e.g., `mean`), while non-numeric
#'    columns take the **first** value.
#' 4. Attributes are declared (if needed) and then values are written to the
#'    network. If attributes already exist and `overwrite_existing = FALSE`,
#'    the function avoids re-declaration but **still sets values** for the
#'    affected actors/edges.
#' 5. A small report is attached to the returned object (as attributes) and,
#'    optionally, saved to disk under `results_dir`.
#'
#' **Normalisers for ID reconciliation**
#' Built-in normaliser steps include: `"identity"`, `"trim"`, `"strip_version"`,
#' `"tolower"`, `"toupper"`, `"rm_dash"` (remove `-` and `_`), `"rm_punct"`
#' (remove punctuation), and `"alnum"` (keep alphanumerics only). You may chain
#' them by supplying a character vector to `actor_key_normalize`; for complete
#' control, pass a custom `actor_key_fun`.
#'
#' **Directedness and canonicalisation**
#' For the purpose of **joining** edge metadata to network edges, endpoints are
#' canonicalised to `(u, v) = (min, max)` when `directed_network = FALSE`; they
#' are left as `(from, to)` when `directed_network = TRUE`.
#'
#' **Outputs and reporting**
#' The modified `ml.network` is returned (invisibly). Two lists are also attached
#' to it as attributes:
#' - `attr(net, "actor_match_report")`: coverage, chosen normaliser, and examples
#'   of unmatched actors (up to `report_n_examples`).
#' - `attr(net, "edge_attach_report")`: edge counts, number of metadata rows
#'   dropped due to unmapped endpoints, and the names of attached edge attributes.
#' If `save_report = TRUE`, both a machine-readable RDS and a human-readable TXT
#' summary are written to `results_dir`.
#'
#' @param net A `multinet::ml.network` object to be enriched.
#' @param nodesMetadata Optional data frame holding actor attributes.
#' @param featureID_col Name of the identifier column in `nodesMetadata` that
#'   matches actors in `net`. Default `"feature_id"`.
#' @param nodeAttrCols Character vector of columns from `nodesMetadata` to attach
#'   as actor attributes. If `NULL` (default) and `nodesMetadata` is provided,
#'   **all** columns except `featureID_col` are attached.
#' @param fillMissingNodes Logical; if `TRUE`, missing actor attributes are
#'   populated with `nodeFillValue`. Default `TRUE`.
#' @param nodeFillValue Scalar (character or numeric) used when
#'   `fillMissingNodes = TRUE`. Numeric columns are filled only if
#'   `nodeFillValue` is numeric.
#' @param actor_key_fun Optional function used to normalise IDs before matching;
#'   if provided, this overrides `actor_key_normalize`.
#' @param actor_key_normalize Character vector specifying a sequence of built-in
#'   normalisers to apply (see Details). Ignored if `actor_key_fun` is supplied.
#' @param auto_detect_keys Logical; if `TRUE`, several candidate normaliser
#'   sequences are tried and the one with highest actor-coverage is used.
#'   Default `TRUE`.
#' @param min_match_warn Numeric in `[0,1]`; if observed coverage is below this
#'   threshold, the function prints examples of unmatched actors. Default `0.7`.
#' @param report_n_examples Integer; maximum number of unmatched actor examples to
#'   include in messages and the saved report. Default `10`.
#' @param map_edges_by_actor_key Logical; if `TRUE`, the same ID key function is
#'   applied to `edgesMetadata` endpoints before matching them to actors in `net`.
#'   Default `TRUE`.
#' @param edgesMetadata Edge metadata as either (i) a data frame containing
#'   endpoints and a layer column, or (ii) a **named list** of per-layer data
#'   frames. In either case, endpoint columns must be present (see `edge_from`,
#'   `edge_to`). Rows that cannot be mapped to actors are dropped (reported).
#' @param edge_from,edge_to Character scalars giving the endpoint column names in
#'   `edgesMetadata`. Defaults `"from"` and `"to"`.
#' @param edge_layer_col For data-frame `edgesMetadata` only: name of the layer
#'   column. If `NULL`, the function tries to auto-detect from
#'   `c("layer","Layer","group","Group")`.
#' @param edgeAttrCols Character vector naming which columns in `edgesMetadata`
#'   to attach as edge attributes. If `NULL`, all columns except endpoints and
#'   (when present) the layer column are used.
#' @param edgeAggregate How to aggregate duplicate metadata rows per `(layer, u, v)`
#'   before attachment. One of `"first"`, `"mean"`, `"median"`, `"sum"`.
#'   Numeric columns are combined using the selected summary; non-numeric columns
#'   always take the first value. Default `"first"`.
#' @param directed_network Logical; whether the underlying network is treated as
#'   directed when reconciling edge metadata (`u = from`, `v = to`) versus
#'   undirected (`u = pmin(from, to)`, `v = pmax(from, to)`). Default `FALSE`.
#' @param fillMissingEdges Logical; if `TRUE`, missing edge attributes are filled
#'   with `edgeFillValue` prior to attachment. Default `FALSE`.
#' @param edgeFillValue Scalar used when `fillMissingEdges = TRUE`. Numeric
#'   columns are filled only if `edgeFillValue` is numeric.
#' @param overwrite_existing Logical; if `FALSE`, existing attributes are not
#'   **re-declared** on the network; however, the function will still **set
#'   values** for the specified actors/edges. Default `TRUE`.
#' @param verbose Logical; print progress messages, coverage summaries, and file
#'   locations. Default `TRUE`.
#' @param results_dir Directory where reports (and optional CSV) are saved.
#'   Default `getOption("mlnet.results_dir", "omicsDNA_results")`.
#' @param save_report Logical; if `TRUE`, save a compact report (RDS + TXT) to
#'   `results_dir`. Default `TRUE`.
#' @param export_edge_join_csv Logical; if `TRUE`, save the final edge-join table
#'   (used to attach attributes) as a CSV under `results_dir` with basename
#'   `edge_join_prefix`. Useful for debugging. Default `FALSE`.
#' @param edge_join_prefix Basename prefix used when `export_edge_join_csv = TRUE`.
#'
#' @return The modified `multinet::ml.network` (returned invisibly). The object
#'   also gains two attributes, `attr(net, "actor_match_report")` and
#'   `attr(net, "edge_attach_report")`, containing concise summaries of the
#'   attachment process.
#'
#' @section Practical notes:
#' - Only **intra-layer** edges are considered when the network stores two layer
#'   columns; cross-layer edges are ignored during attribute attachment.
#' - When mapping edge metadata, rows with endpoints that cannot be matched to
#'   actors in `net` are dropped (the count is reported).
#' - Attribute types are declared as `"numeric"` for numeric vectors and
#'   `"string"` otherwise; this affects how `multinet` stores values internally.
#'
#' @examples
#' \dontrun{
#' # Attach actor attributes (using a simple normaliser pipeline)
#' net <- add_network_attributes(
#'   net,
#'   nodesMetadata       = genes_info,
#'   featureID_col       = "GeneName",
#'   nodeAttrCols        = "GeneType",
#'   actor_key_normalize = c("strip_version", "trim", "tolower"),
#'   auto_detect_keys    = FALSE
#' )
#'
#' # Attach edge attributes from a named list (one data frame per layer)
#' net <- add_network_attributes(
#'   net,
#'   edgesMetadata       = cons_list,
#'   edge_from           = "from",
#'   edge_to             = "to",
#'   edgeAttrCols        = c("n_present", "n_repeats", "prop_present"),
#'   directed_network    = FALSE,
#'   export_edge_join_csv = FALSE
#' )
#' }
#'
#' @seealso
#'   \code{\link{build_multiNet}} for assembling multilayer networks;
#'   \code{\link{edgesFromAdjacency}} and \code{\link{consensusEdges}} for
#'   upstream edge construction and consensus.
#'
#' @importFrom multinet edges_ml actors_ml attributes_ml add_attributes_ml set_values_ml
#' @importFrom utils write.csv head
#' @importFrom stats median
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


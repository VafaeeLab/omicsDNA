
# ------------------------------------------------------------------------------
# 10 - Print a report of attributes attached to a multilayer network
# ------------------------------------------------------------------------------

# A small infix helper used above (so we don't rely on any pkg):
`%||%` <- function(x, y) if (!is.null(x)) x else y

#' Print a report of attributes in a `multinet::ml.network`
#'
#' @description
#' Prints a human‑readable summary of:
#' - **Which attributes were added** (if available from the attachment reports
#'   stored by `add_network_attributes()`), and
#' - **What attributes are currently present** on **actors (nodes)**, **edges**,
#'   and (if supported by your multinet version) **network‑level** attributes.
#'
#' The function is defensive across multinet versions:
#' - It tries `attributes_ml(net, target = "actor"/"edge")`; if unavailable,
#'   it infers attribute names from `actors_ml(net, attributes = TRUE)` and
#'   `edges_ml(net, attributes = TRUE)` by subtracting the canonical ID columns.
#' - For edges, it recognizes many possible column names for endpoints and
#'   layer columns, and reports the **extra columns** as edge attributes.
#'
#' @param net A `multinet::ml.network` object with attributes attached.
#' @param show_examples Logical; also print a few example values for each class
#'   of attributes (actors/edges). Default `TRUE`.
#' @param show_max Integer; maximum number of names/rows to print in each
#'   section. Default `10`.
#' @param verbose Logical; print the report to the console. Default `TRUE`.
#'
#' @return (Invisibly) a list with elements:
#'   - `added` — what was added per the stored reports:
#'       - `actor_attrs` (may be `NULL` if not recorded),
#'       - `edge_attrs` (character vector, if recorded),
#'       - `reports` (the raw `actor_match_report` / `edge_attach_report`).
#'   - `present` — what is currently present on the object:
#'       - `actor` — data.frame with `name` and `type` (or inferred types),
#'       - `edge`  — data.frame with `name` and `type` (or inferred types),
#'       - `network` — data.frame of network‑level attributes (if available).
#'   - `samples` — small tibbles with example rows for actors/edges if requested.
#'
#' @examples
#' \dontrun{
#' # After you call add_network_attributes(...):
#' rep <- report_network_attributes(net, show_examples = TRUE, show_max = 8)
#'
#' # Programmatic access
#' rep$added$edge_attrs
#' rep$present$actor
#' rep$samples$edges
#' }
#' @export
network_attributes <- function(net, show_examples = TRUE, show_max = 10, verbose = TRUE) {
  # ---- small helpers ---------------------------------------------------------
  .pick <- function(cands, nms) {
    z <- cands[cands %in% nms]; if (length(z)) z[1] else NA_character_
  }
  .infer_type <- function(v) if (is.numeric(v)) "numeric" else "string"
  .print_names <- function(x, label, n = show_max) {
    x <- unique(x)
    cat(sprintf("%s (%d): ", label, length(x)))
    if (!length(x)) { cat("<none>\n"); return(invisible(NULL)) }
    cat(paste(utils::head(x, n), collapse = ", "))
    if (length(x) > n) cat(" ...")
    cat("\n")
  }
  .as_df <- function(x) {
    if (isTRUE(is.data.frame(x))) return(x)
    tryCatch(as.data.frame(x, stringsAsFactors = FALSE), error = function(e) NULL)
  }

  # ---- 1) What *was* added (from stored reports, if any) ---------------------
  actor_report <- attr(net, "actor_match_report")
  edge_report  <- attr(net, "edge_attach_report")

  added_actor_attrs <- NULL
  if (is.list(actor_report) && length(actor_report)) {
    # Not all versions stored the attribute names; print coverage/normaliser anyway
    if (verbose) {
      cat("Actor-attribute attachment report:\n")
      cov <- if (!is.null(actor_report$coverage) && is.finite(actor_report$coverage)) actor_report$coverage else NA_real_
      msg <- sprintf("  Coverage: %s   Matched: %s/%s   Normalizer: %s\n",
                     if (is.finite(cov)) sprintf("%.1f%%", 100 * cov) else "NA",
                     actor_report$matched %||% NA_integer_,
                     actor_report$total   %||% NA_integer_,
                     actor_report$normalizer %||% "NA")
      cat(msg)
      if (!is.null(actor_report$unmatched_examples) && length(actor_report$unmatched_examples)) {
        cat("  Unmatched examples: ",
            paste(utils::head(actor_report$unmatched_examples, show_max), collapse = ", "),
            if (length(actor_report$unmatched_examples) > show_max) " ..." else "",
            "\n", sep = "")
      }
    }
    # Try to read stored names if your earlier function added them; else NULL
    if (!is.null(actor_report$attrs_attached)) {
      added_actor_attrs <- as.character(actor_report$attrs_attached)
    }
  } else if (verbose) {
    cat("Actor-attribute attachment report: <none recorded on object>\n")
  }

  added_edge_attrs <- NULL
  if (is.list(edge_report) && length(edge_report)) {
    if (verbose) {
      cat("Edge-attribute attachment report:\n")
      cat(sprintf("  Edges in net: %s | Metadata rows in: %s | Dropped (unmapped): %s\n",
                  edge_report$n_edges_in_net %||% NA_integer_,
                  edge_report$n_meta_rows_in %||% NA_integer_,
                  edge_report$n_meta_dropped %||% NA_integer_))
    }
    if (!is.null(edge_report$attrs_attached)) {
      added_edge_attrs <- as.character(edge_report$attrs_attached)
      if (verbose) .print_names(added_edge_attrs, "  Edge attrs added")
    } else if (verbose) {
      cat("  Edge attrs added: <not recorded by this version>\n")
    }
    if (!is.null(edge_report$edge_join_csv) && verbose) {
      cat("  Edge-join CSV: ", edge_report$edge_join_csv, "\n", sep = "")
    }
  } else if (verbose) {
    cat("Edge-attribute attachment report: <none recorded on object>\n")
  }

  # ---- 2) What is *present now* ----------------------------------------------
  # 2a) Actor attributes
  actor_attrs <- NULL
  ok <- TRUE
  try_actor <- try(multinet::attributes_ml(net, target = "actor"), silent = TRUE)
  if (!inherits(try_actor, "try-error")) actor_attrs <- .as_df(try_actor)
  if (is.null(actor_attrs) || !nrow(actor_attrs)) {
    # Fallback: infer from actors_ml(..., attributes=TRUE)
    act_tbl <- try(multinet::actors_ml(net, attributes = TRUE), silent = TRUE)
    if (!inherits(act_tbl, "try-error") && is.data.frame(act_tbl) && ncol(act_tbl) >= 1) {
      actor_id_col <- names(act_tbl)[1]
      cand_attrs <- setdiff(names(act_tbl), actor_id_col)
      actor_attrs <- data.frame(
        name = cand_attrs,
        type = vapply(act_tbl[cand_attrs], .infer_type, character(1)),
        stringsAsFactors = FALSE
      )
    } else {
      ok <- FALSE
      actor_attrs <- data.frame(name = character(), type = character(), stringsAsFactors = FALSE)
    }
  }

  # 2b) Edge attributes
  edge_attrs <- NULL
  try_edge <- try(multinet::attributes_ml(net, target = "edge"), silent = TRUE)
  if (!inherits(try_edge, "try-error")) edge_attrs <- .as_df(try_edge)
  if (is.null(edge_attrs) || !nrow(edge_attrs)) {
    # Fallback: infer from edges_ml(..., attributes=TRUE)
    ed_tbl <- try(multinet::edges_ml(net, attributes = TRUE), silent = TRUE)
    if (!inherits(ed_tbl, "try-error") && is.data.frame(ed_tbl) && ncol(ed_tbl) >= 1) {
      nms <- names(ed_tbl)
      col_from  <- .pick(c("actor1","i","from","source","u","v1","node1","src"), nms)
      col_to    <- .pick(c("actor2","j","to","target","v","v2","node2","dst"), nms)
      col_l1    <- .pick(c("layer1","from_layer","layer","Layer","l1","L1"), nms)
      col_l2    <- .pick(c("layer2","to_layer","l2","L2"), nms)
      base_cols <- unique(na.omit(c(col_from, col_to, col_l1, col_l2)))
      base_cols <- unique(c(base_cols, "weight"))  # don't treat weight as metadata here
      cand_attrs <- setdiff(nms, base_cols)
      # crude type inference
      edge_attrs <- data.frame(
        name = cand_attrs,
        type = vapply(ed_tbl[cand_attrs], .infer_type, character(1)),
        stringsAsFactors = FALSE
      )
    } else {
      ok <- FALSE
      edge_attrs <- data.frame(name = character(), type = character(), stringsAsFactors = FALSE)
    }
  }

  # 2c) Network-level attributes (if supported)
  network_attrs <- NULL
  try_net <- try(multinet::attributes_ml(net, target = "network"), silent = TRUE)
  if (!inherits(try_net, "try-error")) network_attrs <- .as_df(try_net)
  if (is.null(network_attrs)) {
    # Some versions encode all in a single table; try to detect a 'target' column
    any_table <- try(multinet::attributes_ml(net), silent = TRUE)
    if (!inherits(any_table, "try-error") && is.data.frame(any_table) &&
        "target" %in% names(any_table)) {
      network_attrs <- subset(any_table, target %in% c("network","net","global"))
      network_attrs <- network_attrs[, setdiff(names(network_attrs), "target"), drop = FALSE]
    } else {
      network_attrs <- data.frame(name = character(), type = character(), stringsAsFactors = FALSE)
    }
  }

  # ---- Printing ---------------------------------------------------------------
  if (verbose) {
    cat("\n=== Attribute Report =====================================\n")
    # Added (from reports)
    cat("\n[Added (from stored reports)]\n")
    if (!is.null(added_actor_attrs)) .print_names(added_actor_attrs, "Actor attrs added")
    else cat("Actor attrs added: <not recorded; showing 'present' below>\n")
    if (!is.null(added_edge_attrs))  .print_names(added_edge_attrs,  "Edge attrs added")
    else cat("Edge attrs added:  <not recorded; showing 'present' below>\n")

    # Present now
    cat("\n[Present now]\n")
    .print_names(actor_attrs$name,  "Actor attributes")
    if (nrow(actor_attrs)) print(utils::head(actor_attrs, show_max), row.names = FALSE)

    .print_names(edge_attrs$name,   "Edge attributes")
    if (nrow(edge_attrs)) print(utils::head(edge_attrs, show_max), row.names = FALSE)

    .print_names(network_attrs$name,"Network-level attributes")
    if (nrow(network_attrs)) print(utils::head(network_attrs, show_max), row.names = FALSE)
  }

  # ---- Optional examples ------------------------------------------------------
  samples <- list()
  if (isTRUE(show_examples)) {
    # actor sample table
    atbl <- try(multinet::actors_ml(net, attributes = TRUE), silent = TRUE)
    if (!inherits(atbl, "try-error") && is.data.frame(atbl) && nrow(atbl)) {
      samples$actors <- utils::head(atbl, show_max)
      if (verbose) {
        cat("\n[Actor examples]\n")
        print(samples$actors, row.names = FALSE)
      }
    }

    # edge sample table
    etbl <- try(multinet::edges_ml(net, attributes = TRUE), silent = TRUE)
    if (!inherits(etbl, "try-error") && is.data.frame(etbl) && nrow(etbl)) {
      samples$edges <- utils::head(etbl, show_max)
      if (verbose) {
        cat("\n[Edge examples]\n")
        print(samples$edges, row.names = FALSE)
      }
    }
  }

  invisible(list(
    added = list(
      actor_attrs = added_actor_attrs,
      edge_attrs  = added_edge_attrs,
      reports     = list(actor = actor_report, edge = edge_report)
    ),
    present = list(
      actor   = actor_attrs,
      edge    = edge_attrs,
      network = network_attrs
    ),
    samples = samples
  ))
}



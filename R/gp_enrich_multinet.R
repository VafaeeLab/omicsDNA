

# -----------------------------------------------------------------------------
# 27 - gp_enrich_multinet(): g:Profiler per layer / per community
# -----------------------------------------------------------------------------

#' g:Profiler enrichment for multilayer networks (layer- or community-wise)
#'
#' @title Functional enrichment with g:Profiler for multilayer networks
#'
#' @description
#' For each selected layer (and optionally each community in each layer), this
#' function:
#' \enumerate{
#'   \item extracts gene symbols (from vertex names or community memberships);
#'   \item optionally restricts them to a user-supplied gene set
#'         (\code{restrict_genes});
#'   \item runs \strong{g:Profiler} (\code{\link[gprofiler2]{gost}}) using the
#'         requested sources;
#'   \item saves a \strong{CSV} with results, including a plain-text
#'         \code{intersection_genes} column and simple enrichment metrics;
#'   \item saves a \strong{gost} dot plot of significant terms;
#'   \item saves a \strong{term–gene concept network};
#'   \item saves a \strong{term-only overlap network}.
#' }
#'
#' In \strong{layer scope}, the output structure is conserved:
#' \preformatted{
#'   <results_dir>/<run_name>/
#'     <layer>/
#'       gprof_<layer>.csv
#'       gprof_dot_<layer>.<format>
#'       gprof_cnet_<layer>.<format>
#'       gprof_termnet_<layer>.<format>
#'     SUMMARY.csv
#' }
#'
#' In \strong{community scope}, results are nested inside each layer:
#' \preformatted{
#'   <results_dir>/<run_name>/
#'     <layer>/
#'       cid_<CID>/
#'         gprof_<layer>_cid<CID>.csv
#'         gprof_dot_<layer>_cid<CID>.<format>
#'         gprof_cnet_<layer>_cid<CID>.<format>
#'         gprof_termnet_<layer>_cid<CID>.<format>
#'     SUMMARY.csv
#' }
#'
#' @details
#' \strong{Inputs.}
#' Layers are read from \code{multinet::layers_ml(net)} and converted to
#' \pkg{igraph}s via \code{as.list(net)} (vertex names = gene symbols).
#'
#' \strong{Scopes.}
#' \itemize{
#'   \item \code{enrich_scope = "layer"} (default): one query per layer
#'         (vertex set of the layer, optionally restricted by \code{restrict_genes}).
#'   \item \code{enrich_scope = "community"}: one query per community per layer,
#'         as defined by \code{communities} (\code{actor}, \code{layer}, \code{cid}).
#' }
#'
#' \strong{Communities.} \code{communities} should be a data.frame with at least:
#' \itemize{
#'   \item an \emph{actor} column (gene symbol / actor ID);
#'   \item a \emph{layer} column (matching \code{multinet::layers_ml(net)});
#'   \item a community column. The function normalises this to \code{cid} using:
#'         \code{cid}, \code{com}, or \code{community} (in that order).
#' }
#'
#' \strong{Custom query restriction.}
#' If \code{restrict_genes} is non-\code{NULL}, the query set for each layer
#' (or community) is intersected with this vector before calling g:Profiler.
#' Only communities/layers with at least \code{min_query_size} genes after
#' restriction are tested.
#'
#' \strong{Statistics.}
#' G:Profiler returns adjusted p-values using \code{correction_method}:
#' \code{"g_SCS"}, \code{"fdr"} (Benjamini–Hochberg FDR), or \code{"bonferroni"}.
#' The threshold is given by \code{user_threshold}.
#'
#' Additionally, this function supports \code{correction_method = "none"}:
#' \itemize{
#'   \item g:Profiler is called without significance filtering;
#'   \item raw hypergeometric p-values (\code{p_raw}) are computed from
#'         \code{effective_domain_size}, \code{term_size},
#'         \code{query_size}, and \code{intersection_size};
#'   \item if \code{significant = TRUE}, terms are kept only when
#'         \code{p_raw <= user_threshold};
#'   \item CSVs include an extra \code{p_raw} column, and ranking/plots
#'         use \code{p_raw}.
#' }
#'
#' The result CSV also contains:
#' \itemize{
#'   \item \code{intersection_genes} — semicolon-separated gene symbols;
#'   \item \code{precision = intersection_size / query_size};
#'   \item \code{recall    = intersection_size / term_size}.
#' }
#'
#' @param net A \pkg{multinet} object.
#' @param organism Organism string accepted by g:Profiler (e.g. "hsapiens").
#' @param sources Character vector of g:Profiler sources (e.g. c("GO:BP","REAC")).
#' @param layer_order Optional character vector of layers to process (and order).
#' @param significant Logical; if TRUE, apply a p-value threshold:
#'   for standard corrections this is handled by g:Profiler;
#'   for \code{correction_method = "none"} it is applied to \code{p_raw}.
#' @param user_threshold Numeric significance threshold.
#' @param correction_method One of "g_SCS", "fdr", "bonferroni", or "none".
#'   See Details for behaviour of "none".
#' @param exclude_iea Logical; exclude electronic annotations?
#' @param domain_scope "annotated" or "custom" (with \code{custom_bg}).
#' @param custom_bg Optional background genes (symbols), used when
#'   \code{domain_scope = "custom"}.
#' @param show_terms Number of top terms (by p-value) for dot plot and concept net.
#' @param results_dir Base output directory.
#' @param run_name Optional subfolder name under \code{results_dir}; if NULL a
#'   timestamped name is generated.
#' @param per_layer_subdirs Logical; if TRUE (default) create a subdirectory per
#'   layer under the run directory.
#' @param format Plot format: "png", "pdf", or "jpg".
#' @param width,height,dpi Plot geometry (inches; DPI ignored for PDF).
#' @param show_in_rstudio Logical; print plots to current device.
#' @param verbose Logical; print progress messages.
#' @param ton_show Logical; draw term-only network?
#' @param ton_top_terms Number of top terms to include in term-only network.
#' @param ton_min_overlap Minimum shared genes between terms to draw an edge.
#' @param ton_layout Layout for term-only network ("fr","kk","lgl","mds").
#' @param ton_width,ton_height Size (inches) for term-only network.
#' @param ton_margin_cm Vector c(top,right,bottom,left) in cm for term-only plot.
#' @param communities Optional data.frame describing communities.
#' @param enrich_scope "layer" (default) or "community".
#' @param restrict_genes Optional character vector; if non-NULL, each query
#'   gene set is restricted to its intersection with this vector.
#' @param min_query_size Minimum number of genes in a query to run g:Profiler.
#'
#' @return Invisibly, a list with:
#' \describe{
#'   \item{\code{run_dir}}{Absolute path to the run folder.}
#'   \item{\code{by_layer}}{Layer-keyed list. For layer-scope, each element
#'       contains file paths and \code{n_genes}/\code{n_terms}. For community-
#'       scope, each layer contains a list keyed by community ID (e.g. "cid_C1").}
#'   \item{\code{summary_file}}{Path to run-level \code{SUMMARY.csv}.}
#' }
#'
#' @importFrom gprofiler2 gost gostplot
#' @importFrom ggplot2 ggsave ggplot aes geom_point scale_size_continuous
#'   scale_color_viridis_c labs theme_minimal theme element_text element_blank
#' @importFrom igraph graph_from_data_frame layout_with_fr layout_with_kk
#'   layout_with_lgl layout_with_mds V
#' @importFrom utils write.csv combn
#' @importFrom multinet layers_ml
#' @export
gp_enrich_multinet <- function(
    net,
    organism          = "hsapiens",
    sources           = c("GO:BP","GO:MF","GO:CC","REAC"),
    layer_order       = NULL,
    significant       = TRUE,
    user_threshold    = 0.05,
    correction_method = c("g_SCS","fdr","bonferroni","none"),
    exclude_iea       = FALSE,
    domain_scope      = c("annotated","custom"),
    custom_bg         = NULL,
    show_terms        = 10,
    results_dir       = getOption("mlnet.results_dir","omicsDNA_results"),
    run_name          = NULL,
    per_layer_subdirs = TRUE,
    format            = c("png","pdf","jpg"),
    width             = 9, height = 7, dpi = 300,
    show_in_rstudio   = TRUE,
    verbose           = TRUE,
    # term-only overlap network
    ton_show          = TRUE,
    ton_top_terms     = 10,
    ton_min_overlap   = 1,
    ton_layout        = "fr",
    ton_width         = width,
    ton_height        = height,
    ton_margin_cm     = c(0.5,0.5,0.5,0.5),
    # communities, scope & restriction
    communities       = NULL,
    enrich_scope      = c("layer","community"),
    restrict_genes    = NULL,
    min_query_size    = 3
) {
  # ---- dependencies ----
  if (!requireNamespace("multinet",  quietly = TRUE))
    stop("Please install 'multinet'.")
  if (!requireNamespace("gprofiler2", quietly = TRUE))
    stop("Please install 'gprofiler2': install.packages('gprofiler2')")
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Please install 'ggplot2': install.packages('ggplot2')")
  has_net <- requireNamespace("ggraph", quietly = TRUE) &&
    requireNamespace("igraph", quietly = TRUE)
  if (!has_net && isTRUE(verbose))
    message("Note: 'ggraph' and/or 'igraph' missing — network plots will be skipped.")

  correction_method <- match.arg(correction_method)
  domain_scope      <- match.arg(domain_scope)
  format            <- match.arg(format)
  enrich_scope      <- match.arg(enrich_scope)

  # ---- helpers ---------------------------------------------------------------
  .ensure_dir <- function(d) {
    if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
    invisible(d)
  }
  .stamp   <- function() format(Sys.time(), "%Y-%m-%d_%H%M%S")
  .unit_cm <- function(v) grid::unit(v, "cm")
  .rescale <- function(x, to = c(4,10)) {
    x <- suppressWarnings(as.numeric(x))
    r <- range(x, finite = TRUE)
    if (!is.finite(r[1]) || !is.finite(r[2]) || r[1] == r[2])
      return(rep(mean(to), length(x)))
    (x - r[1])/(r[2] - r[1])*(to[2] - to[1]) + to[1]
  }
  .flatten_df <- function(df) {
    df <- as.data.frame(df, stringsAsFactors = FALSE)
    for (nm in names(df)) if (is.list(df[[nm]]))
      df[[nm]] <- vapply(df[[nm]], function(x) {
        if (is.null(x) || (is.atomic(x) && length(x) == 0L)) return(NA_character_)
        paste(as.character(unlist(x, use.names = FALSE)), collapse = ";")
      }, character(1))
    df
  }
  .extract_symbols <- function(x) {
    if (is.null(x)) return(character(0))
    if (is.list(x)) x <- unlist(x, use.names = FALSE)
    x <- as.character(x)
    if (length(x) == 1L)
      x <- unlist(strsplit(x, "[,;]", perl = TRUE), use.names = FALSE)
    x <- trimws(x)
    x <- sub("\\s*\\(.*\\)$", "", x)  # drop "(EVIDENCE)" etc
    unique(x[nzchar(x)])
  }
  .layout_terms <- function(g, method = "fr") {
    switch(tolower(method),
           fr  = igraph::layout_with_fr(g),
           kk  = igraph::layout_with_kk(g),
           lgl = igraph::layout_with_lgl(g),
           mds = igraph::layout_with_mds(g),
           igraph::layout_with_fr(g))
  }

  # ---- normalise communities (if given) -------------------------------------
  comm <- NULL
  if (!is.null(communities)) {
    comm <- as.data.frame(communities, stringsAsFactors = FALSE)

    # Normalise `cid` column if needed
    if (!"cid" %in% names(comm)) {
      alt <- intersect(c("cid","com","community"), names(comm))[1L]
      if (is.na(alt))
        stop("Could not find cid/com/community column in `communities`.")
      comm$cid <- comm[[alt]]
    }

    # Normalise `actor`
    if (!"actor" %in% names(comm)) {
      alt <- intersect(c("actor","gene","gene_id","symbol"), names(comm))[1L]
      if (is.na(alt))
        stop("Could not find actor/gene/gene_id/symbol column in `communities`.")
      comm$actor <- comm[[alt]]
    }

    # Normalise `layer`
    if (!"layer" %in% names(comm)) {
      alt <- intersect(c("layer","Layer"), names(comm))[1L]
      if (is.na(alt))
        stop("Could not find layer/Layer column in `communities`.")
      comm$layer <- comm[[alt]]
    }

    comm$actor <- as.character(comm$actor)
    comm$layer <- as.character(comm$layer)
    comm$cid   <- as.character(comm$cid)
  }

  if (enrich_scope == "community" && is.null(comm))
    stop("For enrich_scope = 'community', please supply `communities`.")

  # ---- layer selection -------------------------------------------------------
  Ls <- try(multinet::layers_ml(net), silent = TRUE)
  if (inherits(Ls, "try-error") || is.null(Ls) || !length(Ls))
    stop("Could not get layers via multinet::layers_ml(net).")
  Ls <- as.character(Ls)
  layers <- if (is.null(layer_order)) {
    Ls
  } else {
    keep <- intersect(as.character(layer_order), Ls)
    if (!length(keep))
      stop("None of the requested layers are present in 'net'.")
    keep
  }

  # ---- igraph list -----------------------------------------------------------
  glist <- try(as.list(net), silent = TRUE)
  if (inherits(glist, "try-error") || !is.list(glist))
    stop("Could not coerce 'net' to list of igraphs via as.list(net).")
  if (!is.null(names(glist)) && any(names(glist) == "_flat_"))
    glist <- glist[names(glist) != "_flat_"]

  # ---- run folder ------------------------------------------------------------
  .ensure_dir(results_dir)
  if (is.null(run_name) || !nzchar(run_name))
    run_name <- paste0("gprofiler_multinet_", .stamp())
  run_dir <- file.path(results_dir, run_name)
  .ensure_dir(run_dir)
  if (isTRUE(verbose))
    message("Run folder: ", normalizePath(run_dir, winslash = "/", mustWork = FALSE))

  # ---- background logic ------------------------------------------------------
  if (!is.null(custom_bg) && domain_scope != "custom") {
    message("custom_bg provided: switching domain_scope to 'custom'.")
    domain_scope <- "custom"
  }
  if (identical(domain_scope, "custom") && is.null(custom_bg))
    stop("domain_scope='custom' requires a non-NULL 'custom_bg' gene vector.")

  # ---- core worker: run g:Profiler + plots for one query --------------------
  .run_query <- function(query_genes, stub, title_prefix, layer_dir) {
    # basic cleaning + restriction
    query_genes <- unique(as.character(query_genes))
    query_genes <- query_genes[nzchar(query_genes)]

    if (!is.null(restrict_genes)) {
      query_genes <- intersect(query_genes, as.character(restrict_genes))
    }

    n_query <- length(query_genes)
    if (n_query < min_query_size) {
      if (isTRUE(verbose))
        message("  [", title_prefix, "] only ", n_query,
                " genes after restriction (< ", min_query_size, "); skipped.")
      return(list(
        n_genes      = n_query,
        n_terms      = 0L,
        csv_file     = NA_character_,
        dotplot_file = NA_character_,
        cnet_file    = NA_character_,
        termnet_file = NA_character_
      ))
    }

    # choose how to call g:Profiler
    use_raw_p       <- identical(correction_method, "none")
    cm_for_gost     <- if (use_raw_p) "g_SCS" else correction_method
    signif_for_gost <- if (use_raw_p) FALSE   else significant
    thr_for_gost    <- if (use_raw_p) 1       else user_threshold

    # g:Profiler call
    gres <- try(
      gprofiler2::gost(
        query                       = query_genes,
        organism                    = organism,
        ordered_query               = FALSE,
        multi_query                 = FALSE,
        significant                 = signif_for_gost,
        exclude_iea                 = exclude_iea,
        measure_underrepresentation = FALSE,
        evcodes                     = TRUE,
        user_threshold              = thr_for_gost,
        correction_method           = cm_for_gost,
        domain_scope                = domain_scope,
        custom_bg                   = if (!is.null(custom_bg)) custom_bg else NULL,
        sources                     = sources
      ),
      silent = TRUE
    )

    if (inherits(gres, "try-error") ||
        is.null(gres) ||
        is.null(gres$result) ||
        !nrow(gres$result)) {
      if (isTRUE(verbose))
        message("  [", title_prefix, "] no g:Profiler results.")
      return(list(
        n_genes      = n_query,
        n_terms      = 0L,
        csv_file     = NA_character_,
        dotplot_file = NA_character_,
        cnet_file    = NA_character_,
        termnet_file = NA_character_
      ))
    }

    # Flatten + intersection_genes + simple metrics
    inter_chr <- vapply(seq_len(nrow(gres$result)), function(i) {
      s <- .extract_symbols(gres$result$intersection[[i]])
      if (length(s)) paste(s, collapse = ";") else NA_character_
    }, character(1))

    df_flat <- .flatten_df(gres$result)
    df_flat$intersection_genes <- inter_chr

    if (all(c("intersection_size","term_size","query_size") %in% names(df_flat))) {
      qs <- suppressWarnings(as.numeric(df_flat$query_size))
      ts <- suppressWarnings(as.numeric(df_flat$term_size))
      is <- suppressWarnings(as.numeric(df_flat$intersection_size))
      df_flat$precision <- ifelse(qs > 0, is / qs, NA_real_)
      df_flat$recall    <- ifelse(ts > 0, is / ts, NA_real_)
    }

    # ---- p-values: standard vs "none" ---------------------------------------
    if (use_raw_p) {
      req <- c("intersection_size","term_size","query_size","effective_domain_size")
      if (all(req %in% names(df_flat))) {
        k  <- suppressWarnings(as.numeric(df_flat$intersection_size))
        m  <- suppressWarnings(as.numeric(df_flat$term_size))
        q  <- suppressWarnings(as.numeric(df_flat$query_size))
        N  <- suppressWarnings(as.numeric(df_flat$effective_domain_size))
        df_flat$p_raw <- stats::phyper(k - 1, m, N - m, q, lower.tail = FALSE)
      } else {
        df_flat$p_raw <- NA_real_
      }
      p_used <- df_flat$p_raw

      # apply threshold manually if significant = TRUE
      if (isTRUE(significant)) {
        keep <- which(!is.na(p_used) & p_used <= user_threshold)
        if (!length(keep)) {
          if (isTRUE(verbose))
            message("  [", title_prefix,
                    "] no terms pass raw p-value threshold (", user_threshold, ").")
          return(list(
            n_genes      = n_query,
            n_terms      = 0L,
            csv_file     = NA_character_,
            dotplot_file = NA_character_,
            cnet_file    = NA_character_,
            termnet_file = NA_character_
          ))
        }
        df_flat     <- df_flat[keep, , drop = FALSE]
        gres$result <- gres$result[keep, , drop = FALSE]
        p_used      <- p_used[keep]
      }
    } else {
      df_flat$p_raw <- NA_real_
      p_used <- suppressWarnings(as.numeric(df_flat$p_value))
    }

    if (!nrow(gres$result)) {
      if (isTRUE(verbose))
        message("  [", title_prefix, "] (after filtering) no results remain.")
      return(list(
        n_genes      = n_query,
        n_terms      = 0L,
        csv_file     = NA_character_,
        dotplot_file = NA_character_,
        cnet_file    = NA_character_,
        termnet_file = NA_character_
      ))
    }

    # ---- CSV -----------------------------------------------------------------
    csv_file <- file.path(layer_dir, paste0("gprof_", stub, ".csv"))
    utils::write.csv(df_flat, csv_file, row.names = FALSE)
    if (isTRUE(verbose))
      message("  [", title_prefix, "] CSV: ",
              normalizePath(csv_file, winslash = "/", mustWork = FALSE))

    # ---- ordering for plots (always by p_used) ------------------------------
    ord <- order(p_used, na.last = TRUE)
    gresTop <- gres
    gresTop$result <- gres$result[ord, , drop = FALSE]
    if (nrow(gresTop$result) > show_terms)
      gresTop$result <- head(gresTop$result, show_terms)

    # ---- gost dot plot ------------------------------------------------------
    dotplot_file <- NA_character_
    gp <- try(gprofiler2::gostplot(gresTop, capped = TRUE, interactive = FALSE),
              silent = TRUE)
    if (!inherits(gp, "try-error")) {
      dotplot_file <- file.path(layer_dir, paste0("gprof_dot_", stub, ".", format))
      if (isTRUE(show_in_rstudio)) print(gp)
      ggplot2::ggsave(dotplot_file, plot = gp, width = width, height = height,
                      dpi = if (identical(format, "pdf")) NA_integer_ else dpi,
                      bg = "white")
      if (isTRUE(verbose))
        message("  [", title_prefix, "] dot plot: ",
                normalizePath(dotplot_file, winslash = "/", mustWork = FALSE))
    }

    # ---- concept (term–gene) network ---------------------------------------
    cnet_file <- NA_character_
    if (show_terms > 0 && has_net) {
      dft <- gres$result[ord, , drop = FALSE]
      dft <- head(dft, n = min(show_terms, nrow(dft)))
      if (nrow(dft) > 0) {
        ed_rows <- lapply(seq_len(nrow(dft)), function(i) {
          gi <- .extract_symbols(dft$intersection[[i]])
          if (!length(gi)) return(NULL)
          data.frame(term_id = dft$term_id[i], gene = gi, stringsAsFactors = FALSE)
        })
        ed_rows <- Filter(Negate(is.null), ed_rows)
        edges <- if (length(ed_rows)) do.call(rbind, ed_rows) else NULL

        if (!is.null(edges) && nrow(edges)) {
          tn <- unique(dft[, c("term_id","term_name","p_value","intersection_size")])
          tn$lp <- -log10(pmax(as.numeric(tn$p_value), .Machine$double.xmin))
          nodes_term <- data.frame(
            id        = tn$term_id,
            label     = tn$term_name,
            node_type = "term",
            lp        = tn$lp,
            size      = .rescale(tn$intersection_size, c(4,10)),
            stringsAsFactors = FALSE
          )
          genes_unique <- unique(edges$gene)
          nodes_gene <- data.frame(
            id        = genes_unique,
            label     = genes_unique,
            node_type = "gene",
            lp        = NA_real_,
            size      = 3,
            stringsAsFactors = FALSE
          )
          nodes <- rbind(nodes_term, nodes_gene)
          g_bi  <- igraph::graph_from_data_frame(
            d = data.frame(from = edges$term_id, to = edges$gene,
                           stringsAsFactors = FALSE),
            directed = FALSE, vertices = nodes
          )

          cnet_file <- file.path(layer_dir, paste0("gprof_cnet_", stub, ".", format))
          pnet <- ggraph::ggraph(g_bi, layout = "fr") +
            ggraph::geom_edge_link(alpha = 0.25) +
            ggraph::geom_node_point(
              ggplot2::aes(size = size,
                           color = ifelse(nodes$node_type == "term", nodes$lp, NA_real_))
            ) +
            ggplot2::scale_size_continuous(range = c(2, 9), guide = "none") +
            ggplot2::scale_color_viridis_c(na.value = "grey60",
                                           name = expression(-log[10](p))) +
            ggraph::geom_node_text(
              ggplot2::aes(label = label),
              size = 3, color = "black", check_overlap = TRUE
            ) +
            ggplot2::labs(title = paste0("Term–gene concept network: ", title_prefix)) +
            ggplot2::theme_minimal(base_size = 12) +
            ggplot2::theme(legend.position = "right")

          if (isTRUE(show_in_rstudio)) print(pnet)
          ggplot2::ggsave(cnet_file, plot = pnet, width = width, height = height,
                          dpi = if (identical(format, "pdf")) NA_integer_ else dpi,
                          bg = "white")
          if (isTRUE(verbose))
            message("  [", title_prefix, "] concept net: ",
                    normalizePath(cnet_file, winslash = "/", mustWork = FALSE))
        } else if (isTRUE(verbose)) {
          message("  [", title_prefix, "] no term–gene intersections to plot.")
        }
      }
    } else if (show_terms > 0 && !has_net) {
      warning("Skipping concept network for '", title_prefix,
              "' (need ggraph + igraph).")
    }

    # ---- term-only overlap network -----------------------------------------
    termnet_file <- NA_character_
    if (isTRUE(ton_show) && has_net) {
      dfo <- head(gres$result[ord, , drop = FALSE],
                  n = min(ton_top_terms, nrow(gres$result)))
      if (nrow(dfo) > 0) {
        term_list <- lapply(seq_len(nrow(dfo)),
                            function(i) .extract_symbols(dfo$intersection[[i]]))
        names(term_list) <- dfo$term_id
        ids <- dfo$term_id

        e_from <- character(0); e_to <- character(0); e_w <- integer(0)
        if (length(ids) >= 2L) {
          for (pair in utils::combn(ids, 2, simplify = FALSE)) {
            gi <- term_list[[pair[1]]]; gj <- term_list[[pair[2]]]
            if (length(gi) && length(gj)) {
              ov <- length(intersect(gi, gj))
              if (ov >= ton_min_overlap) {
                e_from <- c(e_from, pair[1])
                e_to   <- c(e_to,   pair[2])
                e_w    <- c(e_w,    ov)
              }
            }
          }
        }

        term_nodes <- data.frame(
          id    = dfo$term_id,
          label = dfo$term_name,
          lp    = -log10(pmax(as.numeric(dfo$p_value), .Machine$double.xmin)),
          k     = as.integer(dfo$intersection_size),
          size  = .rescale(dfo$intersection_size, c(4,10)),
          stringsAsFactors = FALSE
        )

        g_terms <- if (length(e_from)) {
          igraph::graph_from_data_frame(
            d = data.frame(from = e_from, to = e_to, weight = e_w,
                           stringsAsFactors = FALSE),
            directed = FALSE, vertices = term_nodes
          )
        } else {
          igraph::graph_from_data_frame(d = NULL, directed = FALSE,
                                        vertices = term_nodes)
        }

        lay <- .layout_terms(g_terms, ton_layout)
        nodes_df <- cbind(term_nodes, data.frame(x = lay[,1], y = lay[,2]))
        edf <- if (length(e_from)) data.frame(
          from = e_from, to = e_to, weight = e_w,
          x1 = nodes_df$x[match(e_from, nodes_df$id)],
          y1 = nodes_df$y[match(e_from, nodes_df$id)],
          x2 = nodes_df$x[match(e_to,   nodes_df$id)],
          y2 = nodes_df$y[match(e_to,   nodes_df$id)]
        ) else data.frame(from=character(0), to=character(0), weight=integer(0),
                          x1=numeric(0), y1=numeric(0),
                          x2=numeric(0), y2=numeric(0))

        p_ton <- ggplot2::ggplot() +
          ggplot2::geom_segment(
            data = edf,
            ggplot2::aes(x = x1, y = y1, xend = x2, yend = y2),
            color = "grey70", alpha = 0.5
          ) +
          ggplot2::geom_point(
            data = nodes_df,
            ggplot2::aes(x = x, y = y, size = size, color = lp)
          ) +
          ggplot2::scale_size_continuous(range = c(3, 10), guide = "none") +
          ggplot2::scale_color_viridis_c(name = expression(-log[10](p))) +
          ggplot2::geom_text(
            data = nodes_df,
            ggplot2::aes(x = x, y = y, label = label),
            size = 4, color = "black", vjust = -0.8
          ) +
          ggplot2::geom_text(
            data = nodes_df,
            ggplot2::aes(x = x, y = y, label = k),
            size = 3.8, color = "grey20", vjust = 1.8
          ) +
          ggplot2::labs(
            title = paste0("Term-only overlap network: ", title_prefix),
            x = "x", y = "y"
          ) +
          ggplot2::theme_minimal(base_size = 12) +
          ggplot2::theme(plot.margin = .unit_cm(ton_margin_cm))

        termnet_file <- file.path(layer_dir, paste0("gprof_termnet_", stub, ".", format))
        if (isTRUE(show_in_rstudio)) print(p_ton)
        ggplot2::ggsave(termnet_file, plot = p_ton,
                        width = ton_width, height = ton_height,
                        dpi = if (identical(format, "pdf")) NA_integer_ else dpi,
                        bg = "white")
        if (isTRUE(verbose))
          message("  [", title_prefix, "] term-only net: ",
                  normalizePath(termnet_file, winslash = "/", mustWork = FALSE))
      }
    } else if (isTRUE(ton_show) && !has_net) {
      warning("Skipping term-only network for '", title_prefix,
              "' (need ggraph + igraph).")
    }

    list(
      n_genes      = n_query,
      n_terms      = nrow(gres$result),
      csv_file     = csv_file,
      dotplot_file = if (!is.na(dotplot_file)) dotplot_file else NA_character_,
      cnet_file    = if (!is.na(cnet_file))    cnet_file    else NA_character_,
      termnet_file = if (!is.na(termnet_file)) termnet_file else NA_character_
    )
  }

  # ---- main loops ------------------------------------------------------------
  by_layer     <- stats::setNames(vector("list", length(layers)), layers)
  summary_rows <- list()

  for (ln in layers) {
    if (isTRUE(verbose)) message("Layer: ", ln)

    layer_dir <- if (isTRUE(per_layer_subdirs)) {
      .ensure_dir(file.path(run_dir, ln))
    } else {
      run_dir
    }

    if (enrich_scope == "layer") {
      g <- glist[[ln]]
      if (!inherits(g, "igraph")) {
        warning("  Skipping '", ln, "' (not an igraph).")
        next
      }
      genes <- unique(as.character(igraph::V(g)$name))
      genes <- genes[nzchar(genes)]
      if (!length(genes)) {
        warning("  Skipping '", ln, "' (no vertex names).")
        next
      }

      res <- .run_query(genes, stub = ln, title_prefix = ln, layer_dir = layer_dir)

      by_layer[[ln]] <- res
      summary_rows[[length(summary_rows) + 1L]] <- data.frame(
        scope        = "layer",
        layer        = ln,
        cid          = NA_character_,
        n_genes      = res$n_genes,
        n_terms      = res$n_terms,
        csv_file     = res$csv_file,
        dotplot_file = res$dotplot_file,
        cnet_file    = res$cnet_file,
        termnet_file = res$termnet_file,
        stringsAsFactors = FALSE
      )

    } else if (enrich_scope == "community") {

      comm_l <- comm[comm$layer == ln, , drop = FALSE]
      if (!nrow(comm_l)) {
        if (isTRUE(verbose))
          message("  No communities for layer ", ln, "; skipped.")
        by_layer[[ln]] <- list()
        next
      }

      cids <- sort(unique(comm_l$cid))
      layer_list <- list()

      for (cid in cids) {
        mod_label <- paste0("cid_", cid)
        title     <- paste0(ln, " / ", cid)

        module_dir <- .ensure_dir(file.path(layer_dir, mod_label))
        gset <- unique(as.character(comm_l$actor[comm_l$cid == cid]))
        gset <- gset[nzchar(gset)]

        res <- .run_query(gset,
                          stub         = paste0(ln, "_", mod_label),
                          title_prefix = title,
                          layer_dir    = module_dir)

        layer_list[[mod_label]] <- res
        summary_rows[[length(summary_rows) + 1L]] <- data.frame(
          scope        = "community",
          layer        = ln,
          cid          = cid,
          n_genes      = res$n_genes,
          n_terms      = res$n_terms,
          csv_file     = res$csv_file,
          dotplot_file = res$dotplot_file,
          cnet_file    = res$cnet_file,
          termnet_file = res$termnet_file,
          stringsAsFactors = FALSE
        )
      }

      by_layer[[ln]] <- layer_list
    }
  }

  # ---- manifest --------------------------------------------------------------
  summary_df <- if (length(summary_rows)) {
    do.call(rbind, summary_rows)
  } else {
    data.frame(
      scope        = character(),
      layer        = character(),
      cid          = character(),
      n_genes      = integer(),
      n_terms      = integer(),
      csv_file     = character(),
      dotplot_file = character(),
      cnet_file    = character(),
      termnet_file = character(),
      stringsAsFactors = FALSE
    )
  }

  summary_file <- file.path(run_dir, "SUMMARY.csv")
  utils::write.csv(summary_df, summary_file, row.names = FALSE)
  if (isTRUE(verbose))
    message("Saved manifest: ",
            normalizePath(summary_file, winslash = "/", mustWork = FALSE))

  invisible(list(
    run_dir      = run_dir,
    by_layer     = by_layer,
    summary_file = summary_file
  ))
}



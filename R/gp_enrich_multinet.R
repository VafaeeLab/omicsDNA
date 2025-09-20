

# -----------------------------------------------------------------------------
# 27 - gp_enrich_multinet(): g:Profiler per layer + CSV + dot plot + 2 networks
# -----------------------------------------------------------------------------

#' g:Profiler enrichment per layer (CSV + gost plot + concept network + term-only network)
#'
#' @title Layer-wise functional enrichment with g:Profiler for multilayer networks,
#' and automatic export of CSV, dot plot, concept (term–gene) network, and term-only overlap network
#'
#' @description
#' For each selected layer of a multilayer network, this function:
#' \enumerate{
#'   \item extracts the layer's gene symbols (from vertex names);
#'   \item runs \strong{g:Profiler} enrichment (\code{\link[gprofiler2]{gost}}) against the requested sources;
#'   \item saves a \strong{gost} dot plot of significant terms;
#'   \item saves a \strong{concept (term–gene) network} showing terms connected to their hit genes
#'         (both term names and gene symbols are labelled; terms are colored by \eqn{-\log_{10}(p)} and
#'         sized by intersection size; genes are grey);
#'   \item saves a \strong{term-only overlap network} (nodes = terms; edges drawn when two terms share
#'         at least \code{ton_min_overlap} genes; \emph{numbers near nodes show each term’s
#'         \code{intersection_size}});
#'   \item writes a CSV of g:Profiler results, augmented with a plain-text
#'         \code{intersection_genes} (semicolon-separated gene symbols).
#' }
#'
#' A per-run manifest (\code{SUMMARY.csv}) lists what was saved for each layer.
#'
#' @details
#' \strong{Inputs.} Layers are read from \code{multinet::layers_ml(net)}. Each layer is converted to an
#' \pkg{igraph} using \code{as.list(net)}; the node names of that graph are interpreted as gene symbols.
#'
#' \strong{Enrichment.} The g:Profiler call uses \code{evcodes = TRUE} so that \code{$intersection} carries
#' actual gene members per term. You can limit to significant results and set the multiple-testing method
#' and threshold via \code{significant}, \code{correction_method}, and \code{user_threshold}.
#' Background choice is controlled with \code{domain_scope} and (optionally) \code{custom_bg}:
#' \itemize{
#'   \item \code{domain_scope = "annotated"} (default): test against the annotation domain.
#'   \item \code{domain_scope = "custom"}: test against a user-supplied background gene set
#'         passed in \code{custom_bg} (required in this mode).
#'   \item If you pass a non-\code{NULL} \code{custom_bg} while \code{domain_scope != "custom"}, the function
#'         automatically switches to \code{domain_scope = "custom"} (with a message).
#' }
#'
#' \strong{Plots.}
#' \itemize{
#'   \item The gost dot plot is produced by \code{\link[gprofiler2]{gostplot}} and is downsampled to the
#'         \code{show_terms} most significant results for legibility.
#'   \item The concept (term–gene) network connects each selected term (top \code{show_terms} by p-value)
#'         to its \code{intersection} genes; terms are colored by \eqn{-\log_{10}(p)} and sized by
#'         \code{intersection_size}. \emph{Both term names and gene symbols are printed.}
#'   \item The term-only overlap network places an edge between two terms when they share at least
#'         \code{ton_min_overlap} intersection genes. Node labels show term names; a second number near
#'         each node is the term’s \code{intersection_size}. Node color encodes \eqn{-\log_{10}(p)}.
#' }
#'
#' \strong{Output structure.} Results are written to:
#' \preformatted{
#'   <results_dir>/<run_name>/
#'     <layer>/
#'       gprof_<layer>.csv
#'       gprof_dot_<layer>.<format>
#'       gprof_cnet_<layer>.<format>
#'       gprof_termnet_<layer>.<format>
#'     SUMMARY.csv          (run-level manifest at the root of <run_name>)
#' }
#' When \code{per_layer_subdirs = FALSE}, the files are written directly under
#' \code{<results_dir>/<run_name>/} with layer-specific filenames.
#'
#' @section Scaling, units and aesthetics:
#' \itemize{
#'   \item Plot \code{width}/\code{height} are in \strong{inches}; \code{dpi} is ignored for PDF.
#'   \item Network node sizes are rescaled to a sensible visual range (terms larger than genes).
#'   \item \code{ton_margin_cm} is in \strong{centimetres} and set as a ggplot plot margin.
#' }
#'
#' @param net A \pkg{multinet} object. Layers are discovered via \code{multinet::layers_ml(net)} and
#'   coerced to \pkg{igraph} objects through \code{as.list(net)}; vertex names are assumed to be gene symbols.
#' @param organism Organism string accepted by g:Profiler (e.g., \code{"hsapiens"}, \code{"mmusculus"}).
#' @param sources Character vector of g:Profiler sources to include (e.g., \code{c("GO:BP","GO:MF","GO:CC","REAC")}).
#' @param layer_order Optional character vector specifying which layers to process and in what order. By default,
#'   all layers from \code{net} are used in their native order.
#' @param significant Logical. If \code{TRUE}, g:Profiler filters to significant terms according to
#'   \code{user_threshold} and \code{correction_method}.
#' @param user_threshold Numeric. Significance threshold (e.g., \code{0.05}).
#' @param correction_method Multiple testing method for g:Profiler. One of \code{"g_SCS"}, \code{"fdr"},
#'   or \code{"bonferroni"}.
#' @param exclude_iea Logical. If \code{TRUE}, exclude electronic annotations (IEA) where applicable.
#' @param domain_scope One of \code{"annotated"} (default) or \code{"custom"}.
#'   Use \code{"custom"} when supplying \code{custom_bg}.
#' @param custom_bg Optional character vector of background genes. Required when \code{domain_scope="custom"}.
#' @param show_terms Integer. Number of top terms (by p-value) to display in the dot plot and the concept network.
#' @param results_dir Base output directory. Created if it does not exist.
#' @param run_name Name of the run subfolder under \code{results_dir}. If \code{NULL}, a name of the form
#'   \code{"gprofiler_multinet_<YYYY-mm-dd_HHMMSS>"} is generated.
#' @param per_layer_subdirs Logical; if \code{TRUE} (default) create a subdirectory per layer.
#' @param format Image format for saved figures: \code{"png"}, \code{"pdf"}, or \code{"jpg"}.
#' @param width,height,dpi Plot geometry for saved figures (inches; DPI ignored for PDF).
#' @param show_in_rstudio Logical; if \code{TRUE}, also \code{print()} plots to the current graphics device.
#' @param verbose Logical; if \code{TRUE}, print progress messages.
#' @param ton_show Logical; if \code{TRUE} (default), draw the term-only overlap network.
#' @param ton_top_terms Integer; number of top terms (by p-value) to include in the term-only network.
#' @param ton_min_overlap Integer; minimum \emph{shared intersection genes} between two terms to draw an edge.
#' @param ton_layout Layout for the term-only network. One of \code{"fr"}, \code{"kk"}, \code{"lgl"}, \code{"mds"}.
#' @param ton_width,ton_height Size (inches) for the term-only network figure.
#' @param ton_margin_cm Numeric vector \code{c(top, right, bottom, left)} in centimetres for term-only plot margins.
#'
#' @return
#' Invisibly returns a list with:
#' \describe{
#'   \item{\code{run_dir}}{Absolute path to the run directory.}
#'   \item{\code{by_layer}}{Named list keyed by layer containing file paths for
#'         \code{csv_file}, \code{dotplot_file}, \code{cnet_file}, and \code{termnet_file}.}
#'   \item{\code{summary_file}}{Absolute path to the run-level \code{SUMMARY.csv} manifest.}
#' }
#'
#' @section Notes and troubleshooting:
#' \itemize{
#'   \item \strong{Custom background:} When \code{domain_scope="custom"}, \code{custom_bg} must be provided;
#'         otherwise the function stops with a clear error.
#'   \item \strong{Very many terms:} Increase \code{show_terms}, figure sizes, or use PDF for higher fidelity.
#'   \item \strong{Networks require} \pkg{ggraph} and \pkg{igraph}. If either is missing, network plots are skipped.
#' }
#'
#' @seealso
#' \code{\link[gprofiler2]{gost}}, \code{\link[gprofiler2]{gostplot}},
#' \code{\link[ggraph]{ggraph}}, \code{\link[igraph]{graph_from_data_frame}}
#'
#' @importFrom gprofiler2 gost gostplot
#' @importFrom ggplot2 ggsave ggplot aes geom_point scale_size_continuous scale_color_viridis_c labs theme_minimal theme element_text element_blank
#' @importFrom igraph graph_from_data_frame layout_with_fr layout_with_kk layout_with_lgl layout_with_mds V
#' @importFrom utils write.csv combn
#' @importFrom multinet layers_ml
#'
#' @examples
#' \dontrun{
#' ## Single comprehensive example (demonstrates all arguments)
#' ## Assume you have a 'multinet' object 'net' with layers E1, E2, M1 whose vertex names are gene symbols.
#'
#' # Prepare a custom background as the union of all vertex names across layers:
#' all_bg <- unique(unlist(lapply(as.list(net), function(g) igraph::V(g)$name)))
#'
#' res <- gp_enrich_multinet(
#'   net                = net,
#'   organism           = "hsapiens",
#'   sources            = c("GO:BP","GO:MF","GO:CC","REAC"),
#'   layer_order        = c("E1","E2","M1"),
#'   significant        = TRUE,
#'   user_threshold     = 0.05,
#'   correction_method  = "g_SCS",                  # or "fdr", "bonferroni"
#'   exclude_iea        = FALSE,
#'   domain_scope       = "custom",                 # <-- use the custom background below
#'   custom_bg          = all_bg,
#'   show_terms         = 12,                       # limit for dot plot + concept network
#'   results_dir        = getOption("mlnet.results_dir","omicsDNA_results"),
#'   run_name           = paste0("gprofiler_run_", format(Sys.time(), "%Y%m%d_%H%M%S")),
#'   per_layer_subdirs  = TRUE,
#'   format             = "png",
#'   width              = 9, height = 7, dpi = 300,
#'   show_in_rstudio    = TRUE,
#'   verbose            = TRUE,
#'   # term-only overlap network:
#'   ton_show           = TRUE,
#'   ton_top_terms      = 10,
#'   ton_min_overlap    = 1,
#'   ton_layout         = "fr",
#'   ton_width          = 9,
#'   ton_height         = 7,
#'   ton_margin_cm      = c(0.8,0.8,0.8,0.8)
#' )
#'
#' # Inspect the run manifest:
#' read.csv(res$summary_file, stringsAsFactors = FALSE)
#' }
#' @export
gp_enrich_multinet <- function(
    net,
    organism          = "hsapiens",
    sources           = c("GO:BP","GO:MF","GO:CC","REAC"),
    layer_order       = NULL,
    significant       = TRUE,
    user_threshold    = 0.05,
    correction_method = c("g_SCS","fdr","bonferroni"),
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
    ton_margin_cm     = c(0.5,0.5,0.5,0.5)
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

  # ---- helpers ----
  .ensure_dir <- function(d) if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  .stamp      <- function() format(Sys.time(), "%Y-%m-%d_%H%M%S")
  .unit_cm    <- function(v) grid::unit(v, "cm")
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
        if (is.null(x) || (is.atomic(x) && length(x)==0L)) return(NA_character_)
        paste(as.character(unlist(x, use.names = FALSE)), collapse = ";")
      }, character(1))
    df
  }
  .extract_symbols <- function(x) {
    if (is.null(x)) return(character(0))
    if (is.list(x)) x <- unlist(x, use.names = FALSE)
    x <- as.character(x)
    if (length(x) == 1L) x <- unlist(strsplit(x, "[,;]", perl = TRUE), use.names = FALSE)
    x <- trimws(x)
    x <- sub("\\s*\\(.*\\)$", "", x)  # strip "(EVIDENCE)" if present
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

  # ---- layer selection ----
  Ls <- try(multinet::layers_ml(net), silent = TRUE)
  if (inherits(Ls, "try-error") || is.null(Ls) || !length(Ls))
    stop("Could not get layers via multinet::layers_ml(net).")
  Ls <- as.character(Ls)
  layers <- if (is.null(layer_order)) Ls else {
    keep <- intersect(as.character(layer_order), Ls)
    if (!length(keep)) stop("None of the requested layers are present in 'net'.")
    keep
  }

  # ---- igraph list ----
  glist <- try(as.list(net), silent = TRUE)
  if (inherits(glist, "try-error") || !is.list(glist))
    stop("Could not coerce 'net' to list of igraphs via as.list(net).")
  if (!is.null(names(glist)) && any(names(glist) == "_flat_"))
    glist <- glist[names(glist) != "_flat_"]

  # ---- run folder ----
  .ensure_dir(results_dir)
  if (is.null(run_name) || !nzchar(run_name))
    run_name <- paste0("gprofiler_multinet_", .stamp())
  run_dir <- file.path(results_dir, run_name)
  .ensure_dir(run_dir)
  if (isTRUE(verbose))
    message("Run folder: ", normalizePath(run_dir, winslash = "/", mustWork = FALSE))

  # ---- background logic ----
  if (!is.null(custom_bg) && domain_scope != "custom") {
    message("custom_bg provided: switching domain_scope to 'custom'.")
    domain_scope <- "custom"
  }
  if (identical(domain_scope, "custom") && is.null(custom_bg))
    stop("domain_scope='custom' requires a non-NULL 'custom_bg' gene vector.")

  by_layer <- stats::setNames(vector("list", length(layers)), layers)
  summary_rows <- list()

  # ---- iterate layers ----
  for (ln in layers) {
    if (isTRUE(verbose)) message("Layer: ", ln)
    layer_dir <- if (isTRUE(per_layer_subdirs)) {
      d <- file.path(run_dir, ln); .ensure_dir(d); d
    } else run_dir

    g <- glist[[ln]]
    if (!inherits(g, "igraph")) { warning("  Skipping '", ln, "' (not an igraph)."); next }
    genes <- unique(as.character(igraph::V(g)$name))
    genes <- genes[nzchar(genes)]
    if (!length(genes)) { warning("  Skipping '", ln, "' (no vertex names)."); next }

    # ---- g:Profiler enrichment ----
    gres <- try(
      gprofiler2::gost(
        query                       = genes,
        organism                    = organism,
        ordered_query               = FALSE,
        multi_query                 = FALSE,
        significant                 = significant,
        exclude_iea                 = exclude_iea,
        measure_underrepresentation = FALSE,
        evcodes                     = TRUE,
        user_threshold              = user_threshold,
        correction_method           = correction_method,
        domain_scope                = domain_scope,
        custom_bg                   = if (!is.null(custom_bg)) custom_bg else NULL,
        sources                     = sources
      ),
      silent = TRUE
    )
    if (inherits(gres, "try-error") || is.null(gres) || is.null(gres$result) || !nrow(gres$result)) {
      if (isTRUE(verbose)) message("  No results for '", ln, "'.")
      by_layer[[ln]] <- list(csv_file=NA, dotplot_file=NA, cnet_file=NA, termnet_file=NA)
      next
    }

    # ---- CSV with intersection_genes ----
    inter_chr <- vapply(seq_len(nrow(gres$result)), function(i) {
      s <- .extract_symbols(gres$result$intersection[[i]])
      if (length(s)) paste(s, collapse=";") else NA_character_
    }, character(1))
    df_flat <- .flatten_df(gres$result)
    df_flat$intersection_genes <- inter_chr
    csv_file <- file.path(layer_dir, paste0("gprof_", ln, ".csv"))
    utils::write.csv(df_flat, csv_file, row.names = FALSE)
    if (isTRUE(verbose))
      message("  Saved CSV: ", normalizePath(csv_file, winslash = "/", mustWork = FALSE))

    # ---- gost dot plot (limited to top show_terms by p-value) ----
    dotplot_file <- NA_character_
    ord <- order(suppressWarnings(as.numeric(gres$result$p_value)), na.last = TRUE)
    gresTop <- gres
    gresTop$result <- gres$result[ord, , drop = FALSE]
    if (nrow(gresTop$result) > show_terms)
      gresTop$result <- head(gresTop$result, show_terms)
    gp <- try(gprofiler2::gostplot(gresTop, capped = TRUE, interactive = FALSE), silent = TRUE)
    if (!inherits(gp, "try-error")) {
      dotplot_file <- file.path(layer_dir, paste0("gprof_dot_", ln, ".", format))
      if (isTRUE(show_in_rstudio)) print(gp)
      ggplot2::ggsave(dotplot_file, plot = gp, width = width, height = height,
                      dpi = if (identical(format, "pdf")) NA_integer_ else dpi, bg = "white")
      if (isTRUE(verbose))
        message("  Saved dot plot: ", normalizePath(dotplot_file, winslash = "/", mustWork = FALSE))
    } else {
      warning("  Could not build gostplot for '", ln, "'.")
    }

    # ---- concept (term–gene) network ----
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
            d = data.frame(from = edges$term_id, to = edges$gene, stringsAsFactors = FALSE),
            directed = FALSE, vertices = nodes
          )

          cnet_file <- file.path(layer_dir, paste0("gprof_cnet_", ln, ".", format))
          pnet <- ggraph::ggraph(g_bi, layout = "fr") +
            ggraph::geom_edge_link(alpha = 0.25) +
            ggraph::geom_node_point(ggplot2::aes(size = size,
                                                 color = ifelse(nodes$node_type == "term", nodes$lp, NA_real_))) +
            ggplot2::scale_size_continuous(range = c(2, 9), guide = "none") +
            ggplot2::scale_color_viridis_c(na.value = "grey60", name = expression(-log[10](p))) +
            ggraph::geom_node_text(ggplot2::aes(label = label),
                                   size = 3, color = "black", check_overlap = TRUE) +
            ggplot2::labs(title = paste0("Term–gene concept network: ", ln)) +
            ggplot2::theme_minimal(base_size = 12) +
            ggplot2::theme(legend.position = "right")

          if (isTRUE(show_in_rstudio)) print(pnet)
          ggplot2::ggsave(cnet_file, plot = pnet, width = width, height = height,
                          dpi = if (identical(format, "pdf")) NA_integer_ else dpi, bg = "white")
          if (isTRUE(verbose))
            message("  Saved concept network: ", normalizePath(cnet_file, winslash = "/", mustWork = FALSE))
        } else if (isTRUE(verbose)) {
          message("  No term–gene intersections to draw for '", ln, "'.")
        }
      }
    } else if (show_terms > 0 && !has_net) {
      warning("  Skipping concept network for '", ln, "' (need ggraph + igraph).")
    }

    # ---- term-only overlap network ----
    termnet_file <- NA_character_
    if (isTRUE(ton_show) && has_net) {
      dfo <- head(gres$result[ord, , drop = FALSE], n = min(ton_top_terms, nrow(gres$result)))
      if (nrow(dfo) > 0) {
        term_list <- lapply(seq_len(nrow(dfo)), function(i) .extract_symbols(dfo$intersection[[i]]))
        names(term_list) <- dfo$term_id
        ids <- dfo$term_id

        e_from <- character(0); e_to <- character(0); e_w <- integer(0)
        if (length(ids) >= 2L) {
          for (pair in utils::combn(ids, 2, simplify = FALSE)) {
            gi <- term_list[[pair[1]]]; gj <- term_list[[pair[2]]]
            if (length(gi) && length(gj)) {
              ov <- length(intersect(gi, gj))
              if (ov >= ton_min_overlap) { e_from <- c(e_from, pair[1]); e_to <- c(e_to, pair[2]); e_w <- c(e_w, ov) }
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
            d = data.frame(from = e_from, to = e_to, weight = e_w, stringsAsFactors = FALSE),
            directed = FALSE, vertices = term_nodes
          )
        } else igraph::graph_from_data_frame(d = NULL, directed = FALSE, vertices = term_nodes)

        lay <- .layout_terms(g_terms, ton_layout)
        nodes_df <- cbind(term_nodes, data.frame(x = lay[,1], y = lay[,2]))
        edf <- if (length(e_from)) data.frame(
          from = e_from, to = e_to, weight = e_w,
          x1 = nodes_df$x[match(e_from, nodes_df$id)],
          y1 = nodes_df$y[match(e_from, nodes_df$id)],
          x2 = nodes_df$x[match(e_to,   nodes_df$id)],
          y2 = nodes_df$y[match(e_to,   nodes_df$id)]
        ) else data.frame(from=character(0), to=character(0), weight=integer(0),
                          x1=numeric(0), y1=numeric(0), x2=numeric(0), y2=numeric(0))

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
            data = nodes_df, ggplot2::aes(x = x, y = y, label = label),
            size = 4, color = "black", vjust = -0.8
          ) +
          ggplot2::geom_text(
            data = nodes_df, ggplot2::aes(x = x, y = y, label = k),
            size = 3.8, color = "grey20", vjust = 1.8
          ) +
          ggplot2::labs(title = paste0("Term-only overlap network: ", ln), x = "x", y = "y") +
          ggplot2::theme_minimal(base_size = 12) +
          ggplot2::theme(plot.margin = .unit_cm(ton_margin_cm))

        termnet_file <- file.path(layer_dir, paste0("gprof_termnet_", ln, ".", format))
        if (isTRUE(show_in_rstudio)) print(p_ton)
        ggplot2::ggsave(termnet_file, plot = p_ton, width = ton_width, height = ton_height,
                        dpi = if (identical(format, "pdf")) NA_integer_ else dpi, bg = "white")
        if (isTRUE(verbose))
          message("  Saved term-only overlap network: ", normalizePath(termnet_file, winslash = "/", mustWork = FALSE))
      }
    } else if (isTRUE(ton_show) && !has_net) {
      warning("  Skipping term-only network for '", ln, "' (need ggraph + igraph).")
    }

    by_layer[[ln]] <- list(
      csv_file     = csv_file,
      dotplot_file = if (!is.na(dotplot_file)) dotplot_file else NA_character_,
      cnet_file    = if (!is.na(cnet_file))    cnet_file    else NA_character_,
      termnet_file = if (!is.na(termnet_file)) termnet_file else NA_character_
    )
    summary_rows[[length(summary_rows)+1L]] <- data.frame(
      layer        = ln,
      n_genes      = length(genes),
      n_terms      = nrow(gres$result),
      csv_file     = csv_file,
      dotplot_file = if (!is.na(dotplot_file)) dotplot_file else NA_character_,
      cnet_file    = if (!is.na(cnet_file))    cnet_file    else NA_character_,
      termnet_file = if (!is.na(termnet_file)) termnet_file else NA_character_,
      stringsAsFactors = FALSE
    )
  }

  # ---- manifest ----
  summary_df <- if (length(summary_rows)) do.call(rbind, summary_rows) else
    data.frame(layer=character(), n_genes=integer(), n_terms=integer(),
               csv_file=character(), dotplot_file=character(),
               cnet_file=character(), termnet_file=character(),
               stringsAsFactors = FALSE)
  summary_file <- file.path(run_dir, "SUMMARY.csv")
  utils::write.csv(summary_df, summary_file, row.names = FALSE)
  message("Saved manifest: ", normalizePath(summary_file, winslash = "/", mustWork = FALSE))

  invisible(list(run_dir = run_dir, by_layer = by_layer, summary_file = summary_file))
}


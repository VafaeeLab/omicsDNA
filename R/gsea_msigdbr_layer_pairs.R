


# -------------------------------------------------------------------------
# 28 - Pairwise GSEA (MSigDB) between consecutive layers (fgsea + msigdbr)
# -------------------------------------------------------------------------

#' Pairwise GSEA (MSigDB) between consecutive network layers
#'
#' @title Pairwise GSEA (MSigDB) across consecutive multilayer network layers with dot, enrichment, cNET (terms+genes), and term‑only overlap plots
#'
#' @description
#' For each **consecutive pair** of layers in \code{layer_order} (e.g., \code{E1 -> E2 -> M1 -> M2}),
#' this function builds a **pairwise per‑gene ranking** from layer‑specific correlation and adjusted
#' p‑value, runs preranked GSEA via \pkg{fgsea}, and writes five outputs per pair:
#' \itemize{
#'   \item a tidy CSV of enriched pathways (ES, NES, pval, padj, size, leadingEdge, etc.);
#'   \item a **dot** (bubble) summary of top up/down pathways by adjusted p‑value;
#'   \item classic **running‑sum enrichment** plots for the best L1‑up and best L2‑up pathways;
#'   \item a **cNET (term–gene) concept network** where \strong{both term names and gene symbols are labelled}
#'         (terms are color‑coded by \eqn{-\log_{10}(adj\;p)} and sized by member count; genes are grey);
#'   \item a **term‑only overlap network** (nodes = terms, edges reflect shared genes; numbers at nodes are
#'         the pathway gene counts used; color encodes \eqn{-\log_{10}(adj\;p)}).
#' }
#'
#' @details
#' \strong{Workflow.}
#' \enumerate{
#'   \item Layers are taken from \code{multinet::layers_ml(net)} (or from \code{layer_order} if supplied);
#'         consecutive pairs \code{(L1, L2)} are formed in that order.
#'   \item Per layer, per‑gene statistics come either from vertex attributes (first match among
#'         \code{cor_attr}, \code{padj_attr}) or from \code{scores_by_layer[[layer]]}.
#'   \item A **pairwise ranking** \eqn{r_g} (for gene \eqn{g}) is computed using \code{pair_rank_by}:
#'         \describe{
#'           \item{\code{"delta_signed_log10_padj"} (default)}{
#'             \eqn{\mathrm{sign}(cor_{L1}) \cdot -\log_{10}(padj_{L1}) - \mathrm{sign}(cor_{L2}) \cdot -\log_{10}(padj_{L2})}
#'           }
#'           \item{\code{"delta_cor"}}{\eqn{cor_{L1} - cor_{L2}}}
#'           \item{\code{"signed_log10_padj_L1"}}{\eqn{\mathrm{sign}(cor_{L1}) \cdot -\log_{10}(padj_{L1})}}
#'           \item{\code{"signed_log10_padj_L2"}}{\eqn{\mathrm{sign}(cor_{L2}) \cdot -\log_{10}(padj_{L2})}}
#'           \item{\code{"cor_L1"}}{\eqn{cor_{L1}}}
#'           \item{\code{"cor_L2"}}{\eqn{cor_{L2}}}
#'         }
#'   \item MSigDB pathways (via \pkg{msigdbr}) are intersected with the ranked universe, size‑filtered
#'         by \code{minSize}/\code{maxSize}, and enriched using:
#'         \itemize{
#'           \item \strong{\code{fgseaMultilevel()}} when \code{nperm = NULL} (recommended);
#'           \item classic \strong{\code{fgsea()}} when \code{nperm} is a positive integer.
#'         }
#'   \item Plots/CSVs are saved into a timestamped run folder; when \code{show_in_rstudio = TRUE} they are
#'         also printed to the active device.
#' }
#'
#' \strong{MSigDB compatibility.} Around msigdbr v10 the API shifted from “category” to
#' “collection[/subcategory]”. This function accepts \code{msigdb_collections} and optional
#' \code{msigdb_subcategories}, and internally adapts to both styles. If a requested subcategory
#' is not available for a collection, a short note is emitted and processing continues.
#'
#' \strong{cNET (term–gene).} Term nodes are colored by \eqn{-\log_{10}(adj\;p)} and sized by the number
#' of connected genes (leading edge or full members per \code{cnet_use}). \emph{Gene symbols are always
#' labelled} (smaller, grey) alongside term labels (larger). Positions are computed via \pkg{ggraph}
#' and labels use \pkg{ggrepel} when available.
#'
#' \strong{Term‑only overlap.} Nodes are the selected top terms (by \code{padj}); edges connect terms
#' sharing at least \code{ton_min_overlap} genes. Numbers near nodes give the gene counts used.
#'
#' \strong{Output structure.} For a pair \code{L1_vs_L2}:
#' \preformatted{
#'   <results_dir>/<run_name>_<YYYY-mm-dd_HHMMSS>/<L1_vs_L2>/
#'     gsea_pair_<L1_vs_L2>.csv
#'     gsea_pair_dot_<L1_vs_L2>.<format>
#'     gsea_pair_enrich_topL1_<L1_vs_L2>.<format>
#'     gsea_pair_enrich_topL2_<L1_vs_L2>.<format>
#'     gsea_pair_cnet_<L1_vs_L2>.<format>
#'     gsea_pair_termnet_<L1_vs_L2>.<format>
#' }
#' Plus a run‑level manifest:
#' \preformatted{
#'   <results_dir>/<run_name>_<YYYY-mm-dd_HHMMSS>/SUMMARY.csv
#' }
#'
#' @section Units, scaling and labels:
#' \itemize{
#'   \item All \code{*_margin_cm} arguments are in **centimetres** and applied as plot margins.
#'   \item Network node sizes are rescaled to a sensible visual range (terms > genes).
#'   \item cNET always prints both term names and gene symbols (genes smaller, grey).
#' }
#'
#' @param net A \pkg{multinet} object. Layers come from \code{multinet::layers_ml(net)} and are coerced
#'   to \pkg{igraph}s via \code{as.list(net)}. Vertex attributes provide \code{cor}/\code{padj} unless
#'   overridden by \code{scores_by_layer}.
#' @param layer_order Character vector giving the traversal order of layers. Pairs are formed as
#'   \code{(L1,L2)}, \code{(L2,L3)}, ... in this order. Default: \code{multinet::layers_ml(net)}.
#' @param scores_by_layer Optional named list \code{layer -> data.frame(gene, cor, padj)} that overrides
#'   vertex attributes for the specified layers.
#' @param cor_attr,padj_attr Character vectors of candidate vertex attribute names for correlation and
#'   adjusted p‑value. The first present in a layer is used (case‑insensitive).
#' @param pair_rank_by Ranking recipe. One of:
#'   \code{"delta_signed_log10_padj"} (default), \code{"delta_cor"},
#'   \code{"signed_log10_padj_L1"}, \code{"signed_log10_padj_L2"}, \code{"cor_L1"}, \code{"cor_L2"}.
#'   See Details for formulas.
#' @param missing_policy \code{"impute"} (default; uses \code{impute_cor=0}, \code{impute_padj=1}) or
#'   \code{"drop"} to remove genes with missing values.
#' @param impute_cor,impute_padj Numeric values used under \code{missing_policy="impute"}.
#' @param case_insensitive Logical; if \code{TRUE} (default) gene symbols are upper‑cased for matching.
#' @param msigdb_species Species string for \pkg{msigdbr} (e.g., \code{"Homo sapiens"}).
#' @param msigdb_collections Character vector of MSigDB collections (e.g., \code{c("H","C2","C5")}).
#' @param msigdb_subcategories Optional vector of subcategories (e.g., \code{c("CP:KEGG","GO:BP")}).
#' @param minSize,maxSize Integer bounds for pathway sizes after intersecting with the ranked universe.
#' @param nperm If \code{NULL} (default), uses \code{fgseaMultilevel()} (recommended). If a positive
#'   integer, uses \code{fgsea()} with \code{nperm} permutations.
#' @param results_dir Base output directory (created if needed).
#' @param run_name Prefix for the run folder (a timestamp \code{"_YYYY-mm-dd_HHMMSS"} is appended).
#' @param top_n Number of top up/down pathways by adjusted p‑value to show in the dot plot.
#' @param format Image format for plots: \code{"png"}, \code{"pdf"} or \code{"jpg"} (DPI ignored for PDF).
#' @param width,height,dpi Default plot geometry (inches) and DPI (when per‑plot sizes not given).
#' @param show_in_rstudio If \code{TRUE}, also prints plots to the current device.
#' @param dotplot_width,dotplot_height,enrich_width,enrich_height Per‑plot sizes (inches) for dot and
#'   enrichment plots. Defaults fall back to \code{width}/\code{height}.
#' @param dotplot_margin_cm,enrich_margin_cm Numeric vectors \code{c(top,right,bottom,left)} in **cm**
#'   specifying plot margins.
#' @param cnet_show If \code{TRUE} (default), build the cNET plot.
#' @param cnet_top_terms Number of top terms (by \code{padj}) to include in cNET.
#' @param cnet_use \code{"leadingEdge"} (default; connect terms to FGSEA leading‑edge genes) or
#'   \code{"members"} (connect all member genes present in the ranked universe).
#' @param cnet_layout Layout name passed to \pkg{ggraph} (e.g., \code{"fr"}, \code{"kk"}, \code{"lgl"}).
#' @param cnet_width,cnet_height Size (inches) for the cNET plot.
#' @param cnet_margin_cm Numeric vector \code{c(top,right,bottom,left)} in **cm** for cNET margins.
#' @param cnet_term_label_size,cnet_gene_label_size Numeric text sizes for term and gene labels in cNET.
#' @param cnet_gene_label_color Color for gene labels in cNET (default \code{"grey25"}).
#' @param ton_show If \code{TRUE} (default), build the term‑only overlap network.
#' @param ton_top_terms Number of top terms (by \code{padj}) to include in the term‑only network.
#' @param ton_min_overlap Minimum number of shared genes to draw a term–term edge.
#' @param ton_layout Layout for the term‑only network: \code{"fr"}, \code{"kk"}, \code{"lgl"}, or \code{"mds"}.
#' @param ton_width,ton_height Size (inches) for the term‑only plot.
#' @param ton_margin_cm Numeric vector \code{c(top,right,bottom,left)} in **cm** for term‑only margins.
#' @param seed Optional integer seed for RNG/reproducible layouts.
#' @param verbose If \code{TRUE} (default), print progress messages.
#'
#' @return
#' A list (returned invisibly) with:
#' \describe{
#'   \item{\code{run_dir}}{Path to the timestamped run folder.}
#'   \item{\code{by_pair}}{Named list by \code{"<L1>_vs_<L2>"} containing file paths:
#'         \code{csv_file}, \code{dotplot_file}, \code{enrich_L1_file}, \code{enrich_L2_file},
#'         \code{cnet_file}, \code{termnet_file}.}
#'   \item{\code{summary_file}}{Path to the manifest \code{SUMMARY.csv}.}
#' }
#'
#' @section Notes and troubleshooting:
#' \itemize{
#'   \item Messages like “no rows for collection/subcategory” mean that subcategory is absent for
#'         your \pkg{msigdbr} build; processing continues with available sets.
#'   \item With \code{fgseaMultilevel()}, rows may contain \code{log2err} and omit \code{nMoreExtreme};
#'         the CSV includes both, filling missing values with \code{NA}.
#'   \item Dense cNETs benefit from increasing \code{cnet_width}/\code{cnet_height} or label sizes.
#' }
#'
#' @seealso \code{\link[fgsea]{fgseaMultilevel}}, \code{\link[fgsea]{fgsea}}, \code{\link[fgsea]{plotEnrichment}},
#'   \code{\link[msigdbr]{msigdbr}}, \code{\link[ggraph]{ggraph}}
#'
#' @importFrom fgsea fgsea fgseaMultilevel plotEnrichment
#' @importFrom msigdbr msigdbr
#' @importFrom ggplot2 ggplot aes geom_point scale_size_continuous scale_color_gradient2 labs theme_minimal theme element_text element_blank ggsave
#' @importFrom igraph graph_from_data_frame layout_with_fr layout_with_kk layout_with_lgl layout_with_mds vertex_attr vertex_attr_names V
#' @importFrom utils write.csv combn
#' @importFrom stats setNames
#'
#' @examples
#' \dontrun{
#' ## Single, fully-specified example (demonstrates every argument)
#' ## Assume 'net' is your multinet with vertex attributes 'cor' and 'padj'
#' ## for each layer; otherwise construct 'scores_by_layer' as below.
#'
#' # Optional override of per-layer stats:
#' # genes <- c("TP53","EGFR","MYC","STAT3","BRCA1","BRCA2")
#' # scores_by_layer <- list(
#' #   E1 = data.frame(gene=genes, cor=rnorm(length(genes)), padj=pmin(1,runif(length(genes)))),
#' #   E2 = data.frame(gene=genes, cor=rnorm(length(genes)), padj=pmin(1,runif(length(genes)))),
#' #   M1 = data.frame(gene=genes, cor=rnorm(length(genes)), padj=pmin(1,runif(length(genes))))
#' # )
#'
#' res <- gsea_msigdbr_layer_pairs(
#'   net                  = net,
#'   layer_order          = c("E1","E2","M1","M2"),
#'   scores_by_layer      = NULL,                          # or use list above
#'   cor_attr             = c("cor","r","rho","correlation"),
#'   padj_attr            = c("padj","FDR","q","adj_p","p_adj","p.adj","p_fdr"),
#'   pair_rank_by         = "delta_signed_log10_padj",
#'   missing_policy       = "impute",                      # or "drop"
#'   impute_cor           = 0,
#'   impute_padj          = 1,
#'   case_insensitive     = TRUE,
#'   msigdb_species       = "Homo sapiens",
#'   msigdb_collections   = c("H","C2","C5"),
#'   msigdb_subcategories = c("CP:KEGG","GO:BP"),
#'   minSize              = 15,
#'   maxSize              = 500,
#'   nperm                = NULL,                          # NULL => fgseaMultilevel(); set e.g. 10000 for fgsea()
#'   results_dir          = getOption("mlnet.results_dir","omicsDNA_results"),
#'   run_name             = "gsea_pairs_demo",
#'   top_n                = 10,
#'   format               = "png",
#'   width                = 9, height = 7, dpi = 300,
#'   show_in_rstudio      = TRUE,
#'   # dot/enrich sizes & margins (cm)
#'   dotplot_width        = 8,  dotplot_height  = 6,
#'   enrich_width         = 8,  enrich_height   = 6,
#'   dotplot_margin_cm    = c(0.8,0.8,0.8,0.8),
#'   enrich_margin_cm     = c(0.8,0.8,0.8,0.8),
#'   # cNET (terms + genes always labelled)
#'   cnet_show            = TRUE,
#'   cnet_top_terms       = 12,
#'   cnet_use             = "leadingEdge",                 # or "members"
#'   cnet_layout          = "fr",
#'   cnet_width           = 10, cnet_height = 8,
#'   cnet_margin_cm       = c(0.8,0.8,0.8,0.8),
#'   cnet_term_label_size = 3.2,
#'   cnet_gene_label_size = 2.4,
#'   cnet_gene_label_color= "grey25",
#'   # term-only overlap network
#'   ton_show             = TRUE,
#'   ton_top_terms        = 12,
#'   ton_min_overlap      = 1,
#'   ton_layout           = "fr",
#'   ton_width            = 10, ton_height = 7,
#'   ton_margin_cm        = c(0.8,0.8,0.8,0.8),
#'   seed                 = 1,
#'   verbose              = TRUE
#' )
#'
#' # See the run manifest:
#' read.csv(res$summary_file, stringsAsFactors = FALSE)
#' }
#' @export
gsea_msigdbr_layer_pairs <- function(
    net,
    layer_order             = NULL,
    scores_by_layer         = NULL,
    cor_attr                = c("cor","r","rho","correlation"),
    padj_attr               = c("padj","fdr","q","adj_p","p_adj","p.adj","p_fdr"),
    pair_rank_by            = c("delta_signed_log10_padj","delta_cor",
                                "signed_log10_padj_L1","signed_log10_padj_L2",
                                "cor_L1","cor_L2"),
    missing_policy          = c("impute","drop"),
    impute_cor              = 0,
    impute_padj             = 1,
    case_insensitive        = TRUE,
    msigdb_species          = "Homo sapiens",
    msigdb_collections      = c("H","C2","C5"),
    msigdb_subcategories    = NULL,
    minSize                 = 15,
    maxSize                 = 500,
    # If NULL -> fgseaMultilevel (recommended). If integer -> classic permutations.
    nperm                   = NULL,
    results_dir             = getOption("mlnet.results_dir","omicsDNA_results"),
    run_name                = "gsea_layerpairs",
    top_n                   = 10,
    format                  = c("png","pdf","jpg"),
    width                   = 9,  height = 7,  dpi = 300,
    show_in_rstudio         = TRUE,
    # sizes & margins
    dotplot_width           = width,   dotplot_height  = height,
    enrich_width            = width,   enrich_height   = height,
    dotplot_margin_cm       = c(0.5,0.5,0.5,0.5),
    enrich_margin_cm        = c(0.5,0.5,0.5,0.5),
    # cNET
    cnet_show               = TRUE,
    cnet_top_terms          = 10,
    cnet_use                = c("leadingEdge","members"),
    cnet_layout             = "fr",
    cnet_width              = width,   cnet_height = height,
    cnet_margin_cm          = c(0.5,0.5,0.5,0.5),
    cnet_term_label_size    = 3.2,
    cnet_gene_label_size    = 2.4,
    cnet_gene_label_color   = "grey25",
    # term-only network
    ton_show                = TRUE,
    ton_top_terms           = 10,
    ton_min_overlap         = 1,
    ton_layout              = "fr",
    ton_width               = width,   ton_height = height,
    ton_margin_cm           = c(0.5,0.5,0.5,0.5),
    seed                    = NULL,
    verbose                 = TRUE
) {
  # ---- dependencies ----
  if (!requireNamespace("multinet", quietly = TRUE))
    stop("Please install 'multinet'.")
  if (!requireNamespace("fgsea", quietly = TRUE))
    stop("Please install 'fgsea' (Bioconductor).")
  if (!requireNamespace("msigdbr", quietly = TRUE))
    stop("Please install 'msigdbr'.")
  if (!requireNamespace("igraph", quietly = TRUE))
    stop("Please install 'igraph'.")
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Please install 'ggplot2'.")
  has_ggraph  <- requireNamespace("ggraph",  quietly = TRUE)
  has_ggrepel <- requireNamespace("ggrepel", quietly = TRUE)

  format         <- match.arg(format)
  pair_rank_by   <- match.arg(pair_rank_by)
  missing_policy <- match.arg(missing_policy)
  cnet_use       <- match.arg(cnet_use)

  # ---- helpers ----
  .ensure_dir <- function(d) if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  .stamp      <- function() format(Sys.time(), "%Y-%m-%d_%H%M%S")
  .normsym    <- function(x) { x <- as.character(x); if (isTRUE(case_insensitive)) toupper(x) else x }
  .pick_attr  <- function(avail, cands) { al <- tolower(avail); for (c in cands) { w <- which(al == tolower(c)); if (length(w)) return(avail[w[1]]) }; NA_character_ }
  .flatten_listcol <- function(x) { if (is.null(x)) return(NA_character_); if (is.list(x)) x <- unlist(x, use.names = FALSE); paste(unique(as.character(x)), collapse = ";") }
  .signed_log10 <- function(p) -log10(pmax(p, .Machine$double.xmin))
  .rescale <- function(x, to = c(4,10)) { x <- suppressWarnings(as.numeric(x)); r <- range(x, finite = TRUE); if (!is.finite(r[1]) || !is.finite(r[2]) || r[1]==r[2]) rep(mean(to), length(x)) else (x-r[1])/(r[2]-r[1])*(to[2]-to[1]) + to[1] }
  .unit_cm <- function(v) grid::unit(v, "cm")
  .layout_terms <- function(g, method = "fr") {
    switch(tolower(method),
           fr  = igraph::layout_with_fr(g),
           kk  = igraph::layout_with_kk(g),
           lgl = igraph::layout_with_lgl(g),
           mds = igraph::layout_with_mds(g),
           igraph::layout_with_fr(g))
  }
  .extract_layer_scores <- function(g, cor_attr, padj_attr) {
    if (!inherits(g, "igraph")) return(NULL)
    vnames <- igraph::V(g)$name; if (is.null(vnames)) return(NULL)
    vnames <- .normsym(vnames)
    vAttrs <- igraph::vertex_attr_names(g)
    cor_nm  <- .pick_attr(vAttrs, cor_attr)
    padj_nm <- .pick_attr(vAttrs, padj_attr)
    cor_vec  <- if (!is.na(cor_nm))  suppressWarnings(as.numeric(igraph::vertex_attr(g, cor_nm)))  else rep(NA_real_, length(vnames))
    padj_vec <- if (!is.na(padj_nm)) suppressWarnings(as.numeric(igraph::vertex_attr(g, padj_nm))) else rep(NA_real_, length(vnames))
    names(cor_vec)  <- vnames; names(padj_vec) <- vnames
    list(cor = cor_vec, padj = padj_vec)
  }

  # --- msigdbr wrapper (>=10 and older) ---
  .get_msigdb_df <- function(species, collections, subcats = NULL, verbose = TRUE) {
    fmls <- names(formals(msigdbr::msigdbr))
    uses_collection <- "collection" %in% fmls
    has_subcat_arg  <- "subcategory" %in% fmls || "subcollection" %in% fmls
    pieces <- list()
    for (coll in collections) {
      if (has_subcat_arg && !is.null(subcats) && length(subcats)) {
        for (sc in subcats) {
          args <- list(species = species)
          if (uses_collection) args$collection <- coll else args$category <- coll
          if ("subcategory"   %in% fmls) args$subcategory   <- sc
          if ("subcollection" %in% fmls) args$subcollection <- sc
          res <- try(do.call(msigdbr::msigdbr, args), silent = TRUE)
          if (!inherits(res, "try-error") && !is.null(res) && nrow(res)) pieces[[length(pieces)+1L]] <- res
          else if (isTRUE(verbose)) message("  (msigdbr) no rows for collection=", coll, " subcategory=", sc)
        }
      } else {
        args <- list(species = species)
        if (uses_collection) args$collection <- coll else args$category <- coll
        res <- try(do.call(msigdbr::msigdbr, args), silent = TRUE)
        if (!inherits(res, "try-error") && !is.null(res) && nrow(res)) pieces[[length(pieces)+1L]] <- res
        else if (isTRUE(verbose)) message("  (msigdbr) no rows for collection=", coll)
      }
    }
    if (!length(pieces)) return(NULL)
    if (requireNamespace("dplyr", quietly = TRUE)) {
      mdf <- dplyr::bind_rows(pieces)
    } else {
      cols <- unique(unlist(lapply(pieces, names)))
      pieces2 <- lapply(pieces, function(df) { miss <- setdiff(cols, names(df)); if (length(miss)) df[miss] <- NA; df[, cols, drop = FALSE] })
      mdf <- do.call(rbind, pieces2)
    }
    if (!is.null(subcats) && length(subcats)) {
      cand_cols <- intersect(c("gs_subcat","subcategory","gs_subcategory","gs_subcollection"), names(mdf))
      if (length(cand_cols)) {
        subcol <- cand_cols[1]
        mdf <- mdf[mdf[[subcol]] %in% subcats, , drop = FALSE]
      } else if (isTRUE(verbose)) message("  (msigdbr) no subcategory column; ignoring filter.")
    }
    mdf
  }

  # ---- layers & pairs ----
  Ls <- try(multinet::layers_ml(net), silent = TRUE)
  if (inherits(Ls, "try-error") || is.null(Ls) || !length(Ls))
    stop("Could not retrieve layers via multinet::layers_ml(net).")
  Ls <- as.character(Ls)
  layers <- if (is.null(layer_order)) Ls else {
    keep <- intersect(as.character(layer_order), Ls)
    if (!length(keep)) stop("None of the requested layers are present in 'net'.")
    keep
  }
  if (length(layers) < 2L) stop("Need at least two layers in 'layer_order' to form pairs.")
  pairs_df <- data.frame(L1 = head(layers, -1L), L2 = tail(layers, -1L), stringsAsFactors = FALSE)

  # ---- gene sets ----
  if (isTRUE(verbose)) message("Fetching MSigDB gene sets (", msigdb_species, "; ", paste(msigdb_collections, collapse = ","), ") ...")
  mdf <- .get_msigdb_df(msigdb_species, msigdb_collections, msigdb_subcategories, verbose = verbose)
  if (is.null(mdf) || !nrow(mdf)) stop("msigdbr returned no rows for the requested species/collections/subcategories.")
  mdf$gs_name     <- as.character(mdf$gs_name)
  mdf$gene_symbol <- .normsym(mdf$gene_symbol)
  pathways_all <- split(mdf$gene_symbol, mdf$gs_name)

  # ---- run folder (timestamped) ----
  .ensure_dir(results_dir)
  run_dir <- file.path(results_dir, paste0(run_name, "_", .stamp()))
  .ensure_dir(run_dir)
  if (isTRUE(verbose)) message("Run folder: ", normalizePath(run_dir, winslash = "/", mustWork = FALSE))

  # ---- per-layer scores ----
  glist <- try(as.list(net), silent = TRUE)
  if (inherits(glist, "try-error") || !is.list(glist))
    stop("Could not coerce 'net' to list of igraphs via as.list(net).")
  if (!is.null(names(glist)) && any(names(glist) == "_flat_"))
    glist <- glist[names(glist) != "_flat_"]

  get_scores <- function(ln) {
    if (!is.null(scores_by_layer) && !is.null(scores_by_layer[[ln]])) {
      df <- scores_by_layer[[ln]]
      if (!all(c("gene","cor","padj") %in% names(df)))
        stop("scores_by_layer[['", ln, "']] must contain columns: gene, cor, padj")
      gsym  <- .normsym(df$gene)
      corv  <- suppressWarnings(as.numeric(df$cor));  names(corv)  <- gsym
      padjv <- suppressWarnings(as.numeric(df$padj)); names(padjv) <- gsym
      return(list(cor = corv, padj = padjv))
    } else .extract_layer_scores(glist[[ln]], cor_attr, padj_attr)
  }

  if (!is.null(seed)) set.seed(as.integer(seed))

  by_pair <- stats::setNames(vector("list", nrow(pairs_df)), paste0(pairs_df$L1, "_vs_", pairs_df$L2))
  summary_rows <- list()

  # ---- loop over pairs ----
  for (i in seq_len(nrow(pairs_df))) {
    L1 <- pairs_df$L1[i]; L2 <- pairs_df$L2[i]
    pair_name <- paste0(L1, "_vs_", L2)
    if (isTRUE(verbose)) message("Pair: ", L1, " vs ", L2)
    pair_dir <- file.path(run_dir, pair_name); .ensure_dir(pair_dir)

    s1 <- get_scores(L1); s2 <- get_scores(L2)
    if (is.null(s1) || is.null(s2)) { warning("  Missing scores for pair ", pair_name, "; skipping."); next }

    # union & missing policy
    genes_all <- unique(c(names(s1$cor), names(s2$cor))); genes_all <- genes_all[nzchar(genes_all)]
    cor1  <- s1$cor[genes_all];  cor2  <- s2$cor[genes_all]
    padj1 <- s1$padj[genes_all]; padj2 <- s2$padj[genes_all]
    if (identical(missing_policy, "impute")) {
      cor1[!is.finite(cor1)]   <- impute_cor
      cor2[!is.finite(cor2)]   <- impute_cor
      padj1[!is.finite(padj1)] <- impute_padj
      padj2[!is.finite(padj2)] <- impute_padj
    } else {
      keep <- is.finite(cor1) & is.finite(cor2) & is.finite(padj1) & is.finite(padj2)
      genes_all <- genes_all[keep]; cor1 <- cor1[keep]; cor2 <- cor2[keep]; padj1 <- padj1[keep]; padj2 <- padj2[keep]
    }

    # ranking
    score <- switch(pair_rank_by,
                    delta_signed_log10_padj = (sign(cor1) * .signed_log10(padj1)) - (sign(cor2) * .signed_log10(padj2)),
                    delta_cor               = cor1 - cor2,
                    signed_log10_padj_L1    = sign(cor1) * .signed_log10(padj1),
                    signed_log10_padj_L2    = sign(cor2) * .signed_log10(padj2),
                    cor_L1                  = cor1,
                    cor_L2                  = cor2
    )
    names(score) <- genes_all
    score <- sort(score[is.finite(score)], decreasing = TRUE)

    # pathways in universe & size bounds
    pathways <- lapply(pathways_all, function(x) intersect(x, names(score)))
    pathways <- pathways[vapply(pathways, function(x) length(x) >= minSize && length(x) <= maxSize, logical(1))]
    if (!length(pathways)) { warning("  No pathways within size bounds for pair ", pair_name, "."); next }

    # fgsea
    if (is.null(nperm)) {
      fg <- fgsea::fgseaMultilevel(pathways = pathways, stats = score,
                                   minSize = minSize, maxSize = maxSize, scoreType = "std")
    } else {
      fg <- fgsea::fgsea(pathways = pathways, stats = score,
                         minSize = minSize, maxSize = maxSize,
                         nperm = nperm, scoreType = "std")
    }
    if (!nrow(fg)) { warning("  fgsea returned no results for pair ", pair_name, "."); next }

    # CSV (robust to fgseaMultilevel missing columns like nMoreExtreme)
    fg$direction <- ifelse(fg$NES >= 0, "L1_up", "L2_up")
    le_str <- vapply(fg$leadingEdge, .flatten_listcol, character(1))
    nme <- if ("nMoreExtreme" %in% names(fg)) fg$nMoreExtreme else rep(NA_integer_, nrow(fg))
    log2err <- if ("log2err" %in% names(fg)) fg$log2err else rep(NA_real_, nrow(fg))

    out_tbl <- data.frame(
      pathway      = fg$pathway,
      size         = fg$size,
      ES           = fg$ES,
      NES          = fg$NES,
      pval         = fg$pval,
      padj         = fg$padj,
      log2err      = log2err,
      nMoreExtreme = nme,
      direction    = fg$direction,
      leadingEdge_genes = le_str,
      stringsAsFactors = FALSE
    )
    csv_file <- file.path(pair_dir, paste0("gsea_pair_", pair_name, ".csv"))
    utils::write.csv(out_tbl, csv_file, row.names = FALSE)
    if (isTRUE(verbose)) message("  Saved CSV: ", normalizePath(csv_file, winslash = "/", mustWork = FALSE))

    # ---- dot plot ----
    fg_up   <- fg[fg$NES > 0, , drop = FALSE]
    fg_down <- fg[fg$NES < 0, , drop = FALSE]
    ord_up   <- order(fg_up$padj, fg_up$pval, -fg_up$NES, na.last = NA)
    ord_down <- order(fg_down$padj, fg_down$pval,  fg_down$NES, na.last = NA)
    sel <- rbind(fg_up[head(ord_up,   top_n), , drop = FALSE],
                 fg_down[head(ord_down, top_n), , drop = FALSE])
    dotplot_file <- NA_character_
    if (nrow(sel)) {
      sel$logpadj <- -log10(pmax(sel$padj, .Machine$double.eps))
      sel$pathway <- factor(sel$pathway, levels = sel$pathway[order(sel$NES, decreasing = TRUE)])
      p <- ggplot2::ggplot(sel, ggplot2::aes(x = NES, y = pathway)) +
        ggplot2::geom_point(ggplot2::aes(size = logpadj, color = NES)) +
        ggplot2::scale_size_continuous(name = expression(-log[10](adj.~p)), range = c(2,9)) +
        ggplot2::scale_color_gradient2(low = "#2c7bb6", mid = "grey70", high = "#d7191c", midpoint = 0) +
        ggplot2::labs(title = paste0("GSEA (", L1, " vs ", L2, ")"),
                      x = "NES (positive = L1)", y = "Pathway") +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(panel.grid = ggplot2::element_blank(),
                       plot.margin = .unit_cm(dotplot_margin_cm),
                       axis.text.y = ggplot2::element_text(size = 9))
      dotplot_file <- file.path(pair_dir, paste0("gsea_pair_dot_", pair_name, ".", format))
      if (isTRUE(show_in_rstudio)) print(p)
      ggplot2::ggsave(dotplot_file, plot = p, width = dotplot_width, height = dotplot_height,
                      dpi = if (identical(format, "pdf")) NA_integer_ else dpi, bg = "white")
      if (isTRUE(verbose)) message("  Saved dot plot: ", normalizePath(dotplot_file, winslash = "/", mustWork = FALSE))
    }

    # ---- classic enrichment plots ----
    enr_up_file <- enr_dn_file <- NA_character_
    if (nrow(fg_up)) {
      top_up <- fg_up[order(fg_up$padj, fg_up$pval, -fg_up$NES), , drop = FALSE][1, ]
      gl_up  <- pathways[[as.character(top_up$pathway)]]
      pe <- fgsea::plotEnrichment(gl_up, score) +
        ggplot2::ggtitle(paste0("Enrichment: ", top_up$pathway, " (", L1, " vs ", L2, ")")) +
        ggplot2::theme(plot.margin = .unit_cm(enrich_margin_cm))
      enr_up_file <- file.path(pair_dir, paste0("gsea_pair_enrich_topL1_", pair_name, ".", format))
      if (isTRUE(show_in_rstudio)) print(pe)
      ggplot2::ggsave(enr_up_file, plot = pe, width = enrich_width, height = enrich_height,
                      dpi = if (identical(format, "pdf")) NA_integer_ else dpi, bg = "white")
      if (isTRUE(verbose)) message("  Saved enrichment plot (L1-up): ", normalizePath(enr_up_file, winslash = "/", mustWork = FALSE))
    }
    if (nrow(fg_down)) {
      top_dn <- fg_down[order(fg_down$padj, fg_down$pval,  fg_down$NES), , drop = FALSE][1, ]
      gl_dn  <- pathways[[as.character(top_dn$pathway)]]
      pe2 <- fgsea::plotEnrichment(gl_dn, score) +
        ggplot2::ggtitle(paste0("Enrichment: ", top_dn$pathway, " (", L1, " vs ", L2, ")")) +
        ggplot2::theme(plot.margin = .unit_cm(enrich_margin_cm))
      enr_dn_file <- file.path(pair_dir, paste0("gsea_pair_enrich_topL2_", pair_name, ".", format))
      if (isTRUE(show_in_rstudio)) print(pe2)
      ggplot2::ggsave(enr_dn_file, plot = pe2, width = enrich_width, height = enrich_height,
                      dpi = if (identical(format, "pdf")) NA_integer_ else dpi, bg = "white")
      if (isTRUE(verbose)) message("  Saved enrichment plot (L2-up): ", normalizePath(enr_dn_file, winslash = "/", mustWork = FALSE))
    }

    # ---- cNET: term–gene concept network (labels: terms + genes) ----
    cnet_file <- NA_character_
    if (isTRUE(cnet_show) && has_ggraph) {
      fg_sel <- fg[order(fg$padj, fg$pval, -abs(fg$NES)), , drop = FALSE]
      fg_sel <- fg_sel[is.finite(fg_sel$padj), , drop = FALSE]
      fg_sel <- head(fg_sel, n = min(cnet_top_terms, nrow(fg_sel)))
      if (nrow(fg_sel)) {
        get_genes_for <- function(term) {
          if (identical(cnet_use, "leadingEdge")) {
            i <- match(term, fg_sel$pathway)
            unique(as.character(unlist(fg_sel$leadingEdge[i], use.names = FALSE)))
          } else {
            unique(as.character(pathways[[term]]))
          }
        }
        edges_list <- lapply(as.character(fg_sel$pathway), function(tn) {
          gi <- intersect(get_genes_for(tn), names(score))
          if (length(gi)) data.frame(term = tn, gene = gi, stringsAsFactors = FALSE) else NULL
        })
        edges_list <- Filter(Negate(is.null), edges_list)
        edges <- if (length(edges_list)) do.call(rbind, edges_list) else NULL

        term_names <- as.character(fg_sel$pathway)
        term_padj  <- as.numeric(fg_sel$padj)
        term_sizes <- if (!is.null(edges) && nrow(edges)) vapply(term_names, function(tn) sum(edges$term == tn), integer(1)) else rep(0L, length(term_names))

        nodes_term <- data.frame(
          id        = term_names,
          label     = term_names,
          node_type = "term",
          lp        = -log10(pmax(term_padj, .Machine$double.xmin)),
          size      = .rescale(term_sizes, to = c(5,11)),
          stringsAsFactors = FALSE
        )
        genes_unique <- if (!is.null(edges) && nrow(edges)) unique(edges$gene) else character(0)
        nodes <- if (length(genes_unique)) {
          rbind(
            nodes_term,
            data.frame(id = genes_unique, label = genes_unique, node_type = "gene",
                       lp = NA_real_, size = 3, stringsAsFactors = FALSE)
          )
        } else nodes_term

        g_bi <- igraph::graph_from_data_frame(
          d = if (!is.null(edges) && nrow(edges))
            data.frame(from = edges$term, to = edges$gene, stringsAsFactors = FALSE)
          else NULL,
          directed = FALSE, vertices = nodes
        )
        lay <- ggraph::create_layout(g_bi, layout = cnet_layout)

        cnet_file <- file.path(pair_dir, paste0("gsea_pair_cnet_", pair_name, ".", format))

        lay_terms <- subset(lay, node_type == "term")
        lay_genes <- subset(lay, node_type == "gene")

        pnet <- ggraph::ggraph(lay) +
          ggraph::geom_edge_link(alpha = 0.25, colour = "grey70") +
          ggplot2::geom_point(data = lay_genes, ggplot2::aes(x = x, y = y, size = size),
                              inherit.aes = FALSE, colour = "grey55", alpha = 0.9) +
          ggplot2::geom_point(data = lay_terms, ggplot2::aes(x = x, y = y, size = size, color = lp),
                              inherit.aes = FALSE) +
          ggplot2::scale_size_continuous(range = c(2, 10), guide = "none") +
          ggplot2::scale_color_viridis_c(name = expression(-log[10](adj.~p))) +
          ggplot2::labs(title = paste0("Term–gene concept network: ", L1, " vs ", L2)) +
          ggplot2::theme_minimal(base_size = 12) +
          ggplot2::theme(legend.position = "right",
                         plot.margin = .unit_cm(cnet_margin_cm))

        # Force labels: terms (bigger) + genes (smaller, grey)
        if (has_ggrepel) {
          if (nrow(lay_terms)) {
            pnet <- pnet + ggrepel::geom_text_repel(
              data = lay_terms,
              ggplot2::aes(x = x, y = y, label = label),
              inherit.aes = FALSE,
              size = cnet_term_label_size,
              box.padding = 0.45, point.padding = 0.35,
              max.overlaps = Inf
            )
          }
          if (nrow(lay_genes)) {
            pnet <- pnet + ggrepel::geom_text_repel(
              data = lay_genes,
              ggplot2::aes(x = x, y = y, label = label),
              inherit.aes = FALSE,
              size = cnet_gene_label_size,
              colour = cnet_gene_label_color,
              box.padding = 0.2, point.padding = 0.1,
              max.overlaps = Inf
            )
          }
        } else {
          if (nrow(lay_terms))
            pnet <- pnet + ggraph::geom_node_text(data = lay_terms, ggplot2::aes(label = label),
                                                  size = cnet_term_label_size, colour = "black", check_overlap = TRUE)
          if (nrow(lay_genes))
            pnet <- pnet + ggraph::geom_node_text(data = lay_genes, ggplot2::aes(label = label),
                                                  size = cnet_gene_label_size, colour = cnet_gene_label_color, check_overlap = TRUE)
        }

        if (isTRUE(show_in_rstudio)) print(pnet)
        ggplot2::ggsave(cnet_file, plot = pnet, width = cnet_width, height = cnet_height,
                        dpi = if (identical(format, "pdf")) NA_integer_ else dpi, bg = "white")
        if (isTRUE(verbose)) message("  Saved concept network: ", normalizePath(cnet_file, winslash = "/", mustWork = FALSE))
      }
    } else if (isTRUE(cnet_show) && !has_ggraph) {
      warning("  Skipping cNET for ", pair_name, " (need 'ggraph').")
    }

    # ---- term-only overlap network ----
    termnet_file <- NA_character_
    if (isTRUE(ton_show) && has_ggraph) {
      fg_ord <- fg[order(fg$padj, fg$pval, -abs(fg$NES)), , drop = FALSE]
      fg_sel <- head(fg_ord, n = min(ton_top_terms, nrow(fg_ord)))
      if (nrow(fg_sel) > 0) {
        get_genes_for <- function(term) {
          if (identical(cnet_use, "leadingEdge")) {
            i <- match(term, fg_sel$pathway)
            unique(as.character(unlist(fg_sel$leadingEdge[i], use.names = FALSE)))
          } else {
            unique(as.character(pathways[[term]]))
          }
        }
        term_ids  <- as.character(fg_sel$pathway)
        term_sets <- lapply(term_ids, function(id) intersect(get_genes_for(id), names(score)))
        names(term_sets) <- term_ids

        e_from <- e_to <- character(0); e_w <- integer(0)
        if (length(term_ids) >= 2L) {
          for (pr in utils::combn(term_ids, 2, simplify = FALSE)) {
            gi <- term_sets[[pr[1]]]; gj <- term_sets[[pr[2]]]
            ov <- length(intersect(gi, gj))
            if (ov >= ton_min_overlap) { e_from <- c(e_from, pr[1]); e_to <- c(e_to, pr[2]); e_w <- c(e_w, ov) }
          }
        }

        k_vec <- vapply(term_ids, function(id) length(term_sets[[id]]), integer(1))
        node_df <- data.frame(
          id    = term_ids,
          label = term_ids,
          lp    = -log10(pmax(as.numeric(fg_sel$padj), .Machine$double.xmin)),
          k     = k_vec,
          size  = .rescale(k_vec, to = c(4,10)),
          stringsAsFactors = FALSE
        )

        g_terms <- if (length(e_from)) {
          igraph::graph_from_data_frame(
            d = data.frame(from = e_from, to = e_to, weight = e_w, stringsAsFactors = FALSE),
            directed = FALSE, vertices = node_df
          )
        } else igraph::graph_from_data_frame(d = NULL, directed = FALSE, vertices = node_df)

        lay <- .layout_terms(g_terms, ton_layout)
        nodes_df <- cbind(node_df, data.frame(x = lay[,1], y = lay[,2]))
        edf <- if (length(e_from)) data.frame(
          from = e_from, to = e_to, weight = e_w,
          x1 = nodes_df$x[match(e_from, nodes_df$id)],
          y1 = nodes_df$y[match(e_from, nodes_df$id)],
          x2 = nodes_df$x[match(e_to,   nodes_df$id)],
          y2 = nodes_df$y[match(e_to,   nodes_df$id)]
        ) else data.frame(from=character(0), to=character(0), weight=integer(0),
                          x1=numeric(0), y1=numeric(0), x2=numeric(0), y2=numeric(0))

        p_ton <- ggplot2::ggplot() +
          ggplot2::geom_segment(data = edf,
                                ggplot2::aes(x = x1, y = y1, xend = x2, yend = y2),
                                color = "grey70", alpha = 0.45) +
          ggplot2::geom_point(data = nodes_df,
                              ggplot2::aes(x = x, y = y, size = size, color = lp)) +
          ggplot2::scale_size_continuous(range = c(3, 10), guide = "none") +
          ggplot2::scale_color_viridis_c(name = expression(-log[10](adj.~p))) +
          ggplot2::geom_text(data = nodes_df, ggplot2::aes(x = x, y = y, label = label),
                             size = 4, color = "black", vjust = -0.8) +
          ggplot2::geom_text(data = nodes_df, ggplot2::aes(x = x, y = y, label = k),
                             size = 3.6, color = "grey20", vjust = 1.8) +
          ggplot2::labs(title = paste0("Term-only overlap network: ", L1, " vs ", L2), x = "x", y = "y") +
          ggplot2::theme_minimal(base_size = 12) +
          ggplot2::theme(plot.margin = .unit_cm(ton_margin_cm))

        termnet_file <- file.path(pair_dir, paste0("gsea_pair_termnet_", pair_name, ".", format))
        if (isTRUE(show_in_rstudio)) print(p_ton)
        ggplot2::ggsave(termnet_file, plot = p_ton, width = ton_width, height = ton_height,
                        dpi = if (identical(format, "pdf")) NA_integer_ else dpi, bg = "white")
        if (isTRUE(verbose)) message("  Saved term-only network: ", normalizePath(termnet_file, winslash = "/", mustWork = FALSE))
      }
    } else if (isTRUE(ton_show) && !has_ggraph) {
      warning("  Skipping term-only network for ", pair_name, " (need 'ggraph').")
    }

    by_pair[[pair_name]] <- list(
      csv_file       = csv_file,
      dotplot_file   = if (!is.na(dotplot_file)) dotplot_file else NA_character_,
      enrich_L1_file = if (!is.na(enr_up_file))  enr_up_file  else NA_character_,
      enrich_L2_file = if (!is.na(enr_dn_file))  enr_dn_file  else NA_character_,
      cnet_file      = if (!is.na(cnet_file))     cnet_file    else NA_character_,
      termnet_file   = if (!is.na(termnet_file))  termnet_file else NA_character_
    )

    summary_rows[[length(summary_rows)+1L]] <- data.frame(
      pair           = pair_name,
      n_ranked_genes = length(score),
      n_pathways     = nrow(fg),
      csv_file       = csv_file,
      dotplot_file   = if (!is.na(dotplot_file)) dotplot_file else NA_character_,
      enrich_L1_file = if (!is.na(enr_up_file))  enr_up_file  else NA_character_,
      enrich_L2_file = if (!is.na(enr_dn_file))  enr_dn_file  else NA_character_,
      cnet_file      = if (!is.na(cnet_file))     cnet_file    else NA_character_,
      termnet_file   = if (!is.na(termnet_file))  termnet_file else NA_character_,
      stringsAsFactors = FALSE
    )
  }

  # ---- manifest ----
  summary_df <- if (length(summary_rows)) do.call(rbind, summary_rows) else
    data.frame(pair=character(), n_ranked_genes=integer(), n_pathways=integer(),
               csv_file=character(), dotplot_file=character(),
               enrich_L1_file=character(), enrich_L2_file=character(),
               cnet_file=character(), termnet_file=character(),
               stringsAsFactors = FALSE)
  summary_file <- file.path(run_dir, "SUMMARY.csv")
  utils::write.csv(summary_df, summary_file, row.names = FALSE)
  message("Saved manifest: ", normalizePath(summary_file, winslash = "/", mustWork = FALSE))

  invisible(list(run_dir = run_dir, by_pair = by_pair, summary_file = summary_file))
}


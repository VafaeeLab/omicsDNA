


# ------------------------------------------------------------
# 26 — GO enrichment: single set OR per‑layer / per‑community multinet
# ------------------------------------------------------------

# Helper: convert gene identifiers to Entrez IDs

#' Convert gene identifiers to Entrez IDs (helper)
#'
#' @description
#' Wraps \code{clusterProfiler::bitr()} to convert input IDs to
#' \strong{ENTREZID}. Returns a conversion table and the unique
#' Entrez vector, with a coverage message.
#'
#' @param genes Character vector of gene identifiers (e.g. gene symbols).
#' @param fromType Input ID type (e.g., \code{"SYMBOL"}, \code{"ENSEMBL"},
#'   \code{"ALIAS"}, \code{"ENTREZID"}).
#' @param OrgDb Annotation package (e.g., \code{org.Hs.eg.db}).
#' @param verbose Logical; print coverage information. Default \code{TRUE}.
#'
#' @return A list with components:
#' \itemize{
#'   \item \code{table} — data.frame with columns \code{fromType} and \code{ENTREZID}.
#'   \item \code{entrez} — character vector of unique Entrez IDs.
#' }
#'
#' @export
convert_to_entrez <- function(genes,
                              fromType = "SYMBOL",
                              OrgDb,
                              verbose = TRUE) {
  if (missing(OrgDb))
    stop("Provide `OrgDb` (e.g., org.Hs.eg.db).")
  if (!requireNamespace("clusterProfiler", quietly = TRUE))
    stop("Install 'clusterProfiler' (BiocManager::install('clusterProfiler')).")

  genes <- unique(as.character(genes))
  if (!length(genes)) {
    return(list(table = data.frame(), entrez = character(0)))
  }

  conv <- try(
    clusterProfiler::bitr(
      genes,
      fromType = fromType,
      toType   = "ENTREZID",
      OrgDb    = OrgDb
    ),
    silent = TRUE
  )

  if (inherits(conv, "try-error") || is.null(conv) || !nrow(conv)) {
    if (isTRUE(verbose))
      message("No IDs mapped from ", fromType, " to ENTREZID.")
    return(list(table = data.frame(), entrez = character(0)))
  }

  conv   <- unique(conv[, c(fromType, "ENTREZID")])
  entrez <- unique(as.character(conv$ENTREZID))

  if (isTRUE(verbose)) {
    cov <- round(100 * length(entrez) / length(genes), 1)
    message("ID conversion: mapped ", length(entrez), "/", length(genes),
            " (", cov, "%).")
  }

  list(table = conv, entrez = entrez)
}





#' GO enrichment with compact plots — single set, per‑layer or per‑community
#'
#' @description
#' Runs GO enrichment (\pkg{clusterProfiler}) and produces barplots (optional),
#' dotplots and cnet plots.
#'
#' Modes:
#' \itemize{
#'   \item \strong{Single‑set} — supply \code{genes} (character vector).
#'   \item \strong{Per‑layer} — supply a \code{multinet::ml.network} in
#'         \code{net}; each layer’s vertex set is treated as a gene list.
#'   \item \strong{Per‑community} — supply \code{net} and a community table
#'         \code{communities} with columns \code{actor}, \code{layer}, and
#'         \code{cid} (from \code{detectCom()}); enrichment is performed for
#'         each community (\code{cid}) in each layer.
#' }
#'
#' Additionally, a \strong{restriction gene set} can be supplied via
#' \code{restrict_genes}. In multinet mode this restricts the query genes:
#' \itemize{
#'   \item \code{enrich_scope = "layer"}: query genes for layer \eqn{L} are
#'         \code{intersect( layer_genes(L), restrict_genes )}.
#'   \item \code{enrich_scope = "community"}: query genes for community
#'         (\eqn{L}, \eqn{cid}) are
#'         \code{intersect( community_genes(L,cid), restrict_genes )}.
#' }
#'
#' @details
#' \strong{ID handling.}
#' Input IDs are converted to \strong{ENTREZID} via
#' \code{\link{convert_to_entrez}} using \code{id_from}. GO enrichment is
#' computed with \code{clusterProfiler::enrichGO()}, and optional redundancy
#' reduction is performed with \code{clusterProfiler::simplify()}.
#'
#' \strong{Background.}
#' \itemize{
#'   \item If \code{universe} is supplied, it is passed directly to
#'         \code{enrichGO()} (assumed to be Entrez IDs).
#'   \item If \code{universe_mode = "union"} and \code{universe = NULL},
#'         the background is the union of mapped genes across layers (using
#'         \code{genes_by_layer} or vertex names per layer).
#'   \item If \code{universe_mode = "none"} and \code{universe = NULL},
#'         the default background of \pkg{clusterProfiler} is used.
#' }
#'
#' \strong{Multinet scopes.}
#' \itemize{
#'   \item \code{enrich_scope = "layer"} (default) — one enrichment per layer.
#'   \item \code{enrich_scope = "community"} — one enrichment per community
#'         (\code{cid}) per layer, using the \code{communities} table.
#' }
#'
#' For each enrichment, the CSV output contains all columns from
#' \code{as.data.frame(enrichResult)}, plus an additional column
#' \code{overlap_genes} which lists the overlapping genes (symbols) between
#' the GO term and the query set (layer, community, or restricted gene set).
#'
#' Plots are optionally shown in the current device (e.g. RStudio Plots) and
#' saved to PDF. In multinet mode you can request a single combined PDF
#' (\code{one_pdf = TRUE}) or one PDF per layer / per community.
#'
#' @param genes Character vector of genes (single‑set mode; ignored when
#'   \code{net} is not \code{NULL}). IDs must be of type \code{id_from}.
#' @param net Optional \code{multinet::ml.network} (multinet mode).
#' @param OrgDb Annotation package (e.g., \code{org.Hs.eg.db}).
#' @param communities Optional data.frame with at least columns \code{actor}
#'   and \code{layer}, plus one of \code{cid}, \code{com}, or \code{community}
#'   (community ID). Required when \code{enrich_scope = "community"}.
#' @param restrict_genes Optional character vector of genes (same ID type as
#'   \code{id_from}). In multinet mode, query genes are intersected with this
#'   set at the chosen scope (layer or community).
#' @param enrich_scope Character; when \code{net} is not \code{NULL}, controls
#'   whether enrichment is per-layer (\code{"layer"}, default) or per-community
#'   (\code{"community"}). Ignored in single-set mode.
#' @param id_from,ont See \code{clusterProfiler::bitr},
#'   \code{clusterProfiler::enrichGO} (defaults: \code{id_from = "SYMBOL"},
#'   \code{ont = "BP"}).
#' @param layer_order Optional character vector giving the order of layers
#'   (multinet mode).
#' @param genes_by_layer Optional named list \code{layer -> character vector}
#'   of genes. If provided, overrides the default “vertex names per layer”.
#'   Only used in multinet mode.
#' @param min_genes Minimum number of mapped genes to run enrichment.
#'   Default \code{3}.
#' @param universe Optional background (Entrez IDs). If \code{universe_mode =
#'   "union"}, the background is the union of all mapped genes across layers.
#' @param universe_mode One of \code{c("none","union")}. Default \code{"none"}.
#' @param pAdjustMethod P-value adjustment method passed to
#'   \code{enrichGO()}. One of \code{"BH"}, \code{"bonferroni"},
#'   \code{"BY"}, \code{"fdr"}, or \code{"none"}. If \code{"none"}, raw
#'   \code{pvalue} is used internally for redundancy reduction.
#' @param pvalueCutoff,qvalueCutoff Significance thresholds passed to
#'   \code{enrichGO()}.
#' @param simplify Logical; apply \code{clusterProfiler::simplify()}?
#'   Default \code{TRUE}.
#' @param simplify_cutoff,simplify_by,simplify_measure Simplification options.
#' @param showCategory Number of categories in dotplot; \code{cnet_showCategory}
#'   controls cnet.
#' @param cnet_showCategory Integer; often fewer than dotplot for legibility.
#' @param cnet_left_margin_pt Left margin (pt) for cnet plot to avoid
#'   clipping of labels.
#' @param include_barplot Logical; include a bar plot in single‑set mode.
#' @param include_barplot_multinet Logical; include a bar plot per layer or
#'   per community in multinet mode.
#' @param out_dir,file_prefix Output directory and file basename.
#' @param save_pdf,save_csv Write PDFs and CSVs? Defaults \code{TRUE}.
#' @param one_pdf Multinet mode: if \code{TRUE}, save one combined PDF; else
#'   one PDF per layer (scope = "layer") or per community (scope = "community").
#' @param show_in_rstudio Print plots to current device (RStudio Plots).
#'   Default \code{TRUE}.
#' @param width,height Figure size for PDFs (inches).
#' @param verbose Verbose messages. Default \code{TRUE}.
#'
#' @return
#' \itemize{
#'   \item \strong{Single‑set:}
#'     \code{list(ego, table, pdf_file, csv_file, mapping)}.
#'   \item \strong{Multinet, \code{enrich_scope = "layer"}:}
#'     \code{list(run_dir, by_layer, combined_pdf, layer_files)}, where
#'     \code{by_layer[[layer]]} contains \code{ego}, the result table and
#'     mapping.
#'   \item \strong{Multinet, \code{enrich_scope = "community"}:}
#'     \code{list(run_dir, by_layer, combined_pdf, module_files)}, where
#'     \code{by_layer[[layer]][["module<cid>"]]} holds the enrichment for that
#'     community.
#' }
#'
#' @examples
#' \dontrun{
#' ## Example objects:
#' ##   - net  : your multilayer network (multinet::ml.network)
#' ##   - comm : detectCom() output with columns actor / layer / cid (or com)
#'
#' ## 1) Single set (vector of symbols) ----------------------------------------
#' go_single <- go_enrichment_report(
#'   genes           = c("RHCG","SPRR1A","SPRR2A","SPRR3","KRT13","KRT4","CRNN"),
#'   OrgDb           = org.Hs.eg.db,
#'   ont             = "BP",
#'   file_prefix     = "demo_GO_BP_single",
#'   show_in_rstudio = TRUE
#' )
#'
#' ## 2) Per-layer enrichment (each layer’s vertex set) ------------------------
#' go_layer <- go_enrichment_report(
#'   net                      = net,
#'   OrgDb                    = org.Hs.eg.db,
#'   id_from                  = "SYMBOL",
#'   ont                      = "BP",
#'   enrich_scope             = "layer",
#'   layer_order              = multinet::layers_ml(net),
#'   universe_mode            = "union",
#'   pAdjustMethod            = "BH",
#'   pvalueCutoff             = 0.05,
#'   include_barplot_multinet = FALSE,
#'   cnet_showCategory        = 8,
#'   cnet_left_margin_pt      = 100,
#'   show_in_rstudio          = TRUE,
#'   save_pdf                 = TRUE,
#'   one_pdf                  = TRUE,    # one combined PDF
#'   file_prefix              = "GO_BP_byLayer",
#'   out_dir                  = getOption("mlnet.results_dir","omicsDNA_results")
#' )
#'
#' ## 3) Per-community enrichment (all genes per community) --------------------
#' go_comm <- go_enrichment_report(
#'   net                      = net,
#'   communities              = comm,    # actor / layer / cid (or com)
#'   OrgDb                    = org.Hs.eg.db,
#'   id_from                  = "SYMBOL",
#'   ont                      = "BP",
#'   enrich_scope             = "community",
#'   layer_order              = multinet::layers_ml(net),
#'   universe_mode            = "union",
#'   pAdjustMethod            = "BH",
#'   pvalueCutoff             = 0.05,
#'   include_barplot_multinet = FALSE,
#'   cnet_showCategory        = 8,
#'   cnet_left_margin_pt      = 100,
#'   show_in_rstudio          = TRUE,
#'   save_pdf                 = TRUE,
#'   one_pdf                  = FALSE,   # one PDF per community
#'   file_prefix              = "GO_BP_byCommunity",
#'   out_dir                  = getOption("mlnet.results_dir","omicsDNA_results")
#' )
#'
#' ## 4) Per-community enrichment restricted to a custom gene set -------------
#' ## e.g. oxidative stress gene list:
#' ## oxidative_stress_genes <- pc_genes[1:100]
#'
#' go_comm_ox <- go_enrichment_report(
#'   net                      = net,
#'   communities              = comm,
#'   restrict_genes           = oxidative_stress_genes,
#'   OrgDb                    = org.Hs.eg.db,
#'   id_from                  = "SYMBOL",
#'   ont                      = "BP",
#'   enrich_scope             = "community",
#'   layer_order              = multinet::layers_ml(net),
#'   universe_mode            = "union",
#'   pAdjustMethod            = "BH",    # or "fdr", "bonferroni", "none"
#'   pvalueCutoff             = 0.05,
#'   include_barplot_multinet = FALSE,
#'   cnet_showCategory        = 8,
#'   cnet_left_margin_pt      = 100,
#'   show_in_rstudio          = TRUE,
#'   save_pdf                 = TRUE,
#'   one_pdf                  = FALSE,
#'   file_prefix              = "customGenes_for_oxidative_stress",
#'   out_dir                  = getOption("mlnet.results_dir","omicsDNA_results")
#' )
#' }
#'
#' @importFrom utils write.csv
#' @importFrom grDevices pdf dev.off
#' @export
go_enrichment_report <- function(genes = NULL,
                                 net   = NULL,
                                 OrgDb,
                                 communities    = NULL,
                                 restrict_genes = NULL,
                                 enrich_scope   = c("layer","community"),
                                 id_from        = c("SYMBOL","ENSEMBL","ALIAS","ENTREZID"),
                                 ont            = c("BP","MF","CC"),
                                 layer_order    = NULL,
                                 genes_by_layer = NULL,
                                 min_genes      = 3,
                                 universe       = NULL,
                                 universe_mode  = c("none","union"),
                                 pAdjustMethod  = c("BH","bonferroni","BY","fdr","none"),
                                 pvalueCutoff   = 0.05,
                                 qvalueCutoff   = 0.2,
                                 simplify       = TRUE,
                                 simplify_cutoff = 0.7,
                                 simplify_by     = "p.adjust",
                                 simplify_measure = "Wang",
                                 showCategory      = 10,
                                 cnet_showCategory = 8,
                                 cnet_left_margin_pt = 80,
                                 include_barplot          = TRUE,
                                 include_barplot_multinet = FALSE,
                                 out_dir        = getOption("mlnet.results_dir","omicsDNA_results"),
                                 file_prefix    = NULL,
                                 save_pdf       = TRUE,
                                 save_csv       = TRUE,
                                 one_pdf        = TRUE,
                                 show_in_rstudio = TRUE,
                                 width           = 8.27,
                                 height          = 11.69,
                                 verbose         = TRUE) {

  enrich_scope   <- match.arg(enrich_scope)
  id_from        <- match.arg(id_from)
  ont            <- match.arg(ont)
  universe_mode  <- match.arg(universe_mode)
  pAdjustMethod  <- match.arg(pAdjustMethod)

  ## *** NEW: if no multiple-testing correction is requested, use raw pvalue
  ##          for redundancy reduction (simplify), instead of p.adjust.
  if (identical(pAdjustMethod, "none") && identical(simplify_by, "p.adjust")) {
    if (isTRUE(verbose))
      message("pAdjustMethod = 'none': using raw 'pvalue' for simplify().")
    simplify_by <- "pvalue"
  }

  if (missing(OrgDb))
    stop("Provide `OrgDb` (e.g., org.Hs.eg.db).")
  if (!requireNamespace("clusterProfiler", quietly = TRUE))
    stop("Install 'clusterProfiler'.")
  has_enrichplot <- requireNamespace("enrichplot", quietly = TRUE)
  has_ggplot2    <- requireNamespace("ggplot2", quietly = TRUE)

  if (!dir.exists(out_dir))
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  stamp <- format(Sys.time(), "%Y-%m-%d_%H%M%S")
  base  <- if (is.null(file_prefix) || !nzchar(file_prefix)) {
    paste0("GO_", ont, "_", stamp)
  } else {
    file_prefix
  }

  run_dir <- file.path(out_dir, paste0("go_enrich_", stamp))
  if (!dir.exists(run_dir))
    dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)

  # normalise restrict_genes, if any
  if (!is.null(restrict_genes)) {
    restrict_genes <- unique(as.character(restrict_genes))
    restrict_genes <- restrict_genes[nzchar(restrict_genes)]
  }

  # helper: add overlap_genes column (symbols) from geneID
  .add_overlap_genes <- function(tbl) {
    tbl <- as.data.frame(tbl)
    if ("geneID" %in% names(tbl)) {
      tbl$overlap_genes <- vapply(
        strsplit(as.character(tbl$geneID), "/"),
        function(x) paste(unique(x[nzchar(x)]), collapse = ";"),
        character(1)
      )
    } else if (!"overlap_genes" %in% names(tbl)) {
      tbl$overlap_genes <- NA_character_
    }
    tbl
  }

  # helper: save list of plots to a PDF (no on.exit, so no null-device issues)
  .save_plots_pdf <- function(filename, plots) {
    plots <- Filter(function(p) !inherits(p, "try-error") && !is.null(p), plots)
    if (!length(plots)) return(invisible(NULL))
    grDevices::pdf(filename, width = width, height = height,
                   onefile = TRUE, paper = "special")
    for (p in plots) print(p)
    grDevices::dev.off()
    invisible(filename)
  }

  # helper: build bar/dot/cnet plots for an enrichResult
  .mk_plots <- function(ego, title_prefix = "", do_bar = FALSE,
                        show_cat = showCategory,
                        cnet_cat = cnet_showCategory) {
    tbl <- as.data.frame(ego)
    mk  <- function(s) if (nzchar(title_prefix)) paste0("[", title_prefix, "] ", s) else s

    p_bar <- if (isTRUE(do_bar)) {
      try(
        clusterProfiler::barplot(
          ego,
          showCategory = show_cat,
          title = mk(paste0("Top GO ", ont, " (n = ", nrow(tbl), " terms)"))
        ),
        silent = TRUE
      )
    } else NULL

    p_dot <- if (isTRUE(has_enrichplot)) {
      try(
        enrichplot::dotplot(
          ego,
          showCategory = show_cat,
          title = mk(paste0("Dotplot: top GO ", ont, " (n = ", nrow(tbl), " terms)"))
        ),
        silent = TRUE
      )
    } else NULL

    p_cnet <- if (isTRUE(has_enrichplot)) {
      tmp <- try(
        enrichplot::cnetplot(
          ego,
          showCategory = cnet_cat,
          foldChange   = NULL,
          circular     = FALSE,
          title        = mk("Concept network (selected terms)")
        ),
        silent = TRUE
      )
      if (!inherits(tmp, "try-error") && isTRUE(has_ggplot2)) {
        tmp <- tmp +
          ggplot2::theme(
            plot.margin = ggplot2::margin(5.5, 5.5, 5.5, cnet_left_margin_pt, "pt")
          ) +
          ggplot2::coord_cartesian(clip = "off")
      }
      tmp
    } else NULL

    list(bar = p_bar, dot = p_dot, cnet = p_cnet)
  }

  # =======================================================================
  # A) SINGLE‑SET MODE (no net)
  # =======================================================================
  if (is.null(net)) {
    if (is.null(genes))
      stop("Provide `genes` (single‑set) or `net` (multinet).")

    q_genes <- unique(as.character(genes))
    q_genes <- q_genes[nzchar(q_genes)]
    if (!is.null(restrict_genes)) {
      q_genes <- intersect(q_genes, restrict_genes)
      if (isTRUE(verbose))
        message("Single-set mode: restricting to ", length(q_genes),
                " genes overlapping `restrict_genes`.")
    }

    map <- convert_to_entrez(
      genes    = q_genes,
      fromType = id_from,
      OrgDb    = OrgDb,
      verbose  = verbose
    )

    if (length(map$entrez) < min_genes) {
      message("Too few mapped genes (", length(map$entrez), " < ", min_genes,
              "); skipping enrichment.")
      return(invisible(list(
        ego      = NULL,
        table    = data.frame(),
        pdf_file = NULL,
        csv_file = NULL,
        mapping  = map$table
      )))
    }

    ego <- clusterProfiler::enrichGO(
      gene          = map$entrez,
      OrgDb         = OrgDb,
      keyType       = "ENTREZID",
      ont           = ont,
      universe      = universe,
      pAdjustMethod = pAdjustMethod,
      pvalueCutoff  = pvalueCutoff,
      qvalueCutoff  = qvalueCutoff,
      readable      = TRUE
    )

    if (is.null(ego) || nrow(as.data.frame(ego)) == 0L) {
      if (isTRUE(verbose))
        message("No enriched GO terms (ont = ", ont, ").")
      return(invisible(list(
        ego      = ego,
        table    = data.frame(),
        pdf_file = NULL,
        csv_file = NULL,
        mapping  = map$table
      )))
    }

    if (isTRUE(simplify)) {
      ego_s <- try(
        clusterProfiler::simplify(
          ego,
          cutoff  = simplify_cutoff,
          by      = simplify_by,
          measure = simplify_measure
        ),
        silent = TRUE
      )
      if (!inherits(ego_s, "try-error") && nrow(as.data.frame(ego_s)) > 0L)
        ego <- ego_s
    }

    tbl      <- .add_overlap_genes(as.data.frame(ego))
    csv_file <- if (isTRUE(save_csv)) file.path(run_dir, paste0(base, ".csv")) else NULL
    if (!is.null(csv_file))
      utils::write.csv(tbl, csv_file, row.names = FALSE)

    plots <- .mk_plots(ego, do_bar = include_barplot)
    if (isTRUE(show_in_rstudio)) {
      if (!inherits(plots$bar,  "try-error") && !is.null(plots$bar))  print(plots$bar)
      if (!inherits(plots$dot,  "try-error") && !is.null(plots$dot))  print(plots$dot)
      if (!inherits(plots$cnet, "try-error") && !is.null(plots$cnet)) print(plots$cnet)
    }

    pdf_file <- NULL
    if (isTRUE(save_pdf)) {
      pdf_file <- file.path(run_dir, paste0(base, ".pdf"))
      .save_plots_pdf(pdf_file, plots)
    }

    return(invisible(list(
      ego      = ego,
      table    = tbl,
      pdf_file = pdf_file,
      csv_file = csv_file,
      mapping  = map$table
    )))
  }

  # =======================================================================
  # B) MULTINET MODE
  # =======================================================================
  if (!requireNamespace("multinet", quietly = TRUE))
    stop("Install 'multinet' for net mode.")
  if (!requireNamespace("igraph", quietly = TRUE))
    stop("Install 'igraph' for net mode.")

  # communities required only for community scope
  if (identical(enrich_scope, "community")) {
    if (is.null(communities))
      stop("For enrich_scope = 'community', please supply `communities`.")
    if (!is.data.frame(communities))
      stop("`communities` must be a data.frame.")

    comm <- communities

    # Normalise `cid` column if needed (your requested snippet)
    if (!"cid" %in% names(comm)) {
      alt <- intersect(c("cid","com","community"), names(comm))[1L]
      if (is.na(alt))
        stop("Could not find cid/com/community column in `comm`.")
      comm$cid <- comm[[alt]]
    }

    if (!"actor" %in% names(comm))
      stop("`communities` must contain column `actor`.")
    if (!"layer" %in% names(comm))
      stop("`communities` must contain column `layer`.")

  } else {
    comm <- NULL
  }

  # layers
  Ls <- as.character(try(multinet::layers_ml(net), silent = TRUE))
  if (inherits(Ls, "try-error") || !length(Ls))
    stop("Could not read layers via multinet::layers_ml(net).")

  if (!is.null(layer_order)) {
    keep <- intersect(as.character(layer_order), Ls)
    if (!length(keep))
      stop("No valid layers after applying `layer_order`.")
    Ls <- keep
  }

  # layer gene sets
  if (is.null(genes_by_layer)) {
    glist <- try(as.list(net), silent = TRUE)
    if (inherits(glist, "try-error") || !is.list(glist))
      stop("Could not coerce `net` to list of igraphs via as.list().")
    if (!is.null(names(glist)) && any(names(glist) == "_flat_"))
      glist <- glist[names(glist) != "_flat_"]
    if (!all(Ls %in% names(glist)))
      stop("Missing igraph layers: ", paste(setdiff(Ls, names(glist)), collapse = ", "))

    genes_by_layer <- lapply(Ls, function(ln) {
      unique(as.character(igraph::V(glist[[ln]])$name))
    })
    names(genes_by_layer) <- Ls
  } else {
    if (is.null(names(genes_by_layer)))
      stop("`genes_by_layer` must be a named list (names = layers).")
    missing_layers <- setdiff(Ls, names(genes_by_layer))
    if (length(missing_layers))
      stop("`genes_by_layer` is missing layers: ", paste(missing_layers, collapse = ", "))
  }

  # background universe (Entrez)
  universe_vec <- universe
  if (identical(universe_mode, "union") && is.null(universe_vec)) {
    all_syms <- unique(unlist(genes_by_layer, use.names = FALSE))
    m_all    <- convert_to_entrez(all_syms, fromType = id_from, OrgDb = OrgDb, verbose = FALSE)
    universe_vec <- unique(m_all$entrez)
    if (isTRUE(verbose))
      message("Using union background with ", length(universe_vec), " mapped genes.")
  }

  # -------------------------------------------------------------------
  # SCOPE = "layer"
  # -------------------------------------------------------------------
  if (identical(enrich_scope, "layer")) {
    by_layer     <- vector("list", length(Ls)); names(by_layer) <- Ls
    to_pdf       <- list()
    combined_pdf <- NULL
    layer_files  <- NULL

    for (ln in Ls) {
      base_genes <- unique(as.character(genes_by_layer[[ln]]))
      base_genes <- base_genes[nzchar(base_genes)]

      q_genes <- base_genes
      if (!is.null(restrict_genes)) {
        q_genes <- intersect(base_genes, restrict_genes)
        if (isTRUE(verbose))
          message("Layer ", ln, ": restricting to ", length(q_genes),
                  " genes overlapping `restrict_genes`.")
      }

      map <- convert_to_entrez(q_genes, fromType = id_from, OrgDb = OrgDb, verbose = FALSE)

      if (length(map$entrez) < min_genes) {
        if (isTRUE(verbose))
          message("Layer ", ln, ": only ", length(map$entrez),
                  " mapped genes (< ", min_genes, "); skipped.")
        by_layer[[ln]] <- list(ego = NULL, table = data.frame(),
                               pdf_file = NULL, csv_file = NULL, mapping = map$table)
        next
      }

      ego <- clusterProfiler::enrichGO(
        gene          = map$entrez,
        OrgDb         = OrgDb,
        keyType       = "ENTREZID",
        ont           = ont,
        universe      = universe_vec,
        pAdjustMethod = pAdjustMethod,
        pvalueCutoff  = pvalueCutoff,
        qvalueCutoff  = qvalueCutoff,
        readable      = TRUE
      )

      if (is.null(ego) || nrow(as.data.frame(ego)) == 0L) {
        if (isTRUE(verbose))
          message("Layer ", ln, ": no enriched GO terms.")
        by_layer[[ln]] <- list(ego = ego, table = data.frame(),
                               pdf_file = NULL, csv_file = NULL, mapping = map$table)
        next
      }

      if (isTRUE(simplify)) {
        ego_s <- try(
          clusterProfiler::simplify(
            ego,
            cutoff  = simplify_cutoff,
            by      = simplify_by,
            measure = simplify_measure
          ),
          silent = TRUE
        )
        if (!inherits(ego_s, "try-error") && nrow(as.data.frame(ego_s)) > 0L)
          ego <- ego_s
      }

      tbl      <- .add_overlap_genes(as.data.frame(ego))
      csv_file <- if (isTRUE(save_csv)) file.path(run_dir, paste0(base, "_", ln, ".csv")) else NULL
      if (!is.null(csv_file))
        utils::write.csv(tbl, csv_file, row.names = FALSE)

      plots <- .mk_plots(ego, title_prefix = ln, do_bar = include_barplot_multinet)

      if (isTRUE(show_in_rstudio)) {
        if (!inherits(plots$bar,  "try-error") && !is.null(plots$bar))  print(plots$bar)
        if (!inherits(plots$dot,  "try-error") && !is.null(plots$dot))  print(plots$dot)
        if (!inherits(plots$cnet, "try-error") && !is.null(plots$cnet)) print(plots$cnet)
      }

      if (!inherits(plots$bar,  "try-error") && !is.null(plots$bar))  to_pdf[[length(to_pdf)+1L]] <- plots$bar
      if (!inherits(plots$dot,  "try-error") && !is.null(plots$dot))  to_pdf[[length(to_pdf)+1L]] <- plots$dot
      if (!inherits(plots$cnet, "try-error") && !is.null(plots$cnet)) to_pdf[[length(to_pdf)+1L]] <- plots$cnet

      by_layer[[ln]] <- list(ego = ego, table = tbl,
                             pdf_file = NULL, csv_file = csv_file, mapping = map$table)
    }

    if (isTRUE(save_pdf)) {
      if (isTRUE(one_pdf)) {
        combined_pdf <- file.path(run_dir, paste0(base, "_byLayer.pdf"))
        .save_plots_pdf(combined_pdf, to_pdf)
        if (isTRUE(verbose))
          message("Saved combined PDF: ", normalizePath(combined_pdf, FALSE))
      } else {
        layer_files <- setNames(character(0), character(0))
        for (ln in names(by_layer)) {
          info <- by_layer[[ln]]
          if (is.null(info$ego) || !nrow(info$table)) next
          plots <- .mk_plots(info$ego, title_prefix = ln, do_bar = include_barplot_multinet)
          fn <- file.path(run_dir, paste0(base, "_", ln, ".pdf"))
          .save_plots_pdf(fn, plots)
          layer_files[ln] <- fn
        }
      }
    }

    return(invisible(list(
      run_dir      = run_dir,
      by_layer     = by_layer,
      combined_pdf = if (isTRUE(save_pdf) && isTRUE(one_pdf)) combined_pdf else NULL,
      layer_files  = if (isTRUE(save_pdf) && !isTRUE(one_pdf)) layer_files else NULL
    )))
  }

  # -------------------------------------------------------------------
  # SCOPE = "community"
  # -------------------------------------------------------------------
  by_layer     <- vector("list", length(Ls)); names(by_layer) <- Ls
  to_pdf       <- list()
  combined_pdf <- NULL
  module_files <- NULL

  for (ln in Ls) {
    comm_l <- comm[comm$layer == ln, , drop = FALSE]
    if (!nrow(comm_l)) {
      if (isTRUE(verbose))
        message("Layer ", ln, ": no communities; skipped.")
      by_layer[[ln]] <- list()
      next
    }

    cids <- sort(unique(comm_l$cid))
    layer_modules <- list()

    for (cid in cids) {
      mod_name <- paste0("module", cid)
      base_genes <- unique(as.character(comm_l$actor[comm_l$cid == cid]))
      base_genes <- base_genes[nzchar(base_genes)]

      q_genes <- base_genes
      if (!is.null(restrict_genes)) {
        q_genes <- intersect(base_genes, restrict_genes)
        if (isTRUE(verbose))
          message("Layer ", ln, " / cid ", cid, ": restricting to ",
                  length(q_genes), " genes overlapping `restrict_genes`.")
      }

      map <- convert_to_entrez(q_genes, fromType = id_from, OrgDb = OrgDb, verbose = FALSE)

      if (length(map$entrez) < min_genes) {
        if (isTRUE(verbose))
          message("Layer ", ln, " / cid ", cid, ": only ", length(map$entrez),
                  " mapped genes (< ", min_genes, "); skipped.")
        layer_modules[[mod_name]] <- list(
          ego = NULL, table = data.frame(),
          pdf_file = NULL, csv_file = NULL, mapping = map$table
        )
        next
      }

      ego <- clusterProfiler::enrichGO(
        gene          = map$entrez,
        OrgDb         = OrgDb,
        keyType       = "ENTREZID",
        ont           = ont,
        universe      = universe_vec,
        pAdjustMethod = pAdjustMethod,
        pvalueCutoff  = pvalueCutoff,
        qvalueCutoff  = qvalueCutoff,
        readable      = TRUE
      )

      if (is.null(ego) || nrow(as.data.frame(ego)) == 0L) {
        if (isTRUE(verbose))
          message("Layer ", ln, " / cid ", cid, ": no enriched GO terms.")
        layer_modules[[mod_name]] <- list(
          ego = ego, table = data.frame(),
          pdf_file = NULL, csv_file = NULL, mapping = map$table
        )
        next
      }

      if (isTRUE(simplify)) {
        ego_s <- try(
          clusterProfiler::simplify(
            ego,
            cutoff  = simplify_cutoff,
            by      = simplify_by,
            measure = simplify_measure
          ),
          silent = TRUE
        )
        if (!inherits(ego_s, "try-error") && nrow(as.data.frame(ego_s)) > 0L)
          ego <- ego_s
      }

      tbl      <- .add_overlap_genes(as.data.frame(ego))
      csv_file <- if (isTRUE(save_csv)) {
        file.path(run_dir, paste0(base, "_", ln, "_module", cid, ".csv"))
      } else NULL
      if (!is.null(csv_file))
        utils::write.csv(tbl, csv_file, row.names = FALSE)

      title_prefix <- paste0(ln, " / cid ", cid)
      plots <- .mk_plots(ego, title_prefix = title_prefix,
                         do_bar = include_barplot_multinet)

      if (isTRUE(show_in_rstudio)) {
        if (!inherits(plots$bar,  "try-error") && !is.null(plots$bar))  print(plots$bar)
        if (!inherits(plots$dot,  "try-error") && !is.null(plots$dot))  print(plots$dot)
        if (!inherits(plots$cnet, "try-error") && !is.null(plots$cnet)) print(plots$cnet)
      }

      if (!inherits(plots$bar,  "try-error") && !is.null(plots$bar))  to_pdf[[length(to_pdf)+1L]] <- plots$bar
      if (!inherits(plots$dot,  "try-error") && !is.null(plots$dot))  to_pdf[[length(to_pdf)+1L]] <- plots$dot
      if (!inherits(plots$cnet, "try-error") && !is.null(plots$cnet)) to_pdf[[length(to_pdf)+1L]] <- plots$cnet

      layer_modules[[mod_name]] <- list(
        ego      = ego,
        table    = tbl,
        pdf_file = NULL,
        csv_file = csv_file,
        mapping  = map$table
      )

      # optional: individual PDF per community
      if (isTRUE(save_pdf) && !isTRUE(one_pdf)) {
        if (is.null(module_files)) module_files <- character(0)
        stub <- paste0(base, "_", ln, "_module", cid)
        fn   <- file.path(run_dir, paste0(stub, ".pdf"))
        .save_plots_pdf(fn, plots)
        module_files[paste0(ln, "_module", cid)] <- fn
      }
    }

    by_layer[[ln]] <- layer_modules
  }

  if (isTRUE(save_pdf) && isTRUE(one_pdf)) {
    combined_pdf <- file.path(run_dir, paste0(base, "_byCommunity.pdf"))
    .save_plots_pdf(combined_pdf, to_pdf)
    if (isTRUE(verbose))
      message("Saved combined community PDF: ",
              normalizePath(combined_pdf, FALSE))
  }

  invisible(list(
    run_dir      = run_dir,
    by_layer     = by_layer,
    combined_pdf = if (isTRUE(save_pdf) && isTRUE(one_pdf)) combined_pdf else NULL,
    module_files = if (isTRUE(save_pdf) && !isTRUE(one_pdf)) module_files else NULL
  ))
}



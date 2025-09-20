# ------------------------------------------------------------
# 26 — GO enrichment: single set OR per‑layer multinet
# ------------------------------------------------------------

#' Convert gene identifiers to Entrez IDs (helper)
#'
#' @description
#' Wraps \code{clusterProfiler::bitr()} to convert input IDs to \strong{ENTREZID}.
#' Returns a conversion table and the unique Entrez vector, with a brief coverage message.
#'
#' @param genes Character vector of gene identifiers.
#' @param fromType Input ID type (e.g., "SYMBOL", "ENSEMBL", "ALIAS", "ENTREZID").
#' @param OrgDb   Annotation package (e.g., \code{org.Hs.eg.db}).
#' @param verbose Logical; print coverage. Default TRUE.
#' @return list(table = data.frame, entrez = character())
#' @export
convert_to_entrez <- function(genes,
                              fromType = "SYMBOL",
                              OrgDb,
                              verbose = TRUE) {
  if (missing(OrgDb)) stop("Provide `OrgDb` (e.g., org.Hs.eg.db).")
  if (!requireNamespace("clusterProfiler", quietly = TRUE))
    stop("Install 'clusterProfiler' (BiocManager::install('clusterProfiler')).")
  genes <- unique(as.character(genes))
  if (!length(genes)) return(list(table = data.frame(), entrez = character(0)))

  conv <- try(clusterProfiler::bitr(genes, fromType = fromType, toType = "ENTREZID", OrgDb = OrgDb), silent = TRUE)
  if (inherits(conv, "try-error") || is.null(conv) || !nrow(conv)) {
    if (isTRUE(verbose)) message("No IDs mapped from ", fromType, " to ENTREZID.")
    return(list(table = data.frame(), entrez = character(0)))
  }
  conv <- unique(conv[, c(fromType, "ENTREZID")])
  entrez <- unique(as.character(conv$ENTREZID))
  if (isTRUE(verbose)) {
    cov <- round(100 * length(entrez) / length(genes), 1)
    message("ID conversion: mapped ", length(entrez), "/", length(genes), " (", cov, "%).")
  }
  list(table = conv, entrez = entrez)
}


#' GO enrichment with compact plots — single set or per‑layer multinet
#'
#' @description
#' Runs GO enrichment (\pkg{clusterProfiler}) and produces a dotplot and a cnet plot.
#' In **single‑set mode**, supply a vector of genes. In **multinet mode**, supply a
#' \code{multinet::ml.network}; the function uses each layer's vertex set (by default)
#' as the genes for that layer. Plots are shown in RStudio (once) and then saved to PDF
#' (one combined file or one file per layer) — **no duplicate pages**.
#'
#' @details
#' - IDs are harmonised to ENTREZ; multiple testing via \code{pAdjustMethod} (default "BH").
#' - Cnet labels can be protected from clipping via \code{cnet_left_margin_pt}.
#' - Layers with too few mapped genes or no significant terms are skipped (with messages).
#'
#' @param genes Character vector of genes (single‑set mode).
#' @param net Optional \code{multinet::ml.network} (per‑layer mode).
#' @param OrgDb Annotation package (e.g., \code{org.Hs.eg.db}).
#' @param id_from,ont See \code{clusterProfiler::bitr}, \code{enrichGO} (default id_from="SYMBOL", ont="BP").
#' @param layer_order Optional vector giving the order of layers (multinet mode).
#' @param genes_by_layer Optional **named list** layer -> character vector of genes.
#'   If provided, overrides the default (all vertex names per layer).
#' @param min_genes Minimum number of mapped genes to run enrichment. Default 3.
#' @param universe Optional background (Entrez IDs). If \code{universe_mode="union"},
#'   the background is the union of all mapped genes across layers.
#' @param universe_mode One of c("none","union"). Default "none".
#' @param pAdjustMethod,pvalueCutoff,qvalueCutoff Enrichment thresholds.
#' @param simplify,simplify_cutoff,simplify_by,simplify_measure Redundancy reduction options.
#' @param showCategory Number of categories in dotplot; \code{cnet_showCategory} for cnet.
#' @param cnet_showCategory Integer, often fewer than dotplot for legibility (default 8).
#' @param cnet_left_margin_pt Left margin (pt) added to cnet plot to avoid label clipping (default 80).
#' @param include_barplot Logical; include a bar plot in single‑set mode. Default TRUE.
#' @param include_barplot_multinet Logical; include a bar plot per layer in multinet mode. Default FALSE.
#' @param out_dir,file_prefix Output directory and file basename.
#' @param save_pdf,save_csv Write PDFs and CSVs? Defaults TRUE.
#' @param one_pdf If TRUE (multinet), use one combined PDF; else one file per layer. Default TRUE.
#' @param show_in_rstudio Print plots to current device (RStudio Plots). Default TRUE.
#' @param width,height Figure size for PDFs (inches).
#' @param verbose Verbose messages. Default TRUE.
#'
#' @return
#' - Single‑set: list(ego, table, pdf_file, csv_file, mapping).
#' - Multinet: list(by_layer = <named list>, combined_pdf, layer_files).
#'
#' @examples
#' \dontrun{
#' ## Single set
#' g <- c("RHCG","SPRR1A","SPRR2A","SPRR3","KRT13","KRT4","CRNN")
#' res <- go_enrichment_report(genes = g, OrgDb = org.Hs.eg.db, ont = "BP",
#'                             file_prefix = "g5l7_1_GO_BP", show_in_rstudio = TRUE)
#'
#' ## Multinet (per layer), using each layer's vertex set as the gene list
#' out <- go_enrichment_report(net = net, OrgDb = org.Hs.eg.db, ont = "BP",
#'                             layer_order = multinet::layers_ml(net),
#'                             universe_mode = "union",
#'                             cnet_left_margin_pt = 100,
#'                             include_barplot_multinet = FALSE,
#'                             one_pdf = TRUE, file_prefix = "GO_BP_byLayer")
#' out$combined_pdf
#' }
#' @importFrom utils write.csv
#' @importFrom grDevices pdf dev.off
#' @export
go_enrichment_report <- function(genes = NULL,
                                 net   = NULL,
                                 OrgDb,
                                 id_from        = c("SYMBOL","ENSEMBL","ALIAS","ENTREZID"),
                                 ont            = c("BP","MF","CC"),
                                 layer_order    = NULL,
                                 genes_by_layer = NULL,
                                 min_genes      = 3,
                                 universe       = NULL,
                                 universe_mode  = c("none","union"),
                                 pAdjustMethod  = c("BH","bonferroni","BY","fdr"),
                                 pvalueCutoff   = 0.05,
                                 qvalueCutoff   = 0.2,
                                 simplify       = TRUE,
                                 simplify_cutoff= 0.7,
                                 simplify_by    = "p.adjust",
                                 simplify_measure = "Wang",
                                 showCategory   = 10,
                                 cnet_showCategory = 8,
                                 cnet_left_margin_pt = 80,
                                 include_barplot          = TRUE,
                                 include_barplot_multinet = FALSE,
                                 out_dir        = getOption("mlnet.results_dir","omicsDNA_results"),
                                 file_prefix    = NULL,
                                 save_pdf       = TRUE,
                                 save_csv       = TRUE,
                                 one_pdf        = TRUE,
                                 show_in_rstudio= TRUE,
                                 width          = 8.27,
                                 height         = 11.69,
                                 verbose        = TRUE) {

  id_from       <- match.arg(id_from)
  ont           <- match.arg(ont)
  universe_mode <- match.arg(universe_mode)
  pAdjustMethod <- match.arg(pAdjustMethod)

  if (missing(OrgDb)) stop("Provide `OrgDb` (e.g., org.Hs.eg.db).")
  if (!requireNamespace("clusterProfiler", quietly = TRUE))
    stop("Install 'clusterProfiler'.")
  has_enrichplot <- requireNamespace("enrichplot", quietly = TRUE)
  has_ggplot2    <- requireNamespace("ggplot2", quietly = TRUE)

  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  stamp <- format(Sys.time(), "%Y-%m-%d_%H%M%S")
  base  <- if (is.null(file_prefix) || !nzchar(file_prefix)) paste0("GO_", ont, "_", stamp) else file_prefix

  # -------- helper to build bar/dot/cnet (returns list(bar, dot, cnet)) --------
  .mk_plots <- function(ego, title_prefix = "", do_bar = FALSE,
                        show_cat = showCategory, cnet_cat = cnet_showCategory) {
    tbl <- as.data.frame(ego)
    mk  <- function(s) if (nzchar(title_prefix)) paste0("[", title_prefix, "] ", s) else s

    p_bar <- if (isTRUE(do_bar)) {
      try(clusterProfiler::barplot(ego, showCategory = show_cat,
                                   title = mk(paste0("Top GO ", ont, " (n=", nrow(tbl), " terms)"))),
          silent = TRUE)
    } else NULL

    p_dot <- if (isTRUE(has_enrichplot)) {
      try(enrichplot::dotplot(ego, showCategory = show_cat,
                              title = mk(paste0("Dotplot: top GO ", ont, " (n=", nrow(tbl), " terms)"))),
          silent = TRUE)
    } else NULL

    p_cnet <- if (isTRUE(has_enrichplot)) {
      tmp <- try(enrichplot::cnetplot(ego, showCategory = cnet_cat, foldChange = NULL,
                                      circular = FALSE, title = mk("Concept network (selected terms)")),
                 silent = TRUE)
      if (!inherits(tmp, "try-error") && isTRUE(has_ggplot2)) {
        tmp <- tmp + ggplot2::theme(
          plot.margin = ggplot2::margin(5.5, 5.5, 5.5, cnet_left_margin_pt, "pt")
        ) + ggplot2::coord_cartesian(clip = "off")
      }
      tmp
    } else NULL

    list(bar = p_bar, dot = p_dot, cnet = p_cnet)
  }

  # ==============================
  # A) SINGLE‑SET MODE
  # ==============================
  if (is.null(net)) {
    if (is.null(genes)) stop("Provide `genes` (single‑set) or `net` (multinet).")

    map <- convert_to_entrez(genes, fromType = id_from, OrgDb = OrgDb, verbose = verbose)
    if (length(map$entrez) < min_genes) {
      message("Too few mapped genes (", length(map$entrez), " < ", min_genes, "); skipping enrichment.")
      return(invisible(list(ego = NULL, table = data.frame(), pdf_file = NULL, csv_file = NULL, mapping = map$table)))
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
    if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
      if (isTRUE(verbose)) message("No enriched GO terms (ont = ", ont, ").")
      return(invisible(list(ego = ego, table = data.frame(), pdf_file = NULL, csv_file = NULL, mapping = map$table)))
    }
    if (isTRUE(simplify)) {
      ego_s <- try(clusterProfiler::simplify(ego, cutoff = simplify_cutoff, by = simplify_by, measure = simplify_measure),
                   silent = TRUE)
      if (!inherits(ego_s, "try-error") && nrow(as.data.frame(ego_s)) > 0) ego <- ego_s
    }
    tbl <- as.data.frame(ego)
    csv_file <- if (isTRUE(save_csv)) file.path(out_dir, paste0(base, ".csv")) else NULL
    if (!is.null(csv_file)) utils::write.csv(tbl, csv_file, row.names = FALSE)

    plots <- .mk_plots(ego, do_bar = include_barplot)
    if (isTRUE(show_in_rstudio)) {
      if (!inherits(plots$bar,  "try-error") && !is.null(plots$bar))  print(plots$bar)
      if (!inherits(plots$dot,  "try-error") && !is.null(plots$dot))  print(plots$dot)
      if (!inherits(plots$cnet, "try-error") && !is.null(plots$cnet)) print(plots$cnet)
    }
    pdf_file <- if (isTRUE(save_pdf)) file.path(out_dir, paste0(base, ".pdf")) else NULL
    if (!is.null(pdf_file)) {
      grDevices::pdf(pdf_file, width = width, height = height, onefile = TRUE, paper = "special")
      on.exit(grDevices::dev.off(), add = TRUE)
      if (!inherits(plots$bar,  "try-error") && !is.null(plots$bar))  print(plots$bar)
      if (!inherits(plots$dot,  "try-error") && !is.null(plots$dot))  print(plots$dot)
      if (!inherits(plots$cnet, "try-error") && !is.null(plots$cnet)) print(plots$cnet)
    }
    return(invisible(list(ego = ego, table = tbl, pdf_file = pdf_file, csv_file = csv_file, mapping = map$table)))
  }

  # ==============================
  # B) MULTINET MODE (per layer)
  # ==============================
  if (!requireNamespace("multinet", quietly = TRUE)) stop("Install 'multinet' for net mode.")
  if (!requireNamespace("igraph",   quietly = TRUE)) stop("Install 'igraph' for net mode.")

  Ls <- as.character(try(multinet::layers_ml(net), silent = TRUE))
  if (inherits(Ls, "try-error") || !length(Ls)) stop("Could not read layers via multinet::layers_ml(net).")
  if (!is.null(layer_order)) {
    keep <- intersect(as.character(layer_order), Ls)
    if (!length(keep)) stop("No valid layers after applying `layer_order`.")
    Ls <- keep
  }

  # Choose gene sets per layer
  if (is.null(genes_by_layer)) {
    glist <- try(as.list(net), silent = TRUE)
    if (inherits(glist, "try-error") || !is.list(glist)) stop("Could not coerce `net` to list of igraphs via as.list().")
    if (!is.null(names(glist)) && any(names(glist) == "_flat_")) glist <- glist[names(glist) != "_flat_"]
    if (!all(Ls %in% names(glist))) stop("Missing igraph layers: ", paste(setdiff(Ls, names(glist)), collapse = ", "))
    genes_by_layer <- lapply(Ls, function(ln) unique(as.character(igraph::V(glist[[ln]])$name)))
    names(genes_by_layer) <- Ls
  } else {
    if (is.null(names(genes_by_layer))) stop("`genes_by_layer` must be a named list (names = layers).")
    missing_layers <- setdiff(Ls, names(genes_by_layer))
    if (length(missing_layers)) stop("`genes_by_layer` is missing layers: ", paste(missing_layers, collapse = ", "))
  }

  # Universe handling
  universe_vec <- universe
  if (identical(universe_mode, "union")) {
    all_syms <- unique(unlist(genes_by_layer, use.names = FALSE))
    m_all    <- convert_to_entrez(all_syms, fromType = id_from, OrgDb = OrgDb, verbose = FALSE)
    universe_vec <- unique(m_all$entrez)
    if (isTRUE(verbose)) message("Using union background with ", length(universe_vec), " mapped genes.")
  }

  by_layer <- vector("list", length(Ls)); names(by_layer) <- Ls
  to_pdf   <- list()

  for (ln in Ls) {
    gset <- genes_by_layer[[ln]]
    map  <- convert_to_entrez(gset, fromType = id_from, OrgDb = OrgDb, verbose = FALSE)
    if (length(map$entrez) < min_genes) {
      if (isTRUE(verbose)) message("Layer ", ln, ": only ", length(map$entrez), " mapped genes (< ", min_genes, "); skipped.")
      by_layer[[ln]] <- list(ego = NULL, table = data.frame(), pdf_file = NULL, csv_file = NULL, mapping = map$table)
      next
    }

    ego <- clusterProfiler::enrichGO(
      gene          = map$entrez, OrgDb = OrgDb, keyType = "ENTREZID", ont = ont,
      universe      = universe_vec, pAdjustMethod = pAdjustMethod,
      pvalueCutoff  = pvalueCutoff, qvalueCutoff = qvalueCutoff, readable = TRUE
    )
    if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
      if (isTRUE(verbose)) message("Layer ", ln, ": no enriched GO terms.")
      by_layer[[ln]] <- list(ego = ego, table = data.frame(), pdf_file = NULL, csv_file = NULL, mapping = map$table)
      next
    }
    if (isTRUE(simplify)) {
      ego_s <- try(clusterProfiler::simplify(ego, cutoff = simplify_cutoff, by = simplify_by, measure = simplify_measure),
                   silent = TRUE)
      if (!inherits(ego_s, "try-error") && nrow(as.data.frame(ego_s)) > 0) ego <- ego_s
    }

    tbl <- as.data.frame(ego)
    csv_file <- if (isTRUE(save_csv)) file.path(out_dir, paste0(base, "_", ln, ".csv")) else NULL
    if (!is.null(csv_file)) utils::write.csv(tbl, csv_file, row.names = FALSE)

    # Build plots for this layer (bar optional in multinet)
    plots <- .mk_plots(ego, title_prefix = ln, do_bar = include_barplot_multinet)

    # Show once on RStudio device (no PDF open yet)
    if (isTRUE(show_in_rstudio)) {
      if (!inherits(plots$bar,  "try-error") && !is.null(plots$bar))  print(plots$bar)
      if (!inherits(plots$dot,  "try-error") && !is.null(plots$dot))  print(plots$dot)
      if (!inherits(plots$cnet, "try-error") && !is.null(plots$cnet)) print(plots$cnet)
    }

    # Queue for saving later (prevents duplicates)
    if (!inherits(plots$bar,  "try-error") && !is.null(plots$bar))  to_pdf[[length(to_pdf)+1L]] <- plots$bar
    if (!inherits(plots$dot,  "try-error") && !is.null(plots$dot))  to_pdf[[length(to_pdf)+1L]] <- plots$dot
    if (!inherits(plots$cnet, "try-error") && !is.null(plots$cnet)) to_pdf[[length(to_pdf)+1L]] <- plots$cnet

    by_layer[[ln]] <- list(ego = ego, table = tbl, pdf_file = NULL, csv_file = csv_file, mapping = map$table)
  }

  # Save: either one combined PDF, or one PDF per layer (no duplicate pages)
  combined_pdf <- NULL; layer_files <- NULL
  if (isTRUE(save_pdf)) {
    if (isTRUE(one_pdf)) {
      combined_pdf <- file.path(out_dir, paste0(base, "_byLayer.pdf"))
      grDevices::pdf(combined_pdf, width = width, height = height, onefile = TRUE, paper = "special")
      on.exit(grDevices::dev.off(), add = TRUE)
      for (p in to_pdf) print(p)
      grDevices::dev.off()
      if (isTRUE(verbose)) message("Saved combined PDF: ", normalizePath(combined_pdf, FALSE))
    } else {
      # one file per layer (re‑render plots layer by layer to avoid mixing)
      layer_files <- setNames(character(0), character(0))
      for (ln in names(by_layer)) {
        info <- by_layer[[ln]]
        if (is.null(info$ego) || !nrow(info$table)) next
        plots <- .mk_plots(info$ego, title_prefix = ln, do_bar = include_barplot_multinet)
        fn <- file.path(out_dir, paste0(base, "_", ln, ".pdf"))
        grDevices::pdf(fn, width = width, height = height, onefile = TRUE, paper = "special")
        on.exit(grDevices::dev.off(), add = TRUE)
        if (!inherits(plots$bar,  "try-error") && !is.null(plots$bar))  print(plots$bar)
        if (!inherits(plots$dot,  "try-error") && !is.null(plots$dot))  print(plots$dot)
        if (!inherits(plots$cnet, "try-error") && !is.null(plots$cnet)) print(plots$cnet)
        grDevices::dev.off()
        layer_files[ln] <- fn
      }
    }
  }

  invisible(list(by_layer = by_layer, combined_pdf = combined_pdf, layer_files = layer_files))
}

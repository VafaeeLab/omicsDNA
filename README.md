omicsDNA
================

- [omicsDNA
  <img src="man/figures/logo.png" align="right" height="110"/>](#omicsdna-)
  - [Installation](#installation)
  - [Function‑by‑function examples](#functionbyfunction-examples)
    - [1) `buildAdjacency()` — group‑wise correlation → adjacency (with
      optional
      resampling)](#1-buildadjacency--groupwise-correlation--adjacency-with-optional-resampling)
    - [1b) `sc_buildAdjacency()` — **single‑cell** Seurat →
      per‑cell‑type
      adjacency](#1b-sc_buildadjacency--singlecell-seurat--percelltype-adjacency)
    - [2) `edgesFromAdjacency()` — matrices (or lists of matrices) →
      edge
      tables](#2-edgesfromadjacency--matrices-or-lists-of-matrices--edge-tables)
    - [3) `consensusEdges()` — across‑repeat edge consensus per layer
      (CSV/XLSX/RDS)](#3-consensusedges--acrossrepeat-edge-consensus-per-layer-csvxlsxrds)
    - [4) `build_multiNet()` — assemble multilayer from per‑layer **edge
      lists**](#4-build_multinet--assemble-multilayer-from-perlayer-edge-lists)
    - [5) `build_multinet_from_graphs()` — assemble multilayer from
      **igraph**
      graphs](#5-build_multinet_from_graphs--assemble-multilayer-from-igraph-graphs)
    - [6) `edges_list_to_graphs()` and 7) `graphs_to_edges_list()` —
      convert back and
      forth](#6-edges_list_to_graphs-and-7-graphs_to_edges_list--convert-back-and-forth)
    - [8) `compare_multinets()` — compare two or more multilayer
      networks](#8-compare_multinets--compare-two-or-more-multilayer-networks)
    - [9) `add_network_attributes()` — attach actor/edge metadata to the
      multilayer
      net](#9-add_network_attributes--attach-actoredge-metadata-to-the-multilayer-net)
    - [10) `network_attributes()` — report what attributes are
      attached](#10-network_attributes--report-what-attributes-are-attached)
    - [11) `detectCom()` — community detection on the
      supra‑graph](#11-detectcom--community-detection-on-the-supragraph)
    - [12) `plotCom()` — community layout + CSV export of node
      positions](#12-plotcom--community-layout--csv-export-of-node-positions)
    - [13) `com_vs_samples()` — sanity: community counts vs sample
      counts](#13-com_vs_samples--sanity-community-counts-vs-sample-counts)
    - [14) `analyze_actor_overlap()` — actor Jaccard overlaps across
      layers](#14-analyze_actor_overlap--actor-jaccard-overlaps-across-layers)
    - [15) `analyze_edge_overlap()` — edge Jaccard overlaps across
      layers](#15-analyze_edge_overlap--edge-jaccard-overlaps-across-layers)
    - [16) `layer_metrics()` — per‑layer size/centrality/path
      summaries](#16-layer_metrics--perlayer-sizecentralitypath-summaries)
    - [17) `annotateCom()` — add feature attributes (e.g., GeneType) to
      community
      rows](#17-annotatecom--add-feature-attributes-eg-genetype-to-community-rows)
    - [18) `sumComFeat()` — summaries by feature type
      (actor/community/layer)](#18-sumcomfeat--summaries-by-feature-type-actorcommunitylayer)
    - [19) `get_FeatureDeg()` — degrees (per layer /
      long)](#19-get_featuredeg--degrees-per-layer--long)
    - [20) `gp_enrich_multinet()` — g:Profiler enrichment **per
      layer**](#20-gp_enrich_multinet--gprofiler-enrichment-per-layer)
    - [21) `go_enrichment_report()` — GO (clusterProfiler) for a vector
      **or**
      per‑layer](#21-go_enrichment_report--go-clusterprofiler-for-a-vector-or-perlayer)
    - [22) `gsea_msigdbr_layer_pairs()` — **GSEA** (MSigDB, fgsea)
      across **consecutive layer
      pairs**](#22-gsea_msigdbr_layer_pairs--gsea-msigdb-fgsea-across-consecutive-layer-pairs)
    - [23) `animate_multiNet()` — interactive HTML
      animation](#23-animate_multinet--interactive-html-animation)
    - [24) `filmstrip_multiNet()` — static grid (one panel per
      layer)](#24-filmstrip_multinet--static-grid-one-panel-per-layer)
    - [25) `plotActivityTimeline()` — Gantt‑style lifespans of
      edges/vertices](#25-plotactivitytimeline--ganttstyle-lifespans-of-edgesvertices)
    - [26) `plot_layer_diff()` — pairwise layer “new / lost / common”
      edges](#26-plot_layer_diff--pairwise-layer-new--lost--common-edges)
    - [27) `grid_layer_diffs()` — multi‑panel grid of **all consecutive
      pair**
      differences](#27-grid_layer_diffs--multipanel-grid-of-all-consecutive-pair-differences)
  - [Reproducibility notes](#reproducibility-notes)
  - [🙋 Feedback & contributions](#raising_hand-feedback--contributions)
  - [📄 License](#page_facing_up-license)

# omicsDNA <img src="man/figures/logo.png" align="right" height="110"/>

*A multilayer network toolkit for grouped omics (and other
high‑dimensional) data.*

`omicsDNA` builds **layered networks** from grouped data (e.g., age
groups, disease stages), stabilises edges via resampling + consensus,
assembles a **multilayer** object, detects and visualises
**communities**, quantifies **layer‑wise structure**, compares
**overlaps**, and produces **static** and **dynamic** visuals with
reproducible, timestamped outputs. (Examples below follow the function
order in the package script.)

> Tip: set a default results directory so all artifacts land in one
> place.

``` r
options(mlnet.results_dir = "~/omicsDNA_results")
```

------------------------------------------------------------------------

## Installation

``` r
# Development version
install.packages("devtools")
devtools::install_github("VafaeeLab/omicsDNA")

# Optional for dynamics/exports
install.packages(c("ndtv","networkDynamic"))        # HTML animations / timelines
install.packages(c("gifski","gganimate","png","av"))# GIF / MP4 exporters
```

------------------------------------------------------------------------

## Function‑by‑function examples

Below, each function has **one example with its most informative
arguments** (as in the package script). Replace object names
(`expression_mat`, `metaData`, `genes_info`, …) with your data.

### 1) `buildAdjacency()` — group‑wise correlation → adjacency (with optional resampling)

``` r
set.seed(1)
adjacency <- buildAdjacency(
  dataMatrix      = expression_mat,         # features × samples
  sample_metadata = metaData,               # has a column "AgeGroup"
  group_col       = "AgeGroup",
  feature_ids     = rownames(expression_mat),
  cor_method      = "spearman",
  corr_threshold  = 0.60,                   # keep strong |rho|
  pval_adjust     = "fdr",
  pval_cutoff     = 0.05,
  resample        = TRUE,                   # bootstrap‑like repeats
  samples_per_group = 20,                   # balanced sampling
  n_repeats       = 50,
  save_rds        = TRUE,                   # auto‑named RDS into results dir
  verbose         = TRUE
)
# Later you can reload via: readRDS(attr(adjacency, "rds_file"))
```

### 1b) `sc_buildAdjacency()` — **single‑cell** Seurat → per‑cell‑type adjacency

``` r
## \donttest{
library(Seurat)
data("pbmc_small")
pbmc_small$celltype <- as.character(Idents(pbmc_small))   # or your own labels

adj_list <- sc_buildAdjacency(
  seurat_obj     = pbmc_small,
  assay          = DefaultAssay(pbmc_small),
  slot           = "data",
  cell_type_col  = "celltype",
  cor_method     = "spearman",
  corr_threshold = 0.6,
  pval_adjust    = "fdr",
  pval_cutoff    = 0.1,
  resample       = FALSE,
  save_rds       = FALSE,
  verbose        = TRUE
)
str(adj_list, max.level = 1)
## }
```

### 2) `edgesFromAdjacency()` — matrices (or lists of matrices) → edge tables

``` r
# Dense matrix
A <- matrix(c(0,0.8,0.8,0), 2, 2, dimnames=list(c("a","b"), c("a","b")))
edgesA <- edgesFromAdjacency(A, directed = FALSE, drop_zeros = TRUE, min_abs = 0)

# Per‑group list, flattened with layer label
L <- list(G1 = A, G2 = A * 0.5)
edgesL <- edgesFromAdjacency(L, flatten = TRUE, id_cols = "layer")

# Nested list (repeat -> group); ragged OK
NL <- list(`1`=list(E1=A, E2=A), `2`=list(E1=A))
edgesNL <- edgesFromAdjacency(NL, flatten = TRUE, id_cols = c("repeat","layer"))
```

### 3) `consensusEdges()` — across‑repeat edge consensus per layer (CSV/XLSX/RDS)

``` r
# Tidy edge table with 'layer' and 'rep' columns
df <- data.frame(
  layer  = rep(c("A","A","B","B"), each = 3),
  rep    = rep(c("r1","r2"), times = 6),
  from   = c("g1","g2","g1","g1","g2","g2","g1","g3","g1","g2","g3","g1"),
  to     = c("g2","g3","g3","g2","g3","g1","g3","g2","g2","g1","g1","g3"),
  weight = rnorm(12)
)

cons_list <- consensusEdges(
  df,
  prop_present = 0.7,            # or min_reps = 3 to override
  summary      = "median",       # across repeats
  as_list      = TRUE,           # list per layer
  write_csv    = TRUE,
  write_xlsx   = TRUE,
  save_to_rds  = TRUE,
  csv_prefix   = "consensus_edges"
)
```

### 4) `build_multiNet()` — assemble multilayer from per‑layer **edge lists**

``` r
net <- build_multiNet(
  edgeListPerLayer     = cons_list,  # from consensusEdges(..., as_list=TRUE)
  nodesMetadata        = genes_info, # optional actor attributes
  featureID_col        = "GeneName",
  nodeAttrCols         = "GeneType",
  directed             = FALSE,
  aggregate_duplicates = "mean",     # if duplicate edges
  save_to_rds          = TRUE,
  verbose              = TRUE
)
```

### 5) `build_multinet_from_graphs()` — assemble multilayer from **igraph** graphs

``` r
# Suppose 'graphs' is a named list of per‑layer igraph objects
net2 <- build_multinet_from_graphs(
  graphs_list          = graphs,
  nodesMetadata        = genes_info,
  featureID_col        = "GeneName",
  nodeAttrCols         = "GeneType",
  save_layers_graphml  = TRUE,       # optional per‑layer GraphML exports
  graphml_prefix       = "layer",
  save_to_rds          = TRUE
)
```

### 6) `edges_list_to_graphs()` and 7) `graphs_to_edges_list()` — convert back and forth

``` r
Gs      <- edges_list_to_graphs(cons_list, directed = FALSE)
edgesDf <- graphs_to_edges_list(Gs)  # tidy data.frame with layer column
```

### 8) `compare_multinets()` — compare two or more multilayer networks

``` r
cmp <- compare_multinets(
  nets         = list(A = net, B = net2),
  by           = "union",      # or "intersect", "strict"
  edge_mode    = "undirected",
  ignore_missing = TRUE,
  write_csv    = TRUE,         # CSVs per layer + manifest
  prefix       = "cmp_AB"
)
```

### 9) `add_network_attributes()` — attach actor/edge metadata to the multilayer net

``` r
net <- add_network_attributes(
  net,
  nodesMetadata = genes_info,
  featureID_col = "GeneName",
  nodeAttrCols  = c("GeneType"),
  edgesMetadata = edges_meta_df,   # or a named list per layer
  edge_from     = "from", edge_to = "to", edge_layer_col = "layer",
  edgeAttrCols  = c("score","evidence"),
  edgeAggregate = "mean",
  fillMissingNodes = TRUE, nodeFillValue = NA,
  auto_detect_keys  = TRUE,        # robust actor ID matching
  save_report       = TRUE,        # saves TXT + RDS summary
  export_edge_join_csv = TRUE, edge_join_prefix = "edge_attach_debug",
  verbose = TRUE
)
```

### 10) `network_attributes()` — report what attributes are attached

``` r
attr_report <- network_attributes(
  net,
  show_in_rstudio = TRUE,
  write_csv       = TRUE,
  file_prefix     = "net_attrs"
)
```

### 11) `detectCom()` — community detection on the supra‑graph

``` r
comm <- detectCom(
  net,
  method      = "louvain",         # "infomap", "clique" (k‑clique union), etc.
  edgeWeight  = "count",           # supra edge weighting scheme
  min.actors  = 15,                # prune tiny communities
  min.layers  = 2,                 # keep those spanning ≥ 2 layers
  seed        = 1,
  write_csv   = TRUE               # saves community CSV
)
```

### 12) `plotCom()` — community layout + CSV export of node positions

``` r
p <- plotCom(
  net, communities = comm,
  layout = "multiforce", gravity = 0.30,   # nice cross‑layer layout
  show_in_rstudio = TRUE,
  save_plot = TRUE, format = "png", width = 10, height = 8, dpi = 300,
  save_df   = TRUE, df_format = "csv", df_prefix = "layout_multiforce"
)
```

### 13) `com_vs_samples()` — sanity: community counts vs sample counts

``` r
cs <- com_vs_samples(
  comm,
  num_samples = table(metaData$AgeGroup),  # vector with per‑layer sample sizes
  show_in_rstudio = TRUE,
  save_plot = TRUE
)
```

### 14) `analyze_actor_overlap()` — actor Jaccard overlaps across layers

``` r
Aov <- analyze_actor_overlap(
  net,
  reorder         = TRUE,   # reorder layers to cluster similarity
  drop_zero_diag  = TRUE,
  palette         = "RdBu",
  show_in_rstudio = TRUE, save_plot = TRUE
)
```

### 15) `analyze_edge_overlap()` — edge Jaccard overlaps across layers

``` r
Eov <- analyze_edge_overlap(
  net,
  reorder         = TRUE,
  drop_zero_diag  = TRUE,
  palette         = "RdBu",
  show_in_rstudio = TRUE, save_plot = TRUE
)
```

### 16) `layer_metrics()` — per‑layer size/centrality/path summaries

``` r
LM <- layer_metrics(
  net,
  layers        = multinet::layers_ml(net),
  directed      = FALSE,
  centralities  = c("degree","betweenness","closeness"),
  paths         = TRUE,             # components, APL, diameter
  write_csv     = TRUE,             # per‑layer + manifest CSVs
  prefix        = "layer_metrics"
)
```

### 17) `annotateCom()` — add feature attributes (e.g., GeneType) to community rows

``` r
commA <- annotateCom(
  comm,
  nodesMetadata = genes_info,
  featureID_col = "GeneName",
  attrCols      = c("GeneType"),
  write_csv     = TRUE
)
```

### 18) `sumComFeat()` — summaries by feature type (actor/community/layer)

``` r
summTF <- sumComFeat(
  commA,
  feature_type  = "TF",        # e.g., "TF", "lncRNA"
  ignore_case   = TRUE,
  normalize_ids = TRUE,
  write_csv     = TRUE,        # multiple CSVs (actor, by‑community, by‑layer)
  prefix        = "TF_summaries"
)
```

### 19) `get_FeatureDeg()` — degrees (per layer / long)

``` r
degTF <- get_FeatureDeg(
  net,
  featureList = subset(commA, GeneType == "TF")$actor,
  layers      = multinet::layers_ml(net),
  long        = TRUE
)
```

### 20) `gp_enrich_multinet()` — g:Profiler enrichment **per layer**

``` r
gp <- gp_enrich_multinet(
  net,
  organism    = "hsapiens",
  sources     = c("GO:BP","REAC"),
  layer_order = multinet::layers_ml(net),
  show_terms  = 8,        # top terms in dot plot
  ton_show    = TRUE,     # build term‑only overlap network
  ton_top_terms   = 12,
  ton_min_overlap = 1,
  ton_layout      = "fr",
  format      = "png",
  width       = 10, height = 7, dpi = 300,
  show_in_rstudio = TRUE
)
```

### 21) `go_enrichment_report()` — GO (clusterProfiler) for a vector **or** per‑layer

``` r
## Single set (vector of symbols)
go1 <- go_enrichment_report(
  genes = c("RHCG","SPRR1A","SPRR2A","SPRR3","KRT13","KRT4","CRNN"),
  OrgDb = org.Hs.eg.db, ont = "BP",
  file_prefix = "demo_GO_BP", show_in_rstudio = TRUE
)

## Per layer (each layer’s vertex set)
go2 <- go_enrichment_report(
  net  = net,
  OrgDb = org.Hs.eg.db, ont = "BP",
  layer_order   = multinet::layers_ml(net),
  universe_mode = "union",
  one_pdf       = TRUE,     # stitch into one PDF
  file_prefix   = "GO_BP_byLayer"
)
```

### 22) `gsea_msigdbr_layer_pairs()` — **GSEA** (MSigDB, fgsea) across **consecutive layer pairs**

``` r
gsea <- gsea_msigdbr_layer_pairs(
  net,
  layer_order     = multinet::layers_ml(net), # defines consecutive pairs
  # scores_by_layer optional; defaults to vertex attrs (e.g., cor/padj)
  collections     = c("H","CP:REACTOME"),     # MSigDB 2023+ style
  top_n           = 10,                       # top up/down dot plots
  format          = "png",
  show_in_rstudio = TRUE
)
```

### 23) `animate_multiNet()` — interactive HTML animation

``` r
anim <- animate_multiNet(
  net,
  communities   = comm,      # data.frame: actor, layer, com
  layout        = "kamadakawai",
  slice.par     = list(start = 0, interval = 1, aggregate.dur = 1),
  displaylabels = TRUE, vertex.cex = 0.9,
  seed = 1
)
# => omicsDNA_results/multiNet_animation_<timestamp>.html
```

### 24) `filmstrip_multiNet()` — static grid (one panel per layer)

``` r
fs <- filmstrip_multiNet(
  net, communities = comm, layout = "kamadakawai",
  ncol = 4, format = "png", width = 12, height = 8, dpi = 300, seed = 1
)
```

### 25) `plotActivityTimeline()` — Gantt‑style lifespans of edges/vertices

``` r
timeline_edges  <- plotActivityTimeline(net, type = "edge",   top_n = 300)
timeline_actors <- plotActivityTimeline(net, type = "vertex", top_n = 100, vertex.col = "tomato")
```

### 26) `plot_layer_diff()` — pairwise layer “new / lost / common” edges

``` r
plot_layer_diff(
  net,
  L1 = "E1", L2 = "E2",
  layout = "fr", edge_alpha = 0.7, vertex_cex = 0.7,
  save_edge_csv = TRUE, file = "diff_E1_vs_E2.png"
)
```

### 27) `grid_layer_diffs()` — multi‑panel grid of **all consecutive pair** differences

``` r
grid_layer_diffs(
  net,
  ncol = 4,
  layout = "fr",
  edge_alpha = 0.7, vertex_cex = 0.7,
  file = "diffs_all_pairs.png", width = 12, height = 10, dpi = 300
)
```

------------------------------------------------------------------------

## Reproducibility notes

Most functions save outputs under
`getOption("mlnet.results_dir","omicsDNA_results")` with timestamped
names; many return a `file` / `files` attribute pointing to what was
written. Set a seed before stochastic steps (resampling, layouts,
communities):

``` r
set.seed(1)
```

------------------------------------------------------------------------

## 🙋 Feedback & contributions

Issues and pull requests are welcome. Please provide a minimal example
and `sessionInfo()` where possible.

------------------------------------------------------------------------

## 📄 License

MIT © VafaeeLab

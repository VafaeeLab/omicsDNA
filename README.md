omicsDNA
================

- [omicsDNA
  <img src="man/figures/logo.png" align="right" height="110"/>](#omicsdna-)
  - [Installation](#installation)
  - [Function‑by‑function examples](#functionbyfunction-examples)
    - [1) `buildAdjacency()` — group‑wise correlation → adjacency (with
      optional
      resampling)](#1-buildadjacency--groupwise-correlation--adjacency-with-optional-resampling)
    - [1.1) `sc_buildAdjacency()` — **single‑cell** Seurat →
      per‑cell‑type
      adjacency](#11-sc_buildadjacency--singlecell-seurat--percelltype-adjacency)
    - [2) `edgesFromAdjacency()` — matrices (or lists of matrices) →
      edge
      tables](#2-edgesfromadjacency--matrices-or-lists-of-matrices--edge-tables)
    - [3) `consensusEdges()` — across‑repeat edge consensus per layer
      (CSV/XLSX/RDS)](#3-consensusedges--acrossrepeat-edge-consensus-per-layer-csvxlsxrds)
    - [4) `build_multiNet()` — assemble multilayer from per‑layer **edge
      lists**](#4-build_multinet--assemble-multilayer-from-perlayer-edge-lists)
    - [5) `edges_list_to_graphs()` and 6) `graphs_to_edges_list()` —
      convert back and
      forth](#5-edges_list_to_graphs-and-6-graphs_to_edges_list--convert-back-and-forth)
    - [7) `build_multinet_from_graphs()` — assemble multilayer from
      **igraph**
      graphs](#7-build_multinet_from_graphs--assemble-multilayer-from-igraph-graphs)
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
    - [13) `edges_vs_samples()` — sanity: community counts vs sample
      counts](#13-edges_vs_samples--sanity-community-counts-vs-sample-counts)
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
      (actor)](#18-sumcomfeat--summaries-by-feature-type-actor)
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
(`expression_mat`, `metaData`, …) with your data.

### 1) `buildAdjacency()` — group‑wise correlation → adjacency (with optional resampling)

``` r

setwd("~/Abir/test_data/")

# Importing the expression matrix with sample IDs as column names and genes as row names.
expression_mat <- read.csv("~/Abir/test_data/norm_DEGs_n_stableGenes_V.csv")
n <- expression_mat$X
expression_mat <- expression_mat[ , -1]
rownames(expression_mat) <- n

# Metadata dataframe must have the "sample column"
metaData <- read_csv("~/Abir/test_data/metaData_V.csv")
metaData <- as.data.frame(metaData)
colnames(metaData)[1] <- "sample"
head(metaData)
                    sample muscle_grp age_group2 age_label       group
1 GTEX-111FC-0826-SM-5GZWO   left_Ven    [65,70)        M4 left_Ven_M4
2 GTEX-111YS-0426-SM-5987O   left_Ven    [65,70)        M4 left_Ven_M4
3 GTEX-1122O-0826-SM-5GICV   left_Ven    [60,65)        M3 left_Ven_M3
4 GTEX-117YX-1126-SM-5H128   left_Ven    [60,65)        M3 left_Ven_M3
5 GTEX-11DXX-0326-SM-5PNWC   left_Ven    [65,70)        M4 left_Ven_M4
6 GTEX-11DXZ-0626-SM-5GU77   left_Ven    [65,70)        M4 left_Ven_M4

head(colnames(expression_mat))
[1] "GTEX.111FC.0826.SM.5GZWO" "GTEX.111YS.0426.SM.5987O" "GTEX.1122O.0826.SM.5GICV" "GTEX.117YX.1126.SM.5H128"
[5] "GTEX.11DXX.0326.SM.5PNWC" "GTEX.11DXZ.0626.SM.5GU77"

# match the names in metadata and column names expression matrix.
metaData$sample <- chartr("-", ".", metaData$sample) 

metaData$age_label %>% table()
# E1 E2 M1 M2 M3 M4 
# 40 43 41 44 46 21

set.seed(1)
adjacency <- buildAdjacency(
  dataMatrix      = expression_mat,        
  sample_metadata = metaData,               
  group_col       = "age_label",
  feature_ids     = rownames(expression_mat),
  cor_method      = "spearman",
  corr_threshold  = 0.35,                   
  pval_adjust     = "none",
  pval_cutoff     = 0.05,
  resample        = TRUE,                   
  samples_per_group = 20,                   
  n_repeats       = 10,
  save_rds        = TRUE,                   
  verbose         = TRUE
)

# Later you can reload via: readRDS(attr(adjacency, "rds_file"))
```

### 1.1) `sc_buildAdjacency()` — **single‑cell** Seurat → per‑cell‑type adjacency

For a comprehensive detailed analysis, refer to following vignette on
human foetal kidney single cell RNA-seq data:

[sc_omicsDNA](https://rstudio.anc-lab.cloud.edu.au/~mpir0002/22_09_2025_sc_omicsDNA_v3.html)

``` r

library(Seurat)
data("pbmc_small")
pbmc_small$celltype <- as.character(Idents(pbmc_small))   # or your own labels

sc_adjacency <- sc_buildAdjacency(
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

```

### 2) `edgesFromAdjacency()` — matrices (or lists of matrices) → edge tables

``` r
# Dense matrix
A <- matrix(c(0,0.8,0.8,0), 2, 2, dimnames=list(c("a","b"), c("a","b")))
edgesA <- edgesFromAdjacency(A, directed = FALSE, drop_zeros = TRUE, min_abs = 0)

# Nested list (repeat -> group); ragged OK
NL <- list(`1`=list(E1=A, E2=A), `2`=list(E1=A))
edgesNL <- edgesFromAdjacency(NL, flatten = FALSE, id_cols = c("repeat","layer"))

edges_list <- edgesFromAdjacency(
  adjacency,
  flatten = FALSE,
  id_cols = c("repeat","group")
)
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


cons_list <- consensusEdges(
  edges_list,
  prop_present = 0.7,
  summary      = "mean",
  as_list      = TRUE
)

```

### 4) `build_multiNet()` — assemble multilayer from per‑layer **edge lists**

``` r

gene_info <- read_csv("gene_info2.csv", skip = 1)
colnames(gene_info) <- c("GeneId" , "GeneName" , "GeneType")
head(gene_info)
# A tibble: 6 × 3
  GeneId            GeneName        GeneType      
  <chr>             <chr>           <chr>         
1 ENSG00000243485.5 MIR1302-2HG     lncRNA        
2 ENSG00000237613.2 FAM138A         lncRNA        
3 ENSG00000186092.7 OR4F5           protein_coding
4 ENSG00000238009.6 ENSG00000238009 lncRNA        
5 ENSG00000239906.1 ENSG00000239906 lncRNA        
6 ENSG00000241860.7 ENSG00000241860 lncRNA    


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

### 5) `edges_list_to_graphs()` and 6) `graphs_to_edges_list()` — convert back and forth

``` r
list_graphs <- edges_list_to_graphs(cons_list, directed = TRUE)

edges_back <- graphs_to_edges_list(
  graphs_list             = list_graphs,
  weight_attr             = "weight",   # use the weight attribute
  default_weight          = 1,          # fill in if missing
  canonicalize_undirected = TRUE,       # put smaller name first for undirected graphs
  verbose                 = TRUE
)  
```

### 7) `build_multinet_from_graphs()` — assemble multilayer from **igraph** graphs

``` r

# Suppose 'list_graphs' is a named list of per‑layer igraph objects
layer_order <- names(cons_list)  # keep the current naming/ordering
net_graphs <- build_multinet_from_graphs(
  graphsPerLayer     = list_graphs,
  layerOrder         = layer_order,    # same ordering
  nodesMetadata      = NULL,           # (match what you passed above)
  featureID_col      = "GeneName",
  nodeAttrCols       = "GeneType",
  save_layers_graphml= TRUE,       # optional per‑layer GraphML exports
  verbose            = FALSE,
  graphml_prefix     = "layer",
  save_to_rds        = TRUE
)
```

Build a second graph to make a comparison between two multi-layer graphs
using the next function

``` r

net_edges <- build_multiNet(
  edgeListPerLayer   = cons_list,
  layerOrder         = layer_order,
  directed           = FALSE,
  nodesMetadata      = NULL,           # or your genes_info + featureID_col, etc.
  aggregate_duplicates = "none",
  save_to_rds        = FALSE,
  verbose            = FALSE
)
```

### 8) `compare_multinets()` — compare two or more multilayer networks

``` r

diff_report <- compare_multinets(
  net_a    = net_edges,
  net_b    = net_graphs,
  tol      = 1e-9,
  show_max = 12,   # show up to 12 items per mismatch category
  verbose  = TRUE  # print a full, readable diff to the console
)

```

**compare_multinets() helps compare common and different edges between
layers of two multilayer graphs, can be used in the downstream
enrichment analysis with aim of comparing two biological contexts across
same time points.**

### 9) `add_network_attributes()` — attach actor/edge metadata to the multilayer net

``` r
net <- add_network_attributes(
  net,
  nodesMetadata = genes_info,
  featureID_col = "GeneName",
  nodeAttrCols  = "GeneType",
  edgesMetadata = cons_list,                          # named list: one df per layer
  edge_from     = "from",
  edge_to       = "to",
  edgeAttrCols  = c("n_present","n_repeats","prop_present"),
  actor_key_normalize = c("strip_version","trim","tolower"),
  auto_detect_keys    = FALSE,
  map_edges_by_actor_key = TRUE,
  directed_network    = FALSE,
  verbose             = TRUE,
  save_report         = TRUE,                         # writes a small report
  export_edge_join_csv = FALSE                        # flip to TRUE if debugging joins
)

# List attributes available (names + types)
multinet::attributes_ml(net)
```

### 10) `network_attributes()` — report what attributes are attached

``` r

# After you enriched the network with add_network_attributes(...):
rep <- network_attributes(
  net,
  show_examples = TRUE,  # print small tables of actors/edges with attrs
  show_max      = 12,    # up to 12 names/rows per section
  verbose       = TRUE
)

rep$added$edge_attrs         # what edge attrs were attached (per report)
rep$present$actor            # current actor-attribute definitions
rep$samples$edges            # first rows of edges with attribute columns
```

### 11) `detectCom()` — community detection on the supra‑graph

``` r

# Generalised Louvain (canonical method), count-based supra-graph weights
comm <- detectCom(
  net,
  method      = "glouvain",   # <- preferred spelling
  edgeWeight  = "count",
  min.actors  = 15,
  min.layers  = 2,
  seed        = 1,
  save_to_rds = TRUE,
  write_csv   = TRUE,
  verbose     = TRUE
)


# Clique-based detection with explicit k; keep only sizeable, multilayer communities
comm_clq <- detectCom(
  net,
  method      = "clique",
  clique.k    = 4,        # cohesion threshold in the supra-graph
  min.actors  = 20,       # post-detection retention threshold
  min.layers  = 2,
  edgeWeight  = "count",
  write_summary_csv = TRUE,
  verbose     = TRUE
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

### 13) `edges_vs_samples()` — sanity: community counts vs sample counts

``` r

ns <- c(E1 = 46, E2 = 22, M1 = 16, M2 = 27, M3 = 23, M4 = 26)

edges_per_layer <- edges_vs_samples(
  communities = NULL,
  num_samples = ns,
  net = net
)
```

### 14) `analyze_actor_overlap()` — actor Jaccard overlaps across layers

``` r

Aov <- analyze_actor_overlap(
  net,
  reorder         = TRUE   # reorder layers to cluster similarity
)
```

### 15) `analyze_edge_overlap()` — edge Jaccard overlaps across layers

``` r

Eov <- analyze_edge_overlap(
  net,
  reorder         = TRUE
)
```

### 16) `layer_metrics()` — per‑layer size/centrality/path summaries

``` r

LM <- layer_metrics(net)

LM$summary[, c("layer","n_nodes","n_edges","density")]
```

### 17) `annotateCom()` — add feature attributes (e.g., GeneType) to community rows

``` r

comm_annot <- annotateCom(
  communities   = comm,
  nodesMetadata = genes_info,
  featureID_col = "GeneName",
  attribute     = "GeneType",
  actor_normalize = c("strip_version","trim","tolower"),
  fill_missing  = "Unknown",
  write_csv     = TRUE,
  output_prefix = "communities_with_GeneType",  # alias kept for convenience
  verbose       = TRUE
)

head(comm_annot)
table(comm_annot$GeneType, useNA = "ifany")
```

### 18) `sumComFeat()` — summaries by feature type (actor)

``` r

# lncRNA
lnc_summary <- sumComFeat(
  communities   = comm_annot,
  feature_col   = "GeneType",
  feature_type  = "lncRNA",        # e.g., "TF", "lncRNA"
  write_csv     = TRUE
)

attr(lnc_summary, "files")  # where the three CSVs were written

# protein-coding
pc_summary <- sumComFeat(comm_annot, feature_col = "GeneType", feature_type = "protein_coding")

# multiple feature types at once
mix_summary <- sumComFeat(comm_annot, feature_col = "GeneType", feature_type = c("lncRNA","TF"), write_csv = TRUE)
```

### 19) `get_FeatureDeg()` — degrees (per layer / long)

``` r
degTF <- get_FeatureDeg(
  net,
  featureList = subset(commA, GeneType == "TF")$actor,
  layers      = multinet::layers_ml(net)
)

# 1) Use any feature set you like (lncRNAs, TFs, enzymes, …)
feat_list  <- unique(subset(comm_annot, GeneType == "protein_coding")$actor)
head(feat_list)
all_layers <- multinet::layers_ml(net)

deg_lnc <- get_FeatureDeg(
  net, featureList = feat_list, layers = all_layers,
  directed = FALSE, mode = "all",
  write_csv   = TRUE,
  output_file = "protein_coding_degree_byLayer.csv"  # saved under getOption("mlnet.results_dir","results")
)

# 2) Protein-coding instead (no file written)
pc_list <- unique(subset(comm_annot, GeneType == "protein_coding")$actor)
deg_pc  <- get_FeatureDeg(net, featureList = pc_list, layers = all_layers)
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

Absolutely — I’ll paste the files inline so you can copy‑paste them directly.

> **Provenance.** The updates below are assembled from your latest `OmicsDNA_v8` script and your previous README content. &#x20;

---

## 📄 README.Rmd

````markdown
---
title: "omicsDNA"
output:
  github_document:
    toc: true
    toc_depth: 3
editor_options:
  chunk_output_type: console
---

<!-- Updated from your OmicsDNA_v8 script + prior README. :contentReference[oaicite:2]{index=2} :contentReference[oaicite:3]{index=3} -->

# omicsDNA <img src="man/figures/logo.png" align="right" height="120"/>

*A multilayer network toolkit for grouped omics (and other high-dimensional) data.*

`omicsDNA` builds **layered networks** from grouped data (e.g., age groups, disease stages), stabilises edges via resampling + consensus, assembles a **multilayer** object, detects and visualises **communities**, quantifies **layer-wise structure**, compares **overlaps**, and produces **static** and **dynamic** visualisations with reproducible, timestamped outputs. Although examples use gene expression, the workflow is **domain-agnostic**. :contentReference[oaicite:4]{index=4}

---

## 🔍 What does `omicsDNA` do?

- Build **correlation networks per layer** with optional resampling and **consensus**;  
- Assemble a **multilayer** network and **annotate** actors/edges;  
- Detect **communities/modules** on a supra-graph; map them back to layers;  
- Quantify **layer structure** (nodes, edges, density, centralities, paths);  
- Compare **actor/edge overlaps** across layers;  
- Summarise **feature classes** (e.g., lncRNA/TF) within/between communities;  
- Generate **dynamic** (HTML/GIF/MP4) and **static** (filmstrip, timelines) visuals. :contentReference[oaicite:5]{index=5}

---

## 🆕 What’s new / expanded in this update

- **Single‑cell support:** `sc_buildAdjacency()` builds per–cell‑type adjacency matrices directly from a **Seurat** object (no pseudobulk). :contentReference[oaicite:6]{index=6}
- **Richer consensus:** `consensusEdges()` now supports `min_reps`, `as_list`, and **XLSX**/CSV/RDS outputs. :contentReference[oaicite:7]{index=7}
- **Community detection:** `detectCom()` gains `"clique"` (k‑clique union), `"abacus"` (if available), and summary CSVs; robust layer mapping & filters (`min.actors`, `min.layers`). :contentReference[oaicite:8]{index=8}
- **Feature summaries:** `sumComFeat()` now handles case‑insensitive matching, actor ID normalisation, and multi‑CSV outputs. :contentReference[oaicite:9]{index=9}
- **Change visualisation:** `plot_layer_diff()` and `grid_layer_diffs()` for pairwise and multi‑panel layer differences (new/lost/common edges). :contentReference[oaicite:10]{index=10}
- **Dynamics:** `plotActivityTimeline()` (Gantt‑style lifespans). `animate_multiNet()` (HTML). Optional MP4/GIF via `animate_multiNet_mp4()`. :contentReference[oaicite:11]{index=11} :contentReference[oaicite:12]{index=12}
- **Enrichment:** `gp_enrich_multinet()` (g:Profiler across layers), `go_enrichment_report()` (GO with clusterProfiler), and `gsea_msigdbr_layer_pairs()` (MSigDB + fgsea across layer pairs). :contentReference[oaicite:13]{index=13} :contentReference[oaicite:14]{index=14} :contentReference[oaicite:15]{index=15}

---

## 📦 Installation

```r
# Development version
install.packages("devtools")
devtools::install_github("VafaeeLab/omicsDNA")

# Optional (dynamic visuals)
install.packages(c("ndtv","networkDynamic"))   # HTML animations/filmstrips
install.packages(c("gifski","gganimate","png","av"))  # GIF/MP4 exporters (optional)
````

Set a default results directory (recommended):

```r
options(mlnet.results_dir = "~/omicsDNA_results")
```

---

## ⚡ Quick start

```r
library(omicsDNA)
set.seed(1)

# Inputs:
#  - expression_mat: numeric matrix (features x samples)
#  - metaData: data.frame with sample->group mapping (e.g., AgeGroup)
#  - genes_info: feature annotations (e.g., GeneName, GeneType)

# 1) Build adjacency per group (optional resampling)
adj <- buildAdjacency(
  dataMatrix      = expression_mat,
  sample_metadata = metaData,
  group_col       = "AgeGroup",
  feature_ids     = rownames(expression_mat),
  cor_method      = "spearman",
  corr_threshold  = 0.6,
  pval_cutoff     = 0.05,
  resample        = TRUE,
  n_repeats       = 50
)

# 2) Edges + consensus (expanded)
edges <- edgesFromAdjacency(adj)
cons  <- consensusEdges(
  edges,
  prop_present = 0.7,
  # new: can also do min_reps, write CSV/XLSX, and save RDS
  write_csv  = TRUE, write_xlsx = TRUE, save_to_rds = TRUE
)

# 3) Assemble multilayer + annotate
net <- build_multiNet(cons)
net <- add_network_attributes(
  net,
  nodesMetadata = genes_info,
  featureID_col = "GeneName",
  nodeAttrCols  = c("GeneType")
)

# 4) Communities (supra-graph) + visual
comm <- detectCom(net, method = "louvain", edgeWeight = "count",
                  min.actors = 15, min.layers = 2, seed = 1)
plotCom(net, comm, layout = "multiforce", gravity = 0.3,
        show_in_rstudio = TRUE, save_plot = TRUE, save_df = TRUE)

# 5) Structure, overlaps, sampling sanity
LM  <- layer_metrics(net)
Aov <- analyze_actor_overlap(net, reorder = TRUE)
Eov <- analyze_edge_overlap(net, reorder = TRUE)
cs  <- com_vs_samples(comm, num_samples = table(metaData$AgeGroup))

# 6) Feature-aware summaries (expanded)
commA  <- annotateCom(comm, genes_info, "GeneName", "GeneType", write_csv = TRUE)
summTF <- sumComFeat(commA, feature_type = "TF", ignore_case = TRUE, write_csv = TRUE)
degTF  <- get_FeatureDeg(net, featureList = subset(commA, GeneType == "TF")$actor)

# 7) Dynamics (pick any)
animate_multiNet(net, communities = comm, seed = 1)             # HTML (ndtv)
filmstrip_multiNet(net, communities = comm, ncol = 4, seed = 1) # PNG/PDF
plotActivityTimeline(net, type = "edge")                         # Gantt lifespans
```

> The steps above mirror the typical flow in the previous README, with additions noted inline.&#x20;

---

## 🧭 Typical analysis flow

1. **Build & stabilise:** `buildAdjacency()` → `edgesFromAdjacency()` → `consensusEdges()`
2. **Assemble & annotate:** `build_multiNet()` → `add_network_attributes()`
3. **Detect & visualise:** `detectCom()` → `plotCom()`
4. **Quantify & compare:** `layer_metrics()`, `analyze_actor_overlap()`, `analyze_edge_overlap()`
5. **Feature‑aware:** `annotateCom()`, `sumComFeat()`, `get_FeatureDeg()`
6. **Communicate change:** `animate_multiNet()`, `filmstrip_multiNet()`, `plotActivityTimeline()`, `plot_layer_diff()`, `grid_layer_diffs()` &#x20;

---

## 🧪 Methods (concise)

* **Per‑layer correlations:** keep edges if `|r| ≥ corr_threshold` & `p ≤ pval_cutoff`. With resampling, repeat on balanced draws.
* **Consensus:** retain edges present in ≥ `prop_present` (or `min_reps`) of repeats; summarise weights across reps; write **CSV/XLSX/RDS** as needed.&#x20;
* **Multilayer assembly:** per‑layer igraphs combined in `multinet::ml.network`.
* **Communities:** supra‑graph aggregated by **count**/**sum**; methods: Louvain, Infomap, **k‑clique** union, **ABACUS** (if available); filter by **min.actors**/**min.layers** with summary CSVs.&#x20;
* **Structure/overlaps:** degree, betweenness, eigenvector, closeness, PageRank, clustering, coreness; Jaccard for actors/edges; path summaries with guards.

---

## 🧩 Function map (current API)

### Build & assemble

| Function                   | Purpose                                                             |
| -------------------------- | ------------------------------------------------------------------- |
| `buildAdjacency()`         | Create per‑group adjacency matrices (optional resampling).          |
| `sc_buildAdjacency()`      | **Single‑cell**: Seurat → per‑cell‑type adjacency (no pseudobulk).  |
| `edgesFromAdjacency()`     | Convert matrices/lists to tidy edge tables.                         |
| `consensusEdges()`         | Keep edges robust across repeats; summarise weights; CSV/XLSX/RDS.  |
| `build_multiNet()`         | Assemble `multinet::ml.network` from per‑layer edges.               |
| `add_network_attributes()` | Attach node/edge attributes (robust ID matching).                   |

### Community detection & reporting

| Function      | Purpose                                                                  |
| ------------- | ------------------------------------------------------------------------ |
| `detectCom()` | Communities on a supra‑graph (Louvain/Infomap/**k‑clique**/**ABACUS**).  |
| `plotCom()`   | Multilayer community plot (+ save long `(actor, layer, com/cid)` table). |

### Structure, overlaps, sampling

| Function                  | Purpose                                              |
| ------------------------- | ---------------------------------------------------- |
| `layer_metrics()`         | Node‑ & layer‑level metrics; saves a run folder.     |
| `analyze_actor_overlap()` | Jaccard overlap of actor sets + heatmap.             |
| `analyze_edge_overlap()`  | Jaccard overlap of edge sets + heatmap.              |
| `com_vs_samples()`        | Relate #communities per layer to #samples per group. |

### Feature-aware summaries

| Function           | Purpose                                                                         |
| ------------------ | ------------------------------------------------------------------------------- |
| `annotateCom()`    | Add a feature attribute (e.g., `GeneType`) to community rows.                   |
| `sumComFeat()`     | Summaries by actor / by community / counts per community × layer (CSV bundle).  |
| `get_FeatureDeg()` | Degrees for any feature list across layers (wide/long).                         |

### Dynamics & comparative visuals

| Function                 | Purpose                                                        |
| ------------------------ | -------------------------------------------------------------- |
| `animate_multiNet()`     | Interactive HTML animation (`ndtv`).                           |
| `animate_multiNet_mp4()` | MP4/GIF via `gganimate`/`gifski`.                              |
| `filmstrip_multiNet()`   | Static grid (one panel per layer).                             |
| `plotActivityTimeline()` | Vertex/edge “lifespans” (Gantt‑style spells).                  |
| `plot_layer_diff()`      | Pairwise layer‑difference plot (new/lost/common edges).        |
| `grid_layer_diffs()`     | Multi‑panel grid of **layer differences** (publication‑ready). |

### Enrichment & functional analysis

| Function                     | Purpose                                                                               |
| ---------------------------- | ------------------------------------------------------------------------------------- |
| `gp_enrich_multinet()`       | g\:Profiler enrichment **per layer** + summary/term‑overlap networks.                 |
| `go_enrichment_report()`     | GO enrichment (single set **or** per‑layer) with clusterProfiler; PDFs & CSVs.        |
| `gsea_msigdbr_layer_pairs()` | **GSEA** across **layer pairs** using MSigDB + fgsea; dotplots, cNET, term networks.  |

---

## 🎬 Dynamics & change: quick recipes

**Interactive animation (HTML; per‑slice community colours)**

```r
anim <- animate_multiNet(
  net,
  communities   = comm,             # data.frame: actor, layer, com/cid
  layout        = "kamadakawai",
  slice.par     = list(start = 0, interval = 1, aggregate.dur = 1),
  displaylabels = TRUE, vertex.cex = 0.9, seed = 1
)
# => omicsDNA_results/multiNet_animation_<timestamp>.html
```

**Static filmstrip (PNG/PDF; one panel per layer)**

```r
fs <- filmstrip_multiNet(
  net, communities = comm, layout = "kamadakawai",
  ncol = 4, format = "png", width = 12, height = 8, dpi = 300, seed = 1
)
```

**Activity timelines (lifespans of edges/vertices)**

```r
timeline_edges  <- plotActivityTimeline(net, type = "edge",   show_in_rstudio = TRUE, save_plot = TRUE)
timeline_actors <- plotActivityTimeline(net, type = "vertex", show_in_rstudio = TRUE, save_plot = TRUE)
```

**Pairwise & grid differences**

```r
plot_layer_diff(net, L1 = "E1", L2 = "E2", layout = "fr",
                edge_alpha = 0.7, vertex_cex = 0.7)

grid_layer_diffs(net, ncol = 4, file = "diffs_all_pairs.png")
```



---

## 🧬 Single‑cell adjacency (Seurat → layers)

```r
# Build one adjacency matrix per cell type directly from a Seurat object
adj_sc <- sc_buildAdjacency(
  seurat_obj,
  cell_type_col       = "cell_type",
  assay               = "RNA",
  slot                = "data",
  feature_ids         = VariableFeatures(seurat_obj),
  max_cells_per_type  = 1000,
  cor_method          = "spearman",
  pval_adjust         = "fdr",
  pval_cutoff         = 0.05,
  corr_threshold      = 0.25
)
# next: edgesFromAdjacency(adj_sc) → consensusEdges(...) → build_multiNet(...)
```

> Computes correlations **across cells** per cell type; keeps edges by combined correlation & p‑value thresholds; returns a list of same‑ordered adjacency matrices.&#x20;

---

## 🧠 Enrichment examples

**g\:Profiler per layer over a multilayer network**

```r
res <- gp_enrich_multinet(
  net,
  organism          = "hsapiens",
  sources           = c("GO:BP","GO:MF","GO:CC","REAC"),
  layer_order       = multinet::layers_ml(net),
  significant       = TRUE,
  user_threshold    = 0.05,
  show_terms        = 12,          # dotplot + concept network cap
  per_layer_subdirs = TRUE,
  format            = "png",
  show_in_rstudio   = TRUE,
  ton_show          = TRUE, ton_top_terms = 10
)
# Inspect the run manifest:
read.csv(res$summary_file, stringsAsFactors = FALSE)
```

**GO (clusterProfiler): single set and per‑layer**

```r
## Single set
genes <- c("RHCG","SPRR1A","SPRR2A","SPRR3","KRT13","KRT4","CRNN")
go1 <- go_enrichment_report(genes = genes, OrgDb = org.Hs.eg.db, ont = "BP",
                            file_prefix = "demo_GO_BP", show_in_rstudio = TRUE)

## Per layer (each layer's vertex set)
go2 <- go_enrichment_report(net = net, OrgDb = org.Hs.eg.db, ont = "BP",
                            layer_order   = multinet::layers_ml(net),
                            universe_mode = "union",
                            one_pdf = TRUE, file_prefix = "GO_BP_byLayer")
go2$combined_pdf
```

**GSEA across layer pairs (MSigDB + fgsea)**

```r
gsea <- gsea_msigdbr_layer_pairs(
  net,
  layer_order = multinet::layers_ml(net),
  # scores_by_layer is optional; by default it will try to extract per-layer scores from net's vertex attrs
  top_n      = 10,
  format     = "png",
  show_in_rstudio = TRUE
)
# Check outputs in options('mlnet.results_dir')
```

&#x20;&#x20;

---

## 📊 Interpretation highlights

* **Consensus edges** → robust wiring; high `prop_present` / `min_reps` ≈ stability.
* **Communities** → size & span quantify breadth and cross‑layer stability.
* **Layer metrics** → density/components/APL/diameter contextualise modules.
* **Overlaps** → actor vs edge Jaccard separates participation vs wiring similarity.
* **Timelines** → stable vs transient vertices/edges; late emergence.&#x20;

---

## 🗂 Saving conventions & reproducibility

* Most functions write to `getOption("mlnet.results_dir","omicsDNA_results")` with **timestamps**; many return a `file`/`files` attribute.
* Set a seed before stochastic steps (resampling, layouts, communities):

```r
set.seed(1)
```



---

## 🙋 Feedback & contributions

Issues and pull requests are welcome. Please provide a minimal example and `sessionInfo()` where possible.

---

## 📄 License

MIT © VafaeeLab

````

---

## 📄 README.txt

```text
omicsDNA — A multilayer network toolkit for grouped omics (and other high-dimensional) data.

Core capabilities:
- Build correlation networks per layer (optional resampling + consensus)
- Assemble a multilayer network; annotate actors/edges
- Detect communities on a supra-graph; map back to layers
- Quantify layer structure; compare actor/edge overlaps
- Summarise feature classes (e.g., lncRNA/TF)
- Static & dynamic visuals (filmstrip, timelines, HTML/GIF/MP4)

Installation
------------
R:
  install.packages("devtools")
  devtools::install_github("VafaeeLab/omicsDNA")
Optional:
  install.packages(c("ndtv","networkDynamic","gifski","gganimate","png","av"))

Tip: options(mlnet.results_dir = "~/omicsDNA_results")

Quick start
-----------
library(omicsDNA); set.seed(1)

# Build per-group adjacency
adj <- buildAdjacency(expression_mat, metaData, group_col="AgeGroup",
                      cor_method="spearman", corr_threshold=0.6, pval_cutoff=0.05,
                      resample=TRUE, n_repeats=50)

# Edges + consensus (expanded outputs)
edges <- edgesFromAdjacency(adj)
cons  <- consensusEdges(edges, prop_present=0.7,
                        write_csv=TRUE, write_xlsx=TRUE, save_to_rds=TRUE)

# Multilayer + annotate
net <- build_multiNet(cons)
net <- add_network_attributes(net, nodesMetadata=genes_info,
                              featureID_col="GeneName", nodeAttrCols="GeneType")

# Communities + plot
comm <- detectCom(net, method="louvain", edgeWeight="count",
                  min.actors=15, min.layers=2, seed=1)
plotCom(net, comm, layout="multiforce", gravity=0.3,
        show_in_rstudio=TRUE, save_plot=TRUE, save_df=TRUE)

# Structure & overlaps
LM <- layer_metrics(net)
Aov <- analyze_actor_overlap(net, reorder=TRUE)
Eov <- analyze_edge_overlap(net, reorder=TRUE)

# Feature-aware summaries (expanded)
commA  <- annotateCom(comm, genes_info, "GeneName", "GeneType", write_csv=TRUE)
summTF <- sumComFeat(commA, feature_type="TF", ignore_case=TRUE, write_csv=TRUE)
degTF  <- get_FeatureDeg(net, featureList=subset(commA, GeneType=="TF")$actor)

# Dynamics
animate_multiNet(net, communities=comm, seed=1)
filmstrip_multiNet(net, communities=comm, ncol=4, seed=1)
plotActivityTimeline(net, type="edge")

Single-cell (Seurat) → per–cell-type adjacency
----------------------------------------------
adj_sc <- sc_buildAdjacency(seurat_obj,
  cell_type_col="cell_type", assay="RNA", slot="data",
  feature_ids=VariableFeatures(seurat_obj),
  max_cells_per_type=1000, cor_method="spearman",
  pval_adjust="fdr", pval_cutoff=0.05, corr_threshold=0.25)

Enrichment & functional analysis
--------------------------------
# g:Profiler per layer
res_gp <- gp_enrich_multinet(net, organism="hsapiens",
           sources=c("GO:BP","GO:MF","GO:CC","REAC"),
           layer_order=multinet::layers_ml(net),
           significant=TRUE, user_threshold=0.05,
           show_terms=12, per_layer_subdirs=TRUE, format="png")
read.csv(res_gp$summary_file, stringsAsFactors=FALSE)

# GO (clusterProfiler)
genes <- c("RHCG","SPRR1A","SPRR2A","SPRR3","KRT13","KRT4","CRNN")
go1 <- go_enrichment_report(genes=genes, OrgDb=org.Hs.eg.db, ont="BP",
                            file_prefix="demo_GO_BP", show_in_rstudio=TRUE)
go2 <- go_enrichment_report(net=net, OrgDb=org.Hs.eg.db, ont="BP",
                            layer_order=multinet::layers_ml(net),
                            universe_mode="union", one_pdf=TRUE,
                            file_prefix="GO_BP_byLayer")

# GSEA across layer pairs (MSigDB + fgsea)
gsea <- gsea_msigdbr_layer_pairs(net,
          layer_order=multinet::layers_ml(net), top_n=10,
          format="png", show_in_rstudio=TRUE)

Change & dynamics visuals
-------------------------
plot_layer_diff(net, L1="E1", L2="E2", layout="fr",
                edge_alpha=0.7, vertex_cex=0.7)

grid_layer_diffs(net, ncol=4, file="diffs_all_pairs.png")

timeline_edges  <- plotActivityTimeline(net, type="edge",   show_in_rstudio=TRUE, save_plot=TRUE)
timeline_actors <- plotActivityTimeline(net, type="vertex", show_in_rstudio=TRUE, save_plot=TRUE)

Interpretation highlights
-------------------------
- Consensus edges → robust wiring; high prop_present / min_reps ≈ stability
- Communities → size & span quantify breadth and cross-layer stability
- Layer metrics → density/components/APL/diameter contextualise modules
- Overlaps → actor vs edge Jaccard separates participation vs wiring similarity
- Timelines → stable vs transient vertices/edges; late emergence

Reproducibility
---------------
- Outputs are timestamped under options('mlnet.results_dir','omicsDNA_results')
- Use set.seed(1) before stochastic steps

License
-------
MIT © VafaeeLab
````

---

If you’d like me to trim or reorganize sections for CRAN/GitHub badges or add a changelog, say the word and I’ll inline those too.

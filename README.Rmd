````markdown
---
output: github_document
editor_options:
  chunk_output_type: inline
---

# omicsDNA <img src="man/figures/logo.png" align="right" height="120"/>

*A multilayer network toolkit for grouped omics (and other high-dimensional) data.*

`omicsDNA` builds **layered networks** from grouped data (e.g., age groups, disease stages), stabilises edges via resampling + consensus, assembles a **multilayer** object, detects and visualises **communities**, quantifies **layer-wise structure**, compares **overlaps**, and produces **static** and **dynamic** visualisations with reproducible, timestamped outputs.

Although examples use gene expression, the workflow is **domain-agnostic**: “actors” can be genes, sensors, taxa, survey items, etc.

---

## 🔍 What does `omicsDNA` do?

- Build **correlation networks per layer** with optional resampling and **consensus**;
- Assemble a **multilayer** network and **annotate** actors/edges;
- Detect **communities/modules** on a supra-graph; map them back to layers;
- Quantify **layer structure** (nodes, edges, density, centralities, paths);
- Compare **actor/edge overlaps** across layers;
- Summarise **feature classes** (e.g., lncRNA/TF) within/between communities;
- Generate **dynamic** (HTML/GIF/MP4) and **static** (filmstrip, timelines) visuals.

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
#  - feature_ids: optional subset of features to include

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

# 2) Edges + consensus
edges  <- edgesFromAdjacency(adj)
cons   <- consensusEdges(edges, prop_present = 0.7)

# 3) Assemble multilayer + annotate
net <- build_multiNet(cons)
net <- add_network_attributes(net, nodesMetadata = genes_info,
                              featureID_col = "GeneName",
                              nodeAttrCols  = c("GeneType"))

# 4) Communities (supra-graph) + visual
comm <- detectCom(net, method = "louvain", edgeWeight = "count",
                  min.actors = 15, min.layers = 2, seed = 1)
plotCom(net, comm, layout = "multiforce", gravity = 0.3,
        show_in_rstudio = TRUE, save_plot = TRUE, save_df = TRUE)

# 5) Structure, overlaps, sampling sanity
LM   <- layer_metrics(net)
Aov  <- analyze_actor_overlap(net, reorder = TRUE)
Eov  <- analyze_edge_overlap(net, reorder = TRUE)
cs   <- com_vs_samples(comm, num_samples = table(metaData$AgeGroup))

# 6) Feature-aware summaries
commA <- annotateCom(comm, genes_info, "GeneName", "GeneType", write_csv = TRUE)
summTF<- sumComFeat(commA, feature_type = "TF", write_csv = TRUE)
degTF <- get_FeatureDeg(net, featureList = subset(commA, GeneType == "TF")$actor)

# 7) Dynamics (pick any)
animate_multiNet(net, communities = comm, seed = 1)           # HTML
filmstrip_multiNet(net, communities = comm, ncol = 4, seed = 1) # PNG/PDF
plotActivityTimeline(net, type = "edge")                        # Gantt lifespans
```

---

## 🧭 Typical analysis flow

1. **Build & stabilise:** `buildAdjacency()` → `edgesFromAdjacency()` → `consensusEdges()`
2. **Assemble & annotate:** `build_multiNet()` → `add_network_attributes()`
3. **Detect & visualise:** `detectCom()` → `plotCom()`
4. **Quantify & compare:** `layer_metrics()`, `analyze_actor_overlap()`, `analyze_edge_overlap()`
5. **Feature-aware:** `annotateCom()`, `sumComFeat()`, `get_FeatureDeg()`
6. **Communicate change:** `animate_multiNet()`, `filmstrip_multiNet()`, `plotActivityTimeline()`

---

## 🧪 Methods (concise)

* **Per-layer correlations:** keep edges if `|r| ≥ corr_threshold` & `p ≤ pval_cutoff`. With resampling, repeat on balanced draws.
* **Consensus:** retain edges present in ≥ `prop_present` of repeats; carry summary weight & coverage.
* **Multilayer assembly:** per-layer igraphs combined in `multinet::ml.network`.
* **Communities:** supra-graph aggregated by **count** or **sum**; methods: Louvain, Infomap, k-clique union; ABACUS (if available).
* **Structure/overlaps:** degree, betweenness, eigenvector, closeness, PageRank, clustering, coreness; Jaccard for actors/edges; path summaries with guards for disconnected graphs.

---

## 🧩 Function map (current API)

### Build & assemble

| Function                   | Purpose                                                    |
| -------------------------- | ---------------------------------------------------------- |
| `buildAdjacency()`         | Create per-group adjacency matrices (optional resampling). |
| `edgesFromAdjacency()`     | Convert matrices/lists to tidy edge tables.                |
| `consensusEdges()`         | Keep edges robust across repeats; summarise weights.       |
| `build_multiNet()`         | Assemble `multinet::ml.network` from per-layer edges.      |
| `add_network_attributes()` | Attach node/edge attributes (robust ID matching).          |

### Community detection & reporting

| Function      | Purpose                                                              |
| ------------- | -------------------------------------------------------------------- |
| `detectCom()` | Communities on a supra-graph (Louvain/Infomap/k-clique/ABACUS).      |
| `plotCom()`   | Multilayer community plot (+ save long `(actor, layer, cid)` table). |

### Structure, overlaps, sampling

| Function                  | Purpose                                              |
| ------------------------- | ---------------------------------------------------- |
| `layer_metrics()`         | Node- & layer-level metrics; saves a run folder.     |
| `analyze_actor_overlap()` | Jaccard overlap of actor sets + heatmap.             |
| `analyze_edge_overlap()`  | Jaccard overlap of edge sets + heatmap.              |
| `com_vs_samples()`        | Relate #communities per layer to #samples per group. |

### Feature-aware summaries

| Function           | Purpose                                                             |
| ------------------ | ------------------------------------------------------------------- |
| `annotateCom()`    | Add a feature attribute (e.g., `GeneType`) to community rows.       |
| `sumComFeat()`     | Summaries by actor, by community, and counts per community × layer. |
| `get_FeatureDeg()` | Degrees for any feature list across layers (wide/long).             |

### Dynamics & comparative visuals

| Function                 | Purpose                                                             |
| ------------------------ | ------------------------------------------------------------------- |
| `animate_multiNet()`     | Interactive HTML animation (one slice per layer; `ndtv`).           |
| `filmstrip_multiNet()`   | Static grid (one panel per layer) with stable layout.               |
| `animate_multiNet_mp4()` | MP4/GIF via `gganimate`/`gifski` (optional fallback).               |
| `plotActivityTimeline()` | Vertex/edge “lifespans” (Gantt-style spells).                       |
| `grid_layer_diffs()`     | Multi-panel grid of **layer differences** (publication-ready).      |
| `plot_layer_diff()`      | Focused **pairwise** layer-difference plot (new/lost/common edges). |

---

## 🎬 Dynamics: quick recipes

**Interactive animation (HTML; per-slice community colours)**

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
# => omicsDNA_results/multiNet_filmstrip_<timestamp>.png
```

**Activity timelines (lifespans of edges/vertices)**

```r
timeline_edges  <- plotActivityTimeline(net, type = "edge",   show_in_rstudio = TRUE, save_plot = TRUE)
timeline_actors <- plotActivityTimeline(net, type = "vertex", show_in_rstudio = TRUE, save_plot = TRUE)
```

**Pairwise layer difference (who’s new, who disappears, what’s common)**

```r
# L1 vs L2 on a shared layout
plot_layer_diff(net, L1 = "E1", L2 = "E2", layout = "fr",
                edge_alpha = 0.7, vertex_cex = 0.7)
```

**Grid of consecutive differences (scan all transitions at once)**

```r
# E1→E2, E2→E3, ... in a facetted grid; saves a single PNG if 'file' is given
grid_layer_diffs(net, ncol = 4, file = "diffs_all_pairs.png")
```

---

## 📊 Interpretation highlights

* **Consensus edges** → robust wiring; high `prop_present` ≈ stability.
* **Communities** → size & span quantify breadth and cross-layer stability.
* **Layer metrics** → density/components/APL/diameter contextualise modules.
* **Overlaps** → actor vs edge Jaccard separates participation vs wiring similarity.
* **Timelines** → stable vs transient vertices/edges; late emergence.

---

## 🗂 Saving conventions & reproducibility

* Most functions write to `getOption("mlnet.results_dir","omicsDNA_results")` with **timestamps**; many return a `file`/`files` attribute.
* Set a seed before stochastic steps (resampling, layout, communities) to reproduce figures and tables:

  ```r
  set.seed(1)
  ```

---

## 🙋 Feedback & contributions

Issues and pull requests are welcome. Please provide a minimal example and `sessionInfo()` where possible.

---

## 📄 License

MIT © VafaeeLab

```
```

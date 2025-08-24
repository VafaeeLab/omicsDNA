
---

output: github\_document
editor\_options:
chunk\_output\_type: inline

---

# omicsDNA <img src="man/figures/logo.png" align="right" height="120"/>

`omicsDNA` is an R package for **multilayer network analysis** of grouped/conditioned data (e.g., age groups, disease stages). It helps you **build per‑group networks**, **stabilise** edges via resampling, **assemble** a multilayer object, **detect communities**, **quantify structure**, and **visualise changes**—with reproducible, timestamped outputs for easy reporting. Although we demonstrate with gene expression, the framework is **domain‑agnostic** (“features” can be genes, proteins, sensors, survey items, etc.).&#x20;

---

## 🔍 What does `omicsDNA` do?

* Build **correlation networks per layer** with optional resampling/consensus.
* Assemble a **multilayer network** and add node/edge annotations.
* Detect **communities/modules** on a supra‑graph; map them back to layers.
* Quantify **layer‑wise structure** and **cross‑layer overlaps** (actors & edges).
* Summarise **feature classes** within/between communities.
* Produce **static and dynamic** visualisations (filmstrip, interactive animation).&#x20;

---

## 📦 Installation

```r
# Install from GitHub (development version)
install.packages("devtools")
devtools::install_github("VafaeeLab/omicsDNA")

# Optional: set a project-wide results folder once per session
options(mlnet.results_dir = "~/omicsDNA_results")
```

> If your RStudio Server has limited system libraries (e.g., for `sf/GDAL`), you can still use all **ndtv/networkDynamic** visualisations; no system‑level GDAL is required for those.

---

## ⚡ Quick start

```r
library(omicsDNA)

# 0) Inputs
#  - expression_mat: numeric matrix (features x samples)
#  - metaData: data.frame with sample->group mapping (e.g., AgeGroup)
#  - genes_info: feature metadata (e.g., GeneName, GeneType)
#  - feature_ids: optional subset of features to include
set.seed(1)

# 1) Build per-group adjacency (optionally with resampling)
adj <- buildAdjacency(
  expression_mat = expression_mat,
  metaData       = metaData,
  group_col      = "AgeGroup",
  feature_ids    = rownames(expression_mat),
  cor_method     = "spearman",
  corr_threshold = 0.6,
  pval_cutoff    = 0.05,
  resample       = TRUE,    # repeat on balanced draws
  n_repeats      = 50
)

# 2) Convert to edges; keep edges robust across repeats
edges      <- edgesFromAdjacency(adj)
edges_cons <- consensusEdges(edges, prop_present = 0.7)

# 3) Assemble multilayer network and attach annotations
net <- build_multiNet(edges_cons)
net <- add_network_attributes(
  net,
  nodesMetadata  = genes_info,
  featureID_col  = "GeneName",
  nodeAttrCols   = c("GeneType")
)

# 4) Communities on a supra-graph, then visualise
comm <- detectCom(net, method = "louvain", edgeWeight = "count",
                  min.layers = 1, min.actors = 10, seed = 1)

plotCom(net, comm, layout = "multiforce", gravity = 0.3,
        show_in_rstudio = TRUE, save_plot = TRUE, save_df = TRUE)

# 5) Structure, overlaps, and summaries
LM   <- layer_metrics(net)                             # per-layer stats
Aov  <- analyze_actor_overlap(net, reorder = TRUE)     # actor Jaccard
Eov  <- analyze_edge_overlap(net,  reorder = TRUE)     # edge Jaccard

comm_annot <- annotateCom(comm, genes_info, "GeneName", "GeneType", write_csv = TRUE)
summ_TF    <- sumComFeat(comm_annot, feature_type = "TF", write_csv = TRUE)

# 6) Dynamics (choose one)
# (a) Interactive HTML animation (ndtv/networkDynamic)
animate_multiNet(net, communities = comm, layout = "kamadakawai", seed = 1)

# (b) Static filmstrip grid (PNG/PDF)
filmstrip_multiNet(net, communities = comm, layout = "kamadakawai",
                   ncol = 4, format = "png", seed = 1)
```

---

## 🧭 Typical analysis flow (and why)

1. **Build & stabilise**: `buildAdjacency()` → `edgesFromAdjacency()` → `consensusEdges()`
   *Focus on reproducible edges; avoid over‑fitting to a single sample split.*&#x20;
2. **Assemble & annotate**: `build_multiNet()` → `add_network_attributes()`
   *Carry forward feature/edge attributes you’ll interpret later.*&#x20;
3. **Detect & map communities**: `detectCom()` → `plotCom()`
   *Find global modules on a supra‑graph; inspect per‑layer realisations.*&#x20;
4. **Quantify structure & similarity**: `layer_metrics()`, overlaps
   *Relate communities to density, fragmentation, and cross‑layer wiring.*&#x20;
5. **Feature‑aware results**: `annotateCom()`, `sumComFeat()`, `get_FeatureDeg()`
   *Connect modules to feature classes (e.g., TF/lncRNA) and centrality.*&#x20;
6. **Communicate change over layers**: `animate_multiNet()` / `filmstrip_multiNet()` / `plotActivityTimeline()`
   *Make dynamics legible for manuscripts and talks.*

---

## 🧪 Methods

* **Correlation networks (per layer):** Pairwise correlations within group; keep edges if `|r| ≥ corr_threshold` & `p ≤ pval_cutoff`. With `resample=TRUE`, repeat on balanced draws and later use **consensus**.
* **Consensus edges:** Retain edges present in ≥ `prop_present` of repeats; summarise weights.
* **Multilayer assembly:** Wrap one igraph per layer in `multinet::ml.network`.
* **Community detection:** Build a **supra‑graph** across chosen layers; parallel edges are aggregated by **count** or **sum**; support **Louvain**, **Infomap**, **k‑clique union**, and **ABACUS** (when available).
* **Overlaps & metrics:** Jaccard for actors/edges; standard centralities; path‑based summaries with guards for disconnected graphs.
* **Dynamics:** `networkDynamic`/`ndtv` for animations, timelines, and filmstrips.&#x20;

---

## 📊 What the results tell you

* **Consensus edges:** Robust wiring; high `prop_present` ≈ stable relationships.
* **Communities:** Size → breadth of co‑variation; span across layers → cross‑layer stability vs. specificity.
* **Layer metrics:** Density, components, APL/diameter contextualise module patterns.
* **Overlaps:** Actor vs. edge Jaccard disentangle participation from wiring similarity.
* **Feature summaries:** Which classes dominate which modules; degree patterns show class‑specific centrality.&#x20;

---

## 🧩 Function map (20 updated functions)

### Build & assemble

| Function                   | Purpose                                                         |
| -------------------------- | --------------------------------------------------------------- |
| `buildAdjacency()`         | Create per‑group adjacency matrices (with optional resampling). |
| `edgesFromAdjacency()`     | Convert matrices/lists to tidy edge tables.                     |
| `consensusEdges()`         | Keep edges robust across repeats; summarise weights.            |
| `build_multiNet()`         | Assemble `multinet::ml.network` from per‑layer edges.           |
| `add_network_attributes()` | Attach node/edge attributes (robust ID matching).               |

### Community detection & reporting

| Function      | Purpose                                                              |
| ------------- | -------------------------------------------------------------------- |
| `detectCom()` | Communities on a supra‑graph (Louvain/Infomap/k‑clique/ABACUS).      |
| `plotCom()`   | Multilayer community plot (+ save long `(actor, layer, cid)` table). |

### Structure, overlaps, sampling

| Function                  | Purpose                                                 |
| ------------------------- | ------------------------------------------------------- |
| `layer_metrics()`         | Node‑level & layer‑level metrics; saves a run folder.   |
| `analyze_actor_overlap()` | Jaccard overlap of actor sets between layers + heatmap. |
| `analyze_edge_overlap()`  | Jaccard overlap of edge sets between layers + heatmap.  |
| `com_vs_samples()`        | Relate #communities per layer to #samples per group.    |

### Feature‑aware summaries

| Function           | Purpose                                                             |
| ------------------ | ------------------------------------------------------------------- |
| `annotateCom()`    | Add a feature attribute (e.g., `GeneType`) to community rows.       |
| `sumComFeat()`     | Summaries by actor, by community, and counts per community × layer. |
| `get_FeatureDeg()` | Degrees for any feature list across layers (wide/long).             |

### Dynamics & visualisation

| Function                 | Purpose                                                     |
| ------------------------ | ----------------------------------------------------------- |
| `animate_multiNet()`     | Interactive HTML animation (one slice per layer; ndtv).     |
| `filmstrip_multiNet()`   | Static grid (one panel per layer) with stable layout.       |
| `animate_multiNet_mp4()` | MP4/GIF via `gganimate`/`gifski` (fallback; optional).      |
| `plotActivityTimeline()` | Vertex/edge “lifespan” timeline (Gantt‑like).               |
| `grid_layer_diffs()`     | Multi‑panel “diffs” grid across layers (publication‑ready). |
| `plot_layer_diff()`      | Side‑by‑side layer comparison with in/out/keep edge sets.   |

> These tables reflect the updated pipeline and nomenclature adopted across the codebase.&#x20;

---

## 🖼️ Dynamics: quick recipes

**Interactive animation (HTML)**

```r
animate_multiNet(
  net,
  communities   = comm,             # data.frame: actor, layer, com/cid
  layout        = "kamadakawai",
  slice.par     = list(start = 0, interval = 1, aggregate.dur = 1),
  displaylabels = TRUE,
  vertex.cex    = 0.9,
  seed          = 1
)
# => omicsDNA_results/multiNet_animation_<timestamp>.html
```

**Filmstrip (PNG/PDF)**

```r
filmstrip_multiNet(
  net,
  communities = comm,
  layout      = "kamadakawai",
  ncol        = 4,                   # e.g., 8 layers -> 2x4
  format      = "png",
  seed        = 1
)
# => omicsDNA_results/multiNet_filmstrip_<timestamp>.png
```

**Activity timeline (lifespans)**

```r
plotActivityTimeline(
  net,
  type     = "vertex",               # or "edge"
  at       = NULL,                   # auto: one slice per layer
  label.cex= 0.8
)
```

---

## ✅ Best practices

* **Reproducibility:** set a seed; leave the default **timestamped filenames**; results land in `getOption("mlnet.results_dir","omicsDNA_results")`.
* **Thresholding:** start around `corr_threshold = 0.6–0.7`, `pval_cutoff = 0.05`; use `consensusEdges()` to guard against instability.
* **Sampling:** check `com_vs_samples()` to see whether community counts track sample sizes.
* **Interpretation:** use **edge overlap** for wiring, **actor overlap** for participation; relate to `layer_metrics()` for structure.&#x20;

---

## 🙋 Feedback & contributions

Issues and pull requests are welcome. Please open them on GitHub and include:

* a short reproducible example (if applicable),
* your `sessionInfo()`,
* and the function names you used.

---

## 📄 License

MIT © VafaeeLab

---

### Citation

If you use `omicsDNA` in published work, please cite the package and describe the multilayer workflow briefly (build → consensus → assemble → detect → quantify → visualise) so readers can reproduce your pipeline.&#x20;

---

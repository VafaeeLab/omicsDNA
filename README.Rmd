
---

```yaml
---
output: github_document
---
```

# omicsDNA <img src="man/figures/logo.png" align="right" height="120"/>

**omicsDNA** is an R package for constructing, stabilizing, and interpreting **multilayer networks** from grouped data (e.g., gene expression by age group or condition). It supports **robust edge selection**, **community detection on a supra‑graph**, **layer‑wise metrics**, **overlap analysis**, **feature‑aware summaries**, and **dynamic visualizations** (filmstrips, HTML animations, GIFs, timelines).

> Although examples refer to omics, the design is **domain‑agnostic**: “features” can be genes, proteins, survey items, sensors, or behaviors.

---

## 🔍 What does omicsDNA do?

* Builds **per‑group adjacency matrices** and converts them to **edge lists**.
* Stabilizes edges via **resampling + consensus**.
* Assembles a **multilayer** network object and **annotates** nodes/edges.
* Detects **communities** (Louvain/Infomap/clique; ABACUS if available).
* Quantifies structure with **layer‑wise node & graph metrics**.
* Compares layers via **actor/edge overlaps (Jaccard)** and **difference plots**.
* Summarizes **feature classes** (e.g., lncRNA/TF) across modules and layers.
* Visualizes dynamics as **filmstrips**, **HTML animations**, **GIFs**, and **activity timelines**.

---

## 📦 Installation

```r
# Install devtools if needed
install.packages("devtools")

# Install from GitHub (main repo)
devtools::install_github("VafaeeLab/omicsDNA")
# or from your fork (example)
# devtools::install_github("mehranpiran/omicsDNA")
```

> Optional visualisation deps (used defensively):
>
> * **ndtv**, **networkDynamic** for HTML animations, filmstrips, and timelines
> * **gganimate**, **gifski** for GIF export (no admin rights needed)
>   The package falls back gracefully if these are absent.

---

## 🚀 Quick Start (minimal)

```r
library(omicsDNA)

# 0) Set a results folder (timestamps are added to outputs)
options(mlnet.results_dir = "omicsDNA_results")

# 1) Build per-group networks from a feature x sample matrix
adj <- buildAdjacency(
  X           = expr_matrix,           # features x samples
  meta        = sample_metadata,       # includes grouping variable
  group_col   = "Group",               # e.g., age group / condition
  feature_ids = rownames(expr_matrix), # or a subset
  method      = "spearman",
  resample    = TRUE, n_repeats = 50,  # optional stability
  min_per_grp = 3
)

# 2) Convert to edge lists & keep robust edges across resamples
edges   <- edgesFromAdjacency(adj)
cons    <- consensusEdges(edges, prop_present = 0.7)  # retained if present in ≥70% repeats

# 3) Build multilayer object and annotate nodes
net <- build_multiNet(cons)
net <- add_network_attributes(net, nodesMetadata = genes_info,
                              featureID_col = "GeneName",
                              nodeAttrCols  = c("GeneType"))

# 4) Detect communities on a supra‑graph and visualize them
comm <- detectCom(net, method = "louvain", edgeWeight = "count",
                  min.actors = 15, min.layers = 2, seed = 1)
plotCom(net, comm, layout = "multiforce", gravity = 0.3,
        save_plot = TRUE, save_df = TRUE)

# 5) Metrics & overlaps
LM   <- layer_metrics(net)                       # run-folder with summary + node metrics
Ajac <- analyze_actor_overlap(net, reorder=TRUE) # heatmap + matrix
Ejac <- analyze_edge_overlap(net,  reorder=TRUE) # heatmap + matrix

# 6) Feature-aware summaries
comm_annot <- annotateCom(comm, genes_info, featureID_col = "GeneName", attribute = "GeneType")
summ_lnc   <- sumComFeat(comm_annot, feature_type = "lncRNA")
deg_lnc    <- get_FeatureDeg(net, featureList = with(comm_annot, actor[GeneType == "lncRNA"]))

# 7) Dynamics
animate_multiNet(net, communities = comm)           # HTML animation
filmstrip_multiNet(net, communities = comm)         # static grid
plotActivityTimeline(net, type = "vertex")          # timeline (lifespan)
# animate_multiNet_gifski(net, communities = comm)  # GIF (if gganimate/gifski available)
```

---

## 🧭 Typical questions you can answer

* Which feature pairs are **consistently** co‑varying within each group?
* Which **communities** (modules) are robust across layers?
* How **similar** are layers by actors and edges?
* Which actors are **central** in a given layer?
* Do layers with **more samples** have more communities?
* How do specific **feature classes** (e.g., lncRNA, TF) distribute across modules/layers?
* Where do layers **differ most** (edges gained/lost)?

---

## 🧩 Pipeline overview

1. **Build & stabilize** → `buildAdjacency()` → `edgesFromAdjacency()` → `consensusEdges()`
2. **Assemble & annotate** → `build_multiNet()` → `add_network_attributes()`
3. **Communities** → `detectCom()` → `plotCom()`
4. **Metrics & overlaps** → `layer_metrics()` → `analyze_actor_overlap()` / `analyze_edge_overlap()`
5. **Feature summaries** → `annotateCom()` → `sumComFeat()` → `get_FeatureDeg()`
6. **Dynamics** → `animate_multiNet()` / `filmstrip_multiNet()` / `plotActivityTimeline()` / `animate_multiNet_gifski()`
7. **Layer differences** → `grid_layer_diffs()` / `plot_layer_diff()`

---

## 📚 Function reference (grouped)

### 1) Build & Stabilize

| Function               | Purpose                                                              |
| ---------------------- | -------------------------------------------------------------------- |
| `buildAdjacency()`     | Compute per‑group correlation networks (optionally with resampling). |
| `edgesFromAdjacency()` | Convert matrices (incl. nested lists) to tidy edge tables.           |
| `consensusEdges()`     | Keep edges robust to resampling (by proportion or count).            |

### 2) Assemble & Annotate

| Function                   | Purpose                                              |
| -------------------------- | ---------------------------------------------------- |
| `build_multiNet()`         | Wrap per‑layer edges into a multilayer object.       |
| `add_network_attributes()` | Attach node/edge attributes with robust ID matching. |

### 3) Communities

| Function      | Purpose                                                                                 |
| ------------- | --------------------------------------------------------------------------------------- |
| `detectCom()` | Community detection on a **supra‑graph** (Louvain/Infomap/clique; ABACUS if available). |
| `plotCom()`   | Plot multilayer communities; also saves `(actor, layer, cid)` table.                    |

### 4) Metrics & Overlaps

| Function                  | Purpose                                                            |
| ------------------------- | ------------------------------------------------------------------ |
| `layer_metrics()`         | Node‑ and layer‑level metrics; saves run folder with CSV/XLSX.     |
| `analyze_actor_overlap()` | Jaccard overlap of actor sets across layers (heatmap + matrix).    |
| `analyze_edge_overlap()`  | Jaccard overlap of edges across layers (heatmap + matrix).         |
| `com_vs_samples()`        | Relate per‑layer community counts to sample sizes (scatter + OLS). |

### 5) Feature‑Focused Analyses

| Function           | Purpose                                                                                |
| ------------------ | -------------------------------------------------------------------------------------- |
| `annotateCom()`    | Merge actor attributes (e.g., GeneType) into community tables.                         |
| `sumComFeat()`     | Summaries restricted to feature types (by actor / by community / counts by CID×layer). |
| `get_FeatureDeg()` | Degrees across layers for any feature list (actors added as isolates if absent).       |

### 6) Dynamics: Quick Recipes

| Function                    | Purpose                                                                      |
| --------------------------- | ---------------------------------------------------------------------------- |
| `animate_multiNet()`        | Interactive **HTML** animation of the multilayer network over layers (ndtv). |
| `filmstrip_multiNet()`      | Static **grid** (“filmstrip”) with one panel per layer (PNG/PDF).            |
| `plotActivityTimeline()`    | **Timeline** (Gantt‑style) of vertex/edge lifespans across layers.           |
| `animate_multiNet_gifski()` | **GIF** export of layer‑by‑layer networks (gganimate + gifski).              |
| `grid_layer_diffs()`        | Matrix of **pairwise layer difference** plots (big‑picture view).            |
| `plot_layer_diff()`         | Focused plot of **edges gained/lost** between **two** layers.                |

---

## 🧠 Methods (concise)

* **Per‑group networks**: correlations (Spearman/Pearson), filters on |r| and *p*.
* **Consensus**: retain edges appearing in ≥ proportion (or count) of resamples; summarize weight.
* **Multilayer object**: one igraph per layer, wrapped into a `multinet` container.
* **Communities**: build a **single undirected supra‑graph**; aggregate duplicates by **count** or **sum**, run Louvain/Infomap or a **k‑clique union**; ABACUS if available.
* **Metrics**: degree, betweenness, eigenvector, closeness, PageRank, clustering coefficient, coreness; layer density, components, average path length, diameter; degree distribution stats (mean, SD, median, min/max, CV, skew).
* **Overlaps**: Jaccard for actors/edges.
* **Dynamics**: layouts stabilized across slices for filmstrip/animation; activity timelines show “lifespans”.

---

## 🧪 Interpreting outputs

* **Consensus edges**: high presence proportion ⇒ robust wiring; good candidates for interpretation.
* **Community size & span**: large CIDs suggest broad co‑variation; multi‑layer span hints at shared programs; layer‑specific CIDs point to condition‑specific processes.
* **Layer metrics**: density/components/APL/diameter contextualize how “connected” each layer is.
* **Overlaps**: actor overlap ≈ participation similarity; edge overlap ≈ wiring similarity.
* **Feature summaries**: reveal which feature classes dominate modules or act as hubs.
* **Dynamics**: filmstrips/animations expose periods of re‑wiring, module splitting/merging; timelines separate **stable** vs **transient** elements.

---

## 🧩 Worked mini‑example (pseudo‑code)

```r
# Build → stabilize
adj   <- buildAdjacency(X=expr_matrix, meta=sample_metadata, group_col="Group",
                        feature_ids=rownames(expr_matrix), method="spearman",
                        resample=TRUE, n_repeats=50)
edges <- edgesFromAdjacency(adj)
cons  <- consensusEdges(edges, prop_present=0.7)

# Assemble → annotate
net <- build_multiNet(cons)
net <- add_network_attributes(net, nodesMetadata=genes_info,
                              featureID_col="GeneName", nodeAttrCols="GeneType")

# Communities → viz
comm <- detectCom(net, method="louvain", edgeWeight="count",
                  min.actors=15, min.layers=2, seed=1)
plotCom(net, comm, layout="multiforce", gravity=0.3)

# Metrics & overlaps
LM   <- layer_metrics(net)
Ajac <- analyze_actor_overlap(net, reorder=TRUE)
Ejac <- analyze_edge_overlap(net,  reorder=TRUE)

# Feature‑aware summaries
comm_annot <- annotateCom(comm, genes_info, featureID_col="GeneName", attribute="GeneType")
summ_TF    <- sumComFeat(comm_annot, feature_type="TF")
deg_TF     <- get_FeatureDeg(net, featureList=with(comm_annot, actor[GeneType=="TF"]))

# Dynamics
animate_multiNet(net, communities=comm, seed=1)
filmstrip_multiNet(net, communities=comm, ncol=4, seed=1)
plotActivityTimeline(net, type="edge")
# animate_multiNet_gifski(net, communities=comm, fps=2) # if gganimate/gifski available
```

---

## 🗂️ Outputs & saving conventions

Most functions write into `getOption("mlnet.results_dir", "omicsDNA_results")` with **timestamped** names.
Some create a **run folder** (e.g., `layer_metrics_<YYYY-mm-dd_HHMMSS>/`), others write a single RDS/CSV/PNG/PDF/HTML/GIF.
Return objects often carry `"file"` or `"files"` attributes with absolute paths.

Set once per session or in `.Rprofile`:

```r
options(mlnet.results_dir = "omicsDNA_results")
```

---

## 💡 Tips & cautions

* **Correlation ≠ causation**; interpret communities as co‑variation modules.
* **Multiple testing**: resampling + consensus helps empirically control false positives.
* **Sample imbalance**: check `com_vs_samples()` to guard against density artifacts.
* **ID hygiene**: joins use robust normalization (strip version, trim, lower‑case).
* **Reproducibility**: set `seed` for `detectCom()`, dynamic plots, and any resampling.

---

## 🔁 Old → new function names (compatibility)

| Old name                         | New name                    |
| -------------------------------- | --------------------------- |
| `generate_adjacency_matrices()`  | `buildAdjacency()`          |
| `build_DNA()`                    | `build_multiNet()`          |
| `add_DNA_attributes()` / similar | `add_network_attributes()`  |
| `plot_filtered_ABACUS_network()` | `detectCom()` + `plotCom()` |

> Older wrappers continue to work where possible, but the **new API** above is recommended.

---

## 🤝 Contributing

Issues and pull requests are welcome!
Please include a minimal reproducible example and session info when reporting bugs.

---

## 📄 License

MIT © VafaeeLab

---

## 📣 Citation

If you use **omicsDNA** in your work, please cite this repository and the specific methods you rely on (e.g., Louvain/Infomap; ndtv/networkDynamic; gganimate/gifski).

```text
VafaeeLab. omicsDNA: Multilayer Network Analysis for Omics (and Beyond).
https://github.com/VafaeeLab/omicsDNA
```

---

### Session info (optional)

```r
sessionInfo()
```

---


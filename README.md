
# omicsDNA <img src="man/figures/logo.png" align="right" height="120"/>

`omicsDNA` is an R package designed for multilayer network analysis of
gene expression data across biological conditions, such as age groups or
disease stages. The package provides a flexible framework to construct,
analyze, and interpret gene networks by integrating community detection,
topological metrics, and long non-coding RNA (lncRNA) profiling.

------------------------------------------------------------------------

## 🔍 What Does omicsDNA Do?

-   Constructs **adjacency matrices** from gene expression profiles.
-   Builds **multilayer gene networks** to reflect condition-specific
    interactions.
-   Applies the **ABACUS** algorithm for community detection across
    layers.
-   Computes **overlap metrics** for actors (genes) and edges between
    layers.
-   Provides tools for **annotating and analyzing lncRNA involvement**
    in communities.
-   Supports extraction of **network metrics**, such as degree
    distribution, skewness, and coefficient of variation.

------------------------------------------------------------------------

## 📦 Installation

You can install the development version of `omicsDNA` using:

``` r
# Install devtools if you haven't already
install.packages("devtools")

# Install from a local directory
devtools::install("path/to/omicsDNA")

# Install via GitHub:
devtools::install_github("VafaeeLab/omicsDNA")
```

## Quick Start

``` r

library(omicsDNA)

# Step 1: Generate adjacency matrices from expression data
adj_list <- generate_adjacency_matrices(
  expression_mat = expr_matrix,
  metaData = sample_metadata,
  genesInfo = gene_annotation,
  DEGs = deg_table,
  num_repetitions = 10,
  num_samples = 30
)

# Step 2: Build multilayer network
net <- create_multilayer_network(adj_list)

# Step 3: Add metadata (e.g., gene type) as actor attributes
net <- add_network_attributes(net, genesInfo = gene_annotation)

# Step 4: Detect communities using ABACUS and plot
plot_filtered_ABACUS_network(net, layer_names = c("Young", "Old"))

# Step 5: Analyze lncRNA involvement in detected communities
annotated_communities <- merge_gene_type_with_communities(ABACUS_communities_net, actors_with_type)
lncRNA_summary <- analyze_lncRNA_communities(annotated_communities)
```

------------------------------------------------------------------------

## 🧪 Core Functionalities

| Function                             | Purpose                                            |
|--------------------------------------|----------------------------------------------------|
| `generate_adjacency_matrices()`      | Create condition-specific adjacency matrices       |
| `create_multilayer_network()`        | Build multilayer gene network from adjacency lists |
| `add_network_attributes()`           | Annotate actors with gene metadata                 |
| `plot_filtered_ABACUS_network()`     | Detect & visualize communities                     |
| `analyze_actor_overlap()`            | Assess gene overlap between layers                 |
| `analyze_edge_overlap()`             | Compare edge similarity across layers              |
| `compute_lncRNA_degree()`            | Calculate lncRNA connectivity per layer            |
| `merge_gene_type_with_communities()` | Join ABACUS output with gene types                 |
| `analyze_lncRNA_communities()`       | Summarize lncRNA distribution across communities   |

------------------------------------------------------------------------

## 🧬 Use Case

This package is ideal for researchers interested in:

-   Comparing gene co-expression networks across age, condition, or
    phenotype
-   Understanding condition-specific community structure in regulatory
    networks
-   Profiling lncRNA roles in systems biology using topological and
    community metrics

------------------------------------------------------------------------

## 📄 License

MIT License © VafaeeLab  
For academic or research use. Please cite appropriately.

------------------------------------------------------------------------

## 🤝 Feedback and Contributions

Contributions are welcome! Please submit issues or pull requests on
GitHub.

\`\`\`

------------------------------------------------------------------------

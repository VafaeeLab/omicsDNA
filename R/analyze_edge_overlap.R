#' Edge overlap across layers via Jaccard similarity (with heatmap):
#'
#' @description
#' Quantify how similar layers are in terms of their **undirected edge sets** by
#' computing pairwise **Jaccard indices**, and (optionally) visualise the result
#' as a heatmap.
#'
#' @details
#' **Computation**
#' - If supported, the function uses
#'   `multinet::layer_comparison_ml(..., method = "jaccard.edges")` and converts
#'   the result to a base R matrix.
#' - Otherwise, it **falls back** to a manual computation. For each layer, it
#'   retrieves the edge table via `multinet::edges_ml(net, layers = L)`, coerces
#'   to a data frame if needed, and heuristically picks the **first two character
#'   columns** as the endpoints. Pairs are canonicalised so that `(u, v)` equals
#'   `(min, max)`; edges are then represented as unique keys and the Jaccard
#'   index \eqn{|E_i \cap E_j| / |E_i \cup E_j|} is computed for each pair of
#'   layers. The diagonal is set to `1`. If both edge sets are empty, the
#'   corresponding entry is `NA`.
#'
#' **Layer selection and reordering**
#' - If `layers` is `NULL`, all layers from `multinet::layers_ml(net)` are used.
#'   If a subset is provided, it is intersected with available layers.
#' - When `reorder = TRUE`, the matrix is re‑indexed by hierarchical clustering
#'   on `1 - Jaccard` (average linkage). `NA` entries are temporarily treated as
#'   `0` for the dendrogram only.
#'
#' **Heatmap**
#' - The Jaccard matrix is converted to long form and plotted with
#'   `ggplot2::geom_tile()`. The fill scale spans `[0, 1]`, `NA`s render in light
#'   grey, aspect ratio is equal, and x‑axis labels are rotated for readability.
#'   Colours can be customised via `palette` (fed to `scale_fill_gradientn()`).
#'
#' @param net A multilayer network object compatible with `multinet`, i.e., one
#'   that responds to `multinet::layers_ml()` and `multinet::edges_ml()`.
#' @param layers Optional character vector naming the layers to include. Default:
#'   all layers returned by `multinet::layers_ml(net)`.
#' @param palette Character vector of **two or more** colours for the heatmap
#'   gradient (passed to `scale_fill_gradientn(colours = ...)`). Default
#'   `c("white", "steelblue")`.
#' @param reorder Logical; if `TRUE`, reorder rows/columns by hierarchical
#'   clustering on `1 - Jaccard`. Default `FALSE`.
#' @param print_plot Logical; if `TRUE`, print the heatmap to the active device.
#'   Default `TRUE`.
#'
#' @return A symmetric numeric **matrix** of Jaccard indices in `[0, 1]` with
#'   row/column names equal to the layer names. The returned object carries an
#'   attribute `"plot"` containing the `ggplot` heatmap object (invisibly
#'   printed when `print_plot = TRUE`).
#'
#' @section Notes and limitations
#' - In the fallback path, the endpoint identification is **heuristic** (first
#'   two character columns). If your edge table contains additional character
#'   columns (e.g., a textual layer label) before the endpoints, re‑order columns
#'   or pre‑process the table to avoid misidentification.
#' - Edges are treated as **undirected**; for directed analyses, build directed
#'   keys explicitly upstream.
#'
#' @examples
#' \dontrun{
#' # Compute and visualise edge overlap across all layers
#' E <- analyze_edge_overlap(net, reorder = TRUE)
#' attr(E, "plot")   # access the ggplot object programmatically
#' }
#'
#' @seealso \code{\link{analyze_actor_overlap}} for actor‑set overlap across layers.
#'
#' @importFrom multinet layers_ml edges_ml layer_comparison_ml
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_gradientn theme_minimal labs ggtitle element_text coord_equal
#' @importFrom stats hclust as.dist
#' @export
analyze_edge_overlap <- function(net,
                                 layers     = NULL,
                                 palette    = c("white", "steelblue"),
                                 reorder    = FALSE,
                                 print_plot = TRUE) {
  ok_layers <- try(multinet::layers_ml(net), silent = TRUE)
  if (inherits(ok_layers, "try-error")) {
    stop("`net` does not respond to multinet::layers_ml(); is it a valid multinet object?")
  }

  if (is.null(layers)) layers <- ok_layers else layers <- intersect(as.character(layers), as.character(ok_layers))
  if (!length(layers)) stop("No layers available.")

  mat <- try(
    as.matrix(multinet::layer_comparison_ml(
      net, layers = layers, method = "jaccard.edges", mode = "all", K = 0
    )),
    silent = TRUE
  )

  if (inherits(mat, "try-error")) {
    get_edge_keys <- function(ly) {
      ed <- multinet::edges_ml(net, layers = ly)
      if (is.null(ed)) return(character())
      if (!is.data.frame(ed)) ed <- try(as.data.frame(ed, stringsAsFactors = FALSE), silent = TRUE)
      if (inherits(ed, "try-error") || !nrow(ed)) return(character())
      char_cols <- which(vapply(ed, function(x) is.character(x) || is.factor(x), logical(1)))
      if (length(char_cols) < 2) return(character())
      a <- as.character(ed[[char_cols[1]]])
      b <- as.character(ed[[char_cols[2]]])
      u <- pmin(a, b); v <- pmax(a, b)
      unique(paste(u, v, sep = "\t"))
    }
    sets <- lapply(layers, get_edge_keys); names(sets) <- layers
    n <- length(layers)
    mat <- matrix(NA_real_, n, n, dimnames = list(layers, layers))
    for (i in seq_len(n)) {
      Ei <- sets[[i]]
      for (j in i:n) {
        Ej <- sets[[j]]
        u  <- union(Ei, Ej)
        val <- if (length(u)) length(intersect(Ei, Ej)) / length(u) else NA_real_
        mat[i, j] <- mat[j, i] <- val
      }
    }
  }

  diag(mat) <- 1

  if (reorder && nrow(mat) > 1L) {
    m2 <- mat; m2[is.na(m2)] <- 0
    ord <- stats::hclust(stats::as.dist(1 - m2), method = "average")$order
    mat <- mat[ord, ord, drop = FALSE]
  }

  df <- as.data.frame(base::as.table(mat), stringsAsFactors = FALSE)
  names(df) <- c("Layer1", "Layer2", "Jaccard")

  plt <- ggplot2::ggplot(df, ggplot2::aes(Layer1, Layer2, fill = Jaccard)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::scale_fill_gradientn(colours = palette, limits = c(0, 1), na.value = "grey90") +
    ggplot2::coord_equal() +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = "Layers", y = "Layers", fill = "Jaccard") +
    ggplot2::ggtitle("Edge overlap between layers") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))

  if (print_plot) print(plt)
  attr(mat, "plot") <- plt
  mat
}

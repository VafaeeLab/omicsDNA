#' Actor (gene) overlap across layers via Jaccard similarity (with heatmap)
#'
#' @description
#' Quantify how similar layers are in terms of their **actor (gene) sets** by
#' computing the pairwise **Jaccard index** for all layer pairs, and (optionally)
#' visualise the resulting similarity matrix as a heatmap.
#'
#' @details
#' **Computation**
#' - If available in your `multinet` build, the function uses
#'   `multinet::layer_comparison_ml(..., method = "jaccard.actors")` to obtain a
#'   Jaccard matrix directly and converts it to base R `matrix` form.
#' - Otherwise, it **falls back** to a manual computation: for each layer, it
#'   extracts the actor set via `multinet::actors_ml(net, layers = L)` (accepting
#'   either a character vector or a data frame with an `"actor"` column), then
#'   for each pair `(A, B)` computes
#'   \deqn{J(A,B) = |A \cap B| \, / \, |A \cup B|.}
#'   The diagonal is set to 1. If both sets are empty, the corresponding entry is
#'   left as `NA` (rendered in a pale colour on the heatmap).
#'
#' **Layer selection and reordering**
#' - If `layers` is `NULL`, all layers from `multinet::layers_ml(net)` are used.
#'   If a subset is supplied, it is intersected with available layers.
#' - When `reorder = TRUE`, the matrix is re‑indexed by hierarchical clustering
#'   on a dissimilarity of `1 - Jaccard` (average linkage). For the sole purpose
#'   of computing the dendrogram, `NA` entries are temporarily treated as `0`
#'   (i.e., no overlap).
#'
#' **Heatmap**
#' - A tidy long table is produced from the Jaccard matrix and drawn with
#'   `ggplot2::geom_tile()`. The fill scale spans `[0, 1]`, `NA` cells are shown
#'   in a light grey, and the plot is square (`coord_equal`). Colours can be
#'   customised via `palette` (passed to `scale_fill_gradientn()`).
#'
#' @param net A multilayer network object compatible with `multinet`, i.e., one
#'   that responds to `multinet::layers_ml()` and `multinet::actors_ml()`.
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
#' - The fallback path assumes `actors_ml()` returns either a character vector of
#'   actor IDs or a data frame with a column named `"actor"`.
#' - If both layers in a pair have zero actors (rare in practice), the Jaccard
#'   index is set to `NA` for that pair.
#'
#' @examples
#' \dontrun{
#' # Compute and visualise actor overlap across all layers
#' A <- analyze_actor_overlap(net, reorder = TRUE)
#' attr(A, "plot")   # access the ggplot object programmatically
#' }
#'
#' @seealso \code{\link{analyze_edge_overlap}} for edge‑set overlap across layers.
#'
#' @importFrom multinet layers_ml actors_ml layer_comparison_ml
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_gradientn theme_minimal labs ggtitle element_text coord_equal
#' @importFrom stats hclust as.dist
#' @export
analyze_actor_overlap <- function(net,
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

  # try built-in Jaccard first
  mat <- try(
    as.matrix(multinet::layer_comparison_ml(
      net, layers = layers, method = "jaccard.actors", mode = "all", K = 0
    )),
    silent = TRUE
  )

  # fallback if not supported
  if (inherits(mat, "try-error")) {
    get_actors <- function(ly) {
      a <- multinet::actors_ml(net, layers = ly)
      if (is.data.frame(a) && "actor" %in% names(a)) unique(as.character(a$actor))
      else if (is.character(a)) unique(a)
      else character()
    }
    sets <- lapply(layers, get_actors); names(sets) <- layers
    n <- length(layers)
    mat <- matrix(NA_real_, n, n, dimnames = list(layers, layers))
    for (i in seq_len(n)) {
      Ai <- sets[[i]]
      for (j in i:n) {
        Aj <- sets[[j]]
        u  <- union(Ai, Aj)
        val <- if (length(u)) length(intersect(Ai, Aj)) / length(u) else NA_real_
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
    ggplot2::ggtitle("Actors (genes) overlap between layers") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))

  if (print_plot) print(plt)
  attr(mat, "plot") <- plt
  mat
}

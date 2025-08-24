#' Relate community counts to sample sizes across layers (with a quick scatter):
#'
#' @description
#' For each layer, count the number of **unique communities** observed in a
#' community assignment table and relate this count to the **number of samples**
#' in that layer. The function returns a named integer vector of per‑layer
#' community counts and, for convenience, prints a scatter plot with a simple
#' OLS trend line (via `geom_smooth(method = "lm")`). The resulting `ggplot`
#' object is attached to the return value as an attribute (`"plot"`).
#'
#' @details
#' **Inputs expected**
#' - `communities` must contain a `layer` column and exactly one community
#'   identifier column among **`cid`** (e.g., ABACUS output) or **`com`**
#'   (e.g., `detectCom()` output). If `com` is present, it is internally
#'   copied to `cid` before counting.
#' - `num_samples` must be a **named numeric vector** whose names are layer
#'   identifiers and whose values are sample counts (non‑negative, finite).
#'
#' **Behaviour**
#' 1. Compute the number of **distinct** community IDs per layer.
#' 2. Strictly validate that the set of layer names in `communities` matches the
#'    names of `num_samples` (the function errors if they differ).
#' 3. Assemble a small data frame and draw a scatter plot of *samples* (y‑axis)
#'    versus *community count* (x‑axis) with a least‑squares trend line.
#'    The plot is printed to the active graphics device and also stored as
#'    an attribute on the return value.
#'
#' **Note on statistics**
#' The plot includes an OLS line for orientation but **does not** compute or
#' annotate model statistics (R², p‑value). If you wish to report these, extract
#' the underlying data from `attr(x, "plot")` and fit a model explicitly.
#'
#' @param communities A data frame with columns:
#'   - `layer`: layer identifier, and
#'   - one of `cid` **or** `com`: community membership per `(actor, layer)` row.
#'     (If `com` is supplied, it is copied into `cid` internally.)
#'   Extra columns are ignored.
#' @param num_samples A **named numeric vector** whose names are layer identifiers
#'   and values are sample counts (non‑negative, finite). The set of names must
#'   match the set of layers in `communities`.
#'
#' @return A **named integer vector** giving the number of unique communities
#'   per layer. The object carries one attribute:
#'   - `"plot"` — the `ggplot` object (scatter with an OLS line).
#'
#' @examples
#' \dontrun{
#' # Example using detectCom() output (has `com` and `layer`)
#' ns <- c(E1 = 46, E2 = 22, M1 = 16, M2 = 27, M3 = 23, M4 = 26, M5 = 18, M6 = 15)
#' counts <- com_vs_samples(comm, ns)   # prints the plot
#' attr(counts, "plot")                 # access the ggplot object
#'
#' # Example using ABACUS-like output (has `cid` and `layer`)
#' # counts2 <- com_vs_samples(ABACUS_communities, ns)
#' }
#'
#' @seealso \code{\link{detectCom}} for obtaining community assignments.
#'
#' @importFrom dplyr group_by summarise n_distinct
#' @importFrom rlang .data
#' @importFrom ggplot2 ggplot aes geom_point geom_smooth labs theme_minimal
#' @importFrom stats setNames
#' @export
com_vs_samples <- function(communities, num_samples) {
  stopifnot(is.data.frame(communities))

  # normalize to 'cid' if needed
  if (!"cid" %in% names(communities)) {
    if ("com" %in% names(communities)) {
      communities$cid <- communities$com
    } else {
      stop("Input must contain either a 'cid' (ABACUS) or 'com' (detectCom) column.")
    }
  }
  if (!"layer" %in% names(communities)) {
    stop("Input must contain a 'layer' column.")
  }

  # counts per layer
  layer_cid_counts <- communities |>
    dplyr::group_by(.data$layer) |>
    dplyr::summarise(unique_cid_count = dplyr::n_distinct(.data$cid), .groups = "drop")

  message("Number of distinct communities per layer:")
  layer_cid_vector <- stats::setNames(layer_cid_counts$unique_cid_count, layer_cid_counts$layer)
  print(layer_cid_vector)

  # validate num_samples
  if (is.null(names(num_samples))) {
    stop("`num_samples` must be a *named* numeric vector.")
  }
  if (!is.numeric(num_samples) || any(!is.finite(num_samples))) {
    stop("`num_samples` must be numeric and finite.")
  }
  if (any(num_samples < 0)) {
    stop("`num_samples` must be non-negative.")
  }

  # check name match (strict)
  layers_comm <- sort(names(layer_cid_vector))
  layers_samp <- sort(names(num_samples))
  if (!identical(layers_comm, layers_samp)) {
    missing_in_samples <- setdiff(layers_comm, layers_samp)
    missing_in_comm    <- setdiff(layers_samp, layers_comm)
    msg <- c()
    if (length(missing_in_samples)) msg <- c(msg, paste0("Missing in `num_samples`: ", paste(missing_in_samples, collapse=", ")))
    if (length(missing_in_comm))    msg <- c(msg, paste0("Missing in communities: ", paste(missing_in_comm, collapse=", ")))
    stop("Group names in `num_samples` do not match the layer names in the community data.\n", paste(msg, collapse="\n"))
  }

  df <- data.frame(
    cids    = as.numeric(layer_cid_vector),
    samples = as.numeric(num_samples[names(layer_cid_vector)]),
    layer   = names(layer_cid_vector),
    stringsAsFactors = FALSE
  )

  plt <- ggplot2::ggplot(df, ggplot2::aes(x = cids, y = samples)) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_smooth(method = "lm", se = FALSE) +
    ggplot2::labs(
      x = "Number of communities (per layer)",
      y = "Number of samples (per group)",
      title = "Communities vs. sample size across layers"
    ) +
    ggplot2::theme_minimal()

  print(plt)

  # return just the vector (as before); the plot is printed, but also
  # attached as an attribute for convenience
  attr(layer_cid_vector, "plot") <- plt
  invisible(layer_cid_vector)
}

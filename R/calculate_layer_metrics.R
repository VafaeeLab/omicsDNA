
#' Calculate Degree-Based Network Metrics Per Layer
#'
#' This function computes a set of degree-based topological metrics
#' for each layer of a multilayer network constructed using the `multinet` package.
#'
#' @param net A multilayer network object of class `ml.network`, typically
#'        created with functions like create_multilayer_network().
#'
#' @return A data frame in wide format with each row representing a layer
#'         and columns showing computed metrics:
#'         - `max.degree`: Maximum node degree in the layer.
#'         - `mean.degree`: Mean node degree.
#'         - `sd.degree`: Standard deviation of node degrees.
#'         - `skewness.degree`: Skewness of the degree distribution.
#'         - `CV.degree`: Coefficient of variation of degrees.
#'
#' @details
#' For each layer, the function uses `layer_summary_ml()` from the `multinet` package
#' to compute the metrics. If a metric fails to compute, `NA` is assigned.
#'
#' @examples
#' \dontrun{
#' metrics_df <- calculate_layer_metrics(net)
#' }
#'
#' @importFrom multinet layers_ml layer_summary_ml
#' @importFrom tidyr pivot_wider
#' @export
calculate_layer_metrics <- function(net) {
  message("Note: Computing network metrics for each layer...")

  metrics <- c("max.degree", "mean.degree", "sd.degree", "skewness.degree", "CV.degree")
  results <- data.frame(Layer = character(), Metric = character(), Value = numeric(), stringsAsFactors = FALSE)

  for (layer in multinet::layers_ml(net)) {
    for (metric in metrics) {
      tryCatch({
        value <- multinet::layer_summary_ml(net, layer = layer, method = metric, mode = "all")
        results <- rbind(results, data.frame(Layer = layer, Metric = metric, Value = value))
      }, error = function(e) {
        warning(sprintf("Error in layer '%s' for metric '%s': %s", layer, metric, e$message))
        results <- rbind(results, data.frame(Layer = layer, Metric = metric, Value = NA))
      })
    }
  }

  results_wide <- tidyr::pivot_wider(results, names_from = Metric, values_from = Value)
  print(results_wide)
  return(results_wide)
}


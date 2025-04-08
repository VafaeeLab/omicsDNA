
#' Compute Connectivity (Degree) of lncRNAs in Each Layer
#'
#' This function calculates the degree (number of connections) for each lncRNA
#' in every layer of a multilayer network.
#'
#' @param net A multilayer network object of class `ml.network`, typically created using
#'        `create_multilayer_network()`.
#' @param genes_list A character vector of lncRNA gene names to analyze.
#' @param layers A character vector of layer names.
#' @param mode Degree mode for multilayer network (default = "all").
#' @param output_file File name to save the wide-format result table. Default = `"lncRNAs-degrees_bylayer.csv"`.
#'
#' @return A data frame in wide format with layers as rows and lncRNAs as columns.
#'
#' @examples
#' \dontrun{
#' lncRNAs <- c("AC011483.1", "DSG1-AS1", "TMEM99")
#' layers <- layers_ml(net)
#' compute_lncRNA_degree(net, genes_list = lncRNAs, layers = layers)
#' }
#'
#' @importFrom multinet degree_ml
#' @importFrom tidyr pivot_wider
#' @importFrom utils write.csv
#' @export
compute_lncRNA_degree <- function(net, genes_list, layers, mode = "all", output_file = "lncRNAs-degrees_bylayer.csv") {
  degree_results <- data.frame(Gene = character(), Layer = character(), Degree = numeric(), stringsAsFactors = FALSE)

  for (gene in genes_list) {
    for (layer in layers) {
      tryCatch({
        degree <- multinet::degree_ml(net, actors = gene, layers = layer, mode = mode)
        degree_results <- rbind(degree_results, data.frame(Gene = gene, Layer = layer, Degree = degree))
      }, error = function(e) {
        warning(sprintf("Warning: Error processing actor '%s' in layer '%s': %s", gene, layer, e$message))
        degree_results <- rbind(degree_results, data.frame(Gene = gene, Layer = layer, Degree = NA))
      })
    }
  }

  results_wide <- tidyr::pivot_wider(degree_results, names_from = Gene, values_from = Degree)

  message("Note: Saved lncRNA degree matrix to: ", output_file)
  utils::write.csv(results_wide, file = output_file, row.names = FALSE)

  return(results_wide)
}

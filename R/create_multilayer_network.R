
#' Create a Multi-layer Network from Selected Edges
#'
#' Constructs a multilayer network using the `multinet` package by converting
#' selected edge lists into individual layer graphs, assigning node-level
#' attributes (e.g., gene types), and adding each layer to the network.
#'
#' @param selected_edges A named list of data frames, where each element represents
#'        the edge list of a specific group (e.g., age group or condition).
#'        Each data frame must contain columns: `from`, `to`, and `weight`.
#' @param genes_info A data frame containing gene annotation information.
#'        It must contain two columns: `GeneName` and `GeneType`.
#'
#' @return An object of class `ml.network` representing the constructed multilayer network.
#'
#' @details
#' - Converts each edge list to an `igraph` object.
#' - Assigns `GeneType` as a node attribute.
#' - Constructs a multilayer network using `multinet`.
#' - Prints the number of genes per layer and returns the multilayer network.
#'
#' @examples
#' \dontrun{
#' net <- create_multilayer_network(selected_edges, genes_info)
#' summary(net)
#' }
#'
#' @importFrom multinet ml_empty add_igraph_layer_ml actors_ml
#' @importFrom igraph graph_from_data_frame V
#' @importFrom dplyr left_join mutate
#' @export
create_multilayer_network <- function(selected_edges, genes_info) {
  graphs_list <- lapply(selected_edges, function(edge_list) {
    igraph::graph_from_data_frame(edge_list, directed = FALSE)
  })

  genes_info <- dplyr::mutate(genes_info, GeneName = as.character(GeneName))

  graphs_list <- lapply(graphs_list, function(g) {
    node_names <- igraph::V(g)$name
    node_data <- dplyr::left_join(data.frame(GeneName = node_names), genes_info, by = "GeneName")
    igraph::V(g)$GeneType <- node_data$GeneType
    return(g)
  })

  layer_names <- sort(names(graphs_list))
  net <- multinet::ml_empty()

  for (name in layer_names) {
    multinet::add_igraph_layer_ml(net, graphs_list[[name]], layer = name)
  }

  num_genes <- sapply(layer_names, function(layer) {
    length(unlist(multinet::actors_ml(net, layers = layer)))
  })

  message("Note: Number of genes per layer:")
  print(data.frame(Layer = layer_names, ActorCount = num_genes))
  return(net)
}



#' Add Node-Level Attributes to a Multilayer Network
#'
#' This function adds gene-level metadata (e.g., GeneType) to actors (nodes)
#' in a multilayer network constructed using the `multinet` package.
#'
#' @param net A multilayer network object of class `ml.network`, created using
#'        functions such as create_multilayer_network().
#' @param genes_info A data frame containing gene metadata. It must include at least
#'        the columns `GeneName` and `GeneType`.
#'
#' @return This function updates the `net` object in-place by assigning the `GeneType`
#'         attribute to all actors. It does not return a value.
#'
#' @details
#' - Extracts actors from the network
#' - Joins metadata
#' - Fills missing values with "Unknown"
#' - Assigns `GeneType` as an actor attribute
#'
#' @examples
#' \dontrun{
#' add_network_attributes(net, genes_info)
#' }
#'
#' @importFrom multinet actors_ml add_attributes_ml set_values_ml get_values_ml
#' @importFrom dplyr left_join mutate
#' @export
add_network_attributes <- function(net, genes_info) {
  actors_in_net <- multinet::actors_ml(net)
  actors_df <- data.frame(GeneName = actors_in_net$actor, stringsAsFactors = FALSE)
  genes_info <- dplyr::mutate(genes_info, GeneName = as.character(GeneName))

  actors_df <- dplyr::left_join(actors_df, genes_info, by = "GeneName")

  missing_types <- sum(is.na(actors_df$GeneType))
  if (missing_types > 0) {
    warning(paste(missing_types, "actors do not have a GeneType. Assigning 'Unknown'."))
    actors_df$GeneType[is.na(actors_df$GeneType)] <- "Unknown"
  } else {
    message("All actors have a GeneType.")
  }

  multinet::add_attributes_ml(net, attributes = "GeneType", target = "actor", type = "string")
  multinet::set_values_ml(net, "GeneType", actors = actors_df$GeneName, values = actors_df$GeneType)

  message("Note: GeneType assignment complete. Preview:")
  print(utils::head(multinet::get_values_ml(net, attribute = "GeneType", actors = actors_in_net)$GeneType))
}

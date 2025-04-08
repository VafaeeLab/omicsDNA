
#' Extract GeneType Attributes from a Multilayer Network
#'
#' This function retrieves the `GeneType` attribute for all actors (genes)
#' present in a multilayer network object created with the `multinet` package.
#'
#' @param net A multilayer network object of class `ml.network`, typically
#'        created using `create_multilayer_network()` and enriched using
#'        `add_network_attributes()`.
#'
#' @return A data frame with two columns:
#' \describe{
#'   \item{GeneName}{The name of each actor (gene) in the network.}
#'   \item{GeneType}{The gene type assigned to each actor.}
#' }
#'
#' @details
#' This function uses `actors_ml()` and `get_values_ml()` to extract actor names
#' and their assigned GeneType attributes.
#'
#' @examples
#' \dontrun{
#' actors_with_type <- extract_gene_types_from_net(net)
#' }
#'
#' @importFrom multinet actors_ml get_values_ml
#' @export
extract_gene_types_from_net <- function(net) {
  actors <- multinet::actors_ml(net)
  gene_types <- multinet::get_values_ml(net, attribute = "GeneType", actors = actors)

  actors_with_type <- data.frame(
    GeneName = actors$actor,
    GeneType = gene_types,
    stringsAsFactors = FALSE
  )

  message("Note: Extracted GeneType information for network actors:")
  print(head(actors_with_type))

  return(actors_with_type)
}

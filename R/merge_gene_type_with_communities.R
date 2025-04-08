
#' Merge GeneType Annotations with ABACUS Community Assignments
#'
#' This function joins gene-level metadata (e.g., GeneType) to the output of
#' community detection performed using the ABACUS algorithm.
#'
#' @param ABACUS_communities_net A data frame produced by `abacus_ml()`, typically
#'        containing columns like `actor`, `cid`, and `layer`.
#' @param actors_with_type A data frame with columns `GeneName` and `GeneType`,
#'        typically created by `extract_gene_types_from_net()`.
#' @param output_prefix A character prefix used for naming the saved `.csv` file.
#'        Default is `"communities_with_GeneType"`.
#'
#' @return A merged data frame with columns: `actor`, `cid`, `layer`, and `GeneType`.
#'         This is also saved as a CSV file.
#'
#' @details
#' This is useful for downstream biological interpretation of network modules.
#'
#' @examples
#' \dontrun{
#' communities <- merge_gene_type_with_communities(ABACUS_communities_net, actors_with_type)
#' }
#'
#' @importFrom dplyr left_join
#' @importFrom utils write.csv
#' @export
merge_gene_type_with_communities <- function(ABACUS_communities_net, actors_with_type, output_prefix = "communities_with_GeneType") {
  colnames(actors_with_type)[colnames(actors_with_type) == "actor"] <- "GeneName"

  communities <- dplyr::left_join(ABACUS_communities_net, actors_with_type, by = c("actor" = "GeneName"))

  output_file <- paste0(output_prefix, "_", Sys.Date(), ".csv")
  utils::write.csv(communities, file = output_file, row.names = FALSE)

  message("Annotated community table saved to: ", output_file)
  return(communities)
}

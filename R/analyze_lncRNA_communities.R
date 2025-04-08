
#' Analyze and Summarize lncRNA Membership in Communities
#'
#' This function explores how long non-coding RNAs (lncRNAs) participate in network communities
#' detected by ABACUS by summarizing their community IDs (CIDs) and the layers in which they occur.
#' It also groups lncRNAs by community.
#'
#' @param communities A data frame of community assignments from ABACUS, typically merged with
#'        gene metadata using `merge_gene_type_with_communities()`. Must contain `actor`, `cid`,
#'        `layer`, and `GeneType` columns.
#' @param lncRNA_output_file CSV file name to save lncRNA-wise summaries. Default: `"lncRNAs_bycid_bylayer.csv"`.
#' @param cid_output_file CSV file name to save CID-wise summaries. Default: `"lncRNAs_grouped-bycid.csv"`.
#'
#' @return A list with two data frames:
#' \describe{
#'   \item{lncRNA_summary}{Summary of CIDs and layers per lncRNA.}
#'   \item{cid_summary}{Summary of lncRNAs grouped by community.}
#' }
#'
#' @examples
#' \dontrun{
#' results <- analyze_lncRNA_communities(communities)
#' }
#'
#' @importFrom dplyr filter group_by summarise
#' @importFrom utils write.csv
#' @export
analyze_lncRNA_communities <- function(communities,
                                       lncRNA_output_file = "lncRNAs_bycid_bylayer.csv",
                                       cid_output_file = "lncRNAs_grouped-bycid.csv") {
  # lncRNA-wise summary
  lncRNA_summary <- communities %>%
    dplyr::filter(GeneType == "lncRNA") %>%
    dplyr::group_by(actor) %>%
    dplyr::summarise(
      CIDs = paste(unique(cid), collapse = ", "),
      Layers = paste(sort(unique(layer)), collapse = ", "),
      .groups = "drop"
    )

  message("Note: Saved lncRNA-to-community summary to: ", lncRNA_output_file)
  utils::write.csv(lncRNA_summary, file = lncRNA_output_file, row.names = FALSE)

  # CID-wise summary
  cid_summary <- communities %>%
    dplyr::filter(GeneType == "lncRNA") %>%
    dplyr::group_by(cid) %>%
    dplyr::summarise(
      lncRNAs = paste(unique(actor), collapse = ", "),
      Layers = paste(unique(layer), collapse = ", "),
      .groups = "drop"
    )

  message("Note: Saved community-to-lncRNA summary to: ", cid_output_file)
  utils::write.csv(cid_summary, file = cid_output_file, row.names = FALSE)

  return(list(lncRNA_summary = lncRNA_summary, cid_summary = cid_summary))
}


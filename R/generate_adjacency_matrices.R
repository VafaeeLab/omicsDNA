
#' Generate Adjacency Matrices Based on Gene Expression and Age Groups
#'
#' This function creates weighted adjacency matrices for each group
#' across multiple repetitions by computing pairwise gene-gene correlations
#' from gene expression data. It filters out low-correlation values based on
#' a specified threshold.
#'
#' @param threshold Numeric. Correlation threshold below which values are set to 0. Default is 0.7.
#' @param num_repetitions Integer. Number of repetitions for sampling and computing adjacency matrices. Default is 50.
#' @param num_samples Integer. Number of samples to randomly select per group for each repetition. Default is 10.
#' @param correlation Character. Type of correlation to compute: `"spearman"` or `"pearson"`. Default is `"spearman"`.
#' @param expression_mat A matrix or data frame of normalized gene expression data. Genes are in rows, samples in columns.
#' @param metaData A data frame containing metadata with at least two columns: `"Groups"` and `"Samples"` (sample names).
#' @param genesInfo A data frame containing information about each gene, with at least two columns; `"GeneName"` and `"GeneType"`.
#' @param DEGs A data frame of differentially expressed genes with at least a `"DEgenes"` column.
#'
#' @return A list of length `num_repetitions`, where each element is a list of adjacency matrices (one per age group).
#' Each adjacency matrix is a weighted, symmetric matrix of gene-gene correlations after thresholding.
#'
#' @examples
#' \dontrun{
#' adjacency_matrices <- generate_adjacency_matrices(
#'   expression_mat = my_expression_matrix,
#'   metaData = my_meta,
#'   genesInfo = my_genes_info,
#'   DEGs = my_deg_list,
#'   threshold = 0.7,
#'   num_repetitions = 10,
#'   num_samples = 5
#' )
#' }
#'
#' @export
generate_adjacency_matrices <- function(threshold = 0.7, num_repetitions = 50, num_samples = 10,
                                        correlation = "spearman", expression_mat,
                                        metaData, genesInfo, DEGs){

  age_groups <- sort(unique(metaData$Groups))
  adjacency_matrices <- list()

  for (i in 1:num_repetitions) {
    repetition_adjacency_matrices <- list()

    for (age_group in age_groups) {
      sub_meta <- metaData[metaData$Groups == age_group, ]
      sub_norm <- expression_mat[rownames(expression_mat) %in% DEGs$DEgenes,
                                 colnames(expression_mat) %in% sub_meta$Samples]
      selected_samples <- sub_norm[, sample(ncol(sub_norm), num_samples)]
      selected_samples_t <- t(selected_samples)

      gsg <- WGCNA::goodSamplesGenes(selected_samples_t, verbose = 0)
      if (!gsg$allOK) {
        selected_samples_t <- selected_samples_t[gsg$goodSamples, gsg$goodGenes]
      }

      spearman_matrix <- cor(selected_samples_t, method = correlation)
      diag(spearman_matrix) <- 0
      weighted_adj_mat <- spearman_matrix
      weighted_adj_mat[abs(spearman_matrix) < threshold] <- 0
      repetition_adjacency_matrices[[age_group]] <- weighted_adj_mat
    }

    adjacency_matrices[[i]] <- repetition_adjacency_matrices
  }

  return(adjacency_matrices)
}



#' Summaries of a clustering.
#'
#' Returns Pielou's J; Gini index and Simpson's D as measures of the inequality
#' in cluster sizes across a clustering.
#'
#' @param clust_df Dataframe with context columns and a column with cluster
#' membership for that context.
#' @param clust_col Name of column in `clust_df` with class membership.
#'
#' @return single row tibble with summary information. `_rescale` cols are
#' `1 - value` so that 1 (higher) values represent even/equal clusters.
#' @export
#'
#' @example inst/examples/make_clusters_ex.R
make_summary_df <- function(clust_df
                            , clust_col = "cluster"
                            ) {

  orig_tab <- table(clust_df[[clust_col]])
  tab <- orig_tab[orig_tab > 0] |> sort()
  prop_tab <- tab / sum(tab)

  tibble::tibble(n_bins = nrow(clust_df)
                 , orig_k = length(orig_tab)
                 , actual_k = length(tab)
                 , pielou_j = (-sum(prop_tab * log2(prop_tab))) / log2(length(tab))
                 , gini = sum((2 * 1:(length(tab)) - length(tab) - 1) * tab) / (length(tab) * sum(tab))
                 , simp_d = sum(prop_tab ^ 2)
                 ) |>
    dplyr::mutate(gini_rescale = 1 - gini
                  , simp_d_rescale = 1 - simp_d
                  )

}

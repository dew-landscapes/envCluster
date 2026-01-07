#' Simple summary of a clustering.
#'
#' @param clust_df Dataframe with context columns and a column with cluster
#' membership for that context.
#' @param clust_col Name of column in `clust_df` with class membership.
#' @param min_sites Desired minimum absolute number of bins in a class.
#'
#' @return single row tibble with summary information
#' @export
#'
#' @example inst/examples/make_clusters_ex.R
make_summary_df <- function(clust_df
                            , clust_col = "cluster"
                            , min_sites = 10
                            ) {

  tab <- table(clust_df[[clust_col]])

  tibble::tibble(min_clust_size = min(tab)
                 , av_clust_size = mean(tab)
                 , max_clust_size = max(tab)
                 , n_min_clusters = sum(tab > min_sites)
                 , n_min_sites = sum(tab[tab > min_sites])
                 , prop_min_clusters = n_min_clusters / length(tab)
                 , prop_min_sites = n_min_sites / nrow(clust_df)
                 )

}

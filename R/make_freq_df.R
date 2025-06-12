#' Proportion of sites and clusters with at least one taxa in common
#'
#' Creates a dataframe of proportion of clusters (and contexts within those
#' clusters) that contain at least one taxa present in greater than `freq_thres`
#' proportion of that clusters contexts. Contexts are usually sites.
#'
#' @param clust_df Dataframe with `context` column(s) and a column with cluster
#' membership for that context. Optional if `clust_col` appears in bio_df.
#' @param bio_df Dataframe containing the site and taxa data in long format.
#' @param context Character. Name(s) of column(s) that define the context.
#' @param clust_col Character. Name of column containing cluster membership.
#' @param taxa_col Character. Name of column containing the taxa names.
#' @param freq_thresh Proportion of sites within a cluster that a taxa should be
#' present
#'
#' @return Single row dataframe of `prop_freq_clusters` and `prop_freq_sites`
#' @export
#'
#' @example inst/examples/make_clusters_ex.R
make_freq_df <- function(clust_df = NULL
                         , bio_df
                         , context
                         , clust_col = "cluster"
                         , taxa_col = "taxa"
                         , freq_thresh = 0.9
                         ) {

  if(!is.null(clust_df)) {

    bio_df <- bio_df %>%
      dplyr::inner_join(clust_df)

  }

  clus_sites <- paste0(clust_col, "_sites")
  clus_recs <- paste0(clust_col, "_records")

  c_sites <- bio_df |>
    dplyr::group_by(dplyr::across(tidyselect::any_of(c(clust_col, context)))) |>
    dplyr::summarise(n_taxa = dplyr::n()) |>
    dplyr::ungroup() |>
    dplyr::count(!!rlang::ensym(clust_col)
                 , name = clus_sites
                 )

  c_recs <- bio_df |>
    dplyr::group_by(dplyr::across(tidyselect::any_of(c(clust_col, taxa_col)))) |>
    dplyr::summarise(!!rlang::ensym(clus_recs) := dplyr::n()) |>
    dplyr::ungroup()

  c_sites |>
    dplyr::left_join(c_recs) |>
    dplyr::mutate(cluster_prop = !!rlang::ensym(clus_recs) / !!rlang::ensym(clus_sites)) |>
    dplyr::group_by(!!rlang::ensym(clust_col), !!rlang::ensym(clus_sites)) |>
    dplyr::summarise(n_freq_taxa = sum(cluster_prop > freq_thresh)) |>
    dplyr::ungroup() %>%
    dplyr::summarise(prop_freq_clusters = sum(n_freq_taxa > 0) / nrow(.)
                     , prop_freq_sites = sum(n_freq_taxa > 0) / sum(!!rlang::ensym(clus_sites))
                     )

}

#' Calculate gap statistic for a clustering
#'
#' @param clust_df Dataframe with context columns and a column with cluster
#' membership for that context.
#' @param clust_col Character. Name of column in `clusters` identifying cluster
#' membership.
#' @param dist_mat Distance matrix (not distance object)
#' @param n_sample Numeric. Number of times to shuffle the membership to
#' recalculate the wss
#'
#' @return Single row dataframe with macro_wss and macro_gap plus a list column
#' of one row per level of `clust_col`, summarising the gap and wss data for
#' each level.
#' @export
#'
#' @example inst/examples/make_clusters_ex.R
make_gap_df <- function(clust_df
                        , dist_mat
                        , clust_col = "cluster"
                        , n_sample = 10
                        ) {

  .clust_col <- clust_col

  tibble::tibble(clusters = list(clust_df)) |>
      dplyr::mutate(wss = purrr::map(clusters
                                     , \(x) envCluster::calc_wss(x
                                                                 , dist_mat = dist_mat
                                                                 , clust_col = !!rlang::ensym(clust_col)
                                                                 )
                                     )
                    ) |>
      dplyr::mutate(macro_wss = purrr::map_dbl(wss, \(x) sum(x$wss))) |>
      dplyr::cross_join(tibble::tibble(run = 1:n_sample)) |>
      dplyr::mutate(wss_for_gap = purrr::map(clusters
                                             , calc_wss
                                             , dist_mat = dist_mat
                                             , clust_col = !!rlang::ensym(.clust_col)
                                             , do_sample = TRUE
                                             )
                    ) |>
      tidyr::unnest(cols = c(wss, wss_for_gap)
                    , names_repair = tidyr::tidyr_legacy
                    ) |>
      dplyr::mutate(gap = wss1 - wss) |>
      dplyr::group_by(!!rlang::ensym(clust_col), wss, macro_wss) %>%
      dplyr::summarise(sample_wss = mean(wss1)
                       , gap_se = sd(gap) / sqrt(dplyr::n())
                       , gap = mean(gap)
                       , samples = n_sample
                       ) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(macro_gap = mean(gap)) |>
      tidyr::nest(gap = -c(tidyselect::matches("macro")))

}


#' Make clusters using mstknnclust::mst.knn
#'
#' @param dist_mat_bio matrix. usually as.matric(distance object)
#' @param clust_col Character. Name to give column representing each of `groups`
#' @param clust_col_out Character. Name to give column representing `groups` in
#' words.
#' @param sites dataframe of sites represented by the distance object
#' @param ... passed to `mstknnclust::mst.knn()` (thus, really just
#' `suggested_k`).
#'
#' @returns List retruned by `mstknnclust::mst.knn()` with an additional
#' `clusters` object as `sites` augmented with their associated cluster.
#' @export
#'
#' @example inst/examples/make_clusters_ex.R
make_knn_clust <- function(dist_mat_bio
                           , clust_col = "clust"
                           , clust_col_out = "cluster"
                           , sites = NULL
                           , ...
                           ) {

  results <- mstknnclust::mst.knn(dist_mat_bio)

  results$clusters <- if(!is.null(sites)) {

    sites |>
      dplyr::left_join(tibble::as_tibble(results$partition) |>
                         dplyr::rename(site_id = object
                                       , !!rlang::ensym(clust_col) := cluster
                                       ) |>
                         dplyr::mutate(site_id = as.integer(site_id)
                                       , !!rlang::ensym(clust_col_out) := envFunc::numbers2words(!!rlang::ensym(clust_col))
                                       , !!rlang::ensym(clust_col_out) := forcats::fct_reorder(!!rlang::ensym(clust_col_out)
                                                                                               , !!rlang::ensym(clust_col)
                                                                                               )
                                       )
                       )

  } else {

    tibble::as_tibble(results$partition) |>
      dplyr::rename(site_id = object
                    , !!rlang::ensym(clust_col) := cluster
                    ) |>
      dplyr::mutate(site_id = as.integer(site_id)
                    , !!rlang::ensym(clust_col_out) := envFunc::numbers2words(!!rlang::ensym(clust_col))
                    , !!rlang::ensym(clust_col_out) := forcats::fct_reorder(!!rlang::ensym(clust_col_out)
                                                                            , !!rlang::ensym(clust_col)
                                                                            )
                    )

  }

  return(results)

}

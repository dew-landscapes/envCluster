#' Apply clustering algorithms
#'
#' @param dend Dendogram
#' @param group_range Numeric vector indicating the range of groups within a
#' clustering.
#' @param clust_col Character. Name to give column representing each of `groups`
#' @param clust_col_out Character. Name to give column representing `groups` in
#' words.
#' @param sites Dataframe used to create `dist_bio` (and, if used,
#' `dist_env`)
#'
#' @return Dataframe with groups and methods columns and a list column of the
#' clustering for that number of groups and method as a tibble with one column
#' (numeric) 'clust_col' indicating group membership.
#' @export
#'
#' @example inst/examples/make_clusters_ex.R
#'

make_clusters <- function(dend
                          , group_range = 2:100
                          , clust_col = "clust"
                          , clust_col_out = "cluster"
                          , sites
                          ) {

  clusters <- sites |>
    dplyr::bind_cols(tibble::as_tibble(stats::cutree(dend, group_range))) %>%
    tidyr::pivot_longer((ncol(sites) + 1):ncol(.)
                        , names_to = "groups"
                        , names_transform = \(x) as.integer(x)
                        , values_to = clust_col
                        ) |>
    tidyr::nest(clusters = -c(groups)) |>
    dplyr::mutate(clusters = purrr::map(clusters
                                        , \(x) x |>
                                          dplyr::mutate(!!rlang::ensym(clust_col_out) := envFunc::numbers2words(!!rlang::ensym(clust_col))
                                                        , !!rlang::ensym(clust_col_out) := forcats::fct_reorder(!!rlang::ensym(clust_col_out)
                                                                                                                , !!rlang::ensym(clust_col)
                                                                                                                )
                                                        )
                                        )
                  )

  return(clusters)

}

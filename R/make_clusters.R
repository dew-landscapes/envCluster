#' Apply clustering algorithms
#'
#' @param dend Dendogram
#' @param group_range Numeric vector indicating the range of groups within a
#' clustering.
#' @param bins Dataframe with bin identifiers used to create `dist_bio` (and, if
#'  used, `dist_env`)
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
                          , bins
                          ) {

  clusts <- stats::cutree(dend, group_range)

  if(length(dim(clusts)) > 1) {

    clusts <- tibble::as_tibble(clusts)

  } else {

    clusts <- tibble::enframe(clusts
                              , name = NULL
                              , value = group_range
                              )

  }

  clusters <- bins |>
    dplyr::bind_cols(clusts) %>%
    tidyr::pivot_longer((ncol(bins) + 1):ncol(.)
                        , names_to = "groups"
                        , names_transform = \(x) as.integer(x)
                        , values_to = "clust"
                        ) |>
    tidyr::nest(clusters = -c(groups)) |>
    dplyr::mutate(clusters = purrr::map(clusters
                                        , \(x) x |>
                                          dplyr::mutate(cluster = envFunc::numbers2words(clust)
                                                        , cluster = forcats::fct_reorder(cluster
                                                                                         , clust
                                                                                         )
                                                        )
                                        )
                  )

  return(clusters)

}

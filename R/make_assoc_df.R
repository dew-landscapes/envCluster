#' Create a dataframe of strength of association between taxa and cluster
#'
#' @param ... Passed to `indicspecies::strassoc()`. For phi, X should be
#' presence/absence. The `func` argument to `indicspecies::strassocc()` is
#' specified within `make_assoc_df` so does not need to be provided (and will
#' be ignored if it is).
#' @param ret Character. Return one of: "prop_only"; "filtered"; or "full"
#' @param rg_thresh Numeric. Threshold for acceptance as an association.
#' @param rg_thresh_count Numeric. How many taxa within a cluster need to meet
#' rg_thresh?
#'
#' @return For "prop_only": a single row tibble with columns: `prop_rg_clusters`
#' (proportion of clusters with at least `rg_thresh_count` taxa meeting
#' `rg_thresh`an association); and `prop_rg_sites`
#' (proportion of sites those clusters). For "filtered": a
#' tibble with any association values below `rg_thresh` filtered. For "full": a
#' tibble with each taxa and its association value with every cluster (`clust`
#' as numeric). If rg_thresh is `NULL` (or `0`), "full" result is returned. Columns `rg`,
#' `b` and `ag` show the matching `indicspecies::strassoc()` result for the
#' corresponding `func` argument (minus the `.` in `func`). Both filtered and
#' full tibbles have the `prop_rg_` columns appended (i.e. the same value
#' repeated for every row).
#' @export
#'
#' @example inst/examples/make_clusters_ex.R
make_assoc_df <- function(...
                          , ret = "prop_only"
                          , rg_thresh = 0.2
                          , rg_thresh_count = 2
                          ) {

  dots_list <- list(...)
  cl <- dots_list$cluster

  # rg -------
  dots_list$func <- "r.g"

  rg <- R.utils::doCall(indicspecies::strassoc
                        , args = dots_list
                        )

  if(!"matrix" %in% class(rg)) {

    rg_l <- rg$lowerCI
    rg_u <- rg$upperCI
    rg <- rg$stat

  }

  rg <- rg |>
    tibble::as_tibble(rownames = "taxa") |>
    tidyr::pivot_longer(tidyselect::any_of(as.character(cl))
                        , values_to = "rg"
                        , names_to = "clust"
                        )

  if(exists("rg_l", inherits = FALSE)) {

    rg <- rg |>
      dplyr::left_join(rg_l |>
                         tibble::as_tibble(rownames = "taxa") |>
                         tidyr::pivot_longer(tidyselect::any_of(as.character(cl))
                                             , names_to = "clust"
                                             , values_to = "lower_ci"
                                             )
                       ) |>
      dplyr::left_join(rg_u |>
                         tibble::as_tibble(rownames = "taxa") |>
                         tidyr::pivot_longer(tidyselect::any_of(as.character(cl))
                                             , names_to = "clust"
                                             , values_to = "upper_ci"
                                             )
                       )

  }

  if(!is.null(rg_thresh)) {

    with_assoc <- rg |>
      dplyr::filter(rg >= rg_thresh) |>
      dplyr::count(clust, taxa) |>
      dplyr::count(clust) |>
      dplyr::filter(n >= rg_thresh_count)

    prop_rg_clusters <- nrow(with_assoc) / length(unique(cl))
    prop_rg_sites <- length(cl[cl %in% with_assoc$clust]) / length(cl)

    if(ret == "prop_only") {

      return(tibble::tibble(prop_rg_clusters = prop_rg_clusters, prop_rg_sites = prop_rg_sites))

    }

    if(ret == "filtered") {

      rg <- rg |>
        dplyr::filter(rg > rg_thresh)

    }

    rg <- rg |>
      dplyr::mutate(prop_rg_clusters = prop_rg_clusters
                    , prop_rg_sites = prop_rg_sites
                    )

  }

  # b --------
  # sensitivity. the proportion of sites in a cluster at which the species occurs
  dots_list$func <- "b"

  b <- R.utils::doCall(indicspecies::strassoc
                       , args = dots_list
                       ) |>
    tibble::as_tibble(rownames = "taxa") |>
    tidyr::pivot_longer(tidyselect::any_of(as.character(cl))
                        , values_to = "b"
                        , names_to = "clust"
                        )

  # ag -------
  # the proportion of the species occurrences that fall in that cluster
  dots_list$func <- "a.g"

  ag <- R.utils::doCall(indicspecies::strassoc
                        , args = dots_list
                        ) |>
    tibble::as_tibble(rownames = "taxa") |>
    tidyr::pivot_longer(tidyselect::any_of(as.character(cl))
                        , values_to = "ag"
                        , names_to = "clust"
                        )

  # result -------
  result <- rg |>
    dplyr::inner_join(b) |>
    dplyr::inner_join(ag) |>
    dplyr::mutate(clust = as.numeric(clust))

  return(result)

}

#' Create a dataframe of strength of association between taxa and cluster
#'
#' @param ... Passed to `indicspecies::assoc()`. For phi, X should be presence/
#' absence
#' @param ret Character. Return one of: "prop_only"; "filtered"; or "full"
#' @param assoc_thresh Numeric. Threshold for acceptance as an association.
#'
#' @return For "prop_only": a single row tibble with prop_ind_clusters
#' (proportion of clusters with an association) and prop_ind_sites
#' (proportion of sites in clusters with an association). For "filtered": a
#' tibble with any association values above assoc_thresh. For "full": a tibble
#' with each taxa and its association value with every cluster (`clust` as
#' numeric). If assoc_thresh is `NULL`, "full" result is returned.
#' @export
#'
#' @example inst/examples/make_clusters_ex.R
make_assoc_df <- function(...
                          , ret = "prop_only"
                          , assoc_thresh = 0.5
                          ){

  dots_list <- list(...)

  result <- R.utils::doCall(indicspecies::strassoc
                            , args = dots_list
                            )

  cl <- dots_list$cluster

  if(!"matrix" %in% class(result)) {

    result_l <- result$lowerCI
    result_u <- result$upperCI
    result <- result$stat

  }

  result <- result |>
    tibble::as_tibble(rownames = "taxa") |>
    tidyr::pivot_longer(tidyselect::any_of(as.character(cl))
                        , values_to = "assoc"
                        , names_to = "clust"
                        )

  if(exists("result_l", inherits = FALSE)) {

    result <- result |>
      dplyr::left_join(result_l |>
                         tibble::as_tibble(rownames = "taxa") |>
                         tidyr::pivot_longer(tidyselect::any_of(as.character(cl))
                                             , names_to = "clust"
                                             , values_to = "lower_ci"
                                             )
                       ) |>
      dplyr::left_join(result_u |>
                         tibble::as_tibble(rownames = "taxa") |>
                         tidyr::pivot_longer(tidyselect::any_of(as.character(cl))
                                             , names_to = "clust"
                                             , values_to = "upper_ci"
                                             )
                       )

  }

  if(!is.null(assoc_thresh)) {

    with_ind <- result |>
      dplyr::filter(assoc > assoc_thresh) |>
      dplyr::count(clust)

    prop_ind_clusters <- nrow(with_ind) / length(unique(cl))
    prop_ind_sites <- length(cl[cl %in% with_ind$clust]) / length(cl)

    if(ret == "prop_only") {

      result <- tibble::tibble(prop_assoc_clusters = prop_ind_clusters
                               , prop_assoc_sites = prop_ind_sites
                               )
    } else if(ret == "filtered") {

      result <- result |>
        dplyr::filter(assoc > assoc_thresh) |>
        dplyr::mutate(prop_assoc_clusters = prop_ind_clusters
                      , prop_assoc_sites = prop_ind_sites
                      )

    } else if(ret == "full") {

      result <- result |>
        plyr::mutate(prop_assoc_clusters = prop_ind_clusters
                     , prop_assoc_sites = prop_ind_sites
                     )

    }

  }

  return(result)

}

#' Create a dataframe of indicator values for a clustering
#'
#' @param clust_df Dataframe with context columns and a column with cluster
#' membership for that context. Optional if `clust_col` appears in bio_wide.
#' @param bio_wide Dataframe containing the site and taxa data in wide format.
#' @param context Character. Name(s) of column(s) that define the context.
#' @param taxas Character. Vector of taxa names (column names in `bio_wide`).
#' Optional if `bio_wide` contains only taxa names and `context` columns.
#' @param clust_col Character. Name of column containing cluster membership.
#' The column `clust_col` in `clust_df` should be a factor.
#' @param p_val_thresh Numeric. Threshold p value at which to accept a taxa as
#' an indicator.
#' @param ... Passed to `labdsv::indval()`
#'
#' @return Dataframe of each taxa and the cluster (clust as numeric, cluster as
#' character) class to which it is most likely an indicator, plus the following
#' values from `labdsv::indval()` output: ind_val, p_val, abu and frq.
#' @export
#'
#' @example inst/examples/make_clusters_ex.R
make_ind_val_df <- function(clust_df = NULL
                            , bio_wide
                            , context
                            , taxas = NULL
                            , clust_col = "cluster"
                            , p_val_thresh = 0.05
                            , ...
                            ){

  # check clust_col is factor #1
  # the order of results is likely to be odd if clust_col is not a factor
  stopifnot("'clust_col' needs to be a factor" = is.factor(clust_df[[clust_col]]))

  if(is.null(taxas)) {

    taxas <- names(bio_wide)[!names(bio_wide) %in% c(context, clust_col)]

  }

  if(!is.null(clust_df)) {

    bio_wide <- bio_wide |>
      dplyr::inner_join(clust_df |>
                          dplyr::distinct(dplyr::across(tidyselect::any_of(c(context, clust_col))))
                        )

  }

  clust_ind <- labdsv::indval(x = bio_wide[,names(bio_wide) %in% c(taxas)]
                              , clustering = as.character(bio_wide[clust_col][[1]])
                              , ...
                              )

  clusts <- levels(bio_wide[clust_col][[1]])

  res <- clust_ind$relabu |>
    tibble::as_tibble(rownames = "taxa") |>
    tidyr::pivot_longer(any_of(clusts)
                        , names_to = clust_col
                        , values_to = "abu"
                        ) |>
    dplyr::inner_join(clust_ind$relfrq |>
                        tibble::as_tibble(rownames = "taxa") |>
                        tidyr::pivot_longer(any_of(clusts)
                                            , names_to = clust_col
                                            , values_to = "frq"
                                            )
                      ) |>
    dplyr::inner_join(clust_ind$indval |>
                        tibble::as_tibble(rownames = "taxa") |>
                        tidyr::pivot_longer(any_of(clusts)
                                            , names_to = clust_col
                                            , values_to = "ind_val"
                                            )
                      ) |>
    dplyr::group_by(taxa) |>
    dplyr::filter(ind_val == max(ind_val)) |>
    dplyr::ungroup() |>
    dplyr::inner_join(clust_ind$pval |>
                        tibble::as_tibble(rownames = "taxa") |>
                        dplyr::rename(p_val = value)
                      ) |>
    dplyr::mutate(!!rlang::ensym(clust_col) := factor(!!rlang::ensym(clust_col), levels = clusts)) |>
    dplyr::select(!!rlang::ensym(clust_col), everything()) |>
    dplyr::arrange(!!rlang::ensym(clust_col), desc(ind_val)) |>
    tidyr::nest(inds = - !!rlang::ensym(clust_col)) %>%
    dplyr::mutate(has_ind = purrr::map_lgl(inds
                                           , \(x) any(x$p_val <= p_val_thresh)
                                           )
                  , prop_ind_clusters = sum(has_ind) / nrow(.)
                  , prop_ind_sites = nrow(. |> dplyr::filter(has_ind) |> dplyr::inner_join(clust_df)) / nrow(clust_df)
                  ) |>
    tidyr::nest(inds = - tidyselect::matches("prop_"))

  return(res)

}

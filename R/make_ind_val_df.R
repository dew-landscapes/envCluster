#' Create a dataframe of indicator values for a clustering
#'
#' @param clust_df Dataframe with context columns and a column with cluster
#' membership for that context. Optional if `clust_col` appears in bio_wide.
#' @param bio_df Dataframe containing the site and taxa data in long format.
#' @param context Character. Name(s) of column(s) that define the context.
#' @param clust_col Character. Name of column containing cluster membership.
#' @param cov_col Character. Name of column containing abundance data (often
#' 'cover' values).
#' @param taxa_col Character. Name of column containing the taxa names.
#' @param ... Passed to `labdsv::indval()`
#'
#' @return Dataframe of each taxa and the cluster (clust as numeric, cluster as
#' character) class to which it is most likely an indicator, plus the following
#' values from `labdsv::indval()` output: ind_val, p_val, abu and frq.
#' @export
#'
#' @examples
make_ind_val_df <- function(clust_df = NULL
                            , bio_df
                            , context
                            , clust_col = "cluster"
                            , cov_col = "use_cover"
                            , taxa_col = "taxa"
                            , ...
                            ){

  taxas <- sort(unique(bio_df$taxa))

  if(is.numeric(clust_col)) clust_col <- names()

  bio_wide <- envCluster::make_wide_df(bio_df
                                       , context = c(context, clust_col)
                                       , num_col = cov_col
                                       , taxa_col = taxa_col
                                       )

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

  clust_ind$relabu %>%
    tibble::as_tibble(rownames = "taxa") %>%
    tidyr::pivot_longer(any_of(clusts)
                        , names_to = clust_col
                        , values_to = "abu"
                        ) %>%
    dplyr::inner_join(clust_ind$relfrq %>%
                        tibble::as_tibble(rownames = "taxa") %>%
                        tidyr::pivot_longer(any_of(clusts)
                                            , names_to = clust_col
                                            , values_to = "frq"
                                            )
                      ) %>%
    dplyr::inner_join(clust_ind$indval %>%
                        tibble::as_tibble(rownames = "taxa") %>%
                        tidyr::pivot_longer(any_of(clusts)
                                            , names_to = clust_col
                                            , values_to = "ind_val"
                                            )
                      ) %>%
    dplyr::group_by(taxa) %>%
    dplyr::filter(ind_val == max(ind_val)) %>%
    dplyr::ungroup() %>%
    dplyr::inner_join(clust_ind$pval %>%
                        tibble::as_tibble(rownames = "taxa") %>%
                        dplyr::rename(p_val = value)
                      ) %>%
    dplyr::mutate(!!rlang::ensym(clust_col) := factor(!!rlang::ensym(clust_col), levels = clusts)) %>%
    dplyr::select(!!rlang::ensym(clust_col), everything()) %>%
    dplyr::arrange(!!rlang::ensym(clust_col), desc(ind_val))

}

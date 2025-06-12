
#' Calculate within cluster sum of squares
#'
#'
#' @param clust_df Dataframe with column of cluster membership.
#' @param dist_obj Distance object from which clusters were generated.
#' @param dist_mat Result from `as.matrix(dist_obj)`, or equivalent. Required if
#' `dist_obj` is not supplied.
#' @param clust_col Character. Name of column in `clust_df` that contains the
#' clusters.
#' @param do_boot Logical. Sample `clust_col` (without replacement) before
#' calculation. As step in gap-statistic calculation.
#'
#' @return Dataframe of within group sum-of-squares for each cluster.
#' @export
#' @keywords internal
#'
#' @examples
calc_wss <- function(clust_df
                    , dist_obj = NULL
                    , dist_mat = NULL
                    , clust_col = "clust"
                    , do_sample = FALSE
                    ) {

  if(isTRUE(is.null(dist_mat))) {

    dist_mat <- as.matrix(dist_obj)

  }

  clust_df %>%
    {if(do_sample) (.) %>% dplyr::mutate(!!rlang::ensym(clust_col) := sample(!!rlang::ensym(clust_col))) else (.)} %>%
    dplyr::mutate(id = dplyr::row_number()) %>%
    dplyr::group_by(!!rlang::ensym(clust_col)) %>%
    tidyr::nest() %>%
    dplyr::ungroup() %>%
    #dplyr::sample_n(2) %>% # TESTING
    dplyr::mutate(wss = purrr::map_dbl(data
                                       , ~sum(dist_mat[.$id,.$id]^2)/(2*nrow(.))
                                       )
                  ) %>%
    dplyr::select(-data) %>%
    dplyr::arrange(!!rlang::ensym(clust_col))

}

#' Create a dataframe of silhouette widths for a clustering
#'
#' @param clust_df Dataframe with context columns and a column with cluster
#' membership for that context.
#' @param dist_obj Distance object with sites in the same order as `clust_df`.
#' @param clust_col Character. Name of column containing cluster membership in
#' `clust_df`.
#'
#' @return `clust_df` with additional columns for silhouette width (as
#' `sil_width`) and neighbouring cluster (as `neighbour`).
#' @export
#'
#' @example inst/examples/make_clusters_ex.R
make_sil_df <- function(clust_df, dist_obj, clust_col = "cluster"){

  clusts <- clust_df[clust_col][[1]]

  if(!is.integer(clusts)) {

    id <- unique(clusts)

    clusts <- match(clusts, id)

  }

  sil_obj <- cluster::silhouette(clusts, dist_obj)

  clust_df %>%
    dplyr::mutate(sil_width = sil_obj[ ,3]
                  , neighbour = sil_obj[ ,2]
                  , macro_sil = mean(sil_width)
                  ) |>
    tidyr::nest(sil = -c(macro_sil))

}

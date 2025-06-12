#' Make a dendogram
#'
#'
#' @param clust_df Dataframe with context columns and a column with cluster
#' membership for that context.
#' @param dist_obj Distance object with sites in the same order as `clust_df`.
#' @param clust_col Character. Name of column containing cluster membership in
#' `clust_df`.
#' @param method Character. Name of method to use to cluster `dist_obj`
#'
#' @return Dendogram object
#' @export
#'
#' @example inst/examples/make_clusters_ex.R
make_dend <- function(clust_df
                      , dist_obj
                      , clust_col = "cluster"
                      , method
                      ) {

  dend <- hclust(dist_obj
                 , method = method
                 ) |>
    as.dendrogram()

  if(! "colour" %in% names(clust_df)) {

    colours <- viridis::viridis(n = length(unique(clust_df[[clust_col]])))

    clust_df <- clust_df |>
      dplyr::mutate(colour = colours[as.integer(as.factor(clust_df[[clust_col]]))])

  }

  k <- length(unique(factor(clust_df$cluster)))

  use_col <- clust_df$colour[order.dendrogram(dend)]
  use_lab <- clust_df$clust[order.dendrogram(dend)]

  use_col_branch <- clust_df |>
    dplyr::slice(order.dendrogram(dend)) |>
    dplyr::distinct(cluster, colour) |>
    dplyr::pull(colour)

  use_lab_branch <- clust_df |>
    dplyr::slice(order.dendrogram(dend)) |>
    dplyr::distinct(clust, colour) |>
    dplyr::pull(clust)

  dend <- dend %>%
    dendextend::set("labels", use_lab) %>%
    dendextend::set("labels_colors", use_col) %>%
    dendextend::set("labels_cex", 0.01) %>%
    dendextend::set("branches_k_color"
                    , value = use_col_branch
                    , k = k
                    ) %>%
    dendextend::colour_branches(k = k
                                , col = use_col_branch
                                , groupLabels = use_lab_branch
                                )

  return(dend)

}

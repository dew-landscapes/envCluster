
#' Mean silhouette width for a cluster
#'
#' @param cluster Name of cluster
#' @param sil_df Output from `make_sil_df()`
#'
#' @return Numeric. Mean of silhouette widths for sites in a cluster.
#' @export
#'
#' @examples
make_cluster_sil_mean <- function(cluster,sil_df) {

  cluster <- as.character(cluster)

  sil_df %>%
    dplyr::filter(cluster == cluster) %>%
    dplyr::pull(sil_width) %>%
    mean

}




#' Explore clusters
#'
#' Explore clusters using silhouette width, within-group sum-of-squares,
#' indicators and widespread taxa.
#'
#' @param clusters_df Dataframe with column 'clusters'.
#' @param cores How many cores to use in `multidplyr::new_cluster()`?
#' @param obj_list List of objects required:
#' \describe{
#'   \item{out_dir}{directory into which to save results.}
#'   \item{dist_flor}{`dist(flor_wide)`.}
#'   \item{dist_flor_mat}{`as.matrix(dist(flor_wide))`.}
#'   \item{dist_env}{`dist(env_wide)`.}
#'   \item{dist_env_mat}{`as.matrix(dist(env_wide))`.}
#'   \item{flor_wide}{Wide version of flor_tidy.}
#'   \item{flor_tidy}{Cleaned data set of taxa observations.}
#'   \item{visit_cols}{context.}
#'   \item{p_thresh}{usually 0.05.}
#'   \item{n_sites}{How many unique contexts in `flor_tidy`.}
#'   \item{most_freq_prop_thresh}{Threshold for proportion of sites in a cluster that must all contain at least one taxa in common.}
#' }
#'
#' @return clusters_df with added columns. Saves each output along the way to `out_dir`.
#' @export
#'
#' @examples
make_clusters_explore <- function(clusters_df
                                  , cores = 1
                                  , obj_list
                                  ) {

  # Variables required come in via obj_list
  purrr::walk2(names(obj_list)
               , obj_list
               , assign
               )

  return_result <- function(df) {

    df %>%
      dplyr::select(-clusters) %>%
      dplyr::collect()

  }

  #-------deal with cl-------

  cl <- multidplyr::new_cluster(cores)


  #---------Start explore-----------

  # load envCluster to each core
  multidplyr::cluster_library(cl
                              , c("envCluster")
                              )

  clusters_use_exp <- clusters_df %>%
    dplyr::select(method, groups, clusters) %>%
    multidplyr::partition(cl)

  explore_res <- list()

  #------silhouette-----

  out_file <- fs::path(out_dir, "clusters_sil.rds")

  if(!file.exists(out_file)) {

    multidplyr::cluster_copy(cl
                             , c("dist_flor")
                             )

    explore_res$clusters_sil <- clusters_use_exp  %>%
      dplyr::mutate(sil = purrr::map(clusters
                                     , make_sil_df
                                     , dist_obj = dist_flor
                                     , clust_col = "clust"
                                     )
                    , macro_sil = purrr::map_dbl(sil
                                                 , ~mean(.$sil_width)
                                                 )
                    ) %>%
      return_result()

    multidplyr::cluster_rm(cl
                           , c("dist_flor")
                           )

    rio::export(explore_res$clusters_sil
                , out_file
                )

  }


  #-------wss--------

  out_file <- fs::path(out_dir, "clusters_wss.rds")

  if(!file.exists(out_file)) {

    dist_flor_mat <- as.matrix(dist_flor)

    multidplyr::cluster_copy(cl
                             , c("dist_flor_mat")
                             )

    explore_res$clusters_wss <- clusters_use_exp %>%
      dplyr::mutate(wss = purrr::map(clusters
                                     , calc_wss
                                     , dist_mat = dist_flor_mat
                                     )
                    , macro_wss = purrr::map_dbl(wss
                                                 , ~sum(.$wss)
                                                 )
                    ) %>%
      return_result()

    multidplyr::cluster_rm(cl
                           , c("dist_flor_mat")
                           )

    rio::export(explore_res$clusters_wss
                , out_file
                )

  }


  #-----silhouette env------

  out_file <- fs::path(out_dir, "clusters_sil_env.rds")

  if(!file.exists(out_file)) {

    multidplyr::cluster_copy(cl
                             , c("dist_env")
                             )

    explore_res$clusters_sil_env <- clusters_use_exp %>%
      dplyr::mutate(sil_env = purrr::map(clusters
                                         , make_sil_df
                                         , dist_obj = dist_env
                                         , clust_col = "clust"
                                         )
                    , macro_sil_env = purrr::map_dbl(sil_env
                                                     , ~mean(.$sil_width)
                                                     )
                    ) %>%
      return_result()

    multidplyr::cluster_rm(cl
                           , c("dist_env")
                           )

    rio::export(explore_res$clusters_sil_env
                , out_file
                )

  }


  #-------wss env------

  out_file <- fs::path(out_dir, "clusters_wss_env.rds")

  if(!file.exists(out_file)) {

    dist_env_mat <- as.matrix(dist_env)

    multidplyr::cluster_copy(cl
                             , c("dist_env_mat")
                             )

    explore_res$clusters_wss_env <- clusters_use_exp %>%
      dplyr::mutate(wss_env = purrr::map(clusters
                                         , calc_wss
                                         , dist_mat = dist_env_mat
                                         )
                    , macro_wss_env = purrr::map_dbl(wss_env
                                                     , ~sum(.$wss)
                                                     )
                    ) %>%
      return_result()

    multidplyr::cluster_rm(cl
                           , c("dist_env_mat")
                           )

    rio::export(explore_res$clusters_wss_env
                , out_file
                )

  }


  #------ind val---------

  out_file <- fs::path(out_dir, "clusters_ind_val.rds")

  if(!file.exists(out_file)) {

    add_to_cluster <- c("flor_wide"
                        , "flor_tidy"
                        , "visit_cols"
                        , "p_thresh"
                        , "n_sites"
                        )

    multidplyr::cluster_copy(cl
                             , add_to_cluster
                             )

    explore_res$clusters_ind_val <- clusters_use_exp %>%
      dplyr::mutate(ind_val = purrr::map(clusters
                                         , make_ind_val_df
                                         , bio_wide = flor_wide
                                         , taxas = unique(flor_tidy$taxa)
                                         , context = visit_cols
                                         )
                    , n_ind_clusters = purrr::map_dbl(ind_val
                                                      , clusters_with_indicator
                                                      , thresh = p_thresh
                                                      )
                    , n_ind_sites = purrr::map2_dbl(ind_val
                                                    , clusters
                                                    , sites_with_indicator
                                                    , thresh = p_thresh
                                                    )
                    , prop_ind_clusters = n_ind_clusters/groups
                    , prop_ind_sites = n_ind_sites/n_sites
                    ) %>%
      return_result()

    multidplyr::cluster_rm(cl
                           , add_to_cluster
                           )

    rio::export(explore_res$clusters_ind_val
                , out_file
                )

  }


  #-------freq-------

  out_file <- fs::path(out_dir, "clusters_freq.rds")

  if(!file.exists(out_file)) {

    add_to_cluster <- c("flor_tidy"
                        , "most_freq_prop_thresh"
                        , "n_sites"
                        )

    multidplyr::cluster_copy(cl
                             , add_to_cluster
                             )

    explore_res$clusters_freq <- clusters_use_exp %>%
      dplyr::mutate(n_freq_clusters = purrr::map_dbl(clusters
                                                     , clusters_with_freq_taxa
                                                     , flor_df = flor_tidy
                                                     , thresh = most_freq_prop_thresh
                                                     )
                    , n_freq_sites = purrr::map_dbl(clusters
                                                    , sites_with_freq_taxa
                                                    , flor_df = flor_tidy
                                                    , thresh = most_freq_prop_thresh
                                                    )
                    , prop_freq_clusters = n_freq_clusters/groups
                    , prop_freq_sites = n_freq_sites/n_sites
                    ) %>%
      return_result()

    multidplyr::cluster_rm(cl
                           , add_to_cluster
                           )

    rio::export(explore_res$clusters_freq
                , out_file
                )

  }


  #-------cleanup-------

  rm(clusters_use_exp)

  multidplyr::cluster_call(cl, rm(list = ls()))

  rm(cl)

  gc()


  #-------Clusters Explore-------

  return(explore_res)

}


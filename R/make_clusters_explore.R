

#' Explore clusters
#'
#' Explore clusters using silhouette width, within-group sum-of-squares,
#' indicators and widespread taxa.
#'
#' @param clusters_df Dataframe with column 'clusters'.
#' @param cores How many cores to use in [multidplyr::new_cluster()]?
#' @param save_results Logical. If true, results will be saved along the way. If
#' the file already exists, that result will not be recreated.
#' @param out_res Directory into which results are saved.
#' @param obj_list List of objects required:
#' \describe{
#'   \item{dist_flor}{`dist(flor_wide)`.}
#'   \item{dist_env}{`dist(env_wide)`.}
#'   \item{flor_wide}{Wide version of flor_tidy.}
#'   \item{flor_tidy}{Cleaned data set of taxa observations.}
#'   \item{visit_cols}{context.}
#'   \item{p_thresh}{usually 0.05.}
#'   \item{n_sites}{How many unique contexts in `flor_tidy`.}
#'   \item{most_freq_prop_thresh}{Threshold for proportion of sites in a cluster
#'    that must all contain at least one taxa in common.}
#' }
#'
#' @return Each output saved along the way to `out_exp`. Outputs are not made if
#' they have already been saved. Created outputs returned in list with elements:
#' \describe{
#'   \item{clusters_sil}{floristics result from [make_sil_df()].}
#'   \item{clusters_wss}{floristics result from [calc_wss()].}
#'   \item{clusters_sil_env}{env result from [make_sil_df()].}
#'   \item{clusters_wss_env}{env result from [calc_wss()].}
#'   \item{clusters_ind_val}{indicator values for floristics from
#'   [make_ind_val_df()], [clusters_with_indicator()] and
#'   [sites_with_indicator()].}
#'   \item{clusters_freq}{floristics results from [clusters_with_freq_taxa()]
#'   and [sites_with_freq_taxa()]}
#' }
#'
#' @export
#'
#' @examples
make_clusters_explore <- function(clusters_df
                                  , cores = 1
                                  , save_results = TRUE
                                  , out_res = if(save_results) tempdir()
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

  out_res <- fs::path(out_dir)

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

  out_file <- fs::path(out_exp, "clusters_sil.rds")

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

    if(save_results) {

      rio::export(explore_res$clusters_sil
                  , out_file
                  )

    }


  }


  #-------wss--------

  out_file <- fs::path(out_exp, "clusters_wss.rds")

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


    if(save_results) {

      rio::export(explore_res$clusters_wss
                  , out_file
                  )

    }


  }


  #-----silhouette env------

  out_file <- fs::path(out_exp, "clusters_sil_env.rds")

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


    if(save_results) {

      rio::export(explore_res$clusters_sil_env
                  , out_file
                  )

    }

  }


  #-------wss env------

  out_file <- fs::path(out_exp, "clusters_wss_env.rds")

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

    if(save_results) {

      rio::export(explore_res$clusters_wss_env
                  , out_file
                  )

    }

  }


  #------ind val---------

  out_file <- fs::path(out_exp, "clusters_ind_val.rds")

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

    if(save_results) {

      rio::export(explore_res$clusters_ind_val
                  , out_file
                  )

    }

  }


  #-------freq-------

  out_file <- fs::path(out_exp, "clusters_freq.rds")

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


    if(save_results) {

      rio::export(explore_res$clusters_freq
                  , out_file
                  )

    }

  }


  #-------cleanup-------

  rm(clusters_use_exp)

  multidplyr::cluster_call(cl, rm(list = ls()))

  rm(cl)

  gc()


  #-------Clusters Explore-------

  return(explore_res)

}


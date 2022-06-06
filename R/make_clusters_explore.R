

#' Explore clusters
#'
#' Explore clusters using silhouette width, within-group sum-of-squares,
#' indicators and widespread taxa.
#'
#' @param clusters_df Dataframe with column 'clusters'.
#' @param n_sample Number of samples to use to calculate gap statistic.
#' @param cores How many cores to use in [multidplyr::new_cluster()]?
#' @param save_results Logical. If true, results will be saved along the way. If
#' the file already exists, that result will not be recreated.
#' @param do_env Logical. Run make env silhouette, wss and gap.
#' @param do_ind Logical. Make
#' @param out_exp Directory into which results are saved.
#' @param exp_type Name to prefix files with and appended to `out_exp`.
#' @param do_gc Run `gc()` after each output to save memory across `cores`.
#' @param obj_list List of objects required:
#' \describe{
#'   \item{dist_bio}{`dist(bio_wide)`.}
#'   \item{dist_env}{`dist(env_wide)`.}
#'   \item{bio_wide}{Wide version of bio_tidy.}
#'   \item{bio_tidy}{Cleaned data set of taxa observations.}
#'   \item{visit_cols}{context.}
#'   \item{p_thresh}{usually 0.05.}
#'   \item{n_sites}{How many unique contexts in `bio_tidy`.}
#'   \item{most_freq_prop_thresh}{Threshold for proportion of sites in a cluster
#'    that must all contain at least one taxa in common.}
#' }
#'
#' @return Each output saved along the way to `out_exp`. Outputs are not made if
#' they have already been saved. Created outputs returned in list with elements:
#' \describe{
#'   \item{clusters_sil}{biota result from [make_sil_df()].}
#'   \item{clusters_wss}{biota result from [calc_wss()].}
#'   \item{clusters_sil_env}{env result from [make_sil_df()].}
#'   \item{clusters_wss_env}{env result from [calc_wss()].}
#'   \item{clusters_ind_val}{indicator values for biota from
#'   [make_ind_val_df()], [clusters_with_indicator()] and
#'   [sites_with_indicator()].}
#'   \item{clusters_freq}{biota results from [clusters_with_freq_taxa()]
#'   and [sites_with_freq_taxa()]}
#' }
#'
#' @export
#'
#' @examples
make_clusters_explore <- function(clusters_df
                                  , n_sample = 10
                                  , cores = 1
                                  , save_results = TRUE
                                  , do_env = TRUE
                                  , do_ind = TRUE
                                  , out_exp = tempdir()
                                  , exp_type = "cluster"
                                  , do_gc = FALSE
                                  , obj_list
                                  ) {

  # create directory
  if(save_results) {

    fs::dir_create(fs::path(out_exp, exp_type))

  }

  # Variables required come in via obj_list
  purrr::walk2(names(obj_list)
               , obj_list
               , assign
               )

  do_env <- if(exists("dist_env") & do_env) TRUE else FALSE

  return_result <- function(df) {

    df %>%
      dplyr::select(-clusters) %>%
      dplyr::collect() %>%
      dplyr::arrange(id)

  }


  #-------deal with func_cl-------

  use_cores <- min(cores, nrow(clusters_df))

  func_cl <- multidplyr::new_cluster(use_cores)


  #---------Start explore-----------

  # load envCluster to each core
  multidplyr::cluster_library(func_cl
                              , c("envCluster")
                              )


  clusters_use_exp <- clusters_df %>%
    dplyr::select(method, groups, clusters) %>%
    dplyr::mutate(id = row_number())

  explore_res <- list()

  #------silhouette-----

  out_file <- fs::path(out_exp,  exp_type, paste0(exp_type, "_sil.rds"))

  if(!file.exists(out_file)) {

    multidplyr::cluster_copy(func_cl
                             , c("dist_bio"
                                 , "exp_type"
                                 )
                             )

    explore_res$sil <- clusters_use_exp %>%
      multidplyr::partition(func_cl) %>%
      dplyr::mutate(sil = purrr::map(clusters
                                     , make_sil_df
                                     , dist_obj = dist_bio
                                     )
                    , macro_sil = purrr::map_dbl(sil
                                                 , ~mean(.$sil_width)
                                                 )
                    ) %>%
      return_result()

    multidplyr::cluster_rm(func_cl
                           , c("dist_bio")
                           )

    if(do_gc) multidplyr::cluster_call(func_cl, gc())

    if(save_results) {

      rio::export(explore_res$sil
                  , out_file
                  )

    }

  }


  #-------wss--------

  out_file <- fs::path(out_exp,  exp_type, paste0(exp_type, "_wss.rds"))

  if(!file.exists(out_file)) {

    dist_bio_mat <- as.matrix(dist_bio)

    multidplyr::cluster_copy(func_cl
                             , c("dist_bio_mat")
                             )

    explore_res$wss <- clusters_use_exp %>%
      multidplyr::partition(func_cl) %>%
      dplyr::mutate(wss = purrr::map(clusters
                                     , calc_wss
                                     , dist_mat = dist_bio_mat
                                     , clust_col = !!rlang::ensym(exp_type)
                                     )
                    , macro_wss = purrr::map_dbl(wss
                                                 , ~sum(.$wss)
                                                 )
                    ) %>%
      return_result()


    explore_res$gap <- tibble::tibble(run = 1:n_sample) %>%
      dplyr::left_join(clusters_use_exp %>%
                         dplyr::left_join(tibble::as_tibble(explore_res$wss))
                       , by = character()
                       )  %>%
      multidplyr::partition(func_cl) %>%
      dplyr::mutate(wss_for_gap = purrr::map(clusters
                                             , calc_wss
                                             , dist_mat = dist_bio_mat
                                             , clust_col = !!rlang::ensym(exp_type)
                                             , do_sample = TRUE
                                             )
                    ) %>%
      return_result() %>%
      tidyr::unnest(cols = c(wss, wss_for_gap)
                    , names_repair = tidyr_legacy
                    ) %>%
      dplyr::mutate(gap = wss1 - wss) %>%
      dplyr::group_by(method, groups, !!rlang::ensym(exp_type)) %>%
      dplyr::summarise(sample_wss = mean(wss1)
                       , gap_se = sd(gap) / sqrt(n())
                       , gap = mean(gap)
                       ) %>%
      dplyr::ungroup() %>%
      tidyr::nest(gap = -c(method, groups)) %>%
      dplyr::mutate(macro_gap = purrr::map_dbl(gap
                                                 , ~sum(.$gap)
                                               )
                    )

    multidplyr::cluster_rm(func_cl
                           , c("dist_bio_mat")
                           )

    if(do_gc) multidplyr::cluster_call(func_cl, gc())


    if(save_results) {

      rio::export(explore_res$wss
                  , out_file
                  )

      rio::export(explore_res$gap
                  , gsub("wss", "gap", out_file)
                  )

    }


  }


  #-----silhouette env------

  if(do_env) {

    out_file <- fs::path(out_exp,  exp_type, paste0(exp_type, "_sil_env.rds"))

    if(!file.exists(out_file)) {

      multidplyr::cluster_copy(func_cl
                               , c("dist_env")
                               )

      explore_res$sil_env <- clusters_use_exp %>%
        multidplyr::partition(func_cl) %>%
        dplyr::mutate(sil_env = purrr::map(clusters
                                           , make_sil_df
                                           , dist_obj = dist_env
                                           )
                      , macro_sil_env = purrr::map_dbl(sil_env
                                                       , ~mean(.$sil_width)
                                                       )
                      ) %>%
        return_result()

      multidplyr::cluster_rm(func_cl
                             , c("dist_env")
                             )

      if(do_gc) multidplyr::cluster_call(func_cl, gc())


      if(save_results) {

        rio::export(explore_res$sil_env
                    , out_file
                    )

      }

    }


    #-------wss env------

    out_file <- fs::path(out_exp,  exp_type, paste0(exp_type, "_wss_env.rds"))

    if(!file.exists(out_file)) {

      dist_env_mat <- as.matrix(dist_env)

      multidplyr::cluster_copy(func_cl
                               , c("dist_env_mat")
                               )

      explore_res$wss_env <- clusters_use_exp %>%
        multidplyr::partition(func_cl) %>%
        dplyr::mutate(wss_env = purrr::map(clusters
                                           , calc_wss
                                           , dist_mat = dist_env_mat
                                           , clust_col = !!rlang::ensym(exp_type)
                                           )
                      , macro_wss_env = purrr::map_dbl(wss_env
                                                       , ~sum(.$wss)
                                                       )
                      ) %>%
        return_result()


      explore_res$gap_env <- tibble::tibble(run = 1:n_sample) %>%
        dplyr::left_join(clusters_use_exp %>%
                           dplyr::left_join(tibble::as_tibble(explore_res$wss_env))
                         , by = character()
                         )  %>%
        multidplyr::partition(func_cl) %>%
        dplyr::mutate(wss_env_for_gap = purrr::map(clusters
                                               , calc_wss
                                               , dist_mat = dist_env_mat
                                               , clust_col = !!rlang::ensym(exp_type)
                                               , do_sample = TRUE
                                               )
                      ) %>%
        return_result() %>%
        tidyr::unnest(cols = c(wss_env, wss_env_for_gap)
                      , names_repair = tidyr_legacy
                      ) %>%
        dplyr::mutate(gap_env = wss1 - wss) %>%
        dplyr::group_by(method, groups, !!rlang::ensym(exp_type)) %>%
        dplyr::summarise(sample_wss_env = mean(wss1)
                         , gap_env_se = sd(gap_env) / sqrt(n())
                         , gap_env = mean(gap_env)
                         ) %>%
        dplyr::ungroup() %>%
        tidyr::nest(gap_env = -c(method, groups)) %>%
        dplyr::mutate(macro_gap_env = purrr::map_dbl(gap_env
                                                 , ~sum(.$gap_env)
                                                 )
                      )

      multidplyr::cluster_rm(func_cl
                             , c("dist_env_mat")
                             )

    }

    if(do_gc) multidplyr::cluster_call(func_cl, gc())

    if(save_results) {

      rio::export(explore_res$wss_env
                  , out_file
                  )

      rio::export(explore_res$gap_env
                  , gsub("wss", "gap", out_file)
                  )

    }

  }


  #------ind val---------

  if(do_ind) {

    out_file <- fs::path(out_exp,  exp_type, paste0(exp_type, "_ind_val.rds"))

    if(!file.exists(out_file)) {

      add_to_cluster <- c("bio_wide"
                          , "bio_tidy"
                          , "visit_cols"
                          , "p_thresh"
                          , "n_sites"
                          )

      multidplyr::cluster_copy(func_cl
                               , add_to_cluster
                               )

      explore_res$ind_val <- clusters_use_exp %>%
        multidplyr::partition(func_cl) %>%
        dplyr::mutate(ind_val = purrr::map(clusters
                                           , make_ind_val_df
                                           , bio_wide = bio_wide
                                           , taxas = unique(bio_tidy$taxa)
                                           , context = visit_cols
                                           , clust_col = exp_type
                                           )
                      , n_ind_clusters = purrr::map_dbl(ind_val
                                                        , clusters_with_indicator
                                                        , thresh = p_thresh
                                                        , clust_col = !!rlang::ensym(exp_type)
                                                        )
                      , n_ind_sites = purrr::map2_dbl(ind_val
                                                      , clusters
                                                      , sites_with_indicator
                                                      , thresh = p_thresh
                                                      , clust_col = !!rlang::ensym(exp_type)
                                                      )
                      , prop_ind_clusters = n_ind_clusters/groups
                      , prop_ind_sites = n_ind_sites/n_sites
                      ) %>%
        return_result()

      multidplyr::cluster_rm(func_cl
                             , add_to_cluster
                             )

      if(do_gc) multidplyr::cluster_call(func_cl, gc())

      if(save_results) {

        rio::export(explore_res$ind_val
                    , out_file
                    )

      }

    }

  }


  #-------freq-------

  out_file <- fs::path(out_exp,  exp_type, paste0(exp_type, "_freq.rds"))

  if(!file.exists(out_file)) {

    add_to_cluster <- c("bio_tidy"
                        , "most_freq_prop_thresh"
                        , "n_sites"
                        )

    multidplyr::cluster_copy(func_cl
                             , add_to_cluster
                             )

    explore_res$freq <- clusters_use_exp %>%
      multidplyr::partition(func_cl) %>%
      dplyr::mutate(n_freq_clusters = purrr::map_dbl(clusters
                                                     , clusters_with_freq_taxa
                                                     , bio_df = bio_tidy
                                                     , thresh = most_freq_prop_thresh
                                                     , clust_col = !!rlang::ensym(exp_type)
                                                     )
                    , n_freq_sites = purrr::map_dbl(clusters
                                                    , sites_with_freq_taxa
                                                    , bio_df = bio_tidy
                                                    , thresh = most_freq_prop_thresh
                                                    , clust_col = !!rlang::ensym(exp_type)
                                                    )
                    , prop_freq_clusters = n_freq_clusters/groups
                    , prop_freq_sites = n_freq_sites/n_sites
                    ) %>%
      return_result()

    multidplyr::cluster_rm(func_cl
                           , add_to_cluster
                           )

    if(do_gc) multidplyr::cluster_call(func_cl, gc())


    if(save_results) {

      rio::export(explore_res$freq
                  , out_file
                  )

    }

  }


  #-------cleanup-------

  rm(clusters_use_exp)

  multidplyr::cluster_call(func_cl, rm(list = ls()))

  rm(func_cl)

  gc()


  #-------Clusters Explore-------

  if(length(explore_res) > 0) {

    names(explore_res) <- paste0(exp_type, "_", names(explore_res))

  }

  return(explore_res)

}


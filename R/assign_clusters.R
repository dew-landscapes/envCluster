#' Assign cluster memberships to sites
#'
#' Biological data is hellinger transformed before distance calculations using
#' `vegan::decostand(bio, method = "hellinger")`.
#'
#' @details
#' Provide a 'method' via `...` to specify one of: `cluster::agnes()` using
#' a value that is valid for the method argument of that function;
#' `genieclust::gclust()` by `method = "gclust"`; `dbscan::hdbscan()`
#' by `method = "hdbscan"`; `stats::kmeans()` by `method = "kmeans"`; or
#' `mclust::Mclust()` by `method = "mclust"`. Note that `stats::kmeans()` and
#' `mclust::Mclust()` are implemented with specific arguments that seem sensible
#' for large-scale floristic clustering.
#'
#' @param bins Dataframe with bin identifiers. For the same sites in the same
#' order as `bio` and, if used, `env`
#' @param bio Biological data as site * species
#' @param env Environmental data as site * environment (in the same order and
#' for the same sites as bio)
#' @param env_weight Numeric. Value >= 0 and <= 1 used as `weights`
#' argument to `analogue::fuse()`.
#' @param ... Passed to: `cluster::agnes()`; `genieclust::gclust()`
#' (use `method = "genie"`); or `dbscan::hdbscan()` (use `method = "dbscan"`).
#'
#' @return A tibble of bins and their associated cluster (as a list column
#' `clusters`) for each group size requested.
#' @export
#'
#' @example
assign_clusters <- function(bins
                            , bio
                            , env = NULL
                            , env_weight = 0.1
                            , group_range = 2:50
                            , ...
                            ) {

  stopifnot(identical(nrow(bins), nrow(bio)))
  if(!is.null(env)) stopifnot(identical(nrow(bio), nrow(env)))

  dots_list <- list(...)

  if(is.null(dots_list$method)) dots_list$method <- "hdbscan"

  # methods ------
  agnes_methods <- c("average", "single", "complete", "ward", "weighted", "flexible", "gaverage")

  stopifnot(dots_list$method %in% c(agnes_methods, "gclust", "hdbscan", "kmeans", "mclust"))

  hell <- vegan::decostand(bio, method = "hellinger")

  # dist based -------
  if(dots_list$method %in% c(agnes_methods, "gclust", "hdbscan")) {

    dbio <- hell |>
      stats::dist()

    denv <- stats::dist(base::scale(env))

    if(env_weight > 0 && env_weight < 1) {

      d <- analogue::fuse(dbio
                          , denv
                          , weights = c(1 - env_weight, env_weight)
                          )

    } else if(env_weight == 0) {

      d <- dbio

    } else if(env_weight == 1) {

      d <- denv

    }

    ## agnes ------
    if(dots_list$method %in% agnes_methods) {

      if(is.null(dots_list$par.method)) dots_list$par.method <- 0.625
      # alpha = 0.625 is used to set beta at -0.25, which is often cited as a good default specifically for ecological community data

      dend <- R.utils::doCall(cluster::agnes
                              , args = dots_list
                              , x = d
                              )

    }

    ## gclust ------
    if(dots_list$method == "gclust") {

      if(!is.null(dots_list$cores)) quitefastmst::omp_set_num_threads(n_threads = dots_list$cores)
      dend <- genieclust::gclust(d)

    }

    ## hdbscan ------
    if(dots_list$method == "hdbscan") {

      dend <- R.utils::doCall(dbscan::hdbscan
                              , x = d
                              , args = dots_list
                              )$hc

    }

    result <- make_clusters(dend
                            , group_range = group_range
                            , bins = bins
                            )

  } else if(dots_list$method %in% c("kmeans", "mclust")) {

    # blended --------
    blend <- function(hell
                      , env
                      , weighting
                      , bio_axes = 30
                      , env_axes = 5
                      ) {

      bio_pca <- irlba::prcomp_irlba(hell, n = bio_axes, center = TRUE, scale. = FALSE)$x
      env_pca <- irlba::prcomp_irlba(base::scale(env), n = env_axes)$x

      tot_var <- function(M) sum(apply(M, 2, var))
      bio_s <- bio_pca / sqrt(tot_var(bio_pca))
      env_s <- env_pca / sqrt(tot_var(env_pca))

      b <- function(w) cbind(sqrt(w) * bio_s, sqrt(1 - w) * env_s)

      blended <- b(weighting)

    }

    blended <- blend(hell
                     , env
                     , weighting = env_weight
                     , bio_axes = 25
                     , env_axes = 5
                     )

    ## kmeans -------
    if(dots_list$method == "kmeans") {

    create_kmeans <- function(blended
                              , use_group
                              ) {

      fit <- stats::kmeans(blended
                           , centers = use_group
                           , nstart = 25
                           , iter.max = 100
                           )

      clusts <- bins |>
        dplyr::bind_cols(tibble::enframe(fit$cluster, name = NULL, value = "clust")) |>
        dplyr::mutate(cluster = envFunc::numbers2words(clust)
                      , cluster = forcats::fct_reorder(cluster
                                                       , clust
                                                       )
                      )

    }

    result <- tibble::enframe(group_range
                              , name = NULL
                              , value = "groups"
                              ) |>
      dplyr::mutate(clusters = purrr::map(groups
                                          , \(x) create_kmeans(blended
                                                               , use_group = x
                                                               )
                                          )
                    )

    }

    ## mclust ------
    if(dots_list$method == "mclust") {

      create_mclust <- function(blended
                                , use_group
                                ) {

        fit <- mclust::Mclust(blended
                              , modelNames = c("EII", "VII", "EEI", "VEI")
                              , G = use_group
                              )

        clusts <- bins |>
          dplyr::bind_cols(tibble::enframe(fit$classification, name = NULL, value = "clust")) |>
          dplyr::mutate(cluster = envFunc::numbers2words(clust)
                        , cluster = forcats::fct_reorder(cluster
                                                         , clust
                                                         )
                        , probability = fit$uncertainty
                        , model_name = fit$modelName
                        )

      }

      if(is.null(dots_list$cores)) dots_list$cores <- 1

      RhpcBLASctl::blas_set_num_threads(dots_list$cores)
      RhpcBLASctl::omp_set_num_threads(dots_list$cores)

      result <- tibble::enframe(group_range
                                , name = NULL
                                , value = "groups"
                                ) |>
        dplyr::mutate(clusters = purrr::map(groups
                                            , \(x) create_mclust(blended
                                                                 , use_group = x
                                                                 )
                                            )
                      )

    }

  }

  return(result)

}

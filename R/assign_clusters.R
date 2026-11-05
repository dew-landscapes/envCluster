#' Assign cluster memberships to sites
#'
#' Biological data is hellinger transformed before distance calculations using
#' `vegan::decostand(bio, method = "hellinger")`
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

  # methods ------
  agnes_methods <- c("average", "single", "complete", "ward", "weighted", "flexible", "gaverage")

  stopifnot(dots_list$method %in% c(agnes_methods, "gclust", "hdbscan"))

  dbio <- vegan::decostand(bio, method = "hellinger") |>
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

  if(dots_list$method %in% agnes_methods) {

    if(is.null(dots_list$par.method)) dots_list$par.method <- 0.625
    # alpha = 0.625 is used to set beta at -0.25, which is often cited as a good default specifically for ecological community data

    dend <- R.utils::doCall(cluster::agnes
                            , args = dots_list
                            , x = d
                            )

  }

  if(dots_list$method == "gclust") {

    if(!is.null(dots_list$cores)) quitefastmst::omp_set_num_threads(n_threads = dots_list$cores)
    dend <- genieclust::gclust(d)

  }

  if(dots_list$method == "hdbscan") {

    dend <- R.utils::doCall(dbscan::hdbscan
                            , x = d
                            , args = dots_list
                            )$hc

  }

  result <- make_clusters(dend
                          , group_range = group_range
                          , clust_col = "clust"
                          , clust_col_out = "cluster"
                          , bins = bins
                          )

  return(result)

}

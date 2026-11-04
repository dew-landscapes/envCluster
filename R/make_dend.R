#' Make a dendogram
#'
#' @param dist_bio Distance object, from, say, `vegan::vegdist()`,
#' `stats::dist()` or `parallelDist::parDist()`.
#' @param dist_env Distance object. Usually environmental distance. Needs to be
#' for the same bins as used to create `dist_bio`. Only needed if `env_weight`
#' is greater than 0.
#' @param env_weight Numeric. Value above 0 and below 1 used as `weights`
#' argument to `analogue::fuse()`.
#' @param ... Passed to `cluster::agnes()`
#'
#' @return An object of class 'hclust'.
#' @export
#'
#' @example inst/examples/make_clusters_ex.R
make_dend <- function(dist_bio
                      , dist_env = NULL
                      , env_weight = 0.1
                      , ...
                      ) {

  dots_list <- list(...)

  if(is.null(dots_list$par.method)) dots_list$par.method = 0.625
  # alpha = 0.625 is used to set beta at -0.25, which is often cited as a good default specifically for ecological community data

  if(grepl("bray", base::attr(dist_bio, "method")) && grepl("ward", dots_list$method)) {

    dist_bio <- sqrt(dist_bio)

  }

  if(env_weight != 0) {

    d <- analogue::fuse(dist_bio
                        , dist_env
                        , weights = c(1 - env_weight, env_weight)
                        )

  } else d <- dist_bio

  dend <- R.utils::doCall(cluster::agnes
                          , dots_list
                          , x = d
                          ) |>
    stats::as.hclust()

  return(dend)

}

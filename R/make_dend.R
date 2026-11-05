#' Make a dendogram
#'
#' Biological data is hellinger transformed before distance calculations using
#' `vegan::decostand(bio, method = "hellinger")`
#'
#' @param bio Biological data as site * species
#' @param env Environmental data as site * environment (in same order for the
#' same sites as bio)
#' @param env_weight Numeric. Value >= 0 and <= 1 used as `weights`
#' argument to `analogue::fuse()`.
#' @param ... Passed to either `cluster::agnes()` or `genieclust::gclust()`. To
#' use `genieclust::gclust()` provide `method = "genie"`
#'
#' @return An object of class 'hclust'.
#' @export
#'
#' @example inst/examples/make_clusters_ex.R
make_dend <- function(bio
                      , env = NULL
                      , env_weight = 0.1
                      , ...
                      ) {

  if(!is.null(env)) stopifnot(identical(nrow(bio), nrow(env)))

  dots_list <- list(...)

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

  if(dots_list$method != "genie") {

    if(is.null(dots_list$par.method)) dots_list$par.method <- 0.625
    # alpha = 0.625 is used to set beta at -0.25, which is often cited as a good default specifically for ecological community data

    dend <- R.utils::doCall(cluster::agnes
                            , dots_list
                            , x = d
                            ) |>
      stats::as.hclust()

  }

  if(dots_list$method == "genie") {

    dend <- genieclust::gclust(d)

  }

  return(dend)

}

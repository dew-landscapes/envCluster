#' Make a dendogram
#'
#' @param method What method to use.
#' @param dist_bio Distance object, from, say, `vegan::vegdist()`,
#' `stats::dist()` or `parallelDist::parDist()`.
#' @param dist_env Distance object. Usually environmental distance. Needs to be
#' for the same sites as used to create `dist_bio`. Only needed if method is
#' "geo".
#' @param geo_alpha Numeric. Value above 0 and below 1 used as `alpha` argument
#' to `ClustGeo::hclustgeo()`.
#'
#' @return An object of class 'hclust'.
#' @export
#'
#' @example inst/examples/make_clusters_ex.R
make_dend <- function(method
                      , dist_bio
                      , dist_env = NULL
                      , geo_alpha = 0.1
                      ) {

  hclust_methods <- c("single", "complete", "average", "mcquitty", "ward.D", "centroid", "median", "ward.D2")
  geo_methods <- "geo"

  dend <- if(any(grepl(method, hclust_methods))) {

    fastcluster::hclust(dist_bio
                        , method
                        )

  } else if(any(grepl(method, geo_methods))) {

    ClustGeo::hclustgeo(dist_bio
                        , dist_env
                        , alpha = geo_alpha
                        )
  } else stop("method (", method, ") should be one of: ", paste0(c(hclust_methods, geo_methods), collapse = ", "))

  return(dend)

}

#' Summaries of a clustering.
#'
#' Returns measures of inequality in cluster sizes. Package `genieclust` is used
#' to generate the gini, bonferroni and devergottini indices which are
#' 'rescaled' (as 1 - value) here so that 1 is perfectly equal.
#'
#' @param cl Vector of cluster membership
#' @param min_abs Numeric. Clusters containing less than `min_abs` bins are
#' deemed 'small'
#' @param min_prop Numeric. Clusters containing less than
#' `min_prop * length(cl) / length(unique(cl))` bins are also deemed 'small'
#'
#' @return single row tibble with summary information. `_rescale` cols are
#' `1 - value` so that 1 (higher) values represent even/equal clusters.
#' @export
#'
#' @example inst/examples/make_clusters_ex.R
make_summary_df <- function(cl
                            , min_abs = 10
                            , min_prop = 0.1
                            ) {

  n  <- length(cl)
  tab_orig <- as.numeric(table(cl))
  tab <- tab_orig[tab_orig > 0]
  k  <- length(tab)
  p  <- tab / n
  H  <- -sum(p * log(p))
  small <- tab < max(min_abs, min_prop * n / k)

  tibble::tibble(orig_k = length(tab_orig)
                 , k = k
                 , prop_top1 = max(p)
                 , prop_top2 = sum(sort(p, decreasing = TRUE)[1:min(2, k)])
                 , n_small = sum(small)
                 , prop_small_clusters = mean(small)
                 , prop_small_bins = sum(tab[small]) / n
                 , gini_rescale = 1 - genieclust::gini_index(cl)
                 , bonf_rescale = 1 - genieclust::bonferroni_index(cl)
                 , deve_rescale = 1 - genieclust::devergottini_index(cl)
                 )

}

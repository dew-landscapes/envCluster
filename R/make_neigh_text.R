#' Make text to describe silhouette.
#'
#' @param this_clust Name of cluster.
#' @param sil_df Output from `make_sil_df()`
#'
#' @return Sentence for use in silhouette section of report.
#' @export
#'
#' @examples
make_neigh_text <- function(this_clust, sil_df) {

  this_clust <- as.character(this_clust)

  df <- sil_df %>%
    dplyr::mutate(sites = nrow(.)) %>%
    dplyr::group_by(neighbour,sites) %>%
    dplyr::summarise(neighbours = n()
                     , silMean = mean(sil_width)
                     ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(prop = neighbours/sites)

  dfPropMax <- max(df$prop,na.rm = TRUE) > 0.25

  df <- df %>%
    {if(dfPropMax) filter(.,prop > 0.25) else top_n(.,3,prop)} %>%
    dplyr::mutate(text = paste0(numbers2words(neighbour)
                                , " ("
                                , round(100*prop,1)
                                , "% of sites with mean silhouette width of "
                                , round(silMean,2)
                                , ")"
    )
    )

  paste0("Sites in cluster "
         , this_clust
         , " were most frequently neighbours to sites in "
         , if(nrow(df) == 1) "cluster " else "clusters "
         , df$text %>% vec_to_sentence()
         ,"."
  )

}

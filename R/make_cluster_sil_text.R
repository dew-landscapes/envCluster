#' Make text to describe silhouette.
#'
#' @param cluster Name of cluster.
#' @param sil_df Output from `make_sil_df()`
#'
#' @return Sentence for use in silhouette section of report.
#' @export
#'
#' @examples
cluster_sil_text <- function(cluster,sil_df) {

  cluster <- as.character(cluster)

  df <- sil_df %>%
    dplyr::add_count(cluster, name = "sites") %>%
    dplyr::group_by(cluster,neighbour,sites) %>%
    dplyr::summarise(neighbours = n()
                     , silMean = mean(sil_width)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(prop = neighbours/sites) %>%
    dplyr::filter(cluster == cluster)

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
         , unique(df$cluster)
         , " were most frequently neighbours to sites in "
         , if(nrow(df) == 1) "cluster " else "clusters "
         , df$text %>% vec_to_sentence()
         ,"."
  )

}

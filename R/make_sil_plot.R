
#' Make a silhouette plot from a silhouette dataframe
#'
#' @param sil_df Silhouette dataframe. See `make_sil_df()`.
#' @param clust_col Character name of column in `sil_df` containing the clusters.
#' @param clust_colours If null (default) `viridis::viridis` is used to colour
#' clusters. Otherwise a vector of colours the same length as levels in `clust_col`.
#' @param include_labels Logical. Include cluster labels on plot? Can get very
#' crowded if there are many clusters.
#'
#' @return ggplot
#' @export
#'
#' @examples
make_sil_plot <- function(sil_df
                          , clust_col = "cluster"
                          , clust_colours = NULL
                          , include_labels = length(unique(sil_df[clust_col][[1]])) < 30
                          ) {

  df <- sil_df %>%
    dplyr::arrange(!!ensym(clust_col),desc(sil_width)) %>%
    dplyr::mutate(row = row_number()
                  , neigh = numbers2words(neighbour)
                  , neigh = factor(neigh,levels=levels(!!ensym(clust_col)))
                  , neigh = if_else(sil_width > 0
                                    , !!ensym(clust_col)
                                    , neigh
                                    )
                  ) %>%
    dplyr::group_by(!!ensym(clust_col)) %>%
    dplyr::mutate(mid = mean(row)) %>%
    dplyr::ungroup()

  clust_labs <- df %>%
    dplyr::count(!!ensym(clust_col),mid)

  mean_sil <- round(mean(sil_df$sil_width),2)

  if(isTRUE(is.null(clust_colours))) {

    clust_colours <- viridis::viridis(n = length(levels(df$cluster)))

  }

  sil_plot <- ggplot() +
    geom_col(data = df
             , aes(row
                   , sil_width
                   , fill=neigh
                   , color=neigh
                   )
             ) +
    coord_flip() +
    labs(subtitle = paste0("Negative silhouette widths are coloured by neighbouring cluster\nMean silhouette width = "
                           , mean_sil
                           , ". n = "
                           , nrow(df)
                           , " patches"
                           )
         , y = "Silhouette width"
         ) +
    scale_x_reverse() +
    theme(axis.text.y=element_blank()
          , axis.title.y=element_blank()
          , axis.ticks.y=element_blank()
          , legend.position = "none"
          ) +
    scale_colour_manual(values = clust_colours) +
    scale_fill_manual(values = clust_colours)

  if(include_labels) {

    sil_plot <- sil_plot +
      ggrepel::geom_label_repel(data = clust_labs
                       , aes(mid
                             , 0
                             , label= !!ensym(clust_col)
                             , fill= !!ensym(clust_col)
                             )
                       , nudge_y = min(df$sil_width)
                       )

    }

  return(sil_plot)

}

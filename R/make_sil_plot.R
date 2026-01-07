
#' Make a silhouette plot from a silhouette dataframe
#'
#' @param sil_df Silhouette dataframe. See `make_sil_df()`.
#' @param clust_col Character name of column in `sil_df` containing the clusters.
#' @param clust_colours If null (default) `viridis::viridis` is used to colour
#' clusters. Otherwise a dataframe with clusters and their associated colours.
#' Needs to be the same length as levels in `clust_col` (or `levs`).
#' @param include_labels Logical. Include cluster labels on plot? Can get very
#' crowded if there are many clusters.
#'
#' @return ggplot object
#' @export
#'
#' @example inst/examples/make_clusters_ex.R
make_sil_plot <- function(sil_df
                          , clust_col = "cluster"
                          , clust_colours = NULL
                          , include_labels = length(unique(sil_df[clust_col][[1]])) < 30
                          ) {

  levs <- if(isTRUE(is.null(clust_colours))) {

    levels(sil_df[clust_col][[1]])

  } else {

    levels(clust_colours[clust_col][[1]])

  }

  df <- sil_df %>%
    dplyr::mutate(!!rlang::ensym(clust_col) := factor(!!rlang::ensym(clust_col)
                                                      , levels = levs
                                                      )
                  ) %>%
    dplyr::arrange(!!rlang::ensym(clust_col), dplyr::desc(sil_width)) %>%
    dplyr::mutate(row = dplyr::row_number()
                  , neigh = envFunc::numbers2words(neighbour)
                  , neigh = factor(neigh
                                   , levels = levs
                                   )
                  , neigh = dplyr::if_else(sil_width > 0
                                           , !!rlang::ensym(clust_col)
                                           , neigh
                                           )
                  ) %>%
    dplyr::group_by(!!rlang::ensym(clust_col)) %>%
    dplyr::mutate(mid = mean(row)) %>%
    dplyr::ungroup()

  mean_sil <- round(mean(sil_df$sil_width), 2)

  if(isTRUE(is.null(clust_colours))) {

    clust_colours <- tibble::tibble(!!rlang::ensym(clust_col) := levs) %>%
      dplyr::mutate(!!rlang::ensym(clust_col) := factor(!!rlang::ensym(clust_col)
                                                        , levels = levs
                                                        )
                    , colour = viridis::viridis(n = nrow(.))
                    )

  } else {

    clust_colours <- clust_colours %>%
      dplyr::mutate(!!rlang::ensym(clust_col) := factor(!!rlang::ensym(clust_col)
                                                        , levels = levs
                                                        )
                    )

  }

  clust_labs <- df %>%
    dplyr::count(!!rlang::ensym(clust_col), mid) %>%
    dplyr::left_join(clust_colours)

  df_plot <- df %>%
    dplyr::left_join(clust_colours
                     , by = c("neigh" = clust_col)
                     )

  sil_plot <- ggplot2::ggplot() +
    ggplot2::geom_col(data = df_plot
                      , ggplot2::aes(row
                                     , sil_width
                                     , fill = colour
                                     , color = colour
                                     )
                      ) +
    ggplot2::coord_flip() +
    ggplot2::labs(subtitle = paste0("Negative silhouette widths are coloured by neighbouring cluster\nMean silhouette width = "
                                    , mean_sil
                                    , ". n = "
                                    , nrow(df)
                                    , " bins"
                                    )
                  , y = "Silhouette width"
                  ) +
    ggplot2::scale_x_reverse() +
    ggplot2::theme(axis.text.y = ggplot2::element_blank()
                   , axis.title.y = ggplot2::element_blank()
                   , axis.ticks.y = ggplot2::element_blank()
                   , legend.position = "none"
                   ) +
    ggplot2::scale_colour_identity() +
    ggplot2::scale_fill_identity()

  if(include_labels) {

    sil_plot <- sil_plot +
      ggrepel::geom_label_repel(data = clust_labs
                                , ggplot2::aes(mid
                                               , 0
                                               , label = !!rlang::ensym(clust_col)
                                               , fill = colour
                                               )
                                , nudge_y = min(df$sil_width)
                                )

    }

  return(sil_plot)

}

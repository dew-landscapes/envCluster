

#' Plot the results from diagnostic_df
#'
#' @param diag_df Dataframe with results from diagnostic_df.
#' @param group_col Character. Name of column with numeric data for x-axis.
#' @param facet_col Character. Optional name of column to facet on (along with label).
#' @param label Character label for the diagnostics (choose another column from diagnosticDF).
#' @param value Charcter. Name of column with values to plot.
#' @param display_all Display all diagnostics or only those used to select best.
#'
#' @return
#' @export
#'
#' @examples
diagnostic_plot <- function(diag_df
                            , group_col = "groups"
                            , facet_col = "method"
                            , label = "diagnostic"
                            , value = "combo"
                            , display_all = FALSE
                            ) {

  df <- diag_df %>%
    {if(display_all) (.) else (.) %>% dplyr::filter(weight)} %>%
    dplyr::mutate(across(where(is.factor),factor)) %>%
    dplyr::filter(!is.na(value))


  facet_formula <- if(isTRUE(!is.null(facet_col))) {

    as.formula(paste0(label,"~",facet_col))

    } else label


  ggplot2::ggplot(df
         ,aes(!!ensym(group_col)
              , scale
              , colour = !!ensym(value)
              , alpha = if(display_all) weight else NULL
              , label = !!ensym(group_col)
              , size = top
              )
         ) +
    ggplot2::geom_point() +
    ggrepel::geom_text_repel(data = df %>%
                               dplyr::filter(best)
                             , size = 2
                             , show.legend = FALSE
                             , box.padding = 1
                             , min.segment.length = 0
                             , colour = "black"
                             ) +
    ggplot2::facet_grid(facet_formula
                        , scales="free_y"
                        , labeller = label_wrap_gen(20,multi_line = TRUE)
                        ) +
    ggplot2::labs(colour = "Combination"
         , alpha = "Diagnostic used" #paste0("Top ",unique(diagnosticdf$topThresh)*100,"%")
         , title = paste0("Labels indicate top ",numbers2words(unique(df$best_thresh))," results")
         , size = paste0("Best ",(unique(df$top_thresh))*100,"%")
         ) +
    ggplot2::scale_colour_viridis_c() +
    ggplot2::scale_x_continuous(breaks = scales::pretty_breaks()) +
    ggplot2::theme(strip.text.y = element_text(angle = 0)
          , strip.text.x = element_text(angle = 90)
          , axis.text.x = element_text(angle = 90)
          )

}



#' Plot the results from diagnostic_df
#'
#' @param diagnosticdf Dataframe with results from diagnostic_df.
#' @param label Character label for the diagnostics (choose another column from diagnosticDF).
#' @param displayall Display all diagnostics or only those used to select best.
#'
#' @return
#' @export
#'
#' @examples
diagnostic_plot <- function(diagnosticdf
                            , label = "diagnostic"
                            , displayall = FALSE
                            ) {

  df <- diagnosticdf %>%
    {if(displayall) (.) else (.) %>% dplyr::filter(weight)} %>%
    dplyr::mutate(across(where(is.factor),factor)) %>%
    dplyr::filter(!is.na(value))

  ggplot2::ggplot(df
         ,aes(groups
              , scale
              , colour = combo
              , alpha = if(displayall) weight else NULL
              , label = groups
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
    ggplot2::facet_grid(as.formula(paste0(label,"~method"))
               , scales="free_y"
               ,  labeller = label_wrap_gen(20,multi_line = TRUE)
               ) +
    ggplot2::labs(colour = "Combination"
         , alpha = "Diagnostic used" #paste0("Top ",unique(diagnosticdf$topThresh)*100,"%")
         , title = paste0("Labels indicate top ",numbers2words(unique(df$bestThresh))," results")
         , size = paste0("Best ",(unique(df$topThresh))*100,"%")
         ) +
    ggplot2::scale_colour_viridis_c() +
    ggplot2::scale_x_continuous(breaks = scales::pretty_breaks()) +
    ggplot2::theme(strip.text.y = element_text(angle = 0)
          , strip.text.x = element_text(angle = 90)
          , axis.text.x = element_text(angle = 90)
          )

}

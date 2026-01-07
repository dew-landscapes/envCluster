#' Make a wide (usually bin * taxa) data frame
#'
#' @param bio_df Dataframe containing the bin and taxa data in long format.
#' @param context Character. Name of columns defining context.
#' @param taxa_col Character. Name of column containing taxa.
#' @param num_col Character. Name of column containing numeric abundance data (
#' usually 'cover' for plants).
#' @param num_col_NA Value to use to fill NA in wide table.
#'
#' @return wide format version of `bio_df`
#' @export
#'
#' @examples
make_wide_df <- function(bio_df
                         , context
                         , taxa_col = "taxa"
                         , num_col = "use_cover"
                         , num_col_NA = 0
                         ) {

  bio_df |>
    dplyr::group_by(dplyr::across(tidyselect::any_of(c(taxa_col, context)))) |>
    dplyr::summarise(value = max(!!rlang::ensym(num_col), na.rm = TRUE)) |>
    dplyr::ungroup() |>
    dplyr::filter(!is.na(value)) |>
    tidyr::pivot_wider(names_from = !!rlang::ensym(taxa_col)
                       , values_fill = num_col_NA
                       )

}

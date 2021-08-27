

#' Use a set of (continuous) columns to choose a good set of candidate rows
#'
#' @param df Dataframe with columns over which to find good candidates.
#' @param df_with_diag Dataframe mapping the name of possible diagnostics to cases (columns) in which to use those diagnostic.
#' @param context Character. Name of columns in df that define context.
#' @param diag_col Character. Name of df_with_diag column to use in this instance.
#' @param summarise_method Character. Name of method to use in summarising if there is more than one row per context.
#' @param group_col Character. Optional. Name of column to filter on (min_groups, max_groups)
#' @param min_groups Numeric. Minimum number of clusters in a clustering.
#' @param max_groups Numeric. Maximum number of clusters in a clustering.
#' @param top_thresh Numeric specifying the proportion of rows considered 'top'.
#' @param best_thresh Numeric specifying the absolute number of rows considered 'best'.
#'
#' @return
#' @export
#'
#' @examples
diagnostic_df <- function(df
                          , df_with_diag = tibble::tibble(diagnostic = "av_clust_size",high_good = TRUE, clust_sum = TRUE)
                          , context = c("method","groups")
                          , diag_col = "clust_sum"
                          , summarise_method = median
                          , group_col = "groups"
                          , min_groups = NULL
                          , max_groups = NULL
                          , top_thresh = 0.25
                          , best_thresh = 5
                          ) {

  if(!isTRUE(is.null(min_groups))) {

    df <- df %>%
      dplyr::filter(!!ensym(group_col) >= min_groups)

  }

  if(!isTRUE(is.null(max_groups))) {

    df <- df %>%
      dplyr::filter(!!ensym(group_col) <= max_groups)

  }

  df_with_diag <- df_with_diag %>%
    dplyr::mutate(diagnostic = fct_inorder(diagnostic))

  df %>%
    dplyr::select(all_of(context),any_of(df_with_diag$diagnostic)) %>%
    #select_if(~ !all(is.na(.))) %>%
    dplyr::group_by(across(all_of(context))) %>%
    dplyr::summarise(across(any_of(df_with_diag$diagnostic),summarise_method)) %>%
    dplyr::ungroup() %>%
    tidyr::pivot_longer(any_of(df_with_diag$diagnostic),names_to = "diagnostic", values_to = "value") %>%
    dplyr::inner_join(df_with_diag) %>%
    dplyr::group_by(across(any_of(names(df_with_diag)))) %>%
    dplyr::mutate(scale = if_else(high_good
                                  ,scales::rescale(value,to=c(0,1))
                                  ,scales::rescale(desc(value),to=c(0,1))
                                  )
                  ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(combo_init = scale*!!ensym(diag_col)) %>%
    dplyr::group_by(across(all_of(context)),across(!!ensym(diag_col))) %>%
    dplyr::mutate(combo = prod(combo_init)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(across(all_of(context))) %>%
    dplyr::mutate(combo = max(combo, na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(combo = scales::rescale(combo,to=c(0,1))) %>%
    dplyr::mutate(top_thresh = top_thresh
                  , best_thresh = best_thresh
                  , top = combo >= quantile(combo,probs = 1-top_thresh,na.rm = TRUE)
                  , top = if_else(is.na(top),FALSE,top)
                  , best = combo >= sort(unique(.$combo),TRUE)[best_thresh]
                  , best = if_else(is.na(best),FALSE,best)
                  , diagnostic = factor(diagnostic, levels = levels(df_with_diag$diagnostic))
                  , weight = !!ensym(diag_col)
                  )

}

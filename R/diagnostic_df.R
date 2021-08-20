

#' Use a set of (continuous) columns to choose a good set of candidate rows
#'
#' @param df Dataframe with columns over which to find good candidates.
#' @param diagnostic_df Dataframe mapping the name of possible diagnostics to cases (columns) in which to use those diagnostic.
#' @param context Character. Name of columns in df that define context.
#' @param diag_col Character. Name of diagnostic_df column to use in this instance.
#' @param summarise_method Character. Name of method to use in summarising if there is more than one row per context.
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
                          , diagnostic_df = tibble::tibble(diagnostic = "av_clust_size",high_good = TRUE, clust_sum = TRUE)
                          , context = c("method","groups")
                          , diag_col = "clust_sum"
                          , summarise_method = median
                          , min_groups = 10
                          , max_groups = Inf
                          , top_thresh = 0.25
                          , best_thresh = 5
                          ) {

  df %>%
    dplyr::filter(groups >= min_groups
                  , groups <= max_groups
                  ) %>%
    dplyr::select(all_of(context),groups,any_of(diagnostics$diagnostic)) %>%
    #select_if(~ !all(is.na(.))) %>%
    dplyr::group_by(across(all_of(context))) %>%
    dplyr::summarise(across(any_of(diagnostics$diagnostic),summarise_method)) %>%
    dplyr::ungroup() %>%
    tidyr::pivot_longer(any_of(diagnostics$diagnostic),names_to = "diagnostic", values_to = "value") %>%
    dplyr::left_join(diagnostics) %>%
    dplyr::group_by(diagnostic,definition) %>%
    dplyr::mutate(scale = if_else(high_good
                                  ,scales::rescale(value,to=c(0,1))
                                  ,scales::rescale(desc(value),to=c(0,1))
                                  )
                  ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(combo_init = scale*!!ensym(diag_col)) %>%
    dplyr::group_by(method,groups) %>%
    dplyr::mutate(combo = mean(combo_init)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(combo = scales::rescale(combo,to=c(0,1))) %>%
    dplyr::mutate(top_thresh = top_thresh
                  , best_thresh = best_thresh
                  , top = combo >= quantile(combo,probs = 1-top_thresh,na.rm = TRUE)
                  , top = if_else(is.na(top),FALSE,top)
                  , best = combo >= sort(unique(.$combo),TRUE)[best_thresh]
                  , best = if_else(is.na(best),FALSE,best)
                  , diagnostic = factor(diagnostic, levels = levels(diagnostics$diagnostic))
                  , weight = !!ensym(diag_col)
                  )

}

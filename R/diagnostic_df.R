

#' Use a set of (continuous) columns to choose a good set of candidate rows
#'
#' @param df Dataframe with columns over which to find good candidates.
#' @param diagnosticdf Dataframe mapping the name of possible diagnostics to cases (columns) in which to use those diagnostic.
#' @param context Character. Name of columns in df that define context.
#' @param diagcol Character. Name of diagnosticdf column to use in this instance.
#' @param summarisemethod Character. Name of method to use in summarising if there is more than one row per context.
#' @param mingroups Numeric. Minimum number of clusters in a clustering.
#' @param maxgroups Numeric. Maximum number of clusters in a clustering.
#' @param topthresh Numeric specifying the proportion of rows considered 'top'.
#' @param bestthresh Numeric specifying the absolute number of rows considered 'best'.
#'
#' @return
#' @export
#'
#' @examples
diagnostic_df <- function(df
                          , diagnosticdf = tibble::tibble(diagnostic = "avClustSize",highGood = TRUE, weightClustSum = TRUE)
                          , context = c("method","groups")
                          , diagcol = "weightClustSum"
                          , summarisemethod = median
                          , mingroups = 10
                          , maxgroups = Inf
                          , topthresh = 0.25
                          , bestthresh = 5
                          ) {

  df %>%
    dplyr::filter(groups >= mingroups
                  , groups <= maxgroups
                  ) %>%
    dplyr::select(all_of(context),groups,any_of(diagnostics$diagnostic)) %>%
    #select_if(~ !all(is.na(.))) %>%
    dplyr::group_by(across(all_of(context))) %>%
    dplyr::summarise(across(any_of(diagnostics$diagnostic),summarisemethod)) %>%
    dplyr::ungroup() %>%
    tidyr::pivot_longer(any_of(diagnostics$diagnostic),names_to = "diagnostic", values_to = "value") %>%
    dplyr::left_join(diagnostics) %>%
    dplyr::group_by(diagnostic,diagDefinition) %>%
    dplyr::mutate(scale = if_else(highGood
                                  ,scales::rescale(value,to=c(0,1))
                                  ,scales::rescale(desc(value),to=c(0,1))
                                  )
                  ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(comboInit = scale*!!ensym(useWeights)) %>%
    dplyr::group_by(method,groups) %>%
    dplyr::mutate(combo = mean(comboInit)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(combo = scales::rescale(combo,to=c(0,1))) %>%
    dplyr::mutate(topthresh = topthresh
                  , bestthresh = bestthresh
                  , top = combo >= quantile(combo,probs = 1-topthresh,na.rm = TRUE)
                  , top = if_else(is.na(top),FALSE,top)
                  , best = combo >= sort(unique(.$combo),TRUE)[bestthresh]
                  , best = if_else(is.na(best),FALSE,best)
                  , diagnostic = factor(diagnostic, levels = levels(diagnostics$diagnostic))
                  , weight = !!ensym(useWeights)
                  )

}

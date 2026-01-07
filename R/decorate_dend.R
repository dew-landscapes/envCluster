#' Make a dendrogram
#'
#'
#' @param clust_df Dataframe with context columns and a column with cluster
#' membership for that context.
#' @param dend Dendrogram
#' @param clust_col Character. Name of column containing cluster membership in
#' `clust_df`.
#' @param max_clust_labels Numeric. Branches will be labelled if there are less
#' than `max_clust_labels` unique values in `clust_col`
#' @param second_group_col Character. Optional. Name of column in `clust_df`
#' containing group membership in which `clust_col` is nested. Used to create
#' subdendrograms and associated lookup.
#' @param site_id Character. Name of column in `clust_df` containing
#' `1:ncol(clust_df)` that matches the order of the sites used when making
#' `dist_obj`.
#' @colour_col Character name of column in `clust_df` containing colours to use
#' for each level of `clust_col` in the dendrogram. Will be generated if not
#' supplied.
#' @label_col Character name of column in `clust_df` containing labels to use
#' on dendrogram leaves.
#'
#' @return List with dendrogram object (as `$dend`). If `meta_clust_col` provided,
#' the list will also include a named list of subdendrograms and (tibble) lookup
#' for the subdendrograms
#' @export
#'
#' @example inst/examples/make_clusters_ex.R
decorate_dend <- function(clust_df
                          , dend
                          , clust_col = "cluster"
                          , max_clust_labels = 30
                          , second_group_col = "ecotype"
                          , site_id = "site_id"
                          , colour_col = "colour"
                          , label_col = site_id
                          ) {

  if(! site_id %in% names(clust_df)) {

    stop("Need a site_id column in clust_df that matches the order of dist_obj")

  }

  clust_df <- clust_df |>
    dplyr::arrange(site_id)

  dend_raw <- dend |>
    as.dendrogram()

  dend <- list()

  if(! colour_col %in% names(clust_df)) {

    colours <- viridis::viridis(n = length(unique(clust_df[[clust_col]])))

    clust_df <- clust_df |>
      dplyr::mutate(!!rlang::ensym(colour_col) := colours[as.integer(as.factor(clust_df[[clust_col]]))])

  }

  k <- length(unique(factor(clust_df$cluster)))

  use_col <- clust_df[[colour_col]][order.dendrogram(dend_raw)]
  use_lab <- clust_df[[label_col]][order.dendrogram(dend_raw)]

  use_col_branch <- clust_df |>
    dplyr::slice(order.dendrogram(dend_raw)) |>
    dplyr::distinct(cluster, !!rlang::ensym(colour_col)) |>
    dplyr::pull(!!rlang::ensym(colour_col))

  use_lab_branch <- clust_df |>
    dplyr::slice(order.dendrogram(dend_raw)) |>
    dplyr::distinct(!!rlang::ensym(clust_col), !!rlang::ensym(colour_col)) |>
    dplyr::pull(!!rlang::ensym(clust_col))

  dend$dend <- dend_raw

  if(length(use_lab_branch) < max_clust_labels) {

    dend$dend <- dend$dend |>
      dendextend::set("labels", use_lab)

  }

  dend$dend <- dend$dend |>
    dendextend::set("labels_colors", use_col) |>
    dendextend::set("branches_k_color"
                    , value = use_col_branch
                    , k = k
                    ) |>
    dendextend::colour_branches(k = k
                                , col = use_col_branch
                                , groupLabels = use_lab_branch
                                )

  if(!is.null(second_group_col)) {

    two_cols <- clust_df |>
      dplyr::count(dplyr::across(tidyselect::any_of(c(clust_col, second_group_col)))) |>
      nrow()

    one_col <- clust_df |>
      dplyr::count(dplyr::across(!!rlang::ensym(clust_col))) |>
      nrow()

    if(one_col == two_cols) {

      dend$dend_list <- dendextend::get_subdendrograms(dend$dend
                                                       , k = length(unique(clust_df[[second_group_col]]))
                                                       )

      dend$lu_dend_list <- tibble::tibble(sub_dend = 1:(length(unique(clust_df[[second_group_col]])))) |>
        dplyr::mutate(site_id = purrr::map(sub_dend
                                           ,\(x) unique(labels(dend$dend_list[[x]]))
                                           )
                      ) |>
        tidyr::unnest(cols = c(site_id)) |>
        dplyr::mutate(site_id = as.numeric(site_id)) |>
        dplyr::left_join(clust_df)

      dend$dend_list <- purrr::set_names(dend$dend_list
                                         , unique(dend$lu_dend_list$ecotype)
                                         )

    } else {

      warning(clust_col
              , " is not nested within "
              , second_group_col
              , ". No sub dendrograms created"
              )

    }

  }

  return(dend)

}

#' Apply clustering algorithms
#'
#' @param bio_df Dataframe. Needs the 'cols' outlined below.
#' @param env_df Dataframe (wide) of environmental variables per sitecol.
#' @param methods_df Dataframe. What methods to use to create clusters.
#' @param context Character. Name of columns in `bio_df` that define the context.
#' @param taxa_col Character. Name of column with taxa.
#' @param num_col Character. Name of column with numeric data indicating
#' relative abundance of taxa.
#' @param env_cols Character. Name of columns in `env_df` that have the
#' @param groups Numeric vector indicating the range of groups within a
#' clustering.
#' @param cores Numeric. Number of cores to use where parallel options exist.
#' @param clust_col_in Character name of column of raw class labels.
#' @param clust_col_out Character name of column of character class labels.
#'
#' @return Dataframe with groups and methods columns and a list column of the
#' clustering for that number of groups and method as a tibble with one column
#' (numeric) 'clust' indicating group membership.
#' @export
#'
#' @examples
  make_clusters <- function(bio_df
                            , env_df = NULL
                            , methods_df = tibble(method = "average")
                            , context
                            , taxa_col = "taxa"
                            , num_col = NULL
                            , env_cols
                            , groups = 2:100
                            , cores = 1
                            , clust_col_in = "clust"
                            , clust_col_out = "cluster"
                            ) {


    .clust_col_in = clust_col_in
    .clust_col_out = clust_col_out

    if(cores > 1) {

      future::plan(multisession
                   , workers = cores
                   )

    } else {

      future::plan(sequential)

    }


    if(isTRUE(is.null(num_col))) df$p = 1

    .bio_df = bio_df
    .context = context
    .taxa_col = taxa_col
    .num_col = num_col

    df_wide <- make_wide_df(bio_df = .bio_df
                             , context = .context
                             )

    assign("flor_wide"
           , df_wide
           , envir = globalenv()
           )

    dat_wide <- df_wide %>%
      dplyr::select(-all_of(context)) %>%
      as.matrix()

    site_names <- df_wide %>%
      dplyr::select(all_of(context))

    dist_flor <- parallelDist::parDist(dat_wide
                                  , method = "bray"
                                  , threads = cores
                                  )

    assign("dist_flor"
           , dist_flor
           , envir = globalenv()
           )

    assign("sq_dist"
           , as.matrix(dist_flor^2)
           , envir = globalenv()
           )


    if(isTRUE(!is.null(env_df))) {

      cells <- site_names %>%
        dplyr::pull(cell)

      dist_env <- parallelDist::parDist(env_df %>%
                                          dplyr::filter(cell %in% cells) %>%
                                          dplyr::select(all_of(env_cols)) %>%
                                          as.matrix()
                                       , method = "euclidean"
                                       , threads = cores
                                       )

      assign("dist_env"
             , dist_env
             , envir = globalenv()
             )

    }

    dend <- methods_df %>%
      dplyr::mutate(dend = furrr::future_map(method
                                             ,~fastcluster::hclust(dist_flor
                                                                   , .
                                                                   )
                                             )
                    )

    clust <- dend %>%
      dplyr::mutate(clusters = furrr::future_map(dend
                                                 , cutree
                                                 , groups
                                                 )
                    , clusters = furrr::future_map(clusters
                                                   , as_tibble
                                                   )
                    ) %>%
      dplyr::select(-dend) %>%
      tidyr::unnest(clusters) %>%
      tidyr::pivot_longer(2:ncol(.)
                          , names_to = "groups"
                          , values_to = clust_col_in
                          ) %>%
      dplyr::mutate(groups = as.integer(groups)) %>%
      tidyr::nest(clusters = c(!!ensym(clust_col_in))) %>%
      dplyr::mutate(clusters = furrr::future_map(clusters
                                                 , make_cluster_df
                                                 , context_df = site_names
                                                 , context = .context
                                                 , clust_col_in = .clust_col_in
                                                 , clust_col_out = .clust_col_out
                                                 )
                    )


  }

#' Clusters with >= 1 taxa with frequency > than threshold
#'
#' @param clust_df Dataframe with column of cluster membership and join column to
#' flordf
#' @param flor_df Dataframe with taxa
#' @param clust_col Character. Name of column in `clust_df` with clusters.
#' @param taxa_col Character. Name of column in flordf with taxa.
#' @param site_col Character. Name of column in clustdf and flordf with 'sites'.
#' @param thresh Numeric. Proportion of sites in a cluster at which >= 1 taxa
#' need to occur.
#'
#' @return Numeric. Number of clusters across clustering that have at least one
#' taxa that occurs in greater than 'thresh' proportion of sites.
#' @export
#'
#' @examples
clusters_with_freq_taxa <- function(clust_df
                                    , flor_df
                                    , clust_col = "cluster"
                                    , taxa_col = "taxa"
                                    , site_col = "cell"
                                    , thresh = 0.9
                                    ){

  clust_df %>%
    dplyr::left_join(flor_df) %>%
    dplyr::group_by(!!ensym(clust_col)) %>%
    dplyr::mutate(sites = dplyr::n_distinct(cell)) %>%
    dplyr::count(!!ensym(taxa_col),!!ensym(clust_col),sites,name = "records") %>%
    dplyr::mutate(prop = records/sites) %>%
    dplyr::filter(prop == max(prop)) %>%
    dplyr::distinct(!!ensym(clust_col),prop) %>%
    dplyr::filter(prop >= thresh) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(!!ensym(clust_col)) %>%
    nrow()

}


# Clusters with >= 1 taxa with frequency > than threshold
#' Title
#'
#' @param clust_df Dataframe with column of cluster membership and join column to
#' flordf
#' @param flor_df Dataframe with taxa
#' @param clust_col Character. Name of column in clustdf with clusters.
#' @param taxa_col Character. Name of column in flordf with taxa.
#' @param site_col Character. Name of column in clustdf and flordf with 'sites'.
#' @param thresh Numeric. Proportion of sites in a cluster at which >= 1 taxa
#' need to occur.
#'
#' @return Numeric. Number of sites across clustering that are in a cluster with
#' at least one taxa that occurs in greater than 'thresh' proportion of sites
#' in that cluster.
#' @export
#'
#' @examples
sites_with_freq_taxa <- function(clust_df
                                 , flor_df
                                 , clust_col = "clust"
                                 , taxa_col = "taxa"
                                 , site_col = "cell"
                                 , thresh = 0.9
                                 ){

  clust_df %>%
    dplyr::left_join(flor_df) %>%
    dplyr::group_by(!!ensym(clust_col)) %>%
    dplyr::mutate(sites = dplyr::n_distinct(!!ensym(site_col))) %>%
    dplyr::count(!!ensym(taxa_col),!!ensym(clust_col),sites,name = "records") %>%
    dplyr::mutate(prop = records/sites) %>%
    dplyr::filter(prop == max(prop)) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(!!ensym(clust_col),prop,sites) %>%
    dplyr::filter(prop >= thresh) %>%
    dplyr::summarise(sites = sum(sites)) %>%
    dplyr::pull(sites)

}



#' Sites with indicator
#'
#' @param ind_val_df Dataframe from make_indval_df
#' @param clust_df Dataframe with column of cluster membership.
#' @param clust_col Name of column in `ind_val_df` with class membership.
#' @param thresh Numeric. p_val threshold for acceptance as indicator taxa.
#'
#' @return Numeric. Number of sites across clustering that are in a cluster with
#' at least one indicator taxa.
#' @export
#'
#' @examples
sites_with_indicator <- function(ind_val_df
                                 , clust_df
                                 , clust_col = "cluster"
                                 , thresh = 0.05
                                 ){

  ind_val_df %>%
    dplyr::distinct(!!ensym(clust_col),p_val) %>%
    dplyr::group_by(!!ensym(clust_col)) %>%
    dplyr::filter(p_val == min(p_val)) %>%
    dplyr::filter(p_val < thresh) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(!!ensym(clust_col)) %>%
    dplyr::inner_join(clust_df) %>%
    nrow()

}

# Clusters with indicator
#' Clusters with indicator
#'
#' @param ind_val_df Dataframe from make_indval_df
#' @param clust_col Name of column in `ind_val_df` with class membership.
#' @param thresh Numeric. p_val threshold for acceptance as indicator taxa.
#'
#' @return Numeric. Number of clusters with at least one indicator taxa.
#' @export
#'
#' @examples
clusters_with_indicator <- function(ind_val_df
                                    , clust_col = "cluster"
                                    , thresh = 0.05
                                    ){

  ind_val_df %>%
    dplyr::select(!!ensym(clust_col)
                  , p_val
                  ) %>%
    dplyr::group_by(!!ensym(clust_col)) %>%
    dplyr::filter(p_val == min(p_val)) %>%
    dplyr::filter(p_val < thresh) %>%
    dplyr::ungroup() %>%
    dplyr::count(!!ensym(clust_col)) %>%
    nrow()

}


#' Simple summary of a clustering.
#'
#' @param df Dataframe containing clusters as a list column.
#' @param clust_col Name of column in `df` with class membership.
#' @param min_sites Desired minimum absolute number of sites in a class.
#'
#' @return single row tibble with summary information
#' @export
#'
#' @examples
  cluster_summarise <- function(df
                                , clust_col = "cluster"
                                , min_sites = 10
                                ) {

    clust <- df %>%
      dplyr::select(all_of(clust_col))

    tab <- table(clust)

    tibble(min_clust_size = min(tab)
           , av_clust_size = mean(tab)
           , max_clust_size = max(tab)
           , n_min_clusters = sum(tab > min_sites)
           , n_min_sites = sum(tab[tab > min_sites])
           , prop_min_clusters = n_min_clusters/length(tab)
           , prop_min_sites = n_min_sites/nrow(clust)
           )

  }




#'  Make a cluster data frame
#'
#' @param raw_clusters Dataframe of cluster membership. Can be result of
#' `make_clusters()`.
#' @param context_df Dataframe with `context` columns in same order as
#' `raw_clusters`.
#' @param context Character. Name of columns in `site_df` defining the context.
#' @param clust_col_in Character. Name of column in `raw_clusters` with classes.
#' @param clust_col_out Character. Name of column to ouput with character class
#' labels.
#'
#' @return Dataframe with 'site','clust' (numeric) and 'cluster' (character)
#' @export
#'
#' @examples
  make_cluster_df <- function(raw_clusters
                              , context_df
                              , context
                              , clust_col_in = "clust"
                              , clust_col_out = "cluster"
                              ) {

    context_df %>%
      dplyr::select(all_of(context)) %>%
      dplyr::bind_cols(raw_clusters %>%
                         dplyr::mutate(!!ensym(clust_col_out) := numbers2words(!!ensym(clust_col_in))
                                       , !!ensym(clust_col_out) := forcats::fct_reorder(!!ensym(clust_col_out)
                                                                                        , !!ensym(clust_col_in)
                                                                                        )
                                       )
                       )

  }




#' Create a dataframe of silhouette widths
#'
#' @param clust_df Dataframe with column of cluster membership.
#' @param dist_obj Distance object with sites in the same order as clustdf.
#' @param clust_col Name or index of column containing cluster membership in
#' `clust_df`.
#'
#' @return Dataframe of silhouette width per site, including neighbouring
#' cluster.
#' @export
#'
#' @examples
make_sil_df <- function(clust_df, dist_obj, clust_col = "clust"){

  clusts <- clust_df[clust_col][[1]]

  sil_obj <- cluster::silhouette(clusts,dist_obj)

  clust_df %>%
    dplyr::bind_cols(tibble(neighbour = sil_obj[,2],sil_width = sil_obj[,3]))

}


#' Calculate within cluster sum of squares
#'
#'
#' @param clust_df Dataframe with column of cluster membership.
#' @param dist_obj Distance object from which clusters were generated.
#' @param clust_col Character. Name of column in clustdf that contains the clusters.
#'
#' @return Dataframe of within group sum-of-squares for each cluster. If it
#' doesn't already exist in global workspace, sqDist is added there as it is
#' needed for calculations.
#' @export
#'
#' @examples
calc_ss <- function(clust_df
                    , dist_obj
                    , clust_col = "clust"
                    ) {

  if(!exists("sq_dist")) assign("sq_dist"
                               , as.matrix(dist_obj^2)
                               , envir = globalenv()
                               )

  clust_df %>%
    dplyr::mutate(id = row_number()) %>%
    dplyr::group_by(!!ensym(clust_col)) %>%
    tidyr::nest() %>%
    dplyr::ungroup() %>%
    #dplyr::sample_n(2) %>% # TESTING
    dplyr::mutate(wss = map_dbl(data
                                , ~sum(sq_dist[.$id,.$id])/(2*nrow(.))
                                )
                  , wss_rank = rank(wss)
                  ) %>%
    dplyr::select(-data)

}



#' Create a dataframe of indicator values
#'
#' @param clust_df Dataframe including a column with cluster membership
#' @param bio_wide Dataframe of sites * taxa. In the same order as clustdf.
#' @param clust_col Column containing cluster membership in clustdf as name or
#' index.
#' @param taxas Character. Names of taxa to use.
#' @param context Character. Names of columns in `clust_df` that define the context.
#'
#' @return Dataframe of each taxa and the cluster (clust as numeric, cluster as
#' character) class to which it is most likely an indicator, plus the following
#' values from labdsv::indval output: ind_val, p_val, abu and frq.
#' @export
#'
#' @examples
make_ind_val_df <- function(clust_df
                            , bio_wide
                            , clust_col = "cluster"
                            , taxas
                            , context
                            ){

  bio_wide <- clust_df %>%
    dplyr::inner_join(bio_wide)

  clust_ind <- labdsv::indval(bio_wide[,names(bio_wide) %in% c(taxas)]
                             , as.character(bio_wide[clust_col][[1]])
                             )

  clusts <- levels(clust_df[clust_col][[1]])


  clust_ind$relabu %>%
    tibble::as_tibble(rownames = "taxa") %>%
    tidyr::pivot_longer(any_of(clusts)
                        , names_to = clust_col
                        , values_to = "abu"
                        ) %>%
    dplyr::inner_join(clust_ind$relfrq %>%
                        tibble::as_tibble(rownames = "taxa") %>%
                        tidyr::pivot_longer(any_of(clusts)
                                            , names_to = clust_col
                                            , values_to = "frq"
                                            )
                      ) %>%
    dplyr::inner_join(clust_ind$indval %>%
                        tibble::as_tibble(rownames = "taxa") %>%
                        tidyr::pivot_longer(any_of(clusts)
                                            , names_to = clust_col
                                            , values_to = "ind_val"
                                            )
                      ) %>%
    dplyr::group_by(taxa) %>%
    dplyr::filter(ind_val == max(ind_val)) %>%
    dplyr::ungroup()  %>%
    dplyr::inner_join(clust_ind$pval %>%
                        tibble::as_tibble(rownames = "taxa") %>%
                        dplyr::rename(p_val = value)
                      ) %>%
    dplyr::mutate(!!ensym(clust_col) := factor(!!ensym(clust_col), levels = clusts)) %>%
    dplyr::select(!!ensym(clust_col),everything()) %>%
    dplyr::arrange(!!ensym(clust_col),desc(ind_val))

}



#' Make a wide (usually site * taxa) data frame
#'
#' @param bio_df Dataframe containing the site and taxa data.
#' @param context Character. Name of columns defining context.
#' @param taxa_col Character. Name of column containing taxa.
#' @param num_col Character. Name of column containing numeric abundance data (
#' usually 'cover' for plants).
#' @param num_col_NA Value to use to fill NA in wide table.
#'
#' @return
#' @export
#'
#' @examples
make_wide_df <- function(bio_df
                         , context
                         , taxa_col = "taxa"
                         , num_col = "use_cover"
                         , num_col_NA = 0
                         ) {

  bio_df %>%
    dplyr::group_by(!!ensym(taxa_col),across(all_of(context))) %>%
    dplyr::summarise(value = max(!!ensym(num_col),na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(!is.na(value)) %>%
    tidyr::pivot_wider(names_from = !!ensym(taxa_col)
                       , values_fill = num_col_NA
                       )

}

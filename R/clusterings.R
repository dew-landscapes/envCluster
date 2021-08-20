#' Apply clustering algorithms
#'
#' @param df Dataframe. Needs the 'cols' outlined below.
#' @param env_df Dataframe (wide) of environmental variables per sitecol.
#' @param methods_df Dataframe. What methods to use to create clusters.
#' @param taxa_col Character. Name of column with taxa.
#' @param site_col Chracter. Name of column with 'site'.
#' @param num_col Character. Name of column with numeric data indicating
#' relative abundance of taxa.
#' @param groups Numeric vector indicating the range of groups within a
#' clustering.
#'
#' @return Dataframe with groups and methods columns and a list column of the
#' clustering for that number of groups and method as a tibble with one column
#' (numeric) 'clust' indicating group membership.
#' @export
#'
#' @examples
  make_clusters <- function(df
                            , env_df = NULL
                            , methods_df = tibble(method = "average")
                            , taxa_col = "taxa"
                            , site_col = "list"
                            , num_col = NULL
                            , groups = 2:100
                            ) {

    if(isTRUE(is.null(num_col))) df$p = 1

    dat_wide <- df %>%
      dplyr::group_by(!!ensym(taxa_col),!!ensym(site_col)) %>%
      dplyr::summarise(value = max(!!ensym(num_col),na.rm = TRUE)) %>%
      dplyr::filter(!is.na(value)) %>%
      tidyr::pivot_wider(names_from = all_of(taxa_col), values_fill = 0) %>%
      dplyr::arrange(!!ensym(site_col)) %>%
      tibble::column_to_rownames(site_col) %>%
      as.matrix()

    assign("dat_wide",dat_wide,envir = globalenv())

    site_names <- rownames(dat_wide)

    dist_flor <- parallelDist::parDist(dat_wide
                                  , method = "bray"
                                  , threads = if(exists("use_cores")) use_cores else 1
                                  )

    assign("dist_flor",dist_flor,envir = globalenv())

    if(isTRUE(!is.null(env_df))) {

      dist_env <- parallelDist::parDist(env_df %>%
                                         dplyr::filter(cell %in% site_names) %>%
                                         as.matrix()
                                       , method = "euclidean"
                                       , threads = if(exists("useCores")) useCores else 1
                                       )

      assign("dist_env",dist_env,envir = globalenv())

    }

    assign("sq_dist",as.matrix(dist_flor^2),pos = .GlobalEnv)

    dend <- methods_df %>%
      dplyr::mutate(dend = map(method
                               ,~fastcluster::hclust(dist_flor, .)
                               )
                    )

    clust <- dend %>%
      dplyr::mutate(clusters = map(dend
                                   , cutree
                                   , groups
                                   )
                    , clusters = map(clusters
                                     , as_tibble
                                     )
                    ) %>%
      dplyr::select(-dend) %>%
      tidyr::unnest(clusters) %>%
      tidyr::pivot_longer(2:ncol(.),names_to = "groups",values_to ="clust") %>%
      dplyr::mutate(groups = as.integer(groups)) %>%
      tidyr::nest(clusters = c(clust))

  }

#' Clusters with >= 1 taxa with frequency > than threshold
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
#' @return Numeric. Number of clusters across clustering that have at least one
#' taxa that occurs in greater than 'thresh' proportion of sites.
#' @export
#'
#' @examples
clusters_with_freq_taxa <- function(clust_df
                                    , flor_df
                                    , clust_col = "clust"
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
#' @param thresh Numeric. p_val threshold for acceptance as indicator taxa.
#'
#' @return Numeric. Number of sites across clustering that are in a cluster with
#' at least one indicator taxa.
#' @export
#'
#' @examples
sites_with_indicator <- function(ind_val_df,clust_df,thresh = 0.05){

  ind_val_df %>%
    dplyr::distinct(cluster,p_val) %>%
    dplyr::group_by(cluster) %>%
    dplyr::filter(p_val == min(p_val)) %>%
    dplyr::filter(p_val < thresh) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(cluster) %>%
    dplyr::inner_join(clust_df) %>%
    nrow()

}

# Clusters with indicator
#' Clusters with indicator
#'
#' @param ind_val_df Dataframe from make_indval_df
#' @param thresh Numeric. p_val threshold for acceptance as indicator taxa.
#'
#' @return Numeric. Number of clusters with at least one indicator taxa.
#' @export
#'
#' @examples
clusters_with_indicator <- function(ind_val_df, thresh = 0.05){

  ind_val_df %>%
    dplyr::select(cluster,p_val) %>%
    dplyr::group_by(cluster) %>%
    dplyr::filter(p_val == min(p_val)) %>%
    dplyr::filter(p_val < thresh) %>%
    dplyr::ungroup() %>%
    dplyr::count(cluster) %>%
    nrow()

}


#' Filter summarised clusters
#'
#' @param df Dataframe containing clusters as a list column.
#' @param clust_col Name of list column containing clusters.
#' @param min_sites Desired minimum absolute number of sites in a cluster.
#'
#' @return Dataframe of filtered clusters.
#' @export
#'
#' @examples
  cluster_summarise <- function(df
                                , clust_col = "clust"
                                , min_sites = 10
                                ) {

    clust <- df %>% dplyr::select(all_of(clust_col))

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
#' make_clusters.
#' @param site_df Dataframe with 'site' ids in same order as rawclusters.
#' @param site_col Character. Name of column in sitedf with sites.
#'
#' @return Dataframe with 'site','clust' (numeric) and 'cluster' (character)
#' @export
#'
#' @examples
  make_cluster_df <- function(raw_clusters,site_df,site_col = "site") {

    site_df %>%
      dplyr::select(!!ensym(site_col)) %>%
      dplyr::bind_cols(raw_clusters %>%
                         dplyr::mutate(cluster = numbers2words(clust)
                                       , cluster = forcats::fct_reorder(cluster,clust)
                                       )
                       )

  }




#' Create a dataframe of silhouette widths
#'
#' @param clust_df Dataframe with column of cluster membership.
#' @param dist_obj Distance object with sites in the same order as clustdf.
#' @param clust_col Name or index of column containing cluster membership in
#' clustdf as name or index.
#'
#' @return Dataframe of silhouette width per site, including neighbouring
#' cluster.
#' @export
#'
#' @examples
make_sil_df <- function(clust_df, dist_obj, clust_col = "clust"){

  clusts <- clustdf[clust_col][[1]]

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
calc_ss <- function(clust_df,dist_obj,clust_col = "clust") {

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
                  ) %>%
    dplyr::select(-data)

}



#' Create a dataframe of indicator values
#'
#' @param clust_df Dataframe including a column with cluster membership
#' @param df_wide Dataframe of sites * taxa. In the same order as clustdf.
#' @param clust_col Column containing cluster membership in clustdf as name or
#' index.
#'
#' @return Dataframe of each taxa and the cluster (clust as numeric, cluster as
#' character) class to which it is most likely an indicator, plus the following
#' values from labdsv::indval output: indval, p_val, abu and frq.
#' @export
#'
#' @examples
make_indval_df <- function(clust_df, df_wide, clust_col = "clust"){

  clust_ind <- labdsv::indval(df_wide
                             , clust_df[clust_col][[1]]
                             )

  tibble(taxa =names(clust_ind$maxcls)
         , clust = clust_ind$maxcls
         , ind_val = clust_ind$indcls
         , p_val = clust_ind$pval
         ) %>%
    dplyr::inner_join(clust_ind$relabu %>%
                        as_tibble(rownames = "taxa") %>%
                        tidyr::gather(clust,abu,names(.)[names(.) %in% unique(clust_df$clust)]) %>%
                        dplyr::mutate(clust = as.numeric(clust)) %>%
                        dplyr::filter(abu > 0)
                      ) %>%
    dplyr::inner_join(clust_ind$relfrq %>%
                        as_tibble(rownames = "taxa") %>%
                        tidyr::gather(clust,frq,names(.)[names(.) %in% unique(clust_df$clust)]) %>%
                        dplyr::mutate(clust = as.numeric(clust)) %>%
                        dplyr::filter(frq > 0)
                      ) %>%
    dplyr::mutate(clust = as.numeric(clust)
                  , cluster = numbers2words(as.numeric(clust))
                  , cluster = fct_reorder(cluster,clust)
                  , taxa = gsub("\\."," ",taxa)
                  )

}


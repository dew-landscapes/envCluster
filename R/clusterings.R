
#' Clusters with >= 1 taxa with frequency > than threshold
#'
#' @param clustdf Dataframe with column of cluster membership and join column to
#' flordf
#' @param flordf Dataframe with taxa
#' @param clustcol Character. Name of column in clustdf with clusters.
#' @param taxacol Character. Name of column in flordf with taxa.
#' @param sitecol Character. Name of column in clustdf and flordf with 'sites'.
#' @param thresh Numeric. Proportion of sites in a cluster at which >= 1 taxa
#' need to occur.
#'
#' @return Numeric. Number of clusters across clustering that have at least one
#' taxa that occurs in greater than 'thresh' proportion of sites.
#' @export
#'
#' @examples
clusters_with_freq_taxa <- function(clustdf
                                    , flordf
                                    , clustcol = "clust"
                                    , taxacol = "taxa"
                                    , sitecol = "cell"
                                    , thresh = 0.9
                                    ){

  clustdf %>%
    dplyr::left_join(flordf) %>%
    dplyr::group_by(!!ensym(clustcol)) %>%
    dplyr::mutate(sites = dplyr::n_distinct(cell)) %>%
    dplyr::count(!!ensym(taxacol),!!ensym(clustcol),sites,name = "records") %>%
    dplyr::mutate(prop = records/sites) %>%
    dplyr::filter(prop == max(prop)) %>%
    dplyr::distinct(!!ensym(clustcol),prop) %>%
    dplyr::filter(prop >= thresh) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(!!ensym(clustcol)) %>%
    nrow()

}


# Clusters with >= 1 taxa with frequency > than threshold
#' Title
#'
#' @param clustdf Dataframe with column of cluster membership and join column to
#' flordf
#' @param flordf Dataframe with taxa
#' @param clustcol Character. Name of column in clustdf with clusters.
#' @param taxacol Character. Name of column in flordf with taxa.
#' @param sitecol Character. Name of column in clustdf and flordf with 'sites'.
#' @param thresh Numeric. Proportion of sites in a cluster at which >= 1 taxa
#' need to occur.
#'
#' @return Numeric. Number of sites across clustering that are in a cluster with
#' at least one taxa that occurs in greater than 'thresh' proportion of sites
#' in that cluster.
#' @export
#'
#' @examples
sites_with_freq_taxa <- function(clustdf
                                 , flordf
                                 , clustcol = "clust"
                                 , taxacol = "taxa"
                                 , sitecol = "cell"
                                 , thresh = 0.9
                                 ){

  clustdf %>%
    dplyr::left_join(flordf) %>%
    dplyr::group_by(!!ensym(clustcol)) %>%
    dplyr::mutate(sites = dplyr::n_distinct(!!ensym(sitecol))) %>%
    dplyr::count(!!ensym(taxacol),!!ensym(clustcol),sites,name = "records") %>%
    dplyr::mutate(prop = records/sites) %>%
    dplyr::filter(prop == max(prop)) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(!!ensym(clustcol),prop,sites) %>%
    dplyr::filter(prop >= thresh) %>%
    dplyr::summarise(sites = sum(sites)) %>%
    dplyr::pull(sites)

}



#' Sites with indicator
#'
#' @param indvaldf Dataframe from make_indval_df
#' @param clustdf Dataframe with column of cluster membership.
#' @param thresh Numeric. pval threshold for acceptance as indicator taxa.
#'
#' @return Numeric. Number of sites across clustering that are in a cluster with
#' at least one indicator taxa.
#' @export
#'
#' @examples
sites_with_indicator <- function(indvaldf,clustdf,thresh = 0.05){

  indvaldf %>%
    dplyr::distinct(cluster,pval) %>%
    dplyr::group_by(cluster) %>%
    dplyr::filter(pval == min(pval)) %>%
    dplyr::filter(pval < thresh) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(cluster) %>%
    dplyr::inner_join(clustdf) %>%
    nrow()

}

# Clusters with indicator
#' Clusters with indicator
#'
#' @param indvaldf Dataframe from make_indval_df
#' @param thresh Numeric. pval threshold for acceptance as indicator taxa.
#'
#' @return Numeric. Number of clusters with at least one indicator taxa.
#' @export
#'
#' @examples
clusters_with_indicator <- function(indvaldf, thresh = 0.05){

  indvaldf %>%
    dplyr::select(cluster,pval) %>%
    dplyr::group_by(cluster) %>%
    dplyr::filter(pval == min(pval)) %>%
    dplyr::filter(pval < thresh) %>%
    dplyr::ungroup() %>%
    dplyr::count(cluster) %>%
    nrow()

}


#' Filter summarised clusters
#'
#' @param df Dataframe containing clusters as a list column.
#' @param clustcol Name of list column containing clusters.
#' @param minsites Desired minimum absolute number of sites in a cluster.
#'
#' @return Dataframe of filtered clusters.
#' @export
#'
#' @examples
  cluster_summarise <- function(df
                                , clustcol = "clust"
                                , minsites = 10
                                ) {

    clust <- df %>% dplyr::select(all_of(clustcol))

    tab <- table(clust)

    tibble(minclustsize = min(tab)
           , avclustsize = mean(tab)
           , maxclustsize = max(tab)
           , nminclusters = sum(tab > minsites)
           , nminsites = sum(tab[tab > minsites])
           , propminclusters = nminclusters/length(tab)
           , propminsites = nminsites/nrow(clust)
           )

  }

#' Apply clustering algorithms
#'
#' @param df Dataframe. Needs the 'cols' outlined below.
#' @param envdf Dataframe (wide) of environmental variables per sitecol.
#' @param methodsdf Dataframe. What methods to use to create clusters.
#' @param taxacol Character. Name of column with taxa.
#' @param sitecol Chracter. Name of column with 'site'.
#' @param numcol Character. Name of column with numeric data indicating
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
                            , envdf = NULL
                            , methodsdf = tibble(method = "average")
                            , taxacol = "taxa"
                            , sitecol = "list"
                            , numcol = NULL
                            , groups = 2:100
                            ) {

    if(isTRUE(is.null(numcol))) df$p = 1

    datwide <- df %>%
      dplyr::group_by(!!ensym(taxacol),!!ensym(sitecol)) %>%
      dplyr::summarise(value = max(!!ensym(numcol),na.rm = TRUE)) %>%
      dplyr::filter(!is.na(value)) %>%
      tidyr::pivot_wider(names_from = all_of(taxacol), values_fill = 0) %>%
      dplyr::arrange(!!ensym(sitecol)) %>%
      tibble::column_to_rownames(sitecol) %>%
      as.matrix()

    assign("datwide",datwide,envir = globalenv())

    siteNames <- rownames(datwide)

    distflor <- parallelDist::parDist(datwide
                                  , method = "bray"
                                  , threads = if(exists("useCores")) useCores else 1
                                  )

    assign("distflor",distflor,envir = globalenv())

    if(isTRUE(!is.null(envdf))) {

      distenv <- parallelDist::parDist(envdf %>%
                                         dplyr::filter(cell %in% siteNames) %>%
                                         as.matrix()
                                       , method = "euclidean"
                                       , threads = if(exists("useCores")) useCores else 1
                                       )

      assign("distenv",distenv,envir = globalenv())

    }

    assign("sqdist",as.matrix(distflor^2),pos = .GlobalEnv)

    dend <- methodsdf %>%
      dplyr::mutate(dend = map(method
                               ,~fastcluster::hclust(distflor, .)
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


#'  Make a cluster data frame
#'
#' @param rawclusters Dataframe of cluster membership. Can be result of
#' make_clusters.
#' @param sitedf Dataframe with 'site' ids in same order as rawclusters.
#' @param sitecol Character. Name of column in sitedf with sites.
#'
#' @return Dataframe with 'site','clust' (numeric) and 'cluster' (character)
#' @export
#'
#' @examples
  make_cluster_df <- function(rawclusters,sitedf,sitecol = "site") {

    sitedf %>%
      dplyr::select(!!ensym(sitecol)) %>%
      dplyr::bind_cols(rawclusters %>%
                         dplyr::mutate(cluster = numbers2words(clust)
                                       , cluster = forcats::fct_reorder(cluster,clust)
                                       )
                       )

  }




#' Create a dataframe of silhouette widths
#'
#' @param clustdf Dataframe with column of cluster membership.
#' @param distobj Distance object with sites in the same order as clustdf.
#' @param clustcol Name or index of column containing cluster membership in
#' clustdf as name or index.
#'
#' @return Dataframe of silhouette width per site, including neighbouring
#' cluster.
#' @export
#'
#' @examples
make_sil_df <- function(clustdf, distobj, clustcol = "clust"){

  clusts <- clustdf[clustcol][[1]]

  silobj <- cluster::silhouette(clusts,distobj)

  clustdf %>%
    dplyr::bind_cols(tibble(neighbour = silobj[,2],sil_width = silobj[,3]))

}


#' Calculate within cluster sum of squares
#'
#'
#' @param clustdf Dataframe with column of cluster membership.
#' @param distobj Distance object from which clusters were generated.
#' @param clustcol Character. Name of column in clustdf that contains the clusters.
#'
#' @return Dataframe of within group sum-of-squares for each cluster. If it
#' doesn't already exist in global workspace, sqDist is added there as it is
#' needed for calculations.
#' @export
#'
#' @examples
calc_ss <- function(clustdf,distobj,clustcol = "clust") {

  if(!exists("sqdist")) assign("sqdist"
                               , as.matrix(distobj^2)
                               , envir = globalenv()
                               )

  clustdf %>%
    dplyr::mutate(id = row_number()) %>%
    dplyr::group_by(!!ensym(clustcol)) %>%
    tidyr::nest() %>%
    dplyr::ungroup() %>%
    #dplyr::sample_n(2) %>% # TESTING
    dplyr::mutate(wss = map_dbl(data
                                , ~sum(sqdist[.$id,.$id])/(2*nrow(.))
                                )
                  ) %>%
    dplyr::select(-data)

}



#' Create a dataframe of indicator values
#'
#' @param clustdf Dataframe including a column with cluster membership
#' @param dfwide Dataframe of sites * taxa. In the same order as clustdf.
#' @param clustcol Column containing cluster membership in clustdf as name or
#' index.
#'
#' @return Dataframe of each taxa and the cluster (clust as numeric, cluster as
#' character) class to which it is most likely an indicator, plus the following
#' values from labdsv::indval output: indval, pval, abu and frq.
#' @export
#'
#' @examples
make_indval_df <- function(clustdf, dfwide, clustcol = "clust"){

  clustind <- labdsv::indval(dfwide
                             , clustdf[clustcol][[1]]
                             )

  tibble(taxa =names(clustind$maxcls)
         , clust = clustind$maxcls
         , indval = clustind$indcls
         , pval = clustind$pval
         ) %>%
    dplyr::inner_join(clustind$relabu %>%
                        as_tibble(rownames = "taxa") %>%
                        tidyr::gather(clust,abu,names(.)[names(.) %in% unique(clustdf$clust)]) %>%
                        dplyr::mutate(clust = as.numeric(clust)) %>%
                        dplyr::filter(abu > 0)
                      ) %>%
    dplyr::inner_join(clustind$relfrq %>%
                        as_tibble(rownames = "taxa") %>%
                        tidyr::gather(clust,frq,names(.)[names(.) %in% unique(clustdf$clust)]) %>%
                        dplyr::mutate(clust = as.numeric(clust)) %>%
                        dplyr::filter(frq > 0)
                      ) %>%
    dplyr::mutate(clust = as.numeric(clust)
                  , cluster = numbers2words(as.numeric(clust))
                  , cluster = fct_reorder(cluster,clust)
                  , taxa = gsub("\\."," ",taxa)
                  )

}


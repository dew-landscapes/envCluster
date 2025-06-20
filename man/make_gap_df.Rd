% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/make_gap_df.R
\name{make_gap_df}
\alias{make_gap_df}
\title{Calculate gap statistic for a clustering}
\usage{
make_gap_df(clust_df, dist_mat, clust_col = "cluster", n_sample = 10)
}
\arguments{
\item{clust_df}{Dataframe with context columns and a column with cluster
membership for that context. Optional if \code{clust_col} appears in bio_df.}

\item{dist_mat}{Distance matrix (not distance object)}

\item{clust_col}{Character. Name of column in \code{clusters} identifying cluster
membership.}

\item{n_sample}{Numeric. Number of times to shuffle the membership to
recalculate the wss}
}
\value{
Single row dataframe with macro_wss and macro_gap plus a list column
of one row per level of \code{clust_col}, summarising the gap and wss data for
each level.
}
\description{
Calculate gap statistic for a clustering
}
\examples{

# data from: https://www.davidzeleny.net/anadat-r/doku.php/en:data:dune
bio_wide <- read.delim('https://raw.githubusercontent.com/zdealveindy/anadat-r/master/data/dune2.spe.txt', row.names = 1) |>
  tibble::as_tibble() |>
  dplyr::mutate(site_id = dplyr::row_number())

bio_long <- bio_wide |>
  tidyr::pivot_longer(cols = - site_id)

sites <- read.delim('https://raw.githubusercontent.com/zdealveindy/anadat-r/master/data/dune2.env.txt', row.names = 1) |>
  tibble::as_tibble() |>
  dplyr::mutate(site_id = dplyr::row_number()
                , Moisture = forcats::fct_relevel(as.character(Moisture), "1", "2", "4", "5")
                )

bio_dist <- vegan::vegdist(bio_wide)

use_methods <- c("mcquitty", "ward.D2", "average", "complete")

# make clusters -----------
clusters <- purrr::map(use_methods
                       , \(x) make_clusters(method = x
                                            , dist_obj = bio_dist
                                            , group_range = 2:(nrow(bio_wide) - 1)
                                            , sites = sites
                                            )
                       ) |>
  purrr::set_names(use_methods) |>
  dplyr::bind_rows(.id = "method")

# summarise clusters --------
clusters_summarise <- clusters |>
  dplyr::mutate(summary = purrr::map(clusters, \(x) make_summary_df(x))) |>
  tidyr::unnest(cols = c(summary))

# frequency of taxa --------
clusters_freq <- clusters |>
  dplyr::mutate(freq = purrr::map(clusters, \(x) make_freq_df(clust_df = x
                                                              , bio_df = bio_long
                                                              , context = "site_id"
                                                              )
                                  )
                ) |>
  dplyr::select(-clusters) |>
  tidyr::unnest(cols = c(freq))

# indicator value ---------
clusters_ind_val <- clusters |>
  dplyr::mutate(ind_val = purrr::map(clusters
                                     , \(x) make_ind_val_df(x
                                                            , bio_wide = bio_wide
                                                            , context = "site_id"
                                                            )
                                     )
                ) |>
  dplyr::select(-clusters) |>
  tidyr::unnest(cols = c(ind_val))

# silhouette width ---------
clusters_sil <- clusters |>
  dplyr::mutate(sil = purrr::map(clusters
                                 , \(x) make_sil_df(x
                                                    , dist_obj = bio_dist
                                                    )
                                 )
                ) |>
  dplyr::select(-clusters) |>
  tidyr::unnest(cols = c(sil))


# gap statistic --------
clusters_gap <- clusters |>
  dplyr::mutate(gap = purrr::map(clusters
                                 , \(x) make_gap_df(x
                                                    , dist_mat = as.matrix(bio_dist)
                                                    , n_sample = 30
                                                    )
                                 )
                ) |>
  dplyr::select(-clusters) |>
  tidyr::unnest(cols = c(gap))

# explore -------

clusters_explore <- clusters_summarise |>
  dplyr::left_join(clusters_freq) |>
  dplyr::left_join(clusters_ind_val) |>
  dplyr::left_join(clusters_sil) |>
  dplyr::left_join(clusters_gap)

# find a good cluster ---------
best <- clusters_explore |>
  ## min_clust_size is larger than one site
  dplyr::filter(min_clust_size > 1) |>
  ## most clusters have an indicator
  dplyr::filter(prop_ind_clusters > 0.8) |>
  ## high macro_sil
  dplyr::filter(macro_sil > quantile(macro_sil, probs = 0.5)) |>
  ## high macro_gap
  dplyr::filter(macro_gap > quantile(macro_gap, probs = 0.5)) |>
  ## low wss
  dplyr::filter(macro_wss < quantile(macro_wss, probs = 0.5))

## dendogram --------
dend <- make_dend(clust_df = best$clusters[[1]]
                  , dist_obj = bio_dist
                  , method = best$method[[1]]
                  , second_group_col = "Moisture"
                  , label_col = "Moisture"
                  )

plot(dend$dend)

## silhouette plot ---------

make_sil_plot(best$sil[[1]])
}

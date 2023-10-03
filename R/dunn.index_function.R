#' Function for the Dunn index
#' (Comparisons' comparability for transitivity evaluation)
#'
#' @description
#'   \code{dunn_index} calculates the Dunn index for a specific linkage method
#'   and number of clusters.
#'
#' @param input An object of 'dist' class. It is a lower off-diagonal matrix
#'   with the dissimilarities of all pairs of comparisons.
#' @param method A character string with values \code{"ward.D"},
#'   \code{"ward.D2"}, \code{"single"}, \code{"complete"}, \code{"average"},
#'   \code{"mcquitty"}, \code{"median"}, or \code{"centroid"} for the
#'   linkage method.
#' @param num_clusters A positive integer for the number of clusters. It
#'   takes values from two to the number of comparisons minus one.
#'
#' @return A scalar value for the Dunn index.
#'
#' @details
#'   \code{dunn_index} is integrated in the function
#'   \code{\link{internal_measures_plot}}. The index ranges from zero to
#'   infinity, with values larger than one indicating good clustering.
#'
#'   \code{dunn_index} uses the \code{\link[stats:cutree]{cutree}}
#'   function (alongside the \code{\link[stats:hclust]{hclust}} function) to cut
#'   the dendrogram for a specified number of clusters via the argument
#'   \code{num_clusters} and linkage method via the argument \code{method}.
#'
#' @author {Loukia M. Spineli}
#'
#' @seealso
#'   \code{\link[stats:cutree]{cutree}},
#'   \code{\link[stats:hclust]{hclust}}, \code{\link{internal_measures_plot}}
#'
#' @references
#' Dunn J. Well-separated clusters and optimal fuzzy partitions.
#' \emph{J Cybern} 1974;\bold{4}(1):95--104.
#'
#' @export
dunn_index <- function (input,
                        method,
                        num_clusters) {


  ## Check the defaults
  input <- if (inherits(input, "dist")) {
    input
  } else {
    stop("'input' must be of class 'dist'", call. = FALSE)
  }

  methods_list <- c("ward.D", "ward.D2", "single", "complete", "average",
                    "mcquitty", "median", "centroid")
  m_list1 <- c("'ward.D', 'ward.D2', 'single', 'complete', 'average'")
  m_list2 <- c("'mcquitty', 'median', 'centroid'")
  method <- if (!is.element(method, methods_list)) {
    stop(paste("'method' must be any of the following:", m_list1, m_list2),
         call. = FALSE)
  } else {
    method
  }

  num_clusters <-
    if (num_clusters > dim(as.matrix(input))[1] - 1 || num_clusters < 2) {
      stop(paste0("'num_clusters' must range from 2 to", " ",
                  dim(as.matrix(input))[1] - 1, "."), call. = FALSE)
  } else {
    num_clusters
  }


  ## Cut dendrogram for two clusters
  cut_dend <- cutree(hclust(input, method = method), k = num_clusters)


  ## The cluster where each comparison belongs
  cluster_res <- data.frame(comp = names(cut_dend), cluster = cut_dend)
  rownames(cluster_res) <- NULL


  ## Turn the dissimilarity matrix into 'matrix' ('dist' object currently)
  diss_mat <- round(as.matrix(input), 4)


  ## For each comparison, sort the comparisons from closest to furthest
  ##  based on their dissimilarity
  close_to_further <-
    do.call(rbind,
            lapply(1:dim(cluster_res)[1],
                   function(x) names(sort(diss_mat[x, ], decreasing = FALSE))))


  ## Similar with above but report the dissimilarities
  close_to_further_dist <-
    do.call(rbind,
            lapply(1:dim(cluster_res)[1],
                   function(x) sort(diss_mat[x, ], decreasing = FALSE)))
  rownames(close_to_further_dist) <- colnames(close_to_further_dist) <- NULL


  ## Replace each comparison with its rank
  for (i in 1:length(cluster_res$comp)) {
    close_to_further[close_to_further == cluster_res[i, 1]] <- cluster_res[i, 2]
  }


  ## Turn into numeric
  close_to_further_new <- matrix(as.numeric(close_to_further),
                                 nrow = dim(close_to_further)[1],
                                 ncol = dim(close_to_further)[2])


  ## Indicate the comparisons that belong to the same cluster
  ## If cluster comprises one data, same_cluster is different from 0.
  same_cluster <-
    ifelse(abs(close_to_further_new - close_to_further_new[, 1]) > 0, NA, 1)
  same_cluster[, 1] <- NA


  ## Get the distances and replace diagonal with NA
  same_cluster_dist <- close_to_further_dist * same_cluster


  ## Largest distance between observations found in the same cluster
  max_intra_cluster0 <- max(same_cluster_dist, na.rm = TRUE)


  ## Apply the correction of 0.0001
  max_intra_cluster <-
    ifelse(max_intra_cluster0 == 0, 0.0001, max_intra_cluster0)


  ## Indicate the comparisons that do *not* belong to the same cluster
  diff_cluster <-
    ifelse(abs(close_to_further_new - close_to_further_new[, 1]) > 0,
           close_to_further_new, NA)
  diff_cluster_dummy <-
    ifelse(abs(close_to_further_new - close_to_further_new[, 1]) > 0, 1, NA)


  ## Get distances
  diff_cluster_dist <- close_to_further_dist * diff_cluster_dummy


  ## Minimum distance between observations not found in the same cluster
  min_inter_cluster <- min(diff_cluster_dist, na.rm = TRUE)


  ## Dunn index
  dunn_value <- min_inter_cluster / max_intra_cluster

  return(dunn_value)
}

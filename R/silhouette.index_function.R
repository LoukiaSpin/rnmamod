#' Function for the silhouette width and silhouette plot
#' (Comparisons' comparability for transitivity evaluation)
#'
#' @description
#'   \code{silhouette_index} calculates the average silhouette width and the
#'   comparison-specific silhouette widths for a specific linkage method and
#'   number of clusters.
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
#' @return
#'   \code{silhouette_index} returns the following list of elements:
#'   \item{silhoutte_comparisons}{A data-frame on the silhouette width of each
#'   comparison.}
#'   \item{silhoutte_width}{A scalar on the average silhouette width across the
#'   comparisons.}
#'
#' @details
#'   \code{silhouette_index} is integrated in the function
#'   \code{\link{internal_measures_plot}}. Silhouette width ranges from minus to
#'   plus one, with values closer to 1 indicating higher compactness and
#'   separation for the selected clustering partition, whilst values closer to
#'   -1 indicate possible misclassification.
#'
#'   \code{silhouette_index} uses the \code{\link[stats:cutree]{cutree}}
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
#' Rousseeuw PJ. Silhouettes: A graphical aid to the interpretation and
#' validation of cluster analysis.
#' \emph{J Comput Appl Math} 1987;\bold{20}:53--65.
#'
#' @export
silhouette_index <- function (input,
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


  ## Cut dendrogram for a defined number of clusters (num_clusters)
  cut_dend <- cutree(hclust(input, method = method), k = num_clusters)


  ## Data-frame with the cluster where each comparison belongs
  cluster_res <- data.frame(comp = names(cut_dend), cluster = cut_dend)
  rownames(cluster_res) <- NULL


  ## Turn the dissimilarity matrix into 'matrix' ('dist' object currently)
  diss_mat <- as.matrix(input)


  ## For each comparison (row), sort the comparisons from closest to furthest
  ##  based on their dissimilarity (column per row)
  close_to_further <-
    do.call(rbind,
            lapply(1:dim(cluster_res)[1],
                   function(x) names(sort(diss_mat[x, ], decreasing = FALSE))))


  ## Similar with above but report the dissimilarities
  close_to_further_dist <-
    do.call(rbind,
            lapply(1:dim(cluster_res)[1],
                   function(x) sort(diss_mat[x, ], decreasing = FALSE)))
  colnames(close_to_further_dist) <- NULL


  ## Replace each dissimilarity with the corresponding cluster of the comparison
  for (i in 1:length(cluster_res$comp)) {
    close_to_further[close_to_further == cluster_res[i, 1]] <- cluster_res[i, 2]
  }


  ## Turn into numeric
  close_to_further_new <- matrix(as.numeric(close_to_further),
                                 nrow = dim(close_to_further)[1],
                                 ncol = dim(close_to_further)[2])


  ## Indicate the comparisons that belong to the same cluster
  ## Remove the first column that refers to comparison with the same comparison
  same_cluster <-
    ifelse(abs(close_to_further_new - close_to_further_new[, 1]) > 0, NA, 1)
  same_cluster[, 1] <- NA


  ## Get the distances
  same_cluster_dist <- close_to_further_dist * same_cluster


  ## Mean distance between each comparison and those found in the same cluster
  alpha0 <- apply(same_cluster_dist, 1, mean, na.rm = TRUE)
  alpha <- ifelse(is.na(alpha0), 0, alpha0)


  ## Indicate the comparisons that do *not* belong to the same cluster
  diff_cluster <-
    ifelse(abs(close_to_further_new - close_to_further_new[, 1]) > 0,
           close_to_further_new, NA)


  ## Dummy indicator
  diff_cluster_dummy <-
    ifelse(abs(close_to_further_new - close_to_further_new[, 1]) > 0, 1, NA)


  ## Get distances
  diff_cluster_dist <- close_to_further_dist * diff_cluster_dummy
  #diff_cluster_dist <- ifelse(is.na(diff_cluster_dist0), 0, diff_cluster_dist0)


  ## Mean distance of each comparison with those found in different clusters
  diff_cluster_mean <-
    lapply(1:dim(diss_mat)[1], function(x)
      aggregate(diff_cluster_dist[x, ], by = list(diff_cluster[x, ]), mean))


  ## Mean distance between each comparison and the nearest neighbouring cluster
  beta <- unlist(lapply(diff_cluster_mean, function(x) min(x[, 2])))


  ## Silhouette value by comparison
  silhoutte_comp0 <- (beta - alpha) / apply(cbind(beta, alpha), 1, max)


  ## For clusters with one observation, si = 0
  silhoutte_comp <- data.frame(cluster_res,
                               ifelse(is.na(alpha0) | beta == alpha,
                                      0,
                                      silhoutte_comp0))
  colnames(silhoutte_comp)[3] <- "silhouette"


  ## Silhouette width
  silhoutte_width <- mean(silhoutte_comp[, 3])

  return(list(silhoutte_comparisons = silhoutte_comp,
              silhoutte_width = silhoutte_width))
}

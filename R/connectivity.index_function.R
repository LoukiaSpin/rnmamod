#' Function for the Connectivity index
#' (Comparisons' comparability for transitivity evaluation)
#'
#' @description
#'   \code{connectivity_index} calculates the connectivity index for a specific
#'   linkage method, number of clusters and neighbouring comparisons.
#'
#' @param input An object of 'dist' class. It is a lower off-diagonal matrix
#'   with the dissimilarities of all pairs of comparisons.
#' @param method A character string with values \code{"ward.D"},
#'   \code{"ward.D2"}, \code{"single"}, \code{"complete"}, \code{"average"},
#'   \code{"mcquitty"}, \code{"median"}, or \code{"centroid"} for the
#'   linkage method.
#' @param num_clusters A positive integer for the number of clusters. It
#'   takes values from two to the number of comparisons minus one.
#' @param num_neighb A positive integer for the number of neighbouring
#'   comparisons. It takes values from two to the number of comparisons minus
#'   one. The default argument equals half the number of comparisons.
#'
#' @return A scalar value for the connectivity index.
#'
#' @details
#'   \code{connectivity_index} is integrated in the function
#'   \code{\link{internal_measures_plot}}. The index ranges from zero to
#'   infinity with smaller values indicating higher connectivity for the
#'   selected number of clusters, which is desired.
#'
#'   \code{connectivity_index} uses the \code{\link[stats:cutree]{cutree}}
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
#' Handl J, Knowles J, Kell DB. Computational cluster validation in post-genomic
#' data analysis. \emph{Biometrics} 2005;\bold{21}(15):3201--120.
#' doi: 10.1093/bioinformatics/bti517
#'
#' @export
connectivity_index <- function (input,
                                method,
                                num_clusters,
                                num_neighb) {


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

  num_neighb <- if (missing(num_neighb)) {
    message(paste("- num_neighb =", round(dim(as.matrix(input))[1] / 2, 0),
                  "was used (default)"))
    round(dim(as.matrix(input))[1] / 2, 0)
  } else if (num_neighb > dim(as.matrix(input))[1] || num_neighb < 2) {
    stop(paste0("'num_neighb' must range from 2 to", " ",
                dim(as.matrix(input))[1], "."), call. = FALSE)
  } else {
    num_neighb
  }


  ## Cut dendrogram for two clusters
  cut_dend <- cutree(hclust(input, method = method), k = num_clusters)


  ## The cluster where each comparison belongs
  cluster_res <- data.frame(comp = names(cut_dend), cluster = cut_dend)
  rownames(cluster_res) <- NULL


  ## Turn the dissimilarity matrix into 'matrix' ('dist' object currently)
  diss_mat <- as.matrix(input)


  ## For each comparison, sort the comparisons from closest to furthest
  ##  based on their dissimilarity
  close_to_further <-
    do.call(rbind,
            lapply(1:dim(cluster_res)[1],
                   function(x) names(sort(diss_mat[x, ], decreasing = FALSE))))


  ## Replace each comparison with its rank
  for (i in 1:length(cluster_res$comp)) {
    close_to_further[close_to_further == cluster_res[i, 1]] <- cluster_res[i, 2]
  }


  ## Turn into numeric
  close_to_further_new <- matrix(as.numeric(close_to_further),
                                 nrow = dim(close_to_further)[1],
                                 ncol = dim(close_to_further)[2])[, -1]


  ## Calculate difference in rank between each comparison and its neighbors
  diff_comp_neigh <-
    ifelse((cluster_res$cluster - close_to_further_new) != 0, 1, 0)

  # Calculate total connectivity
  total_conn <-
    sum(t(t(diff_comp_neigh[, 1:num_neighb]) *
            unlist(lapply(1:num_neighb, function(x) {1/x}))))

  return(total_conn)
}

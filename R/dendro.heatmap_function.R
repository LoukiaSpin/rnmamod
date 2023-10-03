#' Dendrogram with amalgamated heatmap
#' (Comparisons' comparability for transitivity evaluation)
#'
#' @description
#'   \code{dendro_heatmap} creates a dendrogram alongside the heatmap of
#'   dissimilarities among the comparisons for a specific dissimilarity measure,
#'   linkage method and number of clusters.
#'
#' @param input An object of S3 class \code{\link{comp_clustering}}. See 'Value'
#'   in \code{\link{comp_clustering}}.
#'
#' @return
#'   \code{dendro_heatmap} uses the \code{\link[heatmaply:heatmaply]{heatmaply}}
#'   function of the R-package
#'   \href{https://CRAN.R-project.org/package=heatmaply}{heatmaply} to create a
#'   cluster heatmap for a selected linkage method and number of clusters. The
#'   function uses different colours to indicate the clusters directly on the
#'   dendrogram. The comparisons in the have been sorted in decreasing order of
#'   their total dissimilarity.
#'
#'  @details
#'    The function inherits the dissimilarity measure (Euclidean or Canberra),
#'    linkage method and number of optimal clusters by the
#'    \code{\link{comp_clustering}} function.
#'
#'    Remember: when using the \code{\link{comp_clustering}} function, inspect
#'    the results of  three internal measures (connectivity index, silhouette
#'    width, and Dunn index) for a wide range of clusters to decide on the optimal
#'    number of clusters.
#'
#' @author {Loukia M. Spineli}
#'
#' @seealso
#'  \code{\link{comp_clustering}}, \code{\link[heatmaply:heatmaply]{heatmaply}}
#'
#' @export
dendro_heatmap <- function (input) {


  ## Check the defaults
  input <- if (inherits(input, "comp_clustering")) {
    input
  } else {
    stop("'input' must be of class 'comp_clustering'", call. = FALSE)
  }


  ## The dissimilarity matrix (based on the dissimilartity measure used in the
  ## comp_clustering function)
  diss <- input$Dissimilarity_table


  ## 'Optimal' dissimilarity method (based on the cophenetic coefficient)
  optimal_dist <- if (is.null(input$Optimal_dist)) {
    stop("The function can be used *only* for heuristic clustering.", call. = FALSE)
  } else {
    input$Optimal_dist
  }


  ## 'Optimal' linkage method (based on the cophenetic coefficient)
  optimal_link <- if (is.null(input$Optimal_link)) {
    stop("The function can be used *only* for heuristic clustering.", call. = FALSE)
  } else {
    input$Optimal_link
  }


  ## Number of 'optimal' clusters
  optimal_clusters <- length(unique(input$Cluster_color[, 2]))


  ## Function for first letter capital (Source: https://stackoverflow.com/questions/18509527/first-letter-to-upper-case)
  firstup <- function(x) {
    substr(x, 1, 1) <- toupper(substr(x, 1, 1))
    x
  }


  ## Turn first letter capital
  optimal_dist_new <- firstup(optimal_dist)


  ## Create heatmap with dendrogram and coloured clusters
  dendro_heatmap <-
    heatmaply(as.matrix(diss),
              cellnote = round(as.matrix(diss), 2),
              xlab = " ",
              ylab = " ",
              scale = "none",
              Colv = reorder(as.dendrogram(hclust(diss)),
                             input$Total_dissimilarity[, 3], max),
              Rowv = reorder(as.dendrogram(hclust(diss)),
                             input$Total_dissimilarity[, 3], max),
              k_col = optimal_clusters,
              k_row = optimal_clusters,
              scale_fill_gradient_fun =
                scale_fill_gradient2(name = paste(optimal_dist_new,
                                                  "dissimilarity"),
                                     low = "white",
                                     high = "red",
                                     na.value = "grey90"#,
                                     #limit = c(limits_scale[1],
                                     #          limits_scale[2]))
              ))

  return(dendro_heatmap)
}

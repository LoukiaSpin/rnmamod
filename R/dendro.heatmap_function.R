#' Dendrogram with amalgamated heatmap
#' (Comparisons' comparability for transitivity evaluation)
#'
#' @description
#'   \code{dendro_heatmap} creates a dendrogram alongside the heatmap of
#'   Gower dissimilarities among the trials in the network for a specific 
#'   linkage method and number of clusters. T
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
#'   dendrogram. The names of the leaves refer to the trials and corresponding 
#'   pairwise comparison.
#'
#'  @details
#'    The function inherits the linkage method and number of optimal clusters by 
#'    the \code{\link{comp_clustering}} function.
#'
#'    Remember: when using the \code{\link{comp_clustering}} function, inspect
#'    the average silhouette width for a wide range of clusters to decide on the 
#'    optimal number of clusters.
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


  ## The dissimilarity matrix (based on the across-comparison dissimilarities)
  diss <- as.dist(input$Trials_diss_table)


  ## 'Optimal' linkage method (based on the cophenetic coefficient)
  optimal_link <- if (is.null(input$Optimal_link)) {
    stop("The function can be used *only* for hierarchical clustering.", 
         call. = FALSE)
  } else {
    input$Optimal_link
  }


  ## Number of 'optimal' clusters
  optimal_clusters <- length(unique(input$Cluster_comp[, 2]))


  ## Create heatmap with dendrogram and coloured clusters
  dendro_heatmap <-
    heatmaply(as.matrix(diss),
              cellnote = round(as.matrix(diss), 2),
              xlab = " ",
              ylab = " ",
              scale = "none",
              k_col = optimal_clusters,
              k_row = optimal_clusters,
              scale_fill_gradient_fun =
                scale_fill_gradient2(name = " ",
                                     low = "white",
                                     high = "red",
                                     na.value = "grey90",
                                     midpoint = 0.0, 
                                     limits = c(0, 1)))

  return(dendro_heatmap)
}

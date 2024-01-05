#' Dendrogram with amalgamated heatmap
#' (Comparisons' comparability for transitivity evaluation)
#'
#' @description
#'   \code{dendro_heatmap} creates a dendrogram alongside the heatmap of
#'   Gower dissimilarities among the trials in the network for a specific
#'   linkage method and number of clusters.
#'
#' @param input An object of S3 class \code{\link{comp_clustering}}. See 'Value'
#'   in \code{\link{comp_clustering}}.
#' @param label_size A positive integer for the font size of the heatmap
#'   elements. \code{label_size} determines the size argument found in the
#'   geom's aesthetic properties in the R-package
#'   \href{https://CRAN.R-project.org/package=ggplot2}{ggplot2}.
#' @param axis_text_size A positive integer for the font size of row and column
#'   names of the heatmap. \code{axis_text_size} determines the axis.text
#'   argument found in the theme's properties in the R-package
#'   \href{https://CRAN.R-project.org/package=ggplot2}{ggplot2}.
#'
#' @return
#'   \code{dendro_heatmap} uses the \code{\link[heatmaply:heatmaply]{heatmaply}}
#'    function of the R-package
#'   \href{https://CRAN.R-project.org/package=heatmaply}{heatmaply} to create a
#'   cluster heatmap for a selected linkage method and number of clusters. The
#'   function uses different colours to indicate the clusters directly on the
#'   dendrogram, specified using the R-package
#'   \href{https://CRAN.R-project.org/package=dendextend}{dendextend}. The names
#'   of the leaves refer to the trials and corresponding pairwise comparison.
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
dendro_heatmap <- function (input,
                            label_size = 12,
                            axis_text_size = 10) {


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


  ## Prepare a data-frame with the 'study-comparison', 'cluster' and 'silhouette width'
  cluster_comp <- cbind(rownames(input$Trials_diss_table), input$Cluster_comp[, 2:3])
  colnames(cluster_comp)[1] <- "study"


  ## Get the labels after ordering
  ordered_labels <-
    labels(reorder(as.dendrogram(hclust(diss, optimal_link)), cluster_comp[, 2]))


  ## Order the cluster and corresponding colour
  ordered_clust_col <-
    cluster_comp[order(match(cluster_comp[, 1], ordered_labels), decreasing = FALSE), ]

  ## Now colour the branches with the correct colour
  row_dend <-
    reorder(as.dendrogram(hclust(as.dist(input$Trials_diss_table), optimal_link)),
            cluster_comp[, 2]) %>%
    set("branches_k_color",
        value = unique(hue_pal()(optimal_clusters)[ordered_clust_col[, 2]]),
        k = optimal_clusters)


  ## Create heatmap with dendrogram and coloured clusters
  dendro_heatmap <-
    heatmaply(as.matrix(diss),
              cellnote = round(as.matrix(diss), 2),
              seriate = "none",
              xlab = " ",
              ylab = " ",
              scale = "none",
              Rowv = row_dend,
              Colv = row_dend,
              scale_fill_gradient_fun =
                scale_fill_gradient2(name = " ",
                                     low = "white",
                                     high = "red",
                                     na.value = "grey90",
                                     midpoint = 0.0,
                                     limits = c(0, 1)),
              cellnote_size = label_size ,
              fontsize_row = axis_text_size,
              fontsize_col = axis_text_size)

  return(dendro_heatmap)
}

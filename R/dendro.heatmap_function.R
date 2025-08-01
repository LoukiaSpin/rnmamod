#' Dendrogram with amalgamated heatmap
#' (Comparisons' comparability for transitivity evaluation)
#'
#' @description
#'   \code{dendro_heatmap} creates a dendrogram alongside the heatmap of
#'   Gower dissimilarities among the trials in the network for a specific
#'   linkage method and number of clusters (Spineli et al., 2025).
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
#' @references
#' Spineli LM, Papadimitropoulou K, Kalyvas C. Exploring the Transitivity
#' Assumption in Network Meta-Analysis: A Novel Approach and Its Implications.
#' \emph{Stat Med} 2025;\bold{44}(7):e70068.
#' doi: 10.1002/sim.70068.
#'
#' @examples
#' \donttest{
#' # Fictional dataset
#' data_set <- data.frame(Trial_name = as.character(1:7),
#'                       arm1 = c("1", "1", "1", "1", "1", "2", "2"),
#'                       arm2 = c("2", "2", "2", "3", "3", "3", "3"),
#'                       sample = c(140, 145, 150, 40, 45, 75, 80),
#'                       age = c(18, 18, 18, 48, 48, 35, 35),
#'                       blinding = factor(c("yes", "yes", "yes", "no", "no", "no", "no")))
#'
#' # Apply hierarchical clustering (informative = FALSE)
#' hier <- comp_clustering(input = data_set,
#'                         drug_names = c("A", "B", "C"),
#'                         threshold = 0.13,  # General research setting
#'                         informative = FALSE,
#'                         optimal_clusters = 3,
#'                         get_plots = TRUE)
#'
#' # Create the dendrogram with integrated heatmap
#' dendro_heatmap(hier)
#' }
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
              fontsize_col = axis_text_size,
              colorbar_xanchor = "middle",
              colorbar_yanchor = "bottom",
              key.title = "Gower's dissimilarities")

  return(dendro_heatmap)
}

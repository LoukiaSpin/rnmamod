#' Internal measures for cluster validation
#' (Comparisons' comparability for transitivity evaluation)
#'
#' @description
#'   \code{internal_measures_plot} currently prepares the table with the results 
#'   of the average silhouette width for a range of clusters, and visualises the 
#'   results using a profile plot.
#'
#' @param input An object of 'dist' class. It is a lower off-diagonal matrix
#'   with the dissimilarities of all pairs of comparisons.
#' @param optimal_link A character string with values \code{"ward.D"},
#'   \code{"ward.D2"}, \code{"single"}, \code{"complete"}, \code{"average"},
#'   \code{"mcquitty"}, \code{"median"}, or \code{"centroid"} for the optimal
#'   linkage method, corresponding to the highest cophenetic correlation
#'   coefficient value.
#'
#' @return
#'   \code{internal_measures_plot} currently returns the following list of 
#'   elements:
#'   \item{Table_internal_measures}{A data-frame of the average silhouette width 
#'   for a range of 2 to P-1 clusters, with P being the number of comparisons.}
#'   \item{Internal_measures_panel}{A profile plot on the average silhouette 
#'   width for a range of 2 to P-1 clusters, with P being the number of 
#'   comparisons. The candidate optimal number of clusters is indicated with a 
#'   red point directly on the line.}
#'
#' @details
#'   \code{internal_measures_plot} also calls the function
#'   \code{\link{comp_clustering}} to define the argument \code{optimal_link} to
#'   create the silhouette plot for the selected number of clusters.
#'
#'   If the network has three observed comparisons,
#'   \code{internal_measures_plot} will return only the
#'   \code{Table_internal_measures}. This is because with only three observed
#'   comparisons, only two clusters can be considered by the internal measures.
#'
#'   \code{internal_measures_plot} is integrated in the function
#'   \code{\link{comp_clustering}}.
#'
#' @author {Loukia M. Spineli}
#'
#' @seealso
#'  \code{\link{comp_clustering}}, \code{\link{silhouette_index}}
#'
#' @references
#' Handl J, Knowles J, Kell DB. Computational cluster validation in post-genomic
#' data analysis. \emph{Biometrics} 2005;\bold{21}(15):3201--120.
#' doi: 10.1093/bioinformatics/bti517
#'
#' Rousseeuw PJ. Silhouettes: A graphical aid to the interpretation and
#' validation of cluster analysis.
#' \emph{J Comput Appl Math} 1987;\bold{20}:53--65.
#'
#' @export
internal_measures_plot <- function (input,
                                    optimal_link) {


  ## Check the defaults
  input <- if (inherits(input, "dist")) {
    input
  } else {
    stop("'input' must be of class 'dist'", call. = FALSE)
  }


  ## 'Optimal' linkage method (based on the cophenetic coefficient)
  methods_list <- c("ward.D", "ward.D2", "single", "complete", "average",
                    "mcquitty", "median", "centroid")
  m_list1 <- c("'ward.D', 'ward.D2', 'single', 'complete', 'average'")
  m_list2 <- c("'mcquitty', 'median', 'centroid'")
  optimal_link <- if (missing(optimal_link)) {
    stop("The argument 'optimal_link' must be defined", call. = FALSE)
  } else if (!is.element(optimal_link, methods_list)) {
    stop(paste("'optimal_link' must be any of the following:",
               m_list1, m_list2), call. = FALSE)
  } else {
    optimal_link
  }


  ## Create data.frame with the number of clusters and internal measures
  internal_meas_res <- data.frame(clusters = 2:(dim(as.matrix(input))[1] - 1))

  
  ## Obtain the average silhouette width for all combinations
  internal_meas_res$silhouette <-
    mapply(function(x) silhouette_index(input = input,
                                        method = optimal_link,
                                        num_clusters = x)$silhoutte_width,
           2:(dim(as.matrix(input))[1] - 1))


  ## Plots results for Silhouette
  if (dim(internal_meas_res)[1] > 1) {
    plot_silhouette <-
      ggplot(internal_meas_res,
             aes(x = clusters,
                 y = silhouette)) +
      geom_line(linewidth = 1.2) +
      geom_segment(data = subset(internal_meas_res,
                                 silhouette == max(silhouette)),
                   aes(clusters, 0, xend = clusters, yend = silhouette),
                   linewidth = 0.8,
                   linetype = 3) +
      geom_point(data = subset(internal_meas_res, silhouette == max(silhouette)),
                 aes(x = clusters,
                     y = silhouette),
                 colour = "red",
                 size = 5) +
      geom_point(data = subset(internal_meas_res, silhouette != max(silhouette)),
                 aes(x = clusters,
                     y = silhouette),
                 colour = "black",
                 size = 2.5) +
      geom_text(aes(label = sprintf("%0.2f", round(silhouette, 2))),
                hjust = 0.5,
                vjust = -0.7,
                size = 4,
                colour = "blue",
                fontface = "bold") +
      labs(x = "Number of clusters",
           y = "Silhouette width") +
      scale_x_continuous(breaks = seq(1, dim(as.matrix(input))[1], 1),
                         labels = seq(1, dim(as.matrix(input))[1], 1)) +
      theme_classic() +
      theme(axis.title = element_text(size = 14, face = "bold"),
            axis.text = element_text(size = 14))
  }


  ## Collect results
  results <- if (dim(internal_meas_res)[1] > 1) {
    list(Table_internal_measures = round(internal_meas_res, 3),
         Internal_measures_panel = plot_silhouette)
  } else {
    list(Table_internal_measures = round(internal_meas_res, 3))
  }

  return(results)
}

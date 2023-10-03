#' Internal measures for cluster validation
#' (Comparisons' comparability for transitivity evaluation)
#'
#' @description
#'   \code{internal_measures_plot} prepares the table with the results of three
#'   internal measures (connectivity index, silhouette width, and Dunn index)
#'   for a range of clusters, and visualises the results using a profile plot
#'   for each internal measure.
#'
#' @param input An object of 'dist' class. It is a lower off-diagonal matrix
#'   with the dissimilarities of all pairs of comparisons.
#' @param num_neighb A positive integer for the number of neighbouring
#'   comparisons. It takes values from two to the number of comparisons minus
#'   one. The default argument equals half the number of comparisons.
#' @param optimal_link A character string with values \code{"ward.D"},
#'   \code{"ward.D2"}, \code{"single"}, \code{"complete"}, \code{"average"},
#'   \code{"mcquitty"}, \code{"median"}, or \code{"centroid"} for the optimal
#'   linkage method, corresponding to the highest cophenetic correlation
#'   coefficient value.
#'
#' @return
#'   \code{internal_measures_plot} returns the following list of elements:
#'   \item{Table_internal_measures}{A data-frame of the connectivity index,
#'   silhouette width, and Dunn index for a range of 2 to P-1 clusters, with P
#'   being the number of comparisons.}
#'   \item{Internal_measures_panel}{A panel of profile plots on the connectivity
#'   index, silhouette width, and Dunn index for a range of 2 to P-1 clusters,
#'   with P being the number of comparisons. The candidate optimal number of
#'   clusters is indicated with a red point directly on the line.}
#'
#' @details
#'   \code{internal_measures_plot} call the functions
#'   \code{\link{connectivity_index}}, \code{\link{silhouette_index}} and
#'   \code{\link{dunn_index}} to calculate the corresponding internal measures.
#'   \code{internal_measures_plot} also call the function
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
#'  \code{\link{comp_clustering}}, \code{\link{connectivity_index}},
#'  \code{\link{dunn_index}}, \code{\link{silhouette_index}}
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
#' Dunn J. Well-separated clusters and optimal fuzzy partitions.
#' \emph{J Cybern} 1974;\bold{4}(1):95--104.
#'
#' @export
internal_measures_plot <- function (input,
                                    num_neighb,
                                    optimal_link) {


  ## Check the defaults
  # Dataset
  input <- if (inherits(input, "dist")) {
    input
  } else {
    stop("'input' must be of class 'dist'", call. = FALSE)
  }

  # Default to be used in 'connectivity_index'
  num_neighb <- if (missing(num_neighb)) {
    message(paste("num_neighb =", round(dim(as.matrix(input))[1] / 2, 0),
                  "was used (default)"))
    round(dim(as.matrix(input))[1] / 2, 0)
  } else if (num_neighb > dim(as.matrix(input))[1] || num_neighb < 2) {
    stop(paste0("'num_neighb' must range from 2 to", " ",
                dim(as.matrix(input))[1] - 1, "."), call. = FALSE)
  } else {
    num_neighb
  }

  # 'Optimal' linkage method (based on the cophenetic coefficient)
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


  ## Obtain connectivity for all combinations
  internal_meas_res$connectivity <-
    mapply(function(x) connectivity_index(input = input,
                                          method = optimal_link,
                                          num_clusters = x,
                                          num_neighb = num_neighb),
           2:(dim(as.matrix(input))[1] - 1))


  ## Obtain silhouette for all combinations
  internal_meas_res$silhouette <-
    mapply(function(x) silhouette_index(input = input,
                                        method = optimal_link,
                                        num_clusters = x)$silhoutte_width,
           2:(dim(as.matrix(input))[1] - 1))


  ## Obtain Dunn for all combinations
  internal_meas_res$dunn <-
    mapply(function(x) dunn_index(input = input,
                                  method = optimal_link,
                                  num_clusters = x),
           2:(dim(as.matrix(input))[1] - 1))


  ## Plots results for Connectivity
  if (dim(internal_meas_res)[1] > 1) {
    plot_connectivity <-
      ggplot(internal_meas_res,
             aes(x = clusters,
                 y = connectivity)) +
      geom_line(linewidth = 1.2) +
      geom_segment(data = subset(internal_meas_res,
                                 connectivity == min(connectivity) &
                                   connectivity > 0),
                   aes(clusters, 0, xend = clusters, yend = connectivity),
                   linewidth = 0.8,
                   linetype = 3) +
      geom_point(data = subset(internal_meas_res,
                               connectivity == min(connectivity) &
                                 connectivity > 0),
                 aes(x = clusters,
                     y = connectivity),
                 colour = "red",
                 size = 5) +
      geom_point(data = subset(internal_meas_res,
                               connectivity != min(connectivity)),
                 aes(x = clusters,
                     y = connectivity),
                 colour = "black",
                 size = 2.5) +
      geom_text(aes(label = sprintf("%0.2f", round(connectivity, 2))),
                hjust = 0.5,
                vjust = -0.7,
                size = 4,
                colour = "blue",
                fontface = "bold") +
      labs(x = " ",
           y = "Connectivity index") +
      scale_x_continuous(breaks = seq(1, dim(as.matrix(input))[1], 1),
                         labels = seq(1, dim(as.matrix(input))[1], 1)) +
      theme_classic() +
      theme(axis.title = element_text(size = 14, face = "bold"),
            axis.text = element_text(size = 14))
  }


  ## Plots results for Dunn
  if (dim(internal_meas_res)[1] > 1) {
    plot_dunn <-
      ggplot(internal_meas_res,
             aes(x = clusters,
                 y = dunn)) +
      geom_line(linewidth = 1.2) +
      geom_segment(data = subset(internal_meas_res, dunn == max(dunn)),
                   aes(clusters, 0, xend = clusters, yend = dunn),
                   linewidth = 0.8,
                   linetype = 3) +
      geom_point(data = subset(internal_meas_res, dunn == max(dunn)),
                 aes(x = clusters,
                     y = dunn),
                 colour = "red",
                 size = 5) +
      geom_point(data = subset(internal_meas_res, dunn != max(dunn)),
                 aes(x = clusters,
                     y = dunn),
                 colour = "black",
                 size = 2.5) +
      geom_text(aes(label = sprintf("%0.2f", round(dunn, 2))),
                hjust = 0.5,
                vjust = -0.7,
                size = 4,
                colour = "blue",
                fontface = "bold") +
      labs(x = "Number of clusters",
           y = "Dunn index") +
      scale_x_continuous(breaks = seq(1, dim(as.matrix(input))[1], 1),
                         labels = seq(1, dim(as.matrix(input))[1], 1)) +
      theme_classic() +
      theme(axis.title = element_text(size = 14, face = "bold"),
            axis.text = element_text(size = 14))
  }


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


  ## Bring plots for all internal measures together
  if (dim(internal_meas_res)[1] > 1) {
    internal_measures_plot <-
      ggarrange(plot_connectivity, plot_silhouette, plot_dunn,
                nrow = 2, ncol = 2,
                common.legend = TRUE,
                legend = "bottom")
  }


  ## Collect results
  results <- if (dim(internal_meas_res)[1] > 1) {
    list(Table_internal_measures = round(internal_meas_res, 3),
         Internal_measures_panel = internal_measures_plot)
  } else {
    list(Table_internal_measures = round(internal_meas_res, 3))
  }

  return(results)
}

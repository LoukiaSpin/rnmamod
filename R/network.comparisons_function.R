#' Network of pairwise comparisons of interventions
#' (Comparisons' comparability for transitivity evaluation)
#'
#' @description
#'   \code{network_comparisons} offers an alternative visualisation of the
#'   dendrogram with amalgamated heatmap. Presenting the comparisons in a
#'   circular layout, typically seen in networks, facilitates interpretation of
#'   the clustering results in the context of transitivity evaluation. It is
#'   relevant *only* for heuristic clustering.
#'
#' @param data An object of S3 class \code{\link{comp_clustering}}. See 'Value'
#'   in \code{\link{comp_clustering}}.
#' @param node_frame_color A character string for the colour of the node frames.
#'   The default argument is "black". For more details refer to the R-package
#'   \href{https://CRAN.R-project.org/package=igraph}{igraph}.
#' @param node_frame_width A positive integer for the size of the node frame.
#'   The default argument is 1. For more details refer to the R-package
#'   \href{https://CRAN.R-project.org/package=igraph}{igraph}.
#' @param node_label_font A positive integer for the font type of the node
#'   labels. The default argument is 1. For more details refer to the R-package
#'   \href{https://CRAN.R-project.org/package=igraph}{igraph}.
#' @param node_label_cex A positive integer for the size of the node labels.
#'   The default argument is 1. For more details refer to the R-package
#'   \href{https://CRAN.R-project.org/package=igraph}{igraph}.
#' @param node_label_dist A positive integer for the distance between label and
#'   node. The default argument is 0. For more details refer to the R-package
#'   \href{https://CRAN.R-project.org/package=igraph}{igraph}.
#' @param node_label_color A character string for the colour of the node labels.
#'   The default argument is "black". For more details refer to the R-package
#'   \href{https://CRAN.R-project.org/package=igraph}{igraph}.
#' @param edge_level A character string with values \code{"low"}, \code{"low"},
#'   \code{"moderate"}, \code{"high"}, and \code{"very high"} to indicate the
#'   level of dissimilarity between two comparisons to draw.
#'   See 'Details' below.
#' @param edge_lty A positive integer for the line type of the edges. The
#'   default argument is 1. For more details refer to the R-package
#'   \href{https://CRAN.R-project.org/package=igraph}{igraph}.
#' @param edge_curved A value in the range 0 to 1 for the edge curvature. The
#'   default argument is 0. For more details refer to the R-package
#'   \href{https://CRAN.R-project.org/package=igraph}{igraph}.
#' @param ... Further graphical arguments of the
#'   \code{\link[igraph:plot.igraph]{plot.igraph}} function of the
#'   R-package \href{https://CRAN.R-project.org/package=igraph}{igraph}
#'
#' @return
#'   The network of comparisons with different colours for the clustered nodes.
#'
#' @details
#'    The nodes refer to the observed comparisons and their size is proportional
#'    to the corresponding normalised total dissimilarity. The nodes are
#'    coloured with different colours to indicate the cluster they belong.
#'    The clusters and dissimilarity measure are inherited by
#'    \code{\link{comp_clustering}}. The edges refer to the dissimilarities of
#'    pairs of comparisons and their size is proportional to
#'    the corresponding normalised dissimilarity value. Four different colours
#'    are used to indicate the level of dissimilarity between two comparisons:
#'    green, yellow, orange, and red for low, moderate, high and very high
#'    dissimilarity (\eqn{d < 0.25}, \eqn{0.25 \qeq d < 0.5},
#'    \eqn{0.50 \qeq d < 0.75}, and \eqn{d \qeq 0.75}, respectively).
#'
#'    Networks with many comparisons appear cluttered when showing all edges
#'    (i.e., \code{edge_level = "all"}). In this case, it is recommended to draw
#'    the network for the dissimilarity level of interest, for instance,
#'    \code{edge_level = "low"}). When a dissimilarity level does not exist, the
#'    network will be drawn without the corresponding edges.
#'
#' @author {Loukia M. Spineli}
#'
#' @seealso
#'  \code{\link{comp_clustering}},
#'  \code{\link[igraph:plot.igraph]{plot.igraph}}
#'
#' @export
network_comparisons <- function(data,
                                node_frame_color = "black",
                                node_frame_width = 1,
                                node_label_font = 1,
                                node_label_cex = 1,
                                node_label_dist = 0,
                                node_label_color = "black",
                                edge_level,
                                edge_lty = 1,
                                #edge_label_color = "black",
                                #edge_label_font = 1,
                                #edge_label_cex = 1,
                                edge_curved = 0,
                                ...) {


  # Check the defaults
  data <- if (inherits(data, "comp_clustering")) {
    data
  } else {
    stop("'input' must be of class 'comp_clustering'", call. = FALSE)
  }


  # Select the distances for the max cophenetic coefficient
  optimal_dist_link <- subset(data$Table_cophenetic_coefficient,
                              results == max(results))


  # When more distances are proper for the same cophenetic coeff.
  if (length(unique(optimal_dist_link[, 1])) > 1) {
    optimal_dist <- optimal_dist_link[1, 1]
  } else if (length(unique(optimal_dist_link[, 1])) == 1) {
    optimal_dist <- unique(optimal_dist_link[, 1])
  } else if (dim(optimal_dist_link[, 1])[1] == 1) {
    optimal_dist <- optimal_dist_link[1]
  }


  # Number of 'optimal' clusters
  optimal_clusters <- length(unique(data$Cluster_color[, 2]))


  # Set the default for 'edge_level'
  levels <- c("low", "moderate", "high", "very high", "all")
  edge_level <- if (missing(edge_level)) {
    a <- c("'low', 'moderate', 'high', 'very high', 'all'")
    stop(paste("The argument 'edge_level' must be one of the following:", a))
  } else if (!is.element(edge_level, levels)) {
    a <- c("'low', 'moderate', 'high', 'very high', 'all'")
    stop(paste("The argument 'edge_level' must be one of the following:", a))
  } else {
    edge_level
  }


  # Get  comparison names sorted in decreasing order of total dissimilarity
  cluster_comp <-
    levels(reorder(data$Total_dissimilarity[, 1], data$Total_dissimilarity[, 3],
                   decreasing = TRUE))


  # Possible pairwise comparisons
  poss_comp <- combn(cluster_comp, 2)


  # Sort total dissimilarity in decreasing order
  total_diss <- sort(data$Total_dissimilarity[, 3], decreasing = TRUE)


  # Match the comparison names of 'Total_dissimilarity' with those of 'Cluster_color'
  cluster_col <-
    data$Cluster_color[match(cluster_comp, data$Cluster_color[, 1]),]


  # Colour node by cluster
  node_col <- cluster_col[, 3]


  # Prepare plot (igraph package)
  g1 <- igraph::graph(edges = poss_comp, directed = FALSE)


  # Weight each node by the corresponding total distance
  igraph::V(g1)$weight <- (0.40 + total_diss) * 20


  # Possible comparisons among comparisons
  poss_comp_comp <- combn(unique(c(poss_comp)), 2)


  # Sort comparisons of comparisons (edges) by the order of 'poss_comp_comp'
  comp_diss_new <- rep(NA, dim(poss_comp_comp)[2])
  for (i in 1:dim(poss_comp_comp)[2]) {
    comp_diss_new[i] <- as.matrix(data$Dissimilarity_table)[poss_comp_comp[1, i], poss_comp_comp[2, i]]
  }


  # Normalise the dissimilarity matrix of comparisons
  norm_comp_diss <- (comp_diss_new - min(comp_diss_new)) /
    (max(comp_diss_new) - min(comp_diss_new))


  # Weight each edge by the corresponding distance (normalised)
  igraph::E(g1)$weight <- (0.40 + norm_comp_diss) * 10


  # Name the edges according to dissimilarity value
  edge_col_num <-
    as.vector(do.call(rbind,
                      lapply(norm_comp_diss,
                             function(x) if (x <= 0.25) "low"
                             else if (x > 0.25 & x <= 0.50) "moderate"
                             else if (x > 0.50 & x <= 0.75) "high"
                             else "very high")))


  # Assign a colour to the name
  edge_colour <-
    as.vector(do.call(rbind,
                      lapply(edge_col_num,
                             function(x) if (x == "low") "#A6D854"
                             else if (x == "moderate") "#E6AB02"
                             else if (x == "high") "#D95F02"
                             else "#E31A1C")))


  # Adjust the selected level for larger networks (at least 5 comparisons)
  edge_col_select0 <- subset(edge_colour, edge_col_num == edge_level)
  edge_col_select <- ifelse(!is.element(edge_colour, edge_col_select0),
                            NA,
                            edge_col_select0)


  # Colour the edges based on the network size
  if (edge_level != "all") {

    # Adjust the selected level for larger networks (at least 5 comparisons)
    edge_col0 <- subset(edge_colour, edge_col_num == edge_level)
    edge_col <- ifelse(!is.element(edge_colour, edge_col0), NA, edge_col0)
  } else {
    edge_col <- edge_colour
  }


  # Create the igraph
  plot(g1,
       layout = layout_in_circle,
       vertex.color = node_col,
       vertex.frame.color = node_frame_color,
       vertex.frame.width = node_frame_width,
       vertex.shape = "circle",
       vertex.size = igraph::V(g1)$weight,
       vertex.label.font = node_label_font,
       vertex.label.cex	= node_label_cex,
       vertex.label.dist = node_label_dist,
       vertex.label.color = node_label_color,
       edge.color = edge_col,
       edge.width = igraph::E(g1)$weight,
       edge.lty = edge_lty,
       #edge.label = edge_name, #E(g1)$names,
       #edge.label.color = edge_label_color,
       #edge.label.font = edge_label_font,
       #edge.label.cex = edge_label_cex,
       edge.curved = edge_curved,
       ...)
}

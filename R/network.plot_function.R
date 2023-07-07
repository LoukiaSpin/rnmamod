#' Network plot
#'
#' @description Illustrates the network plot for one outcome.
#'
#' @param data A data-frame of a one-trial-per-row format containing arm-level
#'   data of each trial. See 'Format' in \code{\link{run_model}}.
#' @param drug_names A vector of labels with the name of the interventions
#'   (nodes) in the order they appear in the argument \code{data}.
#' @param show_multi Logical to indicate whether to colour the closed-loops
#'   informed by multi-arm trials. The default is \code{show_multi = FALSE}.
#' @param multi_frame A numeric scalar to determine the size of the border
#'   around the closed-loops formed by multi-arm trials. The default is -16.
#'  \code{multi_frame} determines the mark.expand argument found in the
#'   \code{\link[igraph:plot.igraph]{plot.igraph}} function of the
#'   R-package \href{https://CRAN.R-project.org/package=igraph}{igraph}).
#' @param alpha_multi_color A numerical scalar with range from 0 to 1 to
#'   determine the opacity of \code{multi_frame} coloured using the
#'   \code{rainbow} colour palette. The default is 0.1.
#' @param layout The layout specification. The default is
#'   \code{layout = layout_in_circle} to plot the nodes in a circular layout.
#'   For more information refer to the R-package
#'   \href{https://CRAN.R-project.org/package=igraph}{igraph},
#' @param node_color A character or vector of characters (with length equal to
#'   the number of nodes) to indicate the colour of the nodes. The default is
#'   \code{node_color = "tomato"}. \code{node_color} determines the vertex.color
#'   argument found in the \code{\link[igraph:plot.igraph]{plot.igraph}}
#'   function of the R-package
#'   \href{https://CRAN.R-project.org/package=igraph}{igraph}).
#' @param node_frame_color A character to indicate the colour of the frame
#'   around the nodes. The default is \code{node_frame_color = "black"}.
#'   \code{node_frame_color} determines the vertex.frame.color argument found in
#'   the \code{\link[igraph:plot.igraph]{plot.igraph}} function of the
#'   R-package \href{https://CRAN.R-project.org/package=igraph}{igraph}).
#' @param node_frame_width A numerical scalar to indicate the width of the frame
#'   around the nodes. The default is 1. \code{node_frame_width} determines the
#'   vertex.frame.width argument found in the
#'   \code{\link[igraph:plot.igraph]{plot.igraph}} function of the
#'   R-package \href{https://CRAN.R-project.org/package=igraph}{igraph}).
#' @param node_shape A character to indicate the shape of the nodes. The default
#'   is \code{node_shape = "circle"}. \code{node_shape} determines the
#'   vertex.shape argument found in the
#'   \code{\link[igraph:plot.igraph]{plot.igraph}} function of the
#'   R-package \href{https://CRAN.R-project.org/package=igraph}{igraph}).
#' @param node_label_color A character to indicate the color of the node labels.
#'   The default is \code{node_label_color = "black"}. \code{node_label_color}
#'   determines the vertex.label.color argument found in the
#'   \code{\link[igraph:plot.igraph]{plot.igraph}} function of the
#'   R-package \href{https://CRAN.R-project.org/package=igraph}{igraph}).
#' @param node_label_font A numerical scalar to indicate the font of the node
#'   labels. The default is 1. \code{node_label_font} determines the
#'   vertex.label.font argument found in the
#'   \code{\link[igraph:plot.igraph]{plot.igraph}} function of the
#'   R-package \href{https://CRAN.R-project.org/package=igraph}{igraph}).
#' @param node_label_cex A numerical scalar to indicate the font size of the
#'   node labels. The default is 1. \code{node_label_cex} determines the
#'   vertex.label.cex argument found in the
#'   \code{\link[igraph:plot.igraph]{plot.igraph}} function of the
#'   R-package \href{https://CRAN.R-project.org/package=igraph}{igraph}).
#' @param node_label_dist A numerical scale between 0 and 1 to indicate the
#'   position of the node labels relative to the node center. The default is 0,
#'   where the label is centered. \code{node_label_dist} determines the
#'   vertex.label.dist argument found in the
#'   \code{\link[igraph:plot.igraph]{plot.igraph}} function of the
#'   R-package \href{https://CRAN.R-project.org/package=igraph}{igraph}).
#' @param edge_color A character or vector of characters (with length equal to
#'   the number of edges) to indicate the colour of the edges. The default is
#'   \code{edge_color = "grey50"}. \code{edge_color} determines the edge.color
#'   argument found in the \code{\link[igraph:plot.igraph]{plot.igraph}}
#'   function of the R-package
#'   \href{https://CRAN.R-project.org/package=igraph}{igraph}).
#' @param edge_arrow_size A numerical scalar between 0 and 1 to indicate the
#'   arrow size. The default is 0.5. \code{edge_arrow_size} determines the
#'   edge.arrow.size argument found in the
#'   \code{\link[igraph:plot.igraph]{plot.igraph}} function of the R-package
#'   \href{https://CRAN.R-project.org/package=igraph}{igraph}). This argument
#'   work only when \code{direction = FALSE}.
#' @param edge_lty A numerical scalar, discrete with values from 0 to 6 to
#'   indicate the line type of the edges. The default is 1 (solid).
#'   \code{edge_lty} determines the edge.lty argument found in the
#'   \code{\link[igraph:plot.igraph]{plot.igraph}} function of the R-package
#'   \href{https://CRAN.R-project.org/package=igraph}{igraph}).
#' @param edge_label A vector of number of characters with length equal to the
#'   number of edges to present the edge label. The default is
#'   \code{edge_label = NULL} and refers to the number of studies investigating
#'   the corresponding comparisons. \code{edge_label} determines the edge.label
#'   argument found in the \code{\link[igraph:plot.igraph]{plot.igraph}}
#'   function of the R-package
#'   \href{https://CRAN.R-project.org/package=igraph}{igraph}).
#' @param edge_label_color A character to indicate the color of the edge labels.
#'   The default is \code{edge_label_color = "black"}. \code{edge_label_color}
#'   determines the edge.label.color argument found in the
#'   \code{\link[igraph:plot.igraph]{plot.igraph}} function of the
#'   R-package \href{https://CRAN.R-project.org/package=igraph}{igraph}).
#' @param edge_label_font A numerical scalar to indicate the font of the edge
#'   labels. The default is 1. \code{edge_label_font} determines the
#'   edge.label.font argument found in the
#'   \code{\link[igraph:plot.igraph]{plot.igraph}} function of the
#'   R-package \href{https://CRAN.R-project.org/package=igraph}{igraph}).
#' @param edge_label_cex A numerical scalar to indicate the font size of the
#'   edge labels. The default is 2. \code{edge_label_cex} determines the
#'   edge.label.cex argument found in the
#'   \code{\link[igraph:plot.igraph]{plot.igraph}} function of the
#'   R-package \href{https://CRAN.R-project.org/package=igraph}{igraph}).
#' @param edge_curved A numerical scalar with range from 0 to 1 that indicates
#'   the edge curvature. The default is 0 (no curvature). \code{edge_curved}
#'   determines the edge.curved argument found in the
#'   \code{\link[igraph:plot.igraph]{plot.igraph}} function of the
#'   R-package \href{https://CRAN.R-project.org/package=igraph}{igraph}).
#' @param direction Logical to draw (TRUE) or not (FALSE) arrow for each edge
#'   according to each direction. The default is \code{direction = FALSE}. For
#'   more information refer to the R-package
#'   \href{https://CRAN.R-project.org/package=igraph}{igraph}).
#' @param ... Further graphical arguments of the
#'   \code{\link[igraph:plot.igraph]{plot.igraph}} function of the
#'   R-package \href{https://CRAN.R-project.org/package=igraph}{igraph}
#'
#' @return A network plot with coloured closed-loops informed by multi-arm
#'   trials. Each node indicates an intervention and each edge an observed
#'   pairwise comparison.
#'
#' @details The edge thickness is proportional to the number of
#'   trials investigating the corresponding comparison. The node size is
#'   weighted by the total sample size of the corresponding intervention.
#'
#'   The user can control many of the arguments found in the
#'   \code{\link[igraph:plot.igraph]{plot.igraph}} function of the
#'   R-package \href{https://CRAN.R-project.org/package=igraph}{igraph}), such
#'   as the colour of the nodes and edges, the node and edge label size, and so
#'   on.
#'
#' @seealso \code{\link[igraph:plot.igraph]{plot.igraph}},
#'   \code{\link{run_model}}
#'
#' @author {Loukia M. Spineli}
#'
#' @examples
#' data("nma.bottomley2011")
#'
#' # Return the first six trials of the dataset
#' head(nma.bottomley2011)
#'
#' # The names of the interventions in the order they appear in the dataset
#' interv_names <- c("betamethasone dipropionate", "betamethasone valerate",
#'                   "calcipotriol", "calcipotriol plus polytar", "capasal",
#'                   "two-compound formulation gel", "placebo")
#'
#' # Create the network plot
#' netplot(data = nma.bottomley2011,
#'         drug_names = interv_names,
#'         edge_label_cex = 1)
#'
#' @export
netplot <- function(data,
                    drug_names,
                    show_multi = FALSE,
                    multi_frame = -16,
                    alpha_multi_color = 0.1,
                    layout = igraph::layout_in_circle,
                    node_color = "tomato",
                    node_frame_color = "black",
                    node_frame_width = 1,
                    node_shape = "circle",
                    node_label_color = "black",
                    node_label_font = 1,
                    node_label_cex = 1,
                    node_label_dist = 0,
                    edge_color = "grey50",
                    edge_arrow_size = 0.5,
                    edge_lty = 1,
                    edge_label = NULL,
                    edge_label_color = "black",
                    edge_label_font = 1,
                    edge_label_cex = 2,
                    edge_curved = 0,
                    direction = FALSE,
                    ...) {

  # Extract the columns with the intervention id
  treat <- if (dim(data[, startsWith(colnames(data), "t")])[2] == 0) {
    stop("The information on the individual arms is missing", call. = FALSE)
  } else {
    data[, startsWith(colnames(data), "t")]
  }

  # Intervention names
  drug_names <- if (missing(drug_names)) {
    1:length(unique(na.omit(unlist(treat))))
  } else {
    drug_names
  }

  # Extract the columns with the sample size
  sample.size <- data[, startsWith(colnames(data), "n")]

  # Rename interventions in 'treat'
  treat_names <- matrix(drug_names[as.numeric(unlist(treat))],
                        nrow = dim(treat)[1],
                        ncol = dim(treat)[2])

  # Turn wide-format comparisons into long-format comparisons
  # Sort so that smaller treatment id refer to starting node and larger to end node
  pairwise_comb <- matrix(unlist(apply(treat, 1, function(x) apply(combn(na.omit(x), 2), 2, sort))),
                          ncol = 2,
                          byrow = TRUE)

  # Assign the intervention names (if applicable)
  pairwise <- matrix(drug_names[as.numeric(unlist(pairwise_comb))],
                  nrow = dim(pairwise_comb)[1],
                  ncol = 2)

  # Keep unique comparisons and turn into vector (small)
  start_to_end <- c(t(pairwise[!duplicated(pairwise), ]))

  # Prepare plot (igraph package)
  g1 <- igraph::graph(edges = start_to_end, directed = direction)

  # Calculate the total sample size per node (intervention)
  node.size0 <- aggregate(na.omit(unlist(sample.size)),
                          by = list(na.omit(c(treat_names))),
                          sum)

  # Sort by 'start_to_end'
  node.size <- node.size0[match(unique(start_to_end), node.size0[, 1]), 2]

  # Weight each node by the corresponding total sample size
  igraph::V(g1)$weight <- (0.40 + ((node.size - min(node.size)) /
                                     (max(node.size) - min(node.size)) )) * 20

  # Name the nodes
  igraph::V(g1)$names <- V(g1)

  # Calculate the number of trials per edge (pairwise comparison)
  comp0 <-
    data.frame(table(apply(pairwise, 1, function(x) paste(x[2], "vs", x[1]))))
  colnames(comp0) <- c("comparison", "frequency")

  # Sort comparisons by the order of 'start_to_end'
  comp <- comp0[match(apply(matrix(start_to_end,
                                   ncol = 2,
                                   byrow = TRUE),
                            1, function(x) paste(x[2], "vs", x[1])),
                      comp0$comparison), ]

  # Multi-arm trials: Count the number of arms per trial
  na_arms <- apply(treat, 1, function(x) length(na.omit(x)))

  # Multi-arm trials: If there are multi-arm trials ...
  show_multi_arms <- if (max(na_arms) > 2 & show_multi == TRUE) {

    # ... reduce the dataset to multi-arm trials only
    multi_trials <- subset(treat_names, na_arms > 2)

    # ... create a list of nodes id for the corresponding multi-arm trials.
    list.multi <- lapply(as.list(data.frame(t(multi_trials))), function(x) x[!is.na(x)])
    list.multi[!duplicated(list.multi)]
  } else {
    NULL
  }

  # Weight each edge by the corresponding number of trials
  igraph::E(g1)$weight <- comp$frequency

  # Label the edges (default is 'NULL')
  igraph::E(g1)$names <- if (is.null(edge_label)) {
    comp$frequency
  } else {
    edge_label
  }

  # Color edges (sorting the order of 'start_to_end')
  edge_color_new <- if (missing(edge_color)) {
    "grey50"
  } else if (length(edge_color) > 1) {
    edge_color[match(apply(matrix(start_to_end,
                                  ncol = 2,
                                  byrow = TRUE),
                           1, function(x) paste(x[2], "vs", x[1])),
                     edge_color$comparison), 3]
  } else if (length(edge_color) == 1) {
    edge_color
  }

  # Get the network plot
  plot(g1,
       layout = layout,
       vertex.color = node_color,
       vertex.frame.color = node_frame_color,
       vertex.frame.width = node_frame_width,
       vertex.shape = node_shape,
       vertex.size = igraph::V(g1)$weight,
       vertex.label.color = node_label_color,
       vertex.label.font = node_label_font,
       vertex.label.cex	= node_label_cex,
       vertex.label.dist = node_label_dist,
       edge.color = edge_color_new,
       edge.width = igraph::E(g1)$weight,
       edge.arrow.size = edge_arrow_size,
       edge.lty = edge_lty,
       edge.label = igraph::E(g1)$names,
       edge.label.color = edge_label_color,
       edge.label.font = edge_label_font,
       edge.label.cex = edge_label_cex,
       edge.curved = edge_curved,
       mark.groups = show_multi_arms,
       mark.col = rainbow(length(show_multi_arms),
                          alpha = alpha_multi_color),
       mark.shape = 0,
       mark.border = NA,
       mark.expand = multi_frame,
       ...)
}

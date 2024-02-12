find_subnetworks <- function (data) {

  # Extract the columns with the intervention id
  treat <- if (dim(data[, startsWith(colnames(data), "t")])[2] == 0) {
    stop("The information on the individual arms is missing", call. = FALSE)
  } else {
    data[, startsWith(colnames(data), "t")]
  }

  # Turn wide-format comparisons into long-format comparisons
  # Sort so that smaller treatment id refer to starting node and larger to end node
  pairwise <- matrix(unlist(apply(matrix(unlist(treat), ncol = 2, byrow = FALSE), 1, function(x) apply(combn(na.omit(x), 2), 2, sort))),
                     ncol = 2,
                     byrow = TRUE)

  # Keep unique comparisons and turn into vector (small)
  start_to_end <- as.character(c(t(pairwise[!duplicated(pairwise), ])))

  # Prepare plot (igraph package)
  g1 <- igraph::graph(edges = start_to_end, directed = FALSE)

  # Number of subnetworks
  num_subnetworks <- igraph::count_components(g1)

  # Obtain distance matrix
  distance_mat <- igraph::distances(g1)
  colnames(distance_mat) <- sort(unique(as.vector(pairwise)))
  rownames(distance_mat) <- colnames(distance_mat)

  # Return results
  return(list(num_subnetworks = num_subnetworks,
              distance_mat = distance_mat))
}


#' Assign the original name of the interventions to their number id
#'
#' @param drug.names An object of S3 class \code{run.sensitivity}.
#' Indicate all possible comparisons in the network.
#' The interventions in the drug.names should follow align with the id you considered for the function 'run.model'
#' @export
drug.names.alloc <- function(drug.names) {


  ## Number of interventions
  nt <- length(drug.names)



  ## Obtain a matrix with all unique parwise comparisons
  comparison <- matrix(combn(drug.names, 2), nrow = length(combn(drug.names, 2))/2, ncol = 2, byrow = T)



  ## Assign the original name of each intervention to a number id: id-number first column < id-number second column
  comparison.numeric <- matrix(combn(1:nt, 2), nrow = length(combn(1:nt, 2))/2, ncol = 2, byrow = T)



  ## A dataframe with both the original name and the id-number where id-number first column > id-number second column
  comp <- data.frame(1:length(comparison[, 1]), comparison[, 2:1], comparison.numeric[, 2:1])
  colnames(comp) <- c("comparison", "name t1", "name t2", "id t1", "id t2")


  return(comp)

}

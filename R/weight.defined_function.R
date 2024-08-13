#' Preparing the study weights based on Gower's dissimilarities 
#' (Transitivity evaluation)
#' 
#' @description
#' Using the Gower's dissimilarities to prepare the study weights in a network 
#' of interventions to be used by the \code{\link{run_model}} function in the 
#' context of the transitivity evaluation as a sensitivity analysis. 
#' 
#' @param diss_res A list of two elements with the following order:
#'   1) an object of S3 class 'comp_clustering' (\code{\link{comp_clustering}}),
#'   and 2) a character with values \code{"uniform"}, or \code{"rms"} to define 
#'   the usage of the similarities (see 'Details') as being sampled from a 
#'   uniform distribution, or being fixed to their root mean square. 
#'   See 'Details'.
#'   
#' @return A list of the following two elements:
#'   \item{weights}{A vector of study weights if \code{"rms"} has been specified
#'   or a matrix with the minimum and maximum similarity values for each study
#'   if \code{"uniform"} has been specified.}
#'   \item{type}{A character indicating the weight type considered: 
#'   \code{"rms"} or \code{"uniform"}.}
#'   
#' @details
#' The function receives the matrix \code{Trials_diss_table} found in the 
#' results of \code{\link{comp_clustering}}. This matrix contains the Gower's
#' dissimilarities of all study pairs in the network for a specific set of 
#' clinical and methodological characteristics that act as effect modifiers.
#' Gower's dissimilarities take values from 0 to 1 (or 0 to 100%), with 0 and 
#' 1 (or 100%) implying perfect similarity and perfect dissimilarity, 
#' respectively. Hence, subtracting each value from 1 yields the similarities of 
#' all study pairs in the network. 
#' 
#' For each study, the Gower's similarities can be either transformed into the
#' root mean square or restricted to the minimum and maximum values. This 
#' depends on the second element specified in the argument \code{diss_res}. 
#' When the \code{"uniform"} choice has been specified in the second element 
#' of the argument \code{diss_res}, the function checks whether the bounds are
#' the same value. If this is the case for at least one study, the uniform 
#' distribution cannot be defined and, hence, the function returns the root mean 
#' square value (\code{"rms"}) for all studies.
#' 
#' \code{weight_defined} can be used in \code{\link{run_model}} via the argument
#' \code{adjust_wgt}. See 'Details' in \code{\link{run_model}}.
#' 
#' @author {Loukia M. Spineli}
#' 
#' @seealso \code{\link{comp_clustering}}, \code{\link{run_model}}
#' 
#' @references
#' Gower J. General Coefficient of Similarity and Some of Its Properties.
#' \emph{Biometrics} 1971;\bold{27}(4):857--71.
#' doi: 10.2307/2528823
#' 
#' @export
weight_defined <- function (diss_res) {
  
  ## Check default
  if (inherits(diss_res[[1]], "comp_clustering") == FALSE) {
    stop("The first argument must be an object of S3 class 'comp_clustering'", 
         call. = FALSE)
  }
  transf <- if (!is.element(diss_res[[2]], c("rms", "uniform"))) {
    stop("Insert one of the following: 'rms', or 'uniform'.", 
         call. = FALSE)
  } else {
    diss_res[[2]]
  }
  
  ## Consider the matrix of study dissimilarities
  diss <- diss_res[[1]]$Trials_diss_table
  
  ## Copy-paste the lower off-diagonal elements to the corresponding upper
  diss[upper.tri(diss)] <- t(diss)[upper.tri(diss)]
  diss
  
  ## Turn into weights (1 - diss)
  weight <- 1 - diss
  
  ## Turn diagonal into NA
  diag(weight) <- NA
  
  ## Distinguish between two-arm and multi-arm studies
  # Get the unique study ID
  index <- sub("\\(.*", "", gsub('.{3}$', '', rownames(weight)))
  
  # Split dataset by 'index'
  split_multi_arms <- split(weight, factor(index, levels = unique(index)))
  
  ## Transform the weights: sampled from a selected distribution or fixed to rms
  tranform_result <-
    if (transf == "uniform") {    # Uniform distribution
    # Lower bound
    lower <- unlist(lapply(split_multi_arms, function(x) min(x, na.rm = TRUE)))
    
    # upper bound
    upper <- unlist(lapply(split_multi_arms, function(x) max(x, na.rm = TRUE)))
    
    data.frame(lower, upper)
  } else if (transf == "rms") { # Average
    sqrt(unlist(lapply(split_multi_arms, function(x) {sum(na.omit(x)^2) / 
        (length(x) - sum(is.na(x)))})))
  }
  
  
  ## Check whether the uniform distribution is undefined
  with_problem <- if (transf == "uniform") {
    any(apply(tranform_result, 1, diff) == 0)
  } else if (transf == "rms") {
    FALSE
  }  

  
  ## When the uniform distribution is undefined, switch to 'rms'
  proper_tranform <- if (all(transf == "uniform" & with_problem == TRUE) == TRUE) {
    message("Undefined uniform distribution. Root mean square was used instead.")
    sqrt(unlist(lapply(split_multi_arms, function(x) {sum(na.omit(x)^2) / 
        (length(x) - sum(is.na(x)))})))
  } else {
    tranform_result
  }

  
  ## Consider the 'proper' transformation: 'uniform', or 'rms'
  transf_new <- if (with_problem == TRUE) {
    "rms"
  } else {
    transf
  }

  return(list(weights = proper_tranform,
              type = transf_new))
}

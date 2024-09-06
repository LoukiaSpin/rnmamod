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
#'   and 2) a character with values \code{"uniform"}, \code{"fixed"}, or
#'   \code{"index"} to define the usage of the similarities (see 'Details') as
#'   being sampled from a uniform distribution, between-comparisons similarity
#'   (fixed weight) or a ratio of their between-comparisons similarity and total
#'   similarity (percentage index). See 'Details'.
#'
#' @return A list of the following two elements:
#'   \item{weights}{A vector of study weights if \code{"fixed"} or
#'   \code{"index"} has been specified or a matrix with the minimum and maximum
#'   similarity values for each study if \code{"uniform"} has been specified.}
#'   \item{type}{A character indicating the weight type considered:
#'   \code{"uniform"}, \code{"fixed"}, or \code{"index"}.}
#'
#' @details
#' The function receives the matrix \code{Trials_diss_table} found in the
#' results of \code{\link{comp_clustering}}. This matrix contains the Gower's
#' dissimilarities of all study pairs in the network for a specific set of
#' clinical and methodological characteristics that act as effect modifiers.
#' Gower's dissimilarities take values from 0 to 1, with 0 and 1 implying
#' perfect similarity and perfect dissimilarity, respectively. Hence,
#' subtracting each value from 1 yields the similarities of all study pairs in
#' the network.
#'
#' For each study, the Gower's similarities can be either transformed into a
#' fixed value (\code{"fixed"}) or percentage index (\code{"index"}) based on
#' the root mean square of the similarities or restricted to the minimum and
#' maximum similarity values (\code{"uniform"}). This depends on the second
#' element specified in the argument \code{diss_res}.
#' When the \code{"uniform"} choice has been specified in the second element
#' of the argument \code{diss_res}, the function checks whether the bounds are
#' the same value. If this is the case for at least one study, the uniform
#' distribution cannot be defined and, hence, the function returns the
#' between-comparisons similarity values (\code{"fixed"}) for each study.
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
  transf <- if (!is.element(diss_res[[2]], c("fixed", "index", "uniform"))) {
    stop("Insert one of the following: 'fixed', 'index', or 'uniform'.",
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
  index <- gsub("^\\s+|\\s+$", "",
                sub("\\(.*", "", gsub('.{3}$', '', rownames(weight))))

  # Get the comparison for each study
  comp_index <- gsub("^\\s+|\\s+$", "",
                     substring(rownames(weight), nchar(rownames(weight)) - 3))

  # Split 'weight' further by 'rownames(weight)'
  split_comp <- split(weight, factor(rownames(weight), levels = unique(rownames(weight))))

  # Split 'split_comp' further by 'comp_index'
  split_study_comp <- lapply(split_comp, function(x) split(x, comp_index))

  # Set of within-comparison similarities per study
  within_set <-
    lapply(1:length(split_study_comp),
           function(x)
             unlist(split_study_comp[[x]][is.element(names(split_study_comp[[x]]),
                                                     comp_index[x])]))

  # Within-comparison similarity per study
  within_comp <-
    sqrt(unlist(lapply(within_set, function(x) {sum(na.omit(x)^2) /
        (length(x) - sum(is.na(x)))})))

  # Set of between-comparisons similarities per study
  between_set <-
    lapply(1:length(split_study_comp),
           function(x)
             unlist(split_study_comp[[x]][!is.element(names(split_study_comp[[x]]),
                                                      comp_index[x])]))

  # Between-comparisons similarity per study
  between_comp <-
    sqrt(unlist(lapply(between_set, function(x) {sum(x^2) / length(x)})))


  ## Calculate the weights
  tranform_result <-
    if (transf == "uniform") {    # Uniform distribution
    # Lower bound (Return the minimum for multi-arm studies)
    lower <-
      unlist(lapply(split(unlist(lapply(between_set, function(x) min(x))),
                          factor(index, levels = index)), function(x) min(x)))

    # upper bound (Return the minimum for multi-arm studies)
    upper <-
      unlist(lapply(split(unlist(lapply(between_set, function(x) max(x))),
                          factor(index, levels = index)), function(x) min(x)))

    data.frame(lower, upper)
  } else if (transf == "index") { # Return the minimum for multi-arm studies
    unlist(lapply(split(between_comp / (within_comp + between_comp),
                        factor(index, levels = index)), function(x) min(x)))
  } else if (transf == "fixed") { # Return the minimum for multi-arm studies
    unlist(lapply(split(between_comp, factor(index, levels = index)),
                  function(x) min(x)))
  }


  ## Check whether the uniform distribution is undefined
  with_problem <- if (transf == "uniform") {
    any(apply(tranform_result, 1, diff) == 0)
  } else if (is.element(transf, c("fixed", "index"))) {
    FALSE
  }


  ## When the uniform distribution is undefined, switch to 'fixed'
  proper_tranform <-
    if (all(transf == "uniform" & with_problem == TRUE) == TRUE) {
      message("Undefined uniform distribution. Fixed weights were used instead.")
      unlist(lapply(split(between_comp, factor(index, levels = index)),
                    function(x) min(x)))
    } else {
      tranform_result
    }


  ## Consider the 'proper' transformation: 'uniform', 'fixed', or 'index'
  transf_new <- if (with_problem == TRUE) {
    "fixed"
  } else {
    transf
  }

  return(list(weights = proper_tranform,
              type = transf_new))
}

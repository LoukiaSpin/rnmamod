#' Determine the prior distribution for the heterogeneity parameter
#'
#' @description
#'   Generates the prior distribution (weakly informative or empirically-based)
#'   for the heterogeneity parameter.
#'   \code{\link{run_model}} inherits \code{heterogeneity_param_prior} via the
#'   argument \code{heter_prior}.
#'
#' @param measure Character string indicating the effect measure with values
#'   \code{"OR"}, \code{"MD"}, \code{"SMD"}, or \code{"ROM"} for the odds ratio,
#'   mean difference, standardised mean difference and ratio of means,
#'   respectively.
#' @param model Character string indicating the analysis model with values
#'   \code{"RE"}, or \code{"FE"} for the random-effects and fixed-effect model,
#'   respectively. The default argument is \code{"RE"}.
#' @param heter_prior A list of three elements with the following order:
#'   1) a character string indicating the distribution with
#'   (currently available) values \code{"halfnormal"}, \code{"uniform"},
#'   \code{"lognormal"}, or \code{"logt"}; 2) two numeric values that refer to
#'   the parameters of the selected distribution.  For \code{"lognormal"}, and
#'   \code{"logt"} these numbers refer to the mean and precision, respectively.
#'   For \code{"halfnorm"}, these numbers refer to zero and the scale parameter
#'   (equal to 4 or 1 being the corresponding precision of the scale parameter
#'   0.5 or 1). For \code{"uniform"}, these numbers refer to the
#'   minimum and maximum value of the distribution.
#'
#' @return A value to be passed to \code{\link{run_model}}.
#'
#' @details
#'   The names of the (current) prior distributions follow the JAGS syntax.
#'   The mean and precision of \code{"lognormal"} and \code{"logt"} should align
#'   with the values proposed by Turner et al., (2015) and Rhodes et al., (2015)
#'   for the corresponding empirically-based prior distributions when
#'   \code{measure} is \code{"OR"} or \code{"SMD"}, respectively.
#'   The users may refer to Dias et al., (2013) to determine the minimum and
#'   maximum value of the uniform distribution, and to Friede et al., (2017)
#'   to determine the mean and precision of the half-normal distribution.
#'   When \code{model} is \code{"FE"}, \code{heterogeneity_param_prior}
#'   is ignored in \code{\link{run_model}}.
#'
#' @author {Loukia M. Spineli}
#'
#' @seealso \code{\link{run_model}}
#'
#' @references
#' Friede T, Roever C, Wandel S, Neuenschwander B. Meta-analysis of two studies
#' in the presence of heterogeneity with applications in rare diseases.
#' \emph{Biom J} 2017;\bold{59}(4):658--671.
#' \doi{10.1002/bimj.201500236}
#'
#' Turner RM, Jackson D, Wei Y, Thompson SG, Higgins JP. Predictive
#' distributions for between-study heterogeneity and simple methods for their
#' application in Bayesian meta-analysis.
#' \emph{Stat Med} 2015;\bold{34}(6):984--98. \doi{10.1002/sim.6381}
#'
#' Rhodes KM, Turner RM, Higgins JP. Predictive distributions were developed
#' for the extent of heterogeneity in meta-analyses of continuous outcome data.
#' \emph{J Clin Epidemiol} 2015;\bold{68}(1):52--60.
#' \doi{10.1016/j.jclinepi.2014.08.012}
#'
#' Dias S, Sutton AJ, Ades AE, Welton NJ. Evidence synthesis for decision
#' making 2: a generalized linear modeling framework for pairwise and
#' network meta-analysis of randomized controlled trials.
#' \emph{Med Decis Making} 2013;\bold{33}(5):607--617.
#' \doi{10.1177/0272989X12458724}
#'
#' @export
heterogeneity_param_prior <- function(measure, model, heter_prior) {

  # Specifying the prior distribution for the between-trial parameter
  if (model == "RE" &
      missing(heter_prior)) {
    stop("The argument 'heter_prior' needs to be defined", call. = FALSE)
  } else if (model == "FE" &
             missing(heter_prior)) {
    list(NA, NA, NA)
  } else if (model == "FE") {
    #message("The argument 'heter_prior' has been ignored.")
    list(NA, NA, NA)
  } else if (model == "RE" &
             measure == "OR" &
             heter_prior[[1]] != "halfnormal" &
             heter_prior[[1]] != "uniform" &
             heter_prior[[1]] != "lognormal") {
    stop("Insert 'halfnormal', 'uniform', or 'lognormal'", call. = FALSE)
  } else if (model == "RE" &
             measure == "SMD" &
             heter_prior[[1]] != "halfnormal" &
             heter_prior[[1]] != "uniform" &
             heter_prior[[1]] != "logt") {
    stop("Insert 'halfnormal', 'uniform', or 'logt'", call. = FALSE)
  } else if (model == "RE" &
             (measure == "MD" || measure == "ROM") &
             heter_prior[[1]] != "halfnormal" &
             heter_prior[[1]] != "uniform") {
    stop("Insert 'halfnormal', or 'uniform'", call. = FALSE)
  } else if (model == "RE" &
             heter_prior[[1]] == "halfnormal") {
    as.numeric(c(0, heter_prior[[3]], 1))
  } else if (model == "RE" &
             heter_prior[[1]] == "uniform") {
    as.numeric(c(0, heter_prior[[3]], 2))
  } else if (model == "RE" &
             heter_prior[[1]] == "lognormal") {
    as.numeric(c(heter_prior[[2]], heter_prior[[3]], 3))
  } else if (model == "RE" &
             heter_prior[[1]] == "logt") {
    as.numeric(c(heter_prior[[2]], heter_prior[[3]], 4))
  } else if (model == "RE" &
             measure == "OR" &
             heter_prior[[1]] == "lognormal")  {
    as.numeric(c(heter_prior[[2]], heter_prior[[3]], 3))
  } else if (model == "RE" &
             measure != "OR" &
             heter_prior[[1]] == "lognormal") {
    stop("Not the proper prior distribution for continuous outcome",
         call. = FALSE)
  } else if (model == "RE" &
             measure == "SMD" &
             heter_prior[[1]] == "logt") {
    as.numeric(c(heter_prior[[2]], heter_prior[[3]], 3))
  } else if (model == "RE" &
             (measure == "MD" || measure == "ROM") &
             heter_prior[[1]] == "logt") {
    aa <- "Currently, no empirically-based prior distributions for MD and ROM."
    bb <- "Choose a half-normal or a uniform prior distribution, instead"
    stop(paste(aa, bb), call. = FALSE)
  } else if (model == "RE" &
             measure == "OR" &
             heter_prior[[1]] == "logt") {
    stop("This is not the proper prior distribution for binary outcome",
         call. = FALSE)
  }
}

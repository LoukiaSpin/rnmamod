#' Function for the hyper-parameters of the prior distribution of
#' the inconsistency variance (network meta-analysis with random inconsistency
#' effects)
#'
#' @description Calculates the mean and standard deviation of the log-normal
#'   distribution and location-scale t-distribution of the inconsistency
#'   variance in the log-odds ratio and standardised mean difference scales,
#'   respectively, based on corresponding empirical distributions for the
#'   between-study variance proposed by Turner et al. (2015) and
#'   Rhodes et al. (2015). It also return the median value of the inconsistency
#'   standard deviation.
#'
#' @param mean_tau2 Mean value from the empirical prior distribution for the
#'   between-study variance.
#' @param sd_tau2 Standard deviation value from the empirical prior distribution
#'   for the between-study variance.
#' @param mean_scale Positive (non-zero value) as a scaling factor of
#'   \code{mean_tau2}. See Law et al. (2016).
#' @param measure Character string indicating the effect measure. For a binary
#'   outcome, use only \code{"OR"} for the odds ratio. For a continuous outcome,
#'   use only \code{"SMD"} for standardised mean difference.
#'
#' @return A list of three elements: the mean and standard deviation for the
#'   prior distribution for the inconsistency variance, and the median
#'   inconsistency standard deviation according to the selected empirical prior
#'   distribution for the between-study variance.
#'
#' @details Law et al. (2016) suggested using the proposed empirical prior
#'   distributions for between-study variance to construct a prior distribution
#'   for the inconsistency variance. The authors provided the formulas for the
#'   hyper-parameters of the inconsistency variance for a binary outcome
#'   measured in the log odds ratio scale. We extended the idea for a continuous
#'   outcome measured in the standardised mean difference scale. Currently,
#'   the empirical prior distributions for the between-study variance have been
#'   proposed for these effect measures only (Turner et al. (2015),
#'   Rhodes et al. (2015)).
#'
#' @author {Loukia M. Spineli}
#'
#' @references
#' Law M, Jackson D, Turner R, Rhodes K, Viechtbauer W. Two new methods to fit
#' models for network meta-analysis with random inconsistency effects.
#' \emph{BMC Med Res Methodol} 2016;\bold{16}:87. doi: 10.1186/s12874-016-0184-5
#'
#' Rhodes KM, Turner RM, Higgins JP. Predictive distributions were developed for
#' the extent of heterogeneity in meta-analyses of continuous outcome data.
#' \emph{J Clin Epidemiol} 2015;\bold{68}(1):52--60.
#' doi: 10.1016/j.jclinepi.2014.08.012
#'
#' Turner RM, Jackson D, Wei Y, Thompson SG, Higgins JP. Predictive
#' distributions for between-study heterogeneity and simple methods for their
#' application in Bayesian meta-analysis.
#' \emph{Stat Med} 2015;\bold{34}(6):984--98. doi: 10.1002/sim.6381
#'
#' @export
inconsistency_variance_prior <- function (mean_tau2,
                                          sd_tau2,
                                          mean_scale,
                                          measure) {

  ## Default options
  mean_tau2 <- if (missing(mean_tau2)) {
    stop("The argument 'mean_tau2' needs to be defined.", call. = FALSE)
  } else  {
    mean_tau2
  }
  sd_tau2 <- if (missing(sd_tau2)) {
    stop("The argument 'sd_tau2' needs to be defined.", call. = FALSE)
  } else if (sd_tau2 <= 0) {
    stop("The argument 'sd_tau2' must be above 0.", call. = FALSE)
  } else {
    sd_tau2
  }
  mean_scale <- if (missing(mean_scale)) {
    stop("The argument 'mean_scale' needs to be defined.", call. = FALSE)
  } else if (mean_scale <= 0) {
    stop("The argument 'mean_scale' must be above 0.", call. = FALSE)
  } else {
    mean_scale
  }
  measure <- if (missing(measure)) {
    stop("Insert 'OR', or 'SMD'.", call. = FALSE)
  } else if (!is.element(measure, c("OR", "SMD"))) {
    stop("Insert 'OR', or 'SMD'.", call. = FALSE)
  } else {
    measure
  }


  ## Calculate depending on the effect measure (OR or SMD)
  if (measure == "OR") {
    # Mean of log-normal distribution
    M <- mean_scale * exp(mean_tau2 + ((sd_tau2 * sd_tau2) / 2))

    # Variance of log-normal distribution
    V <- (exp(sd_tau2 * sd_tau2) - 1) * exp((2 * mean_tau2) + (sd_tau2 * sd_tau2))

    # Mean of tau2.omega prior (log-normal) distribution
    mean_tau2_omega <- log(M) - 0.5 * log((V / (M * M)) + 1)

    # Standard deviation of tau2.omega prior (log-normal) distribution
    sd_tau2_omega <- sqrt(log((V / (M * M)) + 1))

    # Median tau.omega (in *standard deviation* scale) based on its prior (log-normal) distribution
    median_tau_omega <- sqrt(qlnorm(0.5, mean_tau2_omega, sd_tau2_omega))

  } else {

    # Mean of log tau2.omega prior (t) distribution
    mean_tau2_omega <- log(mean_scale) +  mean_tau2

    # Standard deviation of log tau2.omega prior (t) distribution
    sd_tau2_omega <- sqrt((sd_tau2 * sd_tau2) * (5 / (5 - 2)))

    # Median tau.omega (in *standard deviation* scale) based on its prior (t) distribution
    median_tau_omega <- sqrt(exp((qt(0.5, df = 5) * sd_tau2_omega) + mean_tau2_omega))
  }


  ## Obtain the results
  results <- list(mean_tau2_omega = mean_tau2_omega,
                  sd_tau2_omega = sd_tau2_omega,
                  median_tau_omega = median_tau_omega)

  return(results)
}

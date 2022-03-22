#' Pattern-mixture model with Taylor series for continuous outcome
#'
#' @description Applies the pattern-mixture model under a specific assumption
#'   about the informative missingness parameter in trial-arms with
#'   \bold{continuous} missing participant outcome data and uses the Taylor
#'   series to obtain the effect size and standard error for each trial
#'   (Mavridis et al., 2015).
#'
#' @param data A data-frame in the long arm-based format. Two-arm trials occupy
#'   one row in the data-frame. Multi-arm trials occupy as many rows as the
#'   number of possible comparisons among the interventions. See 'Format' for
#'   the specification of the columns.
#' @param measure Character string indicating the effect measure with values
#'   \code{"MD"}, \code{"SMD"}, or \code{"ROM"} for the mean difference,
#'   standardised mean difference, and ratio of means, respectively.
#' @param mean_value A numeric value for the mean of the normal distribution of
#'   the informative missingness parameter. The same value is considered for all
#'   trial-arms of the dataset. The default argument is 0 and corresponds to the
#'   missing-at-random assumption. For the informative missingness ratio of
#'   means, the mean value is defined in the logarithmic scale.
#' @param var_value A positive non-zero number for the variance of the normal
#'   distribution of the informative missingness parameter. When the
#'   \code{measure} is \code{"MD"}, or \code{"SMD"} the default argument is 1;
#'   when the \code{measure} is \code{"ROM"} the default argument is 0.04. The
#'   same value is considered for all trial-arms of the dataset.
#' @param rho A numeric value in the interval [-1, 1] that indicates the
#'   correlation coefficient between two informative missingness parameters in
#'   a trial. The same value is considered across all trials of the dataset.
#'   The default argument is 0 and corresponds to uncorrelated missingness
#'   parameters.
#'
#' @format The columns of the data-frame in the argument \code{data} refer to
#'   the following ordered elements for a continuous outcome:
#'   \tabular{rl}{
#'   \bold{id} \tab A unique identifier for each trial.\cr
#'   \bold{y1} \tab The observed mean outcome in the first arm of the
#'   comparison.\cr
#'   \bold{y2} \tab The observed mean outcome in the second arm of the
#'   comparison.\cr
#'   \bold{sd1} \tab The observed standard deviation of the outcome in the
#'   first arm of the comparison.\cr
#'   \bold{sd2} \tab The observed standard deviation of the outcome in the
#'   second arm of the comparison.\cr
#'   \bold{m1} \tab The number of missing participants in the first arm of
#'   the comparison.\cr
#'   \bold{m2} \tab The number of missing participants in the second arm of
#'   the comparison.\cr
#'   \bold{n1} \tab The number randomised in the first arm of the
#'   comparison.\cr
#'   \bold{n2} \tab The number randomised in the second arm of the
#'   comparison.\cr
#'   \bold{t1} \tab An identifier for the intervention in the first arm of
#'   the comparison.\cr
#'   \bold{t2} \tab An identifier for the intervention in the second arm of
#'   the comparison.\cr
#'   }
#'
#' @return A data-frame that additionally includes the following elements:
#'   \item{EM}{The effect size adjusted for the missing participants and
#'   obtained using the Taylor series.}
#'   \item{se.EM}{The standard error of the effect size adjusted for the missing
#'   participants and obtained using the Taylor series.}
#'
#' @details The \code{taylor_continuous} function is integrated in the
#'   \code{\link{unrelated_effects_plot}} function. The latter uses the
#'   the \code{\link[netmeta:pairwise]{pairwise}} function from the package
#'   \href{https://CRAN.R-project.org/package=netmeta}{netmeta}
#'   to transform the dataset from the wide arm-based format into the long
#'   arm-based format (see, 'Arguments' for \code{data} in
#'   \code{\link{unrelated_effects_plot}}).
#'
#' @seealso \href{https://CRAN.R-project.org/package=netmeta}{pairwise},
#'   \code{\link{run_model}}, \code{\link{unrelated_effects_plot}}
#'
#' @references
#' Mavridis D, White IR, Higgins JP, Cipriani A, Salanti G. Allowing for
#' uncertainty due to missing continuous outcome data in pairwise and network
#' meta-analysis. \emph{Stat Med} 2015;\bold{34}(5):721--41.
#' \doi{10.1002/sim.6365}
#'
#' @author {Loukia M. Spineli}
#'
#' @export
taylor_continuous <- function(data, measure, mean_value, var_value, rho) {

  # Calculate the probability of observing the outcomes
  # Control arm
  a1 <- (data[, 8] - data[, 6]) / data[, 8]
  # Experimental arm
  a2 <- (data[, 9] - data[, 7]) / data[, 9]

  # Calculate the adjusted-mean for MOD in the randomised sample
  if (measure == "MD" || measure == "SMD") {
    # Control arm
    y_all1 <- data[, 2] + mean_value * (1 - a1)
    # Experimental arm
    y_all2 <- data[, 3] + mean_value * (1 - a2)
  } else if (measure == "ROM") {
    # Control arm
    y_all1 <- data[, 2] * (a1 + exp(mean_value) * (1 - a1))
    # Experimental arm
    y_all2 <- data[, 3] * (a2 + exp(mean_value) * (1 - a2))
  }

  # Estimate the adjusted within-trial effect measures
  if (measure == "MD") {
    # Experimental vs Control
    md <- y_all2 - y_all1
  } else if (measure == "SMD") {
    # Control arm
    nominator1 <- (data[, 8] - data[, 6] - 1) * data[, 4] * data[, 4]
    # Experimental arm
    nominator2 <- (data[, 9] - data[, 7] - 1) * data[, 5] * data[, 5]
    denominator <- (data[, 8] - data[, 6] - 1) + (data[, 9] - data[, 7] - 1)
    sd_pooled <- sqrt((nominator1 + nominator2) / denominator)
    # Experimental vs Control
    smd <- (y_all2 - y_all1) / sd_pooled
  } else if (measure == "ROM") {
    # Experimental vs Control
    lrom <- log(y_all2 / y_all1)
  }

  #####################################################################
  ## USE OF TAYLOR APPROXIMATION FOR THE TOTAL WITHIN-TRIAL VARIANCE ##
  #####################################################################

  # Estimating the variance of the mean effect size based on the observed data
  # Derivative of y_all by y.obs per arm
  if (measure == "MD" || measure == "SMD") {
    alpha1 <- alpha2 <- 1
  } else if (measure == "SMD") {
    alpha1 <- alpha2 <- 1
  } else if (measure == "ROM") {
    alpha1 <- alpha1 + exp(mean_value) * (1 - a1)
    alpha2 <- alpha2 + exp(mean_value) * (1 - a2)
  }

  # Variance of y.obs per arm
  beta1 <- (data[, 4] * data[, 4]) / (data[, 8] - data[, 6])
  beta2 <- (data[, 5] * data[, 5]) / (data[, 9] - data[, 7])

  # Derivative of y_all by prob of MOD (i.e. a) per arm
  if (measure == "MD" || measure == "SMD") {
    c1 <- c2 <- -mean_value
  } else if (measure == "ROM") {
    c1 <- data[, 2] * (1 - exp(mean_value))
    c2 <- data[, 3] * (1 - exp(mean_value))
  }

  # Variance of prob of MOD
  d1 <- (a1 * (1 - a1)) / data[, 8]
  d2 <- (a2 * (1 - a2)) / data[, 9]

  # Derivative of link function for MD, SMD and LROM per arm
  if (measure == "MD") {
    e1 <- e2 <- 1
  } else if (measure == "SMD") {
    e1 <- e2 <- 1 / sd_pooled
  } else if (measure == "ROM") {
    e1 <- 1 / y_all1
    e2 <- 1 / y_all2
  }

  # Variance using the observed cases (Equation 13 in PMID: 17703496)
  v_obs <- (((alpha1 * alpha1 * beta1) + (c1 * c1 * d1)) * e1 * e1) +
    (((alpha2 * alpha2 * beta2) + (c2 * c2 * d2)) * e2 * e2)

  # Estimating the variance of the mean effect size arising from the
  # informative missingness parameter
  # Derivative of y_all by delta per arm
  if (measure == "MD" || measure == "SMD") {
    h1 <- (1 - a1)
    h2 <- (1 - a2)
  } else if (measure == "ROM") {
    h1 <- data[, 2] * exp(mean_value) * (1 - a1)
    h2 <- data[, 3] * exp(mean_value) * (1 - a2)
  }

  # Variance due to informative missingness
  v_delta <- (h1 * h1 * var_value * e1 * e1) + (h2 * h2 * var_value * e2 * e2) -
    2 * rho * h1 * h2 * e1 * e2 * sqrt(var_value) * sqrt(var_value)

  # Variance using the randomised sample
  v_all <- v_obs + v_delta

  # Include trial-specific adjusted MDs, SDMs and SEs in the initial dataset
  if (measure == "MD") {
    final <- data.frame(cbind(data, round(md, 3),  round(sqrt(v_all), 3)))
  } else if (measure == "SMD") {
    final <- data.frame(cbind(data, round(smd, 3),  round(sqrt(v_all), 3)))
  } else if (measure == "ROM") {
    final <- data.frame(cbind(data, round(lrom, 3),  round(sqrt(v_all), 3)))
  }
  colnames(final) <- c("id",
                       "mean1", "mean2",
                       "sd1", "sd2",
                       "m1", "m2",
                       "n1", "n2",
                       "t1", "t2",
                       "EM", "se.EM")

  return(final)
}

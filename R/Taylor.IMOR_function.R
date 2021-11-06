#' Pattern-mixture model with Taylor series for binary outcome
#'
#' @description Applies the pattern-mixture model under a specific assumption
#'   about the informative missingness parameter in trial-arms with
#'   \bold{binary} missing participant outcome data (MOD) and uses the Taylor
#'   series to obtain the effect size and standard error for each trial
#'   (White et al., 2008).
#'
#' @param data A data-frame in the long arm-based format. Two-arm trials occupy
#'   one row in the data-frame. Multi-arm trials occupy as many rows as the
#'   number of possible comparisons among the interventions. See 'Format' for
#'   the specification of the columns.
#' @param mean_value A numeric value for the mean of the normal distribution of
#'   the informative missingness odds ratio in the logarithmic scale. The same
#'   value is considered for all trial-arms of the dataset. The default argument
#'   is 0 and corresponds to the missing-at-random assumption.
#' @param var_value A positive non-zero number for the variance of the normal
#'   distribution of the informative missingness odds ratio in the logarithmic
#'   scale. The default argument is 1.
#' @param rho A numeric value in the interval [-1, 1] that indicates the
#'   correlation coefficient between two missingness parameters in a trial. The
#'   same value is considered across all trials of the dataset. The default
#'   argument is 0 and corresponds to uncorrelated missingness parameters.
#'
#' @format The columns of the data-frame in the argument \code{data} refer to
#'   the following ordered elements for a binary outcome:
#'   \tabular{ll}{
#'    \strong{id} \tab A unique identifier for each trial.\cr
#'    \tab \cr
#'    \strong{r1} \tab The observed  number of events in the first arm of the
#'    comparison.\cr
#'    \tab \cr
#'    \strong{r2} \tab The observed  number of events in the second arm of the
#'    comparison.\cr
#'    \tab \cr
#'    \strong{m1} \tab The number of MOD in the first arm of the comparison.\cr
#'    \tab \cr
#'    \strong{m2} \tab The number of MOD in the second arm of the comparison.\cr
#'    \tab \cr
#'    \strong{n1} \tab The number of participants randomised in the first arm of
#'    the comparison.\cr
#'    \tab \cr
#'    \strong{n2} \tab The number of participants randomised in the second arm
#'    of the comparison.\cr
#'    \tab \cr
#'    \strong{t1} \tab An identified for the intervention in the first arm of
#'    the comparison.\cr
#'    \tab \cr
#'    \strong{t2} \tab An identified for the intervention in the second arm of
#'    the comparison.\cr
#'   }
#'
#' @return A data-frame that additionally includes the following elements:
#'   \tabular{ll}{
#'    \strong{EM} \tab The odds ratio in the logarithmic scale (log OR) adjusted
#'    for MOD and obtained using the Taylor series.\cr
#'    \tab \cr
#'    \strong{se.EM} \tab The standard error of the log OR adjusted for MOD and
#'    obtained using the Taylor series.\cr
#'   }
#'
#' @details The \code{taylor_imor} function is found in the
#'   \code{\link{unrelated_effects_plot}} function. The latter uses the
#'   the \code{\link[netmeta]{pairwise}} function from the package
#'   \href{https://CRAN.R-project.org/package=netmeta}{netmeta}
#'   to transform the dataset from the wide arm-based format
#'   (see, 'Arguments' for \code{data} in
#'   \code{\link{unrelated_effects_plot}}) into the long-arm based
#'   format.
#'
#' @seealso \code{\link[netmeta]{pairwise}}, \code{\link{run_model}},
#'   \code{\link{unrelated_effects_plot}}
#'
#'
#' @references
#' White IR, Higgins JP, Wood AM. Allowing for uncertainty due to missing data
#' in meta-analysis--part 1: two-stage methods.
#' \emph{Stat Med} 2008;\bold{27}(5):711--27. \doi{10.1002/sim.3008}
#'
#' @author {Loukia M. Spineli}
#'
#' @export
taylor_imor <- function(data, mean_value, var_value, rho) {

  for (i in seq_len(length(data[, 1]))) {
    # Add 0.5 continuity correction when there is at least on zero cell
    if (data[i, 2] == 0 || data[i, 3] == 0 ||
       data[i, 6] - data[i, 4] - data[i, 2] == 0 ||
       data[i, 7] - data[i, 5] - data[i, 3] == 0) {
      data[i, 2] <- data[i, 2] + 0.5
      data[i, 3] <- data[i, 3] + 0.5
      data[i, 6] <- data[i, 6] + 1
      data[i, 7] <- data[i, 7] + 1
    } else {
      data[i, 2] <- data[i, 2]
      data[i, 3] <- data[i, 3]
      data[i, 6] <- data[i, 6]
      data[i, 7] <- data[i, 7]
    }
  }

  # Calculate the probability of event among completers
  # Control arm
  p_o1 <- data[, 2] / (data[, 6] - data[, 4])
  # Experimental arm
  p_o2 <- data[, 3] / (data[, 7] - data[, 5])

  # Calculate the probability of missing outcome data in arm
  # Control arm
  a1 <- data[, 4] / data[, 6]
  # Experimental arm
  a2 <- data[, 5] / data[, 7]

  # Calculate the probability of event in randomised sample
  # Control arm
  p_all1 <- (1 - a1) * p_o1 + a1 * ((exp(mean_value) * p_o1) /
                                      (exp(mean_value) * p_o1 + 1 - p_o1))
  # Experimental arm
  p_all2 <- (1 - a2) * p_o2 + a2 * ((exp(mean_value) * p_o2) /
                                      (exp(mean_value) * p_o2 + 1 - p_o2))

  # Estimates the odds ratio in the logarithmic scale
  logor <- log(p_all2 / (1 - p_all2)) - log(p_all1 / (1 - p_all1))

  #################################
  ## USE OF TAYLOR APPROXIMATION ##
  #################################

  # Derivative of p_all by p_o per arm
  # (first term in Equation 14 in PMID: 17703496)
  alpha1 <- 1 - a1 + (a1 * exp(mean_value)) /
    (exp(mean_value) * p_o1 + 1 - p_o1)^2
  alpha2 <- 1 - a2 + (a2 * exp(mean_value)) /
    (exp(mean_value) * p_o2 + 1 - p_o2)^2

  # Variance of p_o per arm (second term in Equation 14 in PMID: 17703496)
  beta1 <- (p_o1 * (1 - p_o1)) / (data[, 6] - data[, 4])
  beta2 <- (p_o2 * (1 - p_o2)) / (data[, 7] - data[, 5])

  # Derivative of p_all by prob of MOD (i.e. a) per arm
  # (third term in Equation 14 in PMID: 17703496)
  c1 <- (p_o1 * (1 - p_o1) * (exp(mean_value) - 1)) /
    (exp(mean_value) * p_o1 + 1 - p_o1)
  c2 <- (p_o2 * (1 - p_o2) * (exp(mean_value) - 1)) /
    (exp(mean_value) * p_o2 + 1 - p_o2)

  # Variance of prob of MOD (i.e. a) per arm
  # (fourth term in Equation 14 in PMID: 17703496)
  d1 <- (a1 * (1 - a1)) / data[, 6]
  d2 <- (a2 * (1 - a2)) / data[, 7]

  # Variance of log odds using delta-method per arm
  e1 <- 1 / (p_all1 * (1 - p_all1))
  e2 <- 1 / (p_all2 * (1 - p_all2))

  # Derivative of p_all by delta per arm
  # (second Equation after Equation (15) in PMID: 17703496)
  h1 <- (a1 * p_o1 * (1 - p_o1) * exp(mean_value)) /
    (exp(mean_value) * p_o1 + 1 - p_o1)^2
  h2 <- (a2 * p_o2 * (1 - p_o2) * exp(mean_value)) /
    (exp(mean_value) * p_o2 + 1 - p_o2)^2

  # Variance using the observed cases (Equation 13 in PMID: 17703496)
  v_obs <- (((alpha1 * alpha1 * beta1) + (c1 * c1 * d1)) * e1 * e1) +
    (((alpha2 * alpha2 * beta2) + (c2 * c2 * d2)) * e2 * e2)

  # Variance due to informative missingness
  # (Equation 16 with correlation in PMID: 17703496)
  v_delta <- (h1 * h1 * var_value * e1 * e1) + (h2 * h2 * var_value * e2 * e2) -
    2 * rho * h1 * h2 * sqrt(var_value) * sqrt(var_value) * e1 * e2

  # Variance using the randomised sample (Equation 10 in PMID: 17703496)
  v_all <- v_obs + v_delta

  # Include trial-specific adjusted logORs and SEs in the initial dataset
  final <- data.frame(cbind(data, round(logor, 3),  round(sqrt(v_all), 3)))
  colnames(final) <- c("id",
                       "r1", "r2",
                       "m1", "m2",
                       "n1", "n2",
                       "t1", "t2",
                       "EM", "se.EM")

  return(final)
}

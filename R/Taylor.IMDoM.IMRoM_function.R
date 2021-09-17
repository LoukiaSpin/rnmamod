#' Pattern-mixture model with Taylor series for continuous outcomes
#'
#' @description The function applies pattern-mixture model under a specific assumption about the missingness parameter in trial-arms with
#'   missing participant outcome data (MOD) and uses the Taylor series to obtain the effect size and standard error for each trial.
#'
#' @param data A data-frame in the long arm-based format containing. Two arm-trials occupy one row in the data-frame. Multi-arm trials, occupy as
#'   many rows as the number of possible comparisons among the interventions. See 'Format' for the specification of the columns.
#' @param measure Character string indicating the effect measure with values \code{"MD"}, \code{"SMD"}, or \code{"ROM"} for the mean difference,
#'   standardised mean difference and ratio of means, respectively.
#' @param mean.value A numeric value for the mean of the normal distribution of the informative missingness parameter. The same value is considered for all trial-arms of the dataset.
#'   The default argument is 0 and corresponds to the missing-at-random assumption.
#'   For the informative missingness ratio of means, the mean value is defined in the logarithmic scale.
#' @param var.value A positive non-zero number for the variance of the normal distribution of the informative missingness parameter. When the \code{measure} is \code{"MD"}, or \code{"SMD"}
#'   the default argument is 1; When the \code{measure} is \code{"ROM"} the default argument is 0.04
#' @param rho A numeric value in the interval [-1, 1] that indicates the correlation coefficient between two missingness parameters in a trial. The same value is considered across all trials of the dataset.
#' The default argument is 0 and corresponds to uncorrelated missingness parameters.
#'
#' @format The columns of the data-frame in the argument \code{data} refer to the following ordered elements for a binary outcome:
#' \tabular{ll}{
#'  \strong{id} \tab A unique identifier for each trial.\cr
#'  \tab \cr
#'  \strong{y1} \tab The observed mean outcome in the first arm of the comparison.\cr
#'  \tab \cr
#'  \strong{y2} \tab The observed mean outcome in the second arm of the comparison.\cr
#'  \tab \cr
#'  \strong{sd1} \tab The observed standard deviation of the outcome in the first arm of the comparison.\cr
#'  \tab \cr
#'  \strong{sd2} \tab The observed standard deviation of the outcome in the second arm of the comparison.\cr
#'  \tab \cr
#'  \strong{m1} \tab The number of MOD in the first arm of the comparison.\cr
#'  \tab \cr
#'  \strong{m2} \tab The number of MOD in the second arm of the comparison.\cr
#'  \tab \cr
#'  \strong{c1} \tab The number of completers in the first arm of the comparison.\cr
#'  \tab \cr
#'  \strong{c2} \tab The number of completers in the second arm of the comparison.\cr
#'  \tab \cr
#'  \strong{t1} \tab An identified for the intervention in the first arm of the comparison.\cr
#'  \tab \cr
#'  \strong{t2} \tab An identified for the intervention in the second arm of the comparison.\cr
#' }
#'
#' @return A data-frame that additionally includes the following elements:
#' \tabular{ll}{
#'  \strong{EM} \tab The effect size adjusted for MOD and obtained using the Taylor series.\cr
#'  \tab \cr
#'  \strong{se.EM} \tab The standard error of the effect size for MOD and obtained using the Taylor series.\cr
#' }
#'
#' @details The \code{Taylor.IMDoM.IMRoM} found is found in the \code{\link[rnmamod]{unrelated.effects.plot}} function. The latter uses the
#'   the \code{\link[netmeta]{pairwise}} function from the package \href{https://cran.r-project.org/web/packages/netmeta/netmeta.pdf}{netmeta}
#'   to transform the dataset from the wide arm-based format (see, 'Arguments' for \code{data} in \code{\link[rnmamod]{unrelated.effects.plot}})
#'   into the long-arm based format.
#'
#' @seealso \code{\link[rnmamod]{run.model}}, \code{\link[rnmamod]{unrelated.effects.plot}}, \code{\link[netmeta]{pairwise}}
#'
#' @references
#' White IR, Higgins JP, Wood AM. Allowing for uncertainty due to missing data in meta-analysis--part 1: two-stage methods.
#' \emph{Stat Med} 2008;\bold{27}(5):711--27. [\doi{10.1002/sim.3008}]
#'
#' @author {Loukia M. Spineli}
#'
#' @export
Taylor.IMDoM.IMRoM <- function(data, measure, mean.value, var.value, rho){


  # Calculate the probability of observing the outcomes in arm 1 (Control) and 2 (Experimental)
  a1 <- data[, 8]/(data[, 8] + data[, 6]); a2 <- data[, 9]/(data[, 9] + data[, 7])


  # Calculate the adjusted-mean for MOD in the randomised sample of arm 1 and 2
  if(measure == "MD" || measure == "SMD"){

    y.all1 <- data[, 2] + mean.value*(1 - a1); y.all2 <- data[, 3] + mean.value*(1 - a2)

  } else if(measure == "ROM"){

    y.all1 <- data[, 2]*(a1 + exp(mean.value)*(1 - a1)); y.all2 <- data[, 3]*(a2 + exp(mean.value)*(1 - a2))

  }



  # Estimate the adjusted within-trial effect measures (Experimental vs Control)
  if(measure == "MD"){

    MD <- y.all2 - y.all1

  } else if(measure == "SMD"){

    nominator1 <- (data[, 8] - 1)*data[, 4]*data[, 4]; nominator2 <- (data[, 9] - 1)*data[, 5]*data[, 5]
    denominator <- (data[, 8] - 1) + (data[, 9] - 1)
    sd.pooled <- sqrt((nominator1 + nominator2)/denominator)

    SMD <- (y.all2 - y.all1)/sd.pooled

  } else if(measure == "ROM"){

    LROM <- log(y.all2/y.all1)

  }


  #####################################################################
  ## USE OF TAYLOR APPROXIMATION FOR THE TOTAL WITHIN-TRIAL VARIANCE ##
  #####################################################################

  ## Estimating the variance of the mean effect size based on the observed data
  # Derivative of y.all by y.obs per arm
  if(measure == "MD" || measure == "SMD"){

    A1 <- A2 <- 1

  } else if(measure == "SMD"){

    A1 <- A2 <- 1

  } else if(measure == "ROM"){

    A1 <-  a1 + exp(mean.value)*(1 - a1);  A2 <-  a2 + exp(mean.value)*(1 - a2)

  }


  # Variance of y.obs per arm
  B1 <- (data[, 4]*data[, 4])/data[, 8]; B2 <- (data[, 5]*data[, 5])/data[, 9]


  # Derivative of y.all by prob of MOD (i.e. a) per arm
  if(measure == "MD" || measure == "SMD"){

    C1 <- C2 <- -mean.value

  } else if(measure == "ROM"){

    C1 <- data[, 2]*(1 - exp(mean.value)); C2 <- data[, 3]*(1 - exp(mean.value))

  }


  # Variance of prob of MOD
  D1 <- (a1*(1 - a1))/(data[, 8] + data[, 6]); D2 <- (a2*(1 - a2))/(data[, 9] + data[, 7])


  # Derivative of link function for MD, SMD and LROM per arm
  if(measure == "MD"){

    E1 <- E2 <- 1

  } else if(measure == "SMD"){

    E1 <- E2 <- 1/sd.pooled

  } else if(measure == "ROM"){

    E1 <- 1/y.all1;  E2 <- 1/y.all2

  }


  # Variance using the observed cases (Equation 13 in PMID: 17703496)
  v.obs <- (A1*A1*B1 + C1*C1*D1)*E1*E1 + (A2*A2*B2 + C2*C2*D2)*E2*E2


  ## Estimating the variance of the mean effect size arising from the informative missingness parameter
  # Derivative of y.all by delta per arm
  if(measure == "MD" || measure == "SMD"){

    H1 <- (1 - a1); H2 <- (1 - a2)

  } else if(measure == "ROM"){

    H1 <- data[, 2]*exp(mean.value)*(1 - a1); H2 <- data[, 3]*exp(mean.value)*(1 - a2)

  }


  # Variance due to informative missingness
  v.delta <-  H1*H1*var.value*E1*E1 + H2*H2*var.value*E2*E2 - 2*rho*H1*H2*E1*E2*sqrt(var.value)*sqrt(var.value)


  # Variance using the randomised sample
  v.all <- v.obs + v.delta


  # Include trial-specific adjusted MDs, SDMs and SEs in the initial dataset
  if (measure == "MD") {
    final <- data.frame(cbind(data, round(MD, 3),  round(sqrt(v.all), 3)))
  } else if (measure == "SMD") {
    final <- data.frame(cbind(data, round(SMD, 3),  round(sqrt(v.all), 3)))
  } else if (measure == "ROM") {
    final <- data.frame(cbind(data, round(LROM, 3),  round(sqrt(v.all), 3)))
  }
  colnames(final) <- c("id", "mean1", "mean2", "sd1", "sd2", "m1", "m2", "c1", "c2", "t1", "t2", "EM", "se.EM")

  return(final)
}
## END

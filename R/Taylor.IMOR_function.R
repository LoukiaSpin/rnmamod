#' Pattern-mixture model with Taylor series for binary outcomes
#'
#' @description The function applies pattern-mixture model under a specific assumption about the missingness parameter in trial-arms with
#'   \bold{binary} missing participant outcome data (MOD) and uses the Taylor series to obtain the effect size and standard error for each trial (White et al., 2008).
#'
#' @param data A data-frame in the long arm-based format. Two arm-trials occupy one row in the data-frame. Multi-arm trials occupy as
#'   many rows as the number of possible comparisons among the interventions. See 'Format' for the specification of the columns.
#' @param mean.value A numeric value for the mean of the normal distribution of the informative missingness odds ratio in the logarithmic scale. The same value is considered for all trial-arms of the dataset.
#'   The default argument is 0 and corresponds to the missing-at-random assumption.
#' @param var.value A positive non-zero number for the variance of the normal distribution of the informative missingness odds ratio in the logarithmic scale. The default argument is 1.
#' @param rho A numeric value in the interval [-1, 1] that indicates the correlation coefficient between two missingness parameters in a trial. The same value is considered across all trials of the dataset.
#'   The default argument is 0 and corresponds to uncorrelated missingness parameters.
#'
#' @format The columns of the data-frame in the argument \code{data} refer to the following ordered elements for a binary outcome:
#' \tabular{ll}{
#'  \strong{id} \tab A unique identifier for each trial.\cr
#'  \tab \cr
#'  \strong{r1} \tab The observed  number of events in the first arm of the comparison.\cr
#'  \tab \cr
#'  \strong{r2} \tab The observed  number of events in the second arm of the comparison.\cr
#'  \tab \cr
#'  \strong{m1} \tab The number of MOD in the first arm of the comparison.\cr
#'  \tab \cr
#'  \strong{m2} \tab The number of MOD in the second arm of the comparison.\cr
#'  \tab \cr
#'  \strong{n1} \tab The number of participants randomised in the first arm of the comparison.\cr
#'  \tab \cr
#'  \strong{n2} \tab The number of participants randomised in the second arm of the comparison.\cr
#'  \tab \cr
#'  \strong{t1} \tab An identified for the intervention in the first arm of the comparison.\cr
#'  \tab \cr
#'  \strong{t2} \tab An identified for the intervention in the second arm of the comparison.\cr
#' }
#'
#' @return A data-frame that additionally includes the following elements:
#' \tabular{ll}{
#'  \strong{EM} \tab The odds ratio in the logarithmic scale (log OR) adjusted for MOD and obtained using the Taylor series.\cr
#'  \tab \cr
#'  \strong{se.EM} \tab The standard error of the log OR adjusted for MOD and obtained using the Taylor series.\cr
#' }
#'
#' @details The \code{Taylor.IMOR} function is found in the \code{\link[rnmamod]{unrelated.effects.plot}} function. The latter uses the
#'   the \code{\link[netmeta]{pairwise}} function from the package \href{https://CRAN.R-project.org/package=netmeta}{netmeta}
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
Taylor.IMOR <- function(data, mean.value, var.value, rho){

  for(i in 1:length(data[, 1])){

    # Add 0.5 continuity correction when there is at least on zero cell
    if(data[i, 2] == 0 || data[i, 3] == 0 || data[i, 6] - data[i, 4] - data[i, 2] == 0 || data[i, 7] - data[i, 5] - data[i, 3] == 0){
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

  # Calculate the probability of event among completers in arm 1 and 2
  p.o1 <- data[, 2]/(data[, 6] - data[, 4]); p.o2 <- data[, 3]/(data[, 7] - data[, 5])

  # Calculate the probability of missing outcome data in arm 1 and 2
  a1 <- data[, 4]/data[, 6]; a2 <- data[, 5]/data[, 7]

  # Calculate the probability of event in randomised sample in arm 1 and 2 via pattern-mixture model
  p.all1 <- (1 - a1)*p.o1 + a1*( (exp(mean.value)*p.o1)/(exp(mean.value)*p.o1 + 1 - p.o1) ); p.all2 <- (1 - a2)*p.o2 + a2*( (exp(mean.value)*p.o2)/(exp(mean.value)*p.o2 + 1 - p.o2) )

  # Estimates the odds ratio in the logarithmic scale
  logOR <- log(p.all1/(1 - p.all1)) - log(p.all2/(1 - p.all2))

  #################################
  ## USE OF TAYLOR APPROXIMATION ##
  #################################

  # Derivative of p.all by p.o per arm (first term in Equation 14 in PMID: 17703496)
  A1 <- 1 - a1 + (a1*exp(mean.value))/(exp(mean.value)*p.o1 + 1 - p.o1)^2; A2 <- 1 - a2 + (a2*exp(mean.value))/(exp(mean.value)*p.o2 + 1 - p.o2)^2

  # Variance of p.o per arm (second term in Equation 14 in PMID: 17703496)
  B1 <- (p.o1*(1 - p.o1))/(data[, 6] - data[, 4]); B2 <- (p.o2*(1 - p.o2))/(data[, 7] - data[, 5])

  # Derivative of p.all by prob of MOD (i.e. a) per arm (third term in Equation 14 in PMID: 17703496)
  C1 <- (p.o1*(1 - p.o1)*(exp(mean.value) - 1))/(exp(mean.value)*p.o1 + 1 - p.o1); C2 <- (p.o2*(1 - p.o2)*(exp(mean.value) - 1))/(exp(mean.value)*p.o2 + 1 - p.o2)

  # Variance of prob of MOD (i.e. a) per arm (fourth term in Equation 14 in PMID: 17703496)
  D1 <- (a1*(1 - a1))/data[, 6]; D2 <- (a2*(1 - a2))/data[, 7]

  # Variance of log odds using delta-method per arm
  E1 <- 1/(p.all1*(1 - p.all1)); E2 <- 1/(p.all2*(1 - p.all2))

  # Derivative of p.all by delta per arm (second Equation after Equation (15) in PMID: 17703496)
  H1 <- (a1*p.o1*(1 - p.o1)*exp(mean.value))/(exp(mean.value)*p.o1 + 1 - p.o1)^2; H2 <- (a2*p.o2*(1 - p.o2)*exp(mean.value))/(exp(mean.value)*p.o2 + 1 - p.o2)^2

  # Variance using the observed cases (Equation 13 in PMID: 17703496)
  v.obs <- (A1*A1*B1 + C1*C1*D1)*E1*E1 + (A2*A2*B2 + C2*C2*D2)*E2*E2

  # Variance due to informative missingness (Equation 16 with correlation in PMID: 17703496)
  v.delta <-  H1*H1*var.value*E1*E1 + H2*H2*var.value*E2*E2 - 2*rho*H1*H2*sqrt(var.value)*sqrt(var.value)*E1*E2

  # Variance using the randomised sample (Equation 10 in PMID: 17703496)
  v.all <- v.obs + v.delta

  # Include trial-specific adjusted logORs and SEs in the initial dataset
  final <- data.frame(cbind(data, round(logOR, 3),  round(sqrt(v.all), 3)))
  colnames(final) <- c("id", "r1", "r2", "m1", "m2", "n1", "n2", "t1", "t2", "EM", "se.EM")

  return(final)
}


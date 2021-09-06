#' Robustness index: Investigate the impact of missing participant outcome data
#'
#' @description This function calculates the robustness index, a novel index that quantifies the overall divergence of the sensitivity analysis results from the primary analysis results
#'   (missing-at-random assumption), and considers objective decision rules to infer the presence of lack of robustness of the primary analysis results when conducting a sensitivity analysis (Spineli et al., 2021).
#'   Currently, \code{robustness.index} is used concerning the impact of missing participant outcome data.
#'
#' @param sens An object of S3 class \code{\link{run.sensitivity}}. See 'Value' in \code{\link{run.sensitivity}}.
#' @param threshold A number indicating the threshold of robustness, that is, the minimally allowed deviation between the primary analysis and re-analysis results. See 'Details' below.
#'
#' @return \code{robustness.index} prints on the R console a message on the threshold of robustness determined by the user in green text.
#' Then, the function returns the following list of elements:
#' \tabular{ll}{
#'  \code{RI} \tab A a numeric scalar or vector on the robustnessindex values. In the case of a pairwise meta-analysis (PMA), \code{RI} is scalar as only one summary effect size is obtained.
#'    In the case of network meta-analysis (NMA), \code{RI} is a vector with length equal to the number of possible pairwise comparisons; one robustness index per possible comparison.\cr
#'  \tab \cr
#'  \code{robust} \tab A character vector on whether the primary analysis results are \emph{robust} or \emph{frail} to the different re-analyses.\cr
#'  \tab \cr
#'  \code{KLD} \tab A vector or matrix on the Kullback-Leibler divergence (KLD) measure in the summary effect size from a subsequent re-analysis to the primary analysis.
#'    In the case of a PMA, \code{KLD} is a vector with length equal to the number of total analyses.
#'    The latter equals the square of the number of scenarios indicated in argument \code{mean.scenarios} of \code{run.sensitivity}. Therefore, one KLD value per analysis.
#'    In the case of NMA, \code{RI} is a matrix with number of rows equal to the number of total analyses and number of columns equal to the number of possible pairwise comparisons;
#'    one KLD value per analysis and possible comparison.\cr
#' }
#'
#' @details The user may consider the values 0.28 and 0.17 in the argument \code{threshold} for binary and continuous outcome data (the default values), respectively, or consider other plausible values.
#'   Spineli et al. (2021) offers a discussion on specifying the \code{threshold} of robustness.
#'
#'   Currently, the primary analysis is considered to be the MAR assumption and it corresponds to the middle of the numbers in the argument \code{mean.scenarios}
#'   of the \code{run.sensitivity} function (see 'Arguments' and 'Details' in \code{run.sensitivity}).
#'
#'   In \code{robust}, the value \code{"robust"} appears when \code{RI} \eqn{<} \code{threshold}); otherwise, the value \code{"frail"} appears.
#'
#'   \code{robustness.index} can be used only for when missing participant outcome data have been extracted for at least one trial. Otherwise, the execution of the function will be stopped and an error message will be printed in the R console.
#'
#' @author {Loukia M. Spineli}
#'
#' @seealso \code{\link{run.sensitivity}}
#'
#' @references
#' Spineli LM, Kalyvas C, Papadimitropoulou K. Quantifying the robustness of primary analysis results: A case study on missing outcome data in pairwise and network meta-analysis. \emph{Res Synth Methods} 2021;\bold{12}(4):475--490. [\doi{10.1002/jrsm.1478}]
#'
#' Kullback S, Leibler RA. On information and sufficiency. \emph{Ann Math Stat} 1951;\bold{22}(1):79--86. [\doi{10.1214/aoms/1177729694}]
#'
#' @examples
#' data("nma.baker2009")
#'
#' \dontrun{
#' # Perform a random-effects network meta-analysis
#' res <- run.model(data = nma.baker2009,
#'                  measure = "OR",
#'                  model = "RE",
#'                  assumption = "IDE-ARM",
#'                  heter.prior = list("halfnormal", 0, 1),
#'                  mean.misspar = c(0, 0),
#'                  var.misspar = 1,
#'                  D = 1,
#'                  n.chains = 3,
#'                  n.iter = 10000,
#'                  n.burnin = 1000,
#'                  n.thin = 1)
#'
#' # Perform the sensitivity analysis (using the 'default' of the argument 'mean.scenarios')
#' res.sens <- run.sensitivity(full = res,
#'                             var.misspar = 1,
#'                             n.chains = 3,
#'                             n.iter = 10000,
#'                             n.burnin = 1000,
#'                             n.thin = 1)
#'
#' # Calculate the robustness index
#' robustness.index(sens = res.sens, threshold = 0.28)
#' }
#'
#' @export
robustness.index <- function(sens, threshold){


  options(warn = -1)

  if (is.list(sens)) {
    ES.mat <- sens[[1]]
    measure <- sens[[2]]
    n.scenar <- sens[[3]]
    primary.scenar <- 1

  } else {
    ES.mat <- sens$EM

    measure <- sens$measure

    n.scenar <- length(sens$scenarios)^2

    if (any(is.na(sens))) {
      stop("Missing participant outcome data have *not* been collected. This function cannot be used.", call. = F)
      return(NA)
    }

    primary.scenar <- median(1:n.scenar)
  }


  nt <- (1 + sqrt(1 + 8*(length(ES.mat[, 1])/sqrt(n.scenar)^2)))/2  # The quadratic formula for the roots of the general quadratic equation



  if (missing(threshold) & is.element(measure, "OR")) {
    threshold <- 0.28
    #message("The value 0.28 was assigned on 'threshold' by default")
    message(cat(paste0("\033[0;", col = 32, "m", txt = "The value 0.28 was assigned on 'threshold' by default", "\033[0m", "\n")))
  } else if (missing(threshold) & is.element(measure, c("MD", "SMD", "ROM"))) {
    threshold <- 0.17
    #message("The value 0.17 was assigned on 'threshold' by default")
    message(cat(paste0("\033[0;", col = 32, "m", txt = "The value 0.17 was assigned on 'threshold' by default", "\033[0m", "\n")))
  } else {
    threshold <- threshold
    #message(paste("The value", threshold, "was assigned on 'threshold' for", effect.measure.name(full$measure)))
    message(cat(paste0("\033[0;", col = 32, "m", txt = paste("The value", threshold, "was assigned on 'threshold' for", effect.measure.name(measure)), "\033[0m", "\n")))
  }




  ## Function for the Kullback-Leibler Divergence (comparing two univariate normal distributions)
  KLD.measure.univ <- function(mean.y, sd.y, mean.x, sd.x){

    # x is the 'truth' (e.g. the MAR assumption)
    KLD.xy <- 0.5*(((sd.x/sd.y)^2) + ((mean.y - mean.x)^2)/(sd.y^2) - 1 + 2*log(sd.y/sd.x))

    return(list(KLD.xy = KLD.xy))
  }



  ## A matrix of effect estimates of MCMC standard deviations (or standard errors) for all possible comparisons under each scenario
  sd.mat <- mean.mat <- matrix(NA, nrow = n.scenar, ncol = (nt*(nt - 1))/2)

  for(i in 1:n.scenar){

    for(j in 1:((nt*(nt - 1))/2)){

      mean.mat[i, j] <- ES.mat[j + ((nt*(nt - 1))/2)*(i - 1), 1]

      sd.mat[i, j] <- ES.mat[j + ((nt*(nt - 1))/2)*(i - 1), 2]

    }
  }

  kldxy <- list()

  RI <- rep(NA, (nt*(nt - 1))/2)

  for(i in 1:(nt*(nt - 1))/2){ ## We are interested in all possible pairwise comparisons of the network

    kldxy[[i]] <- rep(NA, n.scenar)

    for(j in 1:n.scenar){ #for(j in (1:n.scenar)[-primary.scenar])


      ## Returns the KLD of informative scenario j when compared with the primary analysis for comparison i
      kldxy[[i]][j] <- KLD.measure.univ(mean.mat[j, i], sd.mat[j, i], mean.mat[primary.scenar, i], sd.mat[primary.scenar, i])[[1]]

    }
    kldxy[[i]][primary.scenar] <- 0  ## This refers to the primary analysis

    ## Returns the Robustness Index of comparison i across all informative scenarios
    RI[i] <- sqrt(round(t(kldxy[[i]]) %*% kldxy[[i]], 3))
    #RI[i] <- sqrt(round(t(kldxy[[i]][-primary.scenar]) %*% kldxy[[i]][-primary.scenar], 3))

  }

  KLD <- matrix(unlist(kldxy), nrow = (nt*(nt - 1))/2, ncol = n.scenar, byrow = T)


  robust <- ifelse(RI < threshold, "robust", "frail")

  return(list(RI = RI, robust = robust, KLD = KLD, measure = sens$measure, threshold = threshold,
              scenarios = sens$scenarios))
}









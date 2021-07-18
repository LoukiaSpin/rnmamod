#' Sensitivity analysis for aggregate missing outcome participant data
#'
#' @description This function performs sensitivity analysis by applying pairwise meata-analysis (PMA) or network meta-analysis (NMA) for a series of different scenarios about the informative missingness parameter,
#'
#' @param full An object of S3 class \code{\link{run.model}}. See 'Value' in \code{\link{run.model}}.
#' @param assumption Character string indicating the structure of the informative missingness parameter.
#'   Set \code{assumption} equal to one of the following: \code{"HIE-ARM"}, or \code{"IDE-ARM"}.
#'   The default argument is \code{"IDE-ARM"}. The abbreviations \code{"IDE"}, and \code{"HIE"}, stand for identical, and hierarchical, respectively. See 'Details'.
#' @param mean.scenarios A vector with numeric values for the mean of the normal distribution of the informative missingness parameter (see 'Details').
#'   The vector should have a length of 5 or larger \emph{positive odd integer}.
#'   The missing-at-random (MAR) assumption should be the median of the vector, so that the same number of informative scenarios appear before and after the MAR.
#'   The default arguments are c(-log(3), -log(2), log(0.9999), log(2), log(3)) and c(-2, -1, 0, 1, 2) for binary and continuous outcome dara, respectively.
#' @param var.misspar A positive non-zero number for the variance of the normal distribution of the informative missingness parameter. When the \code{measure} is \code{"OR"}, \code{"MD"}, or \code{"SMD"}
#'   the default argument is 1; When the \code{measure} is \code{"ROM"} in \code{\link{run.model}} the default argument is 0.04
#' @param n.chains Integer specifying the number of chains for the MCMC sampling; an argument of the \code{\link[R2jags]{jags}} function of the R-package \href{https://CRAN.R-project.org/package=R2jags}{R2jags}.
#'   The default argument is 2.
#' @param n.iter Integer specifying the number of Markov chains for the MCMC sampling; an argument of the \code{\link[R2jags]{jags}} function of the R-package \href{https://CRAN.R-project.org/package=R2jags}{R2jags}.
#'   The default argument is 10000.
#' @param n.burnin Integer specifying the number of iterations to discard at the beginning of the MCMC sampling; an argument of the \code{\link[R2jags]{jags}} function of the R-package \href{https://CRAN.R-project.org/package=R2jags}{R2jags}.
#'   The default argument is 1000.
#' @param n.thin Integer specifying the thinning rate for the MCMC sampling; an argument of the \code{\link[R2jags]{jags}} function of the R-package \href{https://CRAN.R-project.org/package=R2jags}{R2jags}.
#'   The default argument is 1.
#'
#' @return A list of R2jags outputs on the summaries of the posterior distribution, and the Gelman-Rubin convergence diagnostic (Gelman et al., 1992) of the following monitored parameters for a random-effects PMA:
#' \tabular{ll}{
#'  \code{EM} \tab The estimated summary effect measure (according to the argument \code{measure} in \code{run.model}).\cr
#'  \tab \cr
#'  \code{tau} \tab The between-trial standard deviation. This element does not appear in the case of a fixed-effect PMA.\cr
#' }
#'
#' In a random-effects NMA, \code{EM} refer to all possible pairwise comparisons of interventions in the network. Furthermore, \code{tau} is typically assumed to be common for all observed comparisons in the network.
#'
#' @details The model as specified by the arguments of \code{run.sensitivity} and \code{run.model} (the latter via the argument \code{full}) runs in \code{JAGS} and the progress of the simulation appears in the R console.
#'   The number of times \code{run.sensitivity} is used appears in the R console as a text in red and it equals the number of scenarios specified in argument \code{mean.scenarios} (see 'Examples').
#'   The output of \code{run.sensitivity} is used as an S3 object by other functions of the package function to be processed further and provide an end-user-ready output.
#'
#'   In the case of PMA, \code{EM} and \code{tau} have as many rows as the \emph{square of the number of scenarios} indicated in argument \code{mean.scenarios}.
#'   In the case of NMA, each possible pairwise comparison is estimated as many times as the square of the number of scenarios indicated in argument \code{mean.scenarios}.
#'
#'   The informative missingness parameter is assumed to differ only across the interventions of the dataset. Therefore, the user can specify this parameter to be arm-specific and identical (\code{assumption} = \code{"IDE-ARM}),
#'   or arm-specific and hierarchical (\code{assumption} = \code{"HIE-ARM}) (Spineli et al., 2021).
#'
#'   The number of scenarios in \code{mean.scenarios} should be equal to or more than 5 (a positive odd integer) to allow for an adequate sensitivity analysis.
#'   It is important that the scenario corresponding to the MAR assumption is the middle of the numbers in \code{mean.scenarios}.
#'   Under the inforamtive missingness of mean difference parameter, the MAR assumption is equal to 0.
#'   Under the informative missingness odds ratio parameter and the informative missingness ratio of means parameter, the MAR assumption is equal to 1. Both parameters are analysed in the logarithmic scale. We advise using the value \code{0.999} rather than log(1)
#'   in \code{mean.scenarios}; otherwise, the execution of the fucntion will be stopped and the error 'Invalid parent values' will be printed in the R console.
#'   Currently, there are no empirically-based prior distributions for the informative missingness parameters. The users may refer to White et al. (2008), Mavridis et al. (2015), Turner et al. (2015) and Spineli (2019) to determine \code{mean.scenarios}
#'   for an informative missingness mechanism and select a proper value for \code{var.misspar}.
#'
#'   \code{run.sensitivity} does not contain the arguments \code{data}, \code{measure}, \code{model}, and \code{heter.prior} that are found in \code{run.model}.
#'   This is to prevent misspecifying the Bayesian model as it would make the comparison of the primary analysis (via \code{run.model}) with the re-analyses meaningless.
#'   Instead, these arguments are contained in the argument \code{full} of the function. Therefore, the user needs first to apply \code{run.model}, and then use \code{run.sensitivity} (see, 'Examples').
#'
#'   \code{run.sensitivity} can be used only for when missing participant outcome data have been extracted for at least one trial. Otherwise, the execution of the function will be stopped and an error message will be printed in the R console.
#'
#' @author {Loukia M. Spineli}
#'
#' @seealso \code{\link{run.model}}, \href{https://CRAN.R-project.org/package=R2jags}{R2jags}
#'
#' @references
#' Spineli LM, Kalyvas C, Papadimitropoulou K. Quantifying the robustness of primary analysis results: A case study on missing outcome data in pairwise and network meta-analysis. \emph{Res Synth Methods} 2021;\bold{12}(4):475--490. [\doi{10.1002/jrsm.1478}]
#'
#' Spineli LM. An empirical comparison of Bayesian modelling strategies for missing binary outcome data in network meta-analysis. \emph{BMC Med Res Methodol} 2019;\bold{19}(1):86. [\doi{10.1186/s12874-019-0731-y}]
#'
#' Mavridis D, White IR, Higgins JP, Cipriani A, Salanti G. Allowing for uncertainty due to missing continuous outcome data in pairwise and network meta-analysis. \emph{Stat Med} 2015;\bold{34}(5):721--741. [\doi{10.1002/sim.6365}]
#'
#' Turner NL, Dias S, Ades AE, Welton NJ. A Bayesian framework to account for uncertainty due to missing binary outcome data in pairwise meta-analysis. \emph{Stat Med} 2015;\bold{34}(12):2062--2080. [\doi{10.1002/sim.6475}]
#'
#' White IR, Higgins JP, Wood AM. Allowing for uncertainty due to missing data in meta-analysis--part 1: two-stage methods. \emph{Stat Med} 2008;\bold{27}(5):711--727. [\doi{10.1002/sim.3008}]
#'
#' Gelman, A, Rubin, DB. Inference from iterative simulation using multiple sequences. Stat Sci. 1992;7:457–472.
#'
#' @examples
#' data("nma.liu2013.RData")
#'
#' # Perform the sensitivity analysis (using the 'default' of the argument 'mean.scenarios')
#' run.sensitivity(full = res1,
#'                 assumption = "IDE-ARM",
#'                 var.misspar = 1,
#'                 n.chains = 3,
#'                 n.iter = 10000,
#'                 n.burnin = 1000,
#'                 n.thin = 1)
#'
#' @export
run.sensitivity <- function(full, assumption, mean.scenarios, var.misspar, n.chains, n.iter, n.burnin, n.thin){


  ## Turn off warning when variables in the 'data.jag' are not used
  options(warn = -1)

  data <- full$data
  measure <- full$measure
  model <- full$model
  heter.prior <- full$heter.prior
  D <- full$D


  ## Prepare the dataset for the R2jags
  item <- data.preparation(data, measure)


  if (unique(na.omit(unlist(item$I))) == 0) {
    stop("Missing participant outcome data have *not* been collected. This function cannot be used.", call. = F)
    return(NA)
  }


  ## Default arguments
  assumption <- if (missing(assumption)) {
    "IDE-ARM"
  } else if (!is.element(assumption,  c("IDE-ARM", "HIE-ARM"))) {
    stop("Insert 'IDE-ARM', or 'HIE-ARM'", call. = F)
  } else {
    assumption
  }
  ## Scenarios for missingness mechanism in an intervention (PMID: 30223064)
  mean.scenarios <- if (missing(mean.scenarios) & is.element(measure, c("MD", "SMD"))) {
    message(cat(paste0("\033[0;", col = 32, "m", txt = "The following vector of scenarios was considered by default: c(-2, -1, 0, 1, 2)", "\033[0m", "\n")))
    c(-2, -1, 0, 1, 2)
  } else if (missing(mean.scenarios) & is.element(measure, c("OR", "ROM"))) {
    message(cat(paste0("\033[0;", col = 32, "m", txt = "The following vector of scenarios was considered by default: c(-log(3), -log(2), log(0.9999), log(2), log(3))", "\033[0m", "\n")))
    c(-log(3), -log(2), log(0.9999), log(2), log(3))
  } else if (length(mean.scenarios) < 5) {
    stop("The argument 'mean.scenarios' must have a length of at least 5", call. = F)
  } else {
    mean.scenarios
  }
  var.misspar <- ifelse(missing(var.misspar) & (is.element(measure, c("OR", "MD", "SMD"))), 1, ifelse(missing(var.misspar) & measure == "ROM", 0.2^2, var.misspar))
  n.chains <- ifelse(missing(n.chains), 2, n.chains)
  n.iter <- ifelse(missing(n.iter), 10000, n.iter)
  n.burnin <- ifelse(missing(n.burnin), 1000, n.burnin)
  n.thin <- ifelse(missing(n.thin), 1, n.thin)



  ## A 2x2 matrix of 25 reference-specific scenarios (PMID: 30223064)
  mean.misspar <- as.matrix(cbind(rep(mean.scenarios, each = length(mean.scenarios)), rep(mean.scenarios, length(mean.scenarios)))) # 2nd column refers to the reference intervention (control in MA)


  ## Prepare parameters for JAGS
  jagsfit <- data.jag <- list()


  ## Parameters to save
  param.jags <- if (model == "RE") {
    c("EM", "tau")
  } else {
    c("EM")
  }



  ## Calculate time needed for all models
  for(i in 1:length(mean.misspar[, 1])){

    data.jag[[i]] <- list("m" = item$m,
                          "N" = item$N,
                          "t" = item$t,
                          "na" = item$na,
                          "nt" = item$nt,
                          "ns" = item$ns,
                          "ref" = item$ref,
                          "I" = item$I,
                          "meand.phi" = mean.misspar[i, ],
                          "precd.phi" = 1/var.misspar,
                          "D" = D,
                          "heter.prior" = heter.prior,
                          "eff.mod2" = matrix(0, nrow = item$ns, ncol = max(item$na)),
                          "eff.mod" = rep(0, item$ns))


    if (is.element(measure, c("MD", "SMD", "ROM"))) {
      data.jag[[i]] <- append(data.jag[[i]], list("y.o" = item$y0, "se.o" = item$se0))
    } else if (measure == "OR") {
      data.jag[[i]] <- append(data.jag[[i]], list("r" = item$r))
    }



    message(paste(i, "out of", length(mean.misspar[, 1]), "total scenarios"))
    jagsfit[[i]] <- jags(data = data.jag[[i]],
                         parameters.to.save = param.jags,
                         model.file = textConnection(prepare.model(measure, model, covar.assumption = "NO", assumption)),
                         n.chains = n.chains,
                         n.iter = n.iter,
                         n.burnin = n.burnin,
                         n.thin = n.thin)
  }


  ## Obtain the posterior distribution of the necessary model paramters
  EM <- do.call(rbind,lapply(1:length(mean.misspar[, 1]), function(i) jagsfit[[i]]$BUGSoutput$summary[1:(item$nt*(item$nt - 1)*0.5), c("mean", "sd", "2.5%", "97.5%", "Rhat", "n.eff")]))
  if (model == "RE") {
    tau <- do.call(rbind,lapply(1:length(mean.misspar[, 1]), function(i) jagsfit[[i]]$BUGSoutput$summary["tau", c("50%", "sd", "2.5%", "97.5%", "Rhat", "n.eff")]))
  } else {
    tau <- NA
  }


  ## Return results
  results <- if (model == "RE"){
    list(EM = EM, tau = tau, measure = measure, scenarios = mean.scenarios, D = D, heter = heter.prior)
  } else {
    list(EM = EM, measure = measure, scenarios = mean.scenarios, D = D)
  }

  return(results)

}


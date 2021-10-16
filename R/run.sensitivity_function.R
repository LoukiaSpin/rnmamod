#' Sensitivity analysis for aggregate missing outcome participant data
#'
#' @description This function performs sensitivity analysis by applying pairwise
#'   meta-analysis (PMA) or network meta-analysis (NMA) for a series of
#'   different scenarios about the informative missingness parameter,
#'
#' @param full An object of S3 class \code{\link{run_model}}. See 'Value' in
#'   \code{\link{run_model}}.
#' @param mean_scenarios A vector with numeric values for the mean of the normal
#'    distribution of the informative missingness parameter (see 'Details').
#'   The vector should have a length of 5 or larger \emph{positive odd integer}.
#'   The missing-at-random (MAR) assumption should be the median of the vector,
#'   so that the same number of informative scenarios appear before and after
#'   the MAR. The default arguments are c(-log(3), -log(2), log(0.9999), log(2),
#'   log(3)) and c(-2, -1, 0, 1, 2) for binary and continuous outcome data,
#'   respectively.
#' @param var_misspar A positive non-zero number for the variance of the normal
#'   distribution of the informative missingness parameter. When the
#'   \code{measure} is \code{"OR"}, \code{"MD"}, or \code{"SMD"} the default
#'   argument is 1; When the \code{measure} is \code{"ROM"} in
#'   \code{\link{run_model}} the default argument is 0.04
#' @param n_chains Integer specifying the number of chains for the MCMC
#'   sampling; an argument of the \code{\link[R2jags]{jags}} function of the
#'   R-package \href{https://CRAN.R-project.org/package=R2jags}{R2jags}.
#'   The default argument is 2.
#' @param n_iter Integer specifying the number of Markov chains for the MCMC
#'   sampling; an argument of the \code{\link[R2jags]{jags}} function of the
#'   R-package \href{https://CRAN.R-project.org/package=R2jags}{R2jags}.
#'   The default argument is 10000.
#' @param n_burnin Integer specifying the number of iterations to discard at the
#'   beginning of the MCMC sampling; an argument of the
#'   \code{\link[R2jags]{jags}} function of the R-package
#'   \href{https://CRAN.R-project.org/package=R2jags}{R2jags}.
#'   The default argument is 1000.
#' @param n_thin Integer specifying the thinning rate for the MCMC sampling; an
#'   argument of the \code{\link[R2jags]{jags}} function of the R-package
#'   \href{https://CRAN.R-project.org/package=R2jags}{R2jags}.
#'   The default argument is 1.
#'
#' @return A list of R2jags outputs on the summaries of the posterior
#'   distribution, and the Gelman-Rubin convergence diagnostic
#'   (Gelman et al., 1992) of the following monitored parameters for a
#'   random-effects PMA:
#'   \tabular{ll}{
#'    \code{EM} \tab The estimated summary effect measure (according to the
#'    argument \code{measure} in \code{run_model}).\cr
#'    \tab \cr
#'    \code{tau} \tab The between-trial standard deviation. This element does
#'    not appear in the case of a fixed-effect PMA.\cr
#'   }
#'
#'   In a random-effects NMA, \code{EM} refer to all possible pairwise
#'   comparisons of interventions in the network. Furthermore, \code{tau} is
#'   typically assumed to be common for all observed comparisons in the network.
#'
#' @details The model as specified by the arguments of \code{run_sensitivity}
#'   and \code{run_model} (the latter via the argument \code{full}) runs in
#'   \code{JAGS} and the progress of the simulation appears in the R console.
#'   The number of times \code{run_sensitivity} is used appears in the R console
#'   as a text in red and it equals the number of scenarios specified in
#'   argument \code{mean_scenarios} (see 'Examples'). The output of
#'   \code{run_sensitivity} is used as an S3 object by other functions of the
#'   package function to be processed further and provide an end-user-ready
#'   output.
#'
#'   In the case of PMA, \code{EM} and \code{tau} have as many rows as the
#'   \emph{square of the number of scenarios} indicated in argument
#'   \code{mean_scenarios}. In the case of NMA, each possible pairwise
#'   comparison is estimated as many times as the square of the number of
#'   scenarios indicated in argument \code{mean_scenarios}.
#'
#'   The informative missingness parameter is assumed to differ only across the
#'   interventions of the dataset. Therefore, the user can specify this
#'   parameter to be arm-specific and identical
#'   (\code{assumption} = \code{"IDE-ARM}), or arm-specific and hierarchical
#'   (\code{assumption} = \code{"HIE-ARM}) (Spineli et al., 2021).
#'
#'   The number of scenarios in \code{mean_scenarios} should be equal to or more
#'   than 5 (a positive odd integer) to allow for an adequate sensitivity
#'   analysis. It is important that the scenario corresponding to the MAR
#'   assumption is the middle of the numbers in \code{mean_scenarios}.
#'   Under the informative missingness of mean difference parameter, the MAR
#'   assumption is equal to 0. Under the informative missingness odds ratio
#'   parameter and the informative missingness ratio of means parameter, the MAR
#'   assumption is equal to 1. Both parameters are analysed in the logarithmic
#'   scale. We advise using the value \code{0.999} rather than log(1)
#'   in \code{mean_scenarios}; otherwise, the execution of the function will be
#'   stopped and the error 'Invalid parent values' will be printed in the R
#'   console. Currently, there are no empirically-based prior distributions for
#'   the informative missingness parameters. The users may refer to
#'   White et al. (2008), Mavridis et al. (2015), Turner et al. (2015) and
#'   Spineli (2019) to determine \code{mean_scenarios} for an informative
#'   missingness mechanism and select a proper value for \code{var_misspar}.
#'
#'   \code{run_sensitivity} does not contain the arguments \code{data},
#'   \code{measure}, \code{model}, and \code{heter_prior} that are found in
#'   \code{run_model}. This is to prevent misspecifying the Bayesian model as it
#'   would make the comparison of the primary analysis (via \code{run_model})
#'   with the re-analyses meaningless. Instead, these arguments are contained in
#'   the argument \code{full} of the function. Therefore, the user needs first
#'   to apply \code{run_model}, and then use \code{run_sensitivity}
#'   (see, 'Examples').
#'
#'   \code{run_sensitivity} can be used only for when missing participant
#'   outcome data have been extracted for at least one trial. Otherwise, the
#'   execution of the function will be stopped and an error message will be
#'   printed in the R console.
#'
#' @author {Loukia M. Spineli}
#'
#' @seealso \code{\link{run_model}},
#'   \href{https://CRAN.R-project.org/package=R2jags}{R2jags}
#'
#' @references
#' Spineli LM, Kalyvas C, Papadimitropoulou K. Quantifying the robustness of
#' primary analysis results: A case study on missing outcome data in pairwise
#' and network meta-analysis.
#' \emph{Res Synth Methods} 2021;\bold{12}(4):475--490.
#' [\doi{10.1002/jrsm.1478}]
#'
#' Spineli LM. An empirical comparison of Bayesian modelling strategies for
#' missing binary outcome data in network meta-analysis.
#' \emph{BMC Med Res Methodol} 2019;\bold{19}(1):86.
#' [\doi{10.1186/s12874-019-0731-y}]
#'
#' Mavridis D, White IR, Higgins JP, Cipriani A, Salanti G. Allowing for
#' uncertainty due to missing continuous outcome data in pairwise and network
#' meta-analysis. \emph{Stat Med} 2015;\bold{34}(5):721--741.
#' [\doi{10.1002/sim.6365}]
#'
#' Turner NL, Dias S, Ades AE, Welton NJ. A Bayesian framework to account for
#' uncertainty due to missing binary outcome data in pairwise meta-analysis.
#' \emph{Stat Med} 2015;\bold{34}(12):2062--2080. [\doi{10.1002/sim.6475}]
#'
#' White IR, Higgins JP, Wood AM. Allowing for uncertainty due to missing data
#' in meta-analysis--part 1: two-stage methods.
#' \emph{Stat Med} 2008;\bold{27}(5):711--727. [\doi{10.1002/sim.3008}]
#'
#' Gelman, A, Rubin, DB. Inference from iterative simulation using multiple
#' sequences. Stat Sci. 1992;7:457–472.
#'
#' @examples
#' data("nma.baker2009")
#'
#' \dontrun{
#' # Perform a random-effects network meta-analysis
#' res <- run_model(data = nma.baker2009,
#'                  measure = "OR",
#'                  model = "RE",
#'                  assumption = "IDE-ARM",
#'                  heter_prior = list("halfnormal", 0, 1),
#'                  mean_misspar = c(0, 0),
#'                  var_misspar = 1,
#'                  D = 1,
#'                  n_chains = 3,
#'                  n_iter = 10000,
#'                  n_burnin = 1000,
#'                  n_thin = 1)
#'
#' # Perform the sensitivity analysis (missing-at-random assumption)
#' run_sensitivity(full = res,
#'                 var_misspar = 1,
#'                 n_chains = 3,
#'                 n_iter = 10000,
#'                 n_burnin = 1000,
#'                 n_thin = 1)
#' }
#' @export
run_sensitivity <- function(full,
                            mean_scenarios,
                            var_misspar,
                            n_chains,
                            n_iter,
                            n_burnin,
                            n_thin) {

  # Turn off warning when variables in the 'data_jag' are not used
  options(warn = -1)

  data <- full$data
  measure <- full$measure
  assumption <- full$assumption
  model <- full$model
  heter_prior <- full$heter_prior
  D <- full$D

  # Prepare the dataset for the R2jags
  item <- data_preparation(data, measure)

  if (unique(na.omit(unlist(item$I))) == 0) {
    stop("Missing participant outcome data have *not* been collected.
         This function cannot be used.", call. = F)
    return(NA)
  }

  # Scenarios for missingness mechanism in an intervention (PMID: 30223064)
  mean_scenarios <- if (missing(mean_scenarios) &
                        is.element(measure, c("MD", "SMD"))) {
    message(cat(paste0("\033[0;",
                       col = 32,
                       "m",
                       txt = "The following vector of scenarios was considered
                       by default: c(-2, -1, 0, 1, 2)", "\033[0m", "\n")))
    c(-2, -1, 0, 1, 2)
  } else if (missing(mean_scenarios) & is.element(measure, c("OR", "ROM"))) {
    message(cat(paste0("\033[0;",
                       col = 32,
                       "m",
                       txt = "The following vector of scenarios was considered
                       by default: c(-log(3), -log(2), log(0.9999), log(2),
                       log(3))", "\033[0m", "\n")))
    c(-log(3), -log(2), log(0.9999), log(2), log(3))
  } else if (length(mean_scenarios) < 5) {
    stop("The argument 'mean_scenarios' must have a length of at least 5",
         call. = F)
  } else {
    mean_scenarios
  }
  var_misspar <- ifelse(missing(var_misspar) &
                          (is.element(measure, c("OR", "MD", "SMD"))), 1,
                        ifelse(missing(var_misspar) & measure == "ROM", 0.2^2,
                               var_misspar))
  n_chains <- ifelse(missing(n_chains), 2, n_chains)
  n_iter <- ifelse(missing(n_iter), 10000, n_iter)
  n_burnin <- ifelse(missing(n_burnin), 1000, n_burnin)
  n_thin <- ifelse(missing(n_thin), 1, n_thin)

  # A 2x2 matrix of 25 reference-specific scenarios (PMID: 30223064)
  mean_misspar <- as.matrix(cbind(rep(mean_scenarios,
                                      each = length(mean_scenarios)),
                                  rep(mean_scenarios, length(mean_scenarios))))

  # Prepare parameters for JAGS
  jagsfit <- data_jag <- list()

  ## Parameters to save
  param_jags <- if (model == "RE") {
    c("EM", "tau")
  } else {
    c("EM")
  }

  # Calculate time needed for all models
  for (i in seq_len(length(mean_misspar[, 1]))) {
    data_jag[[i]] <- list("m" = item$m,
                          "N" = item$N,
                          "t" = item$t,
                          "na" = item$na,
                          "nt" = item$nt,
                          "ns" = item$ns,
                          "ref" = item$ref,
                          "I" = item$I,
                          "meand.phi" = mean_misspar[i, ],
                          "precd.phi" = 1 / var_misspar,
                          "D" = D,
                          "heter.prior" = heter_prior,
                          "eff.mod2" = matrix(0,
                                              nrow = item$ns,
                                              ncol = max(item$na)),
                          "eff.mod" = rep(0, item$ns))

    if (is.element(measure, c("MD", "SMD", "ROM"))) {
      data_jag[[i]] <- append(data_jag[[i]],
                              list("y.o" = item$y0, "se.o" = item$se0))
    } else if (measure == "OR") {
      data_jag[[i]] <- append(data_jag[[i]], list("r" = item$r))
    }

    message(paste(i, "out of", length(mean_misspar[, 1]), "total scenarios"))
    jagsfit[[i]] <- jags(data = data_jag[[i]],
                         parameters.to.save = param_jags,
                         model.file =
                           textConnection(prepare_model(measure,
                                                        model,
                                                        covar_assumption = "NO",
                                                        assumption)),
                         n.chains = n_chains,
                         n.iter = n_iter,
                         n.burnin = n_burnin,
                         n.thin = n_thin)
  }

  # Obtain the posterior distribution of the necessary model paramters
  EM <- do.call(rbind,
                lapply(
                  seq_len(length(mean_misspar[, 1])),
                  function(i) jagsfit[[i]]$BUGSoutput$summary[
                    1:(item$nt * (item$nt - 1) * 0.5),
                    c("mean", "sd", "2.5%", "97.5%", "Rhat", "n.eff")]))
  if (model == "RE") {
    tau <- do.call(rbind,
                   lapply(
                     seq_len(length(mean_misspar[, 1])),
                     function(i) jagsfit[[i]]$BUGSoutput$summary[
                       "tau", c("50%", "sd", "2.5%", "97.5%", "Rhat", "n.eff"
                                )]))
  } else {
    tau <- NA
  }

  # Return results
  results <- if (model == "RE") {
    list(EM = EM,
         tau = tau,
         measure = measure,
         scenarios = mean_scenarios,
         D = D,
         heter = heter_prior)
  } else {
    list(EM = EM,
         measure = measure,
         scenarios = mean_scenarios,
         D = D)
  }

  return(results)
}

#' Perform sensitivity meta-analysis for an aggregate binary or continuous outcome with missing participant data
#'
#' @param data A data-frame of a one-trial-per-row format containing arm-level data of each trial. This format is widely used for BUGS models. See 'Format' for the specification of the columns.
#' @param measure Character string indicating the effect measure with values \code{"OR"}, \code{"MD"}, \code{"SMD"}, or \code{"ROM"}.
#' @param assumption Character string indicating the structure of the informative missingness parameter. Set \code{assumption} equal to one of the following:  \code{"HIE-ARM"}, or \code{"IDE-ARM"}.
#' @param mean.misspar A positive non-zero number for the mean of the normal distribution of the informative missingness parameter.
#' @param var.misspar A positive non-zero number for the variance of the normal distribution of the informative missingness parameter.
#' @param D A binary number for the direction of the outcome. Set \code{D = 1} for a positive outcome and \code{D = 0} for a negative outcome.
#' @param n.chains Integer specifying the number of chains for the MCMC sampling; an argument of the \code{\link[R2jags]{jags}} function.
#' @param n.iter Integer specifying the number of Markov chains for the MCMC sampling; an argument of the \code{\link[R2jags]{jags}} function.
#' @param n.burnin Integer specifying the number of iterations to discard at the beginning of the MCMC sampling; an argument of the \code{\link[R2jags]{jags}} function.
#' @param n.thin Integer specifying the thinning rate for the MCMC sampling; an argument of the \code{\link[R2jags]{jags}} function.
#'
#' @return An R2jags output on the summaries of the posterior distribution, and the Gelman–Rubin convergence diagnostic of the following parameters:
#' \describe{
#'  \item{\code{EM}}{The effect estimate of all possible comparisons of interventions.}
#'  \item{\code{SUCRA}}{The surface under the cumulative ranking curve for each intervention.}
#'  \item{\code{phi}}{The informative missingness parameter.}
#'  \item{\code{delta}}{The underlying trial-specific effect estimate. For a multi-arm trial, we estimate \emph{T-1} trial-specific effect estimates, where \emph{T} is the number of interventions in the trial.}
#'  \item{\code{tau}}{The between-trial standard deviation assumed to be common for all observed comparisons.}
#' }
#'
#' @format The columns of the data frame \code{data} refer to the following ordered elements for a continuous outcome:
#' \describe{
#'  \item{\strong{t}}{An intervention identifier.}
#'  \item{\strong{y}}{The observed mean value of the outcome.}
#'  \item{\strong{sd}}{The observed standard deviation of the outcome.}
#'  \item{\strong{m}}{The number of missing outcome data.}
#'  \item{\strong{c}}{The number of participants completing the assigned intervention.}
#'  \item{\strong{na}}{The number of compared interventions.}
#' }
#' Apart from \strong{na}, all other elements appear in \code{data} as many times as the maximum number of interventions compared in a trial. See, 'Example'.
#'
#' @seealso \code{\link{R2jags}}
#'
#' @references
#' Gelman, A, Rubin, DB. Inference from iterative simulation using multiple sequences. Stat Sci. 1992;7:457–472.
#'
#' \dontshow{load("./data/NMA Dataset Continuous.RData")}
#' @examples
#' ### Show the data (one-trial-per-row format)
#' (data <- as.data.frame(one.stage.dataset.NMA[[3]]))
#'
#' ### Run a random-effects network meta-analysis with consistency equations for the standardised mean difference
#' ### assuming missing at random for identical, common informative missingness difference of means.
#' run.sensitivity(data = data, measure = "SMD", assumption = "IDE-COMMON", mean.misspar = 0, var.misspar = 1, D = 0, n.chains = 3, n.iter = 10000, n.burnin = 1000, n.thin = 1)
#'
#' @export
run.sensitivity <- function(data, measure, model, assumption, heter.prior, mean.misspar, var.misspar, D, n.chains, n.iter, n.burnin, n.thin){


  ## Prepare the dataset for the R2jags
  item <- data.preparation(data, measure)


  ## Default arguments
  model <- if (missing(model)) {
    "RE"
  } else if (!is.element(model, c("RE", "FE"))) {
    stop("Insert 'RE', or 'FE'")
  } else {
    model
  }
  assumption <- if (missing(assumption)) {
    "IDE-ARM"
  } else if (!is.element(assumption,  c("IDE-ARM", "IDE-TRIAL", "IDE-COMMON", "HIE-ARM", "HIE-TRIAL", "HIE-COMMON", "IND-CORR", "IND-UNCORR"))) {
    stop("Insert 'IDE-ARM', 'IDE-TRIAL', 'IDE-COMMON', 'HIE-ARM', 'HIE-TRIAL', 'HIE-COMMON', 'IND-CORR', or 'IND-UNCORR'")
  } else {
    assumption
  }
  mean.misspar <- missingness.param.prior(assumption, mean.misspar)
  heter.prior <- heterogeneity.param.prior(measure, model, heter.prior)
  var.misspar <- ifelse(missing(var.misspar) & (is.element(measure, c("OR", "MD", "SMD"))), 1, ifelse(missing(var.misspar) & measure == "ROM", 0.2^2, var.misspar))
  n.chains <- ifelse(missing(n.chains), 2, n.chains)
  n.iter <- ifelse(missing(n.iter), 10000, n.iter)
  n.burnin <- ifelse(missing(n.burnin), 1000, n.burnin)
  n.thin <- ifelse(missing(n.thin), 1, n.thin)


  ## Scenarios for missingness mechanism in an intervention (PMID: 30223064)
  if (is.element(measure, c("MD", "SMD", "ROM"))) {
    scenarios <- c(-2, -1, 0, 1, 2)
  } else {
    scenarios <- c(-log(3), -log(2), log(0.9999), log(2), log(3))
   }


  ## A 2x2 matrix of 25 reference-specific scenarios (PMID: 30223064)
  mean.misspar <- as.matrix(cbind(rep(scenarios, each = 5), rep(scenarios, 5))) # 2nd column refers to the reference intervention (control in MA)


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
                         model.file = textConnection(prepare.model(measure, model, assumption)),
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
    list(EM = EM, tau = tau, measure = measure)
  } else {
    list(EM = EM, measure = measure)
  }

  return(results)

}


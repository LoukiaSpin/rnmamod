#' Markov Chain Monte Carlo Diagnostics
#'
#' @param data A data-frame of a one-trial-per-row format containing arm-level data of each trial. This format is widely used for BUGS models. See 'Format' for the specification of the columns.
#' @param measure Character string indicating the effect measure with values \code{"OR"}, \code{"MD"}, \code{"SMD"}, or \code{"ROM"}.
#' @param assumption Character string indicating the structure of the informative missingness parameter. Set \code{assumption} equal to one of the following: \code{"HIE-COMMON"}, \code{"HIE-TRIAL"}, \code{"HIE-ARM"}, \code{"IDE-COMMON"}, \code{"IDE-TRIAL"}, \code{"IDE-ARM"}, \code{"IND-CORR"}, or \code{"IND-UNCORR"}.
#' @param heter.prior A vector of length equal to two with the following values: \code{rep(1, 2)}, \code{rep(2, 2)}, and \code{rep(3, 2)} refers to half-normal distribution with variance 1 or 0.5, and uniform distribution with interval [0, 5], respectively,
#' for the between-trial standard deviation. To indicate an empirically-based prior distribution for the between-trial variance, the first and second values of the vector should be the mean and precision
#' of the selected prior distribution. The empirically-based prior distribution for the between-trial variance is applicable only when \code{"OR"} or \code{"SMD"} is considered.
#' @param mean.misspar A positive non-zero number for the mean of the normal distribution of the informative missingness parameter.
#' @param var.misspar A positive non-zero number for the variance of the normal distribution of the informative missingness parameter.
#' @param D A binary number for the direction of the outcome. Set \code{D = 1} for a positive outcome and \code{D = 0} for a negative outcome.
#' @param n.chains Integer specifying the number of chains for the MCMC sampling; an argument of the \code{\link[R2jags]{jags}} function.
#' @param n.iter Integer specifying the number of Markov chains for the MCMC sampling; an argument of the \code{\link[R2jags]{jags}} function.
#' @param n.burnin Integer specifying the number of iterations to discard at the beginning of the MCMC sampling; an argument of the \code{\link[R2jags]{jags}} function.
#' @param n.thin Integer specifying the thinning rate for the MCMC sampling; an argument of the \code{\link[R2jags]{jags}} function.
#'
#' @return A panel of autocorrelation plots where the rows correspond to the chains and the columns correspond to the monitor parameters (maximum three).
#' Additionally, it uses the \code{\link[mcmcplots]{mcmcplot}} function to create an HTML file with a panel of diagnostic plots (trace, density, and autocorrelation) for each monitored parameter.
#'
#' @format See, function \code{nma.continuous.full.model}.
#'
#' @seealso \code{\link{mcmcplots}}
#'
#' @references
#' Gelman, A, Rubin, DB. Inference from iterative simulation using multiple sequences. Stat Sci. 1992;7:457–472.
#'
#' \dontshow{load("netmodr/data/One-stage model_NMA Dataset.RData")}
#' @examples
#' ### Obtain the diagnostic plots and check convergence for all monitored parameters using the R.hat
#' mcmc.diagnostics(par = c("tau2", "EM[3,1]", "EM[3,2]"), net = res1)
#'
#' @export
mcmc.diagnostics <- function(net, par){


  options(warn = -1)


  par <- if (missing(par)) {
    stop("The argument 'par' needs to be defined", call. = F)
  } else {
    par
  }


  jagsfit <- net$jagsfit

  ## Turn results into a data-frame to select model parameters (using 'dplyr')
  getResults <- as.data.frame(t(jagsfit$BUGSoutput$summary))

  # Effect size of all unique pairwise comparisons
  EM <- t(getResults %>% dplyr::select(starts_with("EM[")))

  # Predictive effects of all unique pairwise comparisons
  EM.pred <- t(getResults %>% dplyr::select(starts_with("EM.pred[")))


  # SUrface under the Cumulative RAnking curve values
  SUCRA <- t(getResults %>% dplyr::select(starts_with("SUCRA")))

  # Within-trial effects size (multi-arm trials with T interventions provide T-1 such effect sizes)
  delta <- t(getResults %>% dplyr::select(starts_with("delta") & !ends_with(",1]")))


  # Ranking probability of each intervention for every rank
  effectiveness <- t(getResults %>% dplyr::select(starts_with("effectiveness")))

  # Between-trial standard deviation
  tau <- t(getResults %>% dplyr::select(starts_with("tau")))


  # Estimated missingness parameter
  phi <- t(getResults %>% dplyr::select(starts_with("phi") | starts_with("mean.phi") | starts_with("mean.phi[") | starts_with("phi[")))

  # Regression coefficient for comparisons with the reference intervention
  beta <- t(getResults %>% dplyr::select(starts_with("beta[")))



  ## Turn 'R2jags' object into 'mcmc.plot' object
  jagsfit.mcmc <- as.mcmc(jagsfit)



  ## A panel of autocorrelation plots for each chain and every monitored parameter
  n.chains <- res1$jagsfit$BUGSoutput$n.chains
  autocorrelation <- par(mfrow = c(3, n.chains))
  for (i in 1:n.chains) {
    autplot1(jagsfit.mcmc[, par[1]], chain = i, main = paste(par[1], "-", "chain", i))
    autplot1(jagsfit.mcmc[, par[2]], chain = i, main = paste(par[2], "-","chain", i))
    autplot1(jagsfit.mcmc[, par[3]], chain = i, main = paste(par[3], "-","chain", i))
  }



  ## An HTML file with a panel of diagnostic plots per monitored paraemter
  mcmcplot <- mcmcplot(jagsfit.mcmc, parms = par)



  ## Keep results on the maximum Rhat for the selected monitored model parameters
  if(is.null(dim(phi))){

    R.hat.max <- append(R.hat.max, c(max(EM[, 8]), max(EM.pred[, 8]), max(delta[, 8]), max(tau[8]), max(SUCRA[, 8]), max(effectiveness[, 8]), phi[8], max(beta[, 8])))
    R.hat.max[c(2:4, 8)] <- ifelse(is.infinite(R.hat.max[c(2:4, 8)]), NA, R.hat.max[c(2:4, 8)])

  } else {

    R.hat.max <- c(max(EM[, 8]), max(EM.pred[, 8]), max(delta[, 8]), max(tau[8]), max(SUCRA[, 8]), max(effectiveness[, 8]), max(phi[, 8]), max(beta[, 8]))
    R.hat.max[c(2:4, 8)] <- ifelse(is.infinite(R.hat.max[c(2:4, 8)]), NA, R.hat.max[c(2:4, 8)])

  }



  ## Indicate whether each model parameter achieved or failed to achieve convergence
  conv <- rep(NA, length(R.hat.max))
  for(i in 1:length(R.hat.max)) {

    #conv[i] <- ifelse(!is.null(R.hat.max[i]) & R.hat.max[i] < 1.1, "achieved", ifelse(!is.null(R.hat.max[i]) & R.hat.max[i] >= 1.1,"failed", "Not applicable"))
    conv[i] <- ifelse(is.na(R.hat.max[i]), "Not applicable", ifelse(!is.na(R.hat.max[i]) & R.hat.max[i] < 1.1,"achieved", "failed"))

  }



  ## A data-frame with results on convergence for all monitored parameters using the Rhat
  convergence <- data.frame(R.hat.max, conv)
  rownames(convergence) <- c("EM", "Pred", "delta", "tau", "SUCRA", "effectiveness", "phi", "beta")
  colnames(convergence) <- c("R.hat max", "convergence status")


  return(list(convergence = convergence))
}


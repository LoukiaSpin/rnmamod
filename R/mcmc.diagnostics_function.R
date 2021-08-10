#' Markov Chain Monte Carlo Diagnostics
#'
#' @description This function evaluates whether convergence has been achieved for the monitored parameters of the Bayesian models.
#'   The Gelman-Rubin convergence diagnostic and relevant diagnostic plots are applied for that purpose.
#'
#' @param net An object of S3 class \code{\link{run.model}} and \code{\link{run.metareg}}. See 'Value' in \code{\link{run.model}}.
#' @param par A vector of three character strings that refer to three monitored parameters in \code{jagsfit} which is an object of S3 class \code{\link{run.model}} and \code{\link{run.metareg}}.
#'   These three selected parameters will be considered in the diagnostic plots (see 'Value').
#'
#' @return This function returns a data-frame that contains the Gelman-Rubin convergence diagnostic, R-hat, and convergence status of the following monitored parameters:
#' \tabular{ll}{
#'  \code{EM} \tab The estimated summary effect measure.\cr
#'  \tab \cr
#'  \code{EM.pred} \tab The predicted summary effect measure.\cr
#'  \tab \cr
#'  \code{delta} \tab The estimated trial-specific effect measure.\cr
#'  \tab \cr
#'  \code{effectiveneness} \tab The ranking probability of each intervention for every rank.\cr
#'  \tab \cr
#'  \code{phi} \tab The informative missingness parameter.\cr
#'  \tab \cr
#'  \code{beta} \tab The regression coefficient.\cr
#' }
#' \code{mcmc.diagnostics} also uses the \code{\link[mcmcplots]{mcmcplot}} function to create an HTML file with a panel of diagnostic plots (trace, density, and autocorrelation) for each monitored parameter.
#'
#' @details For each monitored parameter, \code{mcmc.diagnostics} considers the maximum R-hat and compares it with the threshold 1.1: convergence is achieved for the monitored parameter, when the maximum R-hat
#'   is below that threshold; otherwise, the Markov Chain Monte Carlo algorithm has not converged for that parameter. If the monitored parameter is a vector with the posterior results, there is only one R-hat.
#'   If the monitored parameter is a matrix of the posterior results, there are as many R-hats as the number of rows for that parameter.
#'
#' @author {Loukia M. Spineli}
#'
#' @seealso \code{\link[mcmcplots]{mcmcplot}}, \code{\link{run.model}}, \code{\link{run.metareg}}
#'
#' @references
#' Gelman, A, Rubin, DB. Inference from iterative simulation using multiple sequences. \emph{Stat Sci} 1992;\bold{7}:457--472. [\doi{10.1214/ss/1177011136}]
#'
#' @examples
#' data("nma.liu2013")
#'
#' \dontrun{
#' # Perform a random-effects network meta-analysis
#' res <- run.model(data = nma.liu2013,
#'                  measure = "OR",
#'                  model = "RE",
#'                  assumption = "IDE-ARM",
#'                  heter.prior = list("halfnormal", 0, 1),
#'                  mean.misspar = 0,
#'                  var.misspar = 1,
#'                  D = 1,
#'                  n.chains = 3,
#'                  n.iter = 10000,
#'                  n.burnin = 1000,
#'                  n.thin = 1)
#'
#' # Obtain the diagnostic plots and check convergence for all monitored parameters using the R.hat
#' mcmc.diagnostics(net = res, par = c("tau", "EM[2,1]", "EM.pred[2,1]"))
#' }
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

  item <- data.preparation(net$data, net$measure)

  ## Turn results into a data-frame to select model parameters (using 'dplyr')
  getResults <- as.data.frame(t(jagsfit$BUGSoutput$summary))

  # Effect size of all unique pairwise comparisons
  EM <- t(getResults %>% select(starts_with("EM[")))

  # Predictive effects of all unique pairwise comparisons
  EM.pred <- t(getResults %>% select(starts_with("EM.pred[")))


  # SUrface under the Cumulative RAnking curve values
  SUCRA <- t(getResults %>% select(starts_with("SUCRA")))


  # Within-trial effects size (multi-arm trials with T interventions provide T-1 such effect sizes)
  delta <- t(getResults %>% select(starts_with("delta") & !ends_with(",1]")))


  # Ranking probability of each intervention for every rank
  effectiveness <- t(getResults %>% select(starts_with("effectiveness")))

  # Between-trial standard deviation
  tau <- t(getResults %>% select(starts_with("tau")))


  # Estimated missingness parameter
  #phi <- net$phi
  phi <- if (length(unique(unlist(item$m))) > 2) {
    t(getResults %>% select(starts_with("phi") | starts_with("mean.phi") | starts_with("mean.phi[") | starts_with("phi[")))
  } else {
    NA
  }

  # Regression coefficient
  beta <- t(getResults %>% select(starts_with("beta[") | starts_with("beta")))


  ## Turn 'R2jags' object into 'mcmc' object
  jagsfit.mcmc <- as.mcmc(jagsfit)


  ## A panel of autocorrelation plots for each chain and every monitored parameter
  #n.chains <- res1$jagsfit$BUGSoutput$n.chains
  #autocorrelation <- par(mfrow = c(3, n.chains))
  #for (i in 1:n.chains) {
  #  autplot1(jagsfit.mcmc[, par[1]], chain = i, main = paste(par[1], "-", "chain", i))
  #  autplot1(jagsfit.mcmc[, par[2]], chain = i, main = paste(par[2], "-","chain", i))
  #  autplot1(jagsfit.mcmc[, par[3]], chain = i, main = paste(par[3], "-","chain", i))
  #}


  ## An HTML file with a panel of diagnostic plots per monitored parameter
  mcmcplot <- mcmcplot(jagsfit.mcmc, parms = par)


  ## Keep results on the maximum Rhat for the selected monitored model parameters
  phi.R.hat.max <- if (is.null(dim(phi))) {
    phi[8]
  } else {
    max(phi[, 8])
  }

  beta.R.hat.max <- if (is.null(dim(beta))) {
    beta[8]
  } else {
    max(beta[, 8])
  }


  R.hat.max <- c(max(EM[, 8]), max(EM.pred[, 8]), max(delta[, 8]), tau[8], max(SUCRA[, 8]), max(effectiveness[, 8]), phi.R.hat.max, beta.R.hat.max)
  for (i in 1:length(R.hat.max)) {
    R.hat.max[i] <- ifelse(is.infinite(R.hat.max[i]), NA, R.hat.max[i])
  }



  #if (is.null(dim(phi))) {
  #  R.hat.max <- c(max(EM[, 8]), max(EM.pred[, 8]), max(delta[, 8]), tau[8], max(SUCRA[, 8]), max(effectiveness[, 8]), max(phi[8]), max(beta[, 8]))
    #R.hat.max[c(2:4, 8)] <- ifelse(is.infinite(R.hat.max[c(2:4, 8)]), NA, R.hat.max[c(2:4, 8)])
    #R.hat.max <- ifelse(is.infinite(R.hat.max), NA, R.hat.max)
  #  for (i in 1:length(R.hat.max)) {
  #    R.hat.max[i] <- ifelse(is.infinite(R.hat.max[i]), NA, R.hat.max[i])
  #  }
  #} else {
  #  R.hat.max <- c(max(EM[, 8]), max(EM.pred[, 8]), max(delta[, 8]), tau[8], max(SUCRA[, 8]), max(effectiveness[, 8]), max(phi[, 8]), max(beta[, 8]))
  #  #R.hat.max[c(2:4, 8)] <- ifelse(is.infinite(R.hat.max[c(2:4, 8)]), NA, R.hat.max[c(2:4, 8)])
  #  #R.hat.max <- ifelse(is.infinite(R.hat.max), NA, R.hat.max)
  #  for (i in 1:length(R.hat.max)) {
  #    R.hat.max[i] <- ifelse(is.infinite(R.hat.max[i]), NA, R.hat.max[i])
  #  }
  #}



  ## Indicate whether each model parameter achieved or failed to achieve convergence
  conv <- rep(NA, length(R.hat.max))
  for (i in 1:length(R.hat.max)) {
    conv[i] <- ifelse(is.na(R.hat.max[i]), "Not applicable", ifelse(!is.na(R.hat.max[i]) & R.hat.max[i] < 1.1, "achieved", "failed"))
  }


  ## A data-frame with results on convergence for all monitored parameters using the Rhat
  convergence <- data.frame(R.hat.max, conv)
  rownames(convergence) <- c("EM", "Pred", "delta", "tau", "SUCRA", "effectiveness", "phi", "beta")
  colnames(convergence) <- c("R.hat max", "convergence status")


  return(list(convergence = convergence))
}


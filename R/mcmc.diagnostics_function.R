#' Markov Chain Monte Carlo Diagnostics
#'
#' @description This function evaluates whether convergence has been achieved for the monitored parameters of the Bayesian models.
#'   The Gelman-Rubin convergence diagnostic and relevant diagnostic plots are applied for that purpose.
#'
#' @param net An object of S3 class \code{\link{run.model}}, \code{\link{run.series.meta}}, \code{\link{run.nodesplit}}, \code{\link{run.UME}},
#'   \code{\link{run.sensitivity}} and \code{\link{run.metareg}}. See 'Value' in the functions above.
#' @param par A vector of three character strings that refer to three monitored parameters in \code{jagsfit} which is an object of S3 class
#'   \code{\link{run.model}}, \code{\link{run.UME}} and \code{\link{run.metareg}}.
#'   These three selected parameters will be considered in the diagnostic plots (see 'Value').
#'   This argument will be ignored for objects of S3 class \code{\link{run.series.meta}}, \code{\link{run.nodesplit}}, or \code{\link{run.sensitivity}}.
#'
#' @return This function returns a data-frame that contains the Gelman-Rubin convergence diagnostic, R-hat, and convergence status of the following monitored parameters:
#' \tabular{ll}{
#'  \code{EM} \tab The estimated summary effect measure.\cr
#'  \tab \cr
#'  \code{EM.pred} \tab The predicted summary effect measure.\cr
#'  \tab \cr
#'  \code{delta} \tab The estimated trial-specific effect measure.\cr
#'  \tab \cr
#'  \code{tau} \tab The between-trial standard deviation.\cr
#'  \tab \cr
#'  \code{direct} \tab The direct estimate of the split node (see 'Value' \code{\link{run.nodesplit}}).\cr
#'  \tab \cr
#'  \code{indirect} \tab The indirect estimate of the split node (see 'Value' \code{\link{run.nodesplit}}).\cr
#'  \tab \cr
#'  \code{IF} \tab The inconsistency factor of the split node (see 'Value' \code{\link{run.nodesplit}}).\cr
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
#' @seealso \code{\link[mcmcplots]{mcmcplot}}, \code{\link{run.model}}, \code{\link{run.series.meta}}, \code{\link{run.nodesplit}}, \code{\link{run.UME}}, \code{\link{run.UME}}, \code{\link{run.metareg}}
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
#'                  mean.misspar = c(0, 0),
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

  par <- if (!is.null(net$jagsfit) & missing(par)) {
    stop("The argument 'par' needs to be defined", call. = F)
  } else if (is.null(net$jagsfit)) {
    message(cat(paste0("\033[0;", col = 32, "m", txt = "The argument 'par' is ignored as it is used only for the functions 'run.model', 'run.ume' and 'run.metareg'", "\033[0m", "\n")))
    NULL
  } else {
    par
  }


  if (!is.null(net$jagsfit)) {

    jagsfit <- net$jagsfit

    ## Turn results into a data-frame to select model parameters (using 'dplyr')
    getResults <- as.data.frame(t(jagsfit$BUGSoutput$summary))

    # Effect size of all unique pairwise comparisons
    EM0 <- t(getResults %>% select(starts_with("EM[")))
    EM <- max(EM0[, 8])

    # Predictive effects of all unique pairwise comparisons
    EM.pred <- t(getResults %>% select(starts_with("EM.pred[")))

    # Within-trial effects size (multi-arm trials with T interventions provide T-1 such effect sizes)
    delta <- t(getResults %>% select(starts_with("delta") & !ends_with(",1]")))

    # Between-trial standard deviation
    tau0 <- t(getResults %>% select(starts_with("tau")))
    tau <- tau0[8]

    # Direct estimate from split nodes
    direct <- NULL

    # Indirect estimate from split nodes
    indirect <- NULL

    # Inconsistency factor estimate from split nodes
    diff <- NULL

    item <- data.preparation(net$data, net$measure)

    # Estimated missingness parameter
    phi <- if (length(unique(unlist(item$m))) > 2) {
      t(getResults %>% select(starts_with("phi") | starts_with("mean.phi") | starts_with("mean.phi[") | starts_with("phi[")))
    } else {
      NA
    }

    # Regression coefficient
    beta <- t(getResults %>% select(starts_with("beta[") | starts_with("beta")))

  } else {

    # Effect size of pairwise comparisons with at least two trials
    EM <- if(length(net$EM[1, ]) == 11) {
      #... for each pairwise comparison with at least two trials
      max(net$EM[, 10])
    } else {
      # ... for split node
      max(net$EM[, 5])
    }

    # Predictive effects of all unique pairwise comparisons
    EM.pred <- NULL

    # Within-trial effects size (multi-arm trials with T interventions provide T-1 such effect sizes)
    delta <- NULL

    # Between-trial standard deviation ...
    tau <- if (length(net$tau[1, ]) == 11) {
      #... for each pairwise comparison with at least two trials
      max(net$tau[, 10])
    } else if (length(net$tau[1, ]) == 8) {
      # ... for split node
      max(net$tau[, 7])
    } else if (length(net$tau[1, ]) == 6) {
      # ... sensitivity analysis to different scenarios about the missingness parameter
      max(net$tau[, 5])
    }

    # Direct estimate from split nodes
    direct <- net$direct[, 7]

    # Indirect estimate from split nodes
    indirect <- net$indirect[, 7]

    # Inconsistency factor estimate from split nodes
    diff <- net$diff[, 7]

    # Estimated missingness parameter
    phi <- Inf

    # Regression coefficient
    beta <- Inf
  }


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


  if (!is.null(net$jagsfit)) {
    ## Turn 'R2jags' object into 'mcmc' object
    jagsfit.mcmc <- as.mcmc(jagsfit)

    ## An HTML file with a panel of diagnostic plots per monitored parameter
    mcmcplot <- mcmcplot(jagsfit.mcmc, parms = par)
  }


  R.hat.max <- c(EM, max(EM.pred[, 8]), max(delta[, 8]), tau, max(direct), max(indirect), max(diff), phi.R.hat.max, beta.R.hat.max)
  for (i in 1:length(R.hat.max)) {
    R.hat.max[i] <- ifelse(is.infinite(R.hat.max[i]), NA, R.hat.max[i])
  }


  ## Indicate whether each model parameter achieved or failed to achieve convergence
  conv <- rep(NA, length(R.hat.max))
  for (i in 1:length(R.hat.max)) {
    conv[i] <- ifelse(is.na(R.hat.max[i]), "Not applicable", ifelse(!is.na(R.hat.max[i]) & R.hat.max[i] < 1.1, "achieved", "failed"))
  }


  ## A data-frame with results on convergence for all monitored parameters using the Rhat
  convergence <- data.frame(R.hat.max, conv)
  rownames(convergence) <- c("Effect estimates (EM)",
                             "Predictions (EM.pred)",
                             "Within-trial estimates (delta)",
                             "Between-trial standard deviation (tau)",
                             "Direct effects (node-splitting; direct)",
                             "Indirect effect(s) (node-splitting; indirect)",
                             "Inconsistency factor(s) (node-splitting; IF)",
                             "Informative missingness parameter(s) (phi)",
                             "Regression coefficient(s) (beta)")
  colnames(convergence) <- c("R.hat max", "convergence status")


  return(list(convergence = convergence))
}


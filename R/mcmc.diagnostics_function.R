#' Markov Chain Monte Carlo diagnostics
#'
#' @description Evaluates whether convergence has been achieved for the
#'   monitored parameters of the Bayesian models. The Gelman-Rubin convergence
#'   diagnostic and relevant diagnostic plots are applied.
#'
#' @param net An object of S3 class \code{\link{run_model}},
#'   \code{\link{run_series_meta}}, \code{\link{run_nodesplit}},
#'   \code{\link{run_ume}}, \code{\link{run_sensitivity}} and
#'   \code{\link{run_metareg}}. See 'Value' in the functions above.
#' @param par A vector of three character strings that refer to three monitored
#'   parameters in \code{jagsfit} which is an object of S3 class
#'   \code{\link{run_model}}, \code{\link{run_ume}} and
#'   \code{\link{run_metareg}}. These three selected parameters will be
#'   considered in the diagnostic plots (see 'Value'). This argument will be
#'   ignored for objects of S3 class \code{\link{run_series_meta}},
#'   \code{\link{run_nodesplit}}, or \code{\link{run_sensitivity}}.
#'
#' @return \code{mcmc_diagnostics} returns a data-frame that contains the
#'   Gelman-Rubin convergence diagnostic, R-hat, and convergence status of the
#'   following monitored parameters:
#'   \tabular{ll}{
#'    \code{EM} \tab The estimated summary effect measure.\cr
#'    \tab \cr
#'    \code{EM_pred} \tab The predicted summary effect measure.\cr
#'    \tab \cr
#'    \code{delta} \tab The estimated trial-specific effect measure.\cr
#'    \tab \cr
#'    \code{tau} \tab The between-trial standard deviation.\cr
#'    \tab \cr
#'    \code{direct} \tab The direct estimate of the split node
#'    (see 'Value' in \code{\link{run_nodesplit}}).\cr
#'    \tab \cr
#'    \code{indirect} \tab The indirect estimate of the split node
#'    (see 'Value' in \code{\link{run_nodesplit}}).\cr
#'    \tab \cr
#'    \code{IF} \tab The inconsistency factor of the split node
#'    (see 'Value' in \code{\link{run_nodesplit}}).\cr
#'    \tab \cr
#'    \code{phi} \tab The informative missingness parameter.\cr
#'    \tab \cr
#'    \code{beta} \tab The regression coefficient.\cr
#'   }
#'   \code{mcmc_diagnostics} also uses the \code{\link[mcmcplots]{mcmcplot}}
#'   function to create an HTML file with a panel of diagnostic plots
#'   (trace, density, and autocorrelation) for each monitored parameter.
#'
#' @details For each monitored parameter, \code{mcmc_diagnostics} considers the
#'   maximum R-hat and compares it with the threshold 1.1: convergence is
#'   achieved for the monitored parameter, when the maximum R-hat is below that
#'   threshold; otherwise, the Markov Chain Monte Carlo algorithm has not
#'   converged for that parameter. If the monitored parameter is a vector with
#'   the posterior results, there is only one R-hat. If the monitored parameter
#'   is a matrix of the posterior results, there are as many R-hats as the
#'   number of rows for that parameter.
#'
#' @author {Loukia M. Spineli}
#'
#' @seealso \code{\link[mcmcplots]{mcmcplot}}, \code{\link{run_metareg}},
#'   \code{\link{run_model}}, \code{\link{run_nodesplit}},
#'   \code{\link{run_sensitivity}}, \code{\link{run_series_meta}},
#'   \code{\link{run_ume}}
#'
#' @references
#' Gelman, A, Rubin, DB. Inference from iterative simulation using multiple
#' sequences. \emph{Stat Sci} 1992;\bold{7}:457--472.
#'
#' @examples
#' data("nma.liu2013")
#'
#' \dontrun{
#' # Perform a random-effects network meta-analysis
#' res <- run_model(data = nma.liu2013,
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
#' # Obtain the diagnostic plots and check convergence based on R-hat
#' mcmc_diagnostics(net = res,
#'                  par = c("tau", "EM[2,1]", "EM.pred[2,1]"))
#' }
#' @export
mcmc_diagnostics <- function(net, par) {

  options(warn = -1)

  par <- if (!is.null(net$jagsfit) & missing(par)) {
    stop("The argument 'par' needs to be defined", call. = FALSE)
  } else if (is.null(net$jagsfit)) {
    aa <- "The argument 'par' is ignored as it is used only"
    bb <- "for the functions 'run_model', 'run_ume' and 'run_metareg'"
    message(cat(paste0("\033[0;", col = 32, "m", aa, bb, "\033[0m", "\n")))
    NULL
  } else {
    par
  }

  if (!is.null(net$jagsfit)) {
    jagsfit <- net$jagsfit
    # Turn results into a data-frame to select model parameters (using 'dplyr')
    get_results <- as.data.frame(t(jagsfit$BUGSoutput$summary))

    # Effect size of all unique pairwise comparisons
    EM0 <- t(get_results %>% select(starts_with("EM[")))
    EM <- max(EM0[, 8])

    # Predictive effects of all unique pairwise comparisons
    EM_pred <- t(get_results %>% select(starts_with("EM_pred[")))

    # Within-trial effects size
    delta <- t(get_results %>% select(starts_with("delta") & !ends_with(",1]")))

    # Between-trial standard deviation
    tau0 <- t(get_results %>% select(starts_with("tau")))
    tau <- tau0[8]

    # Direct estimate from split nodes
    direct <- NULL

    # Indirect estimate from split nodes
    indirect <- NULL

    # Inconsistency factor estimate from split nodes
    diff <- NULL

    item <- data_preparation(net$data, net$measure)

    # Estimated missingness parameter
    phi <- if (length(unique(unlist(item$m))) > 2) {
      t(get_results %>% select(starts_with("phi") |
                                starts_with("mean.phi") |
                                starts_with("mean.phi[") |
                                starts_with("phi[")))
    } else {
      NA
    }

    # Regression coefficient
    beta <- t(get_results %>% select(starts_with("beta[") |
                                       starts_with("beta")))
  } else {
    # Effect size of pairwise comparisons with at least two trials
    EM <- if (length(net$EM[1, ]) == 11) {
      #... for each pairwise comparison with at least two trials
      max(net$EM[, 10])
    } else {
      # ... for split node
      max(net$EM[, 5])
    }

    # Predictive effects of all unique pairwise comparisons
    EM_pred <- NULL

    # Within-trial effects size
    delta <- NULL

    # Between-trial standard deviation ...
    tau <- if (length(net$tau[1, ]) == 11) {
      #... for each pairwise comparison with at least two trials
      max(net$tau[, 10])
    } else if (length(net$tau[1, ]) == 8) {
      # ... for split node
      max(net$tau[, 7])
    } else if (length(net$tau[1, ]) == 6) {
      # ... sensitivity analysis to different scenarios on missingness parameter
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

  phi_r_hat_max <- if (is.null(dim(phi))) {
    phi[8]
  } else {
    max(phi[, 8])
  }

  beta_r_hat_max <- if (is.null(dim(beta))) {
    beta[8]
  } else {
    max(beta[, 8])
  }

  if (!is.null(net$jagsfit)) {
    # Turn 'R2jags' object into 'mcmc' object
    jagsfit_mcmc <- as.mcmc(jagsfit)

    # An HTML file with a panel of diagnostic plots per monitored parameter
    mcmcplot <- mcmcplot(jagsfit_mcmc, parms = par)
  }

  r_hat_max <- c(EM,
                 max(EM_pred[, 8]),
                 max(delta[, 8]),
                 tau,
                 max(direct),
                 max(indirect),
                 max(diff),
                 phi_r_hat_max,
                 beta_r_hat_max)
  for (i in seq_len(length(r_hat_max))) {
    r_hat_max[i] <- ifelse(is.infinite(r_hat_max[i]), NA, r_hat_max[i])
  }

  # Indicate whether each model parameter achieved or failed to converge
  conv <- rep(NA, length(r_hat_max))
  for (i in seq_len(length(r_hat_max))) {
    conv[i] <- ifelse(is.na(r_hat_max[i]), "Not applicable",
                      ifelse(!is.na(r_hat_max[i]) & r_hat_max[i] < 1.1,
                             "achieved", "failed"))
  }

  # A data-frame on convergence for all monitored parameters using the Rhat
  convergence <- data.frame(r_hat_max, conv)
  rownames(convergence) <- c("Effect estimates (EM)",
                             "Predictions (EM_pred)",
                             "Within-trial estimates (delta)",
                             "Between-trial standard deviation (tau)",
                             "Direct effects (node-splitting; direct)",
                             "Indirect effect(s) (node-splitting; indirect)",
                             "Inconsistency factor(s) (node-splitting; diff)",
                             "Informative missingness parameter(s) (phi)",
                             "Regression coefficient(s) (beta)")
  colnames(convergence) <- c("R.hat max", "convergence status")

  return(list(convergence = convergence))
}

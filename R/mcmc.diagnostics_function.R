#' Markov Chain Monte Carlo diagnostics
#'
#' @description Evaluates whether convergence has been achieved for the
#'   monitored parameters of the Bayesian models. The Gelman-Rubin convergence
#'   diagnostic, the Markov Chain Monte Carl (MCMC) error and relevant
#'   diagnostic plots are applied.
#'
#' @param net An object of S3 class \code{\link{run_metareg}},
#'   \code{\link{run_model}}, \code{\link{run_nodesplit}},
#'   \code{\link{run_sensitivity}}, \code{\link{run_series_meta}}, and
#'   \code{\link{run_ume}}. See 'Value' in the functions above.
#' @param par A vector of three character strings that refer to three monitored
#'   parameters in \code{jagsfit} which is an object of S3 class
#'   \code{\link{run_metareg}}, \code{\link{run_model}}, and
#'   \code{\link{run_ume}}. These three selected parameters will be
#'   considered in the diagnostic plots (see 'Value'). This argument will be
#'   ignored for objects of S3 class \code{\link{run_nodesplit}},
#'   \code{\link{run_sensitivity}}, and \code{\link{run_series_meta}}.
#'
#' @return \code{mcmc_diagnostics} returns a data-frame that contains the
#'   Gelman-Rubin convergence diagnostic, \strong{R-hat}, the MCMC error, and
#'   the convergence status of the following monitored parameters:
#'   \code{EM} {The estimated summary effect measure.}
#'   \code{EM_pred} {The predicted summary effect measure.}
#'   \code{delta} {The estimated trial-specific effect measure.}
#'   \code{tau} {The between-trial standard deviation.}
#'   \code{direct} {The direct estimate of the split node (see 'Value' in
#'   \code{\link{run_nodesplit}}).}
#'   \code{indirect} {The indirect estimate of the split node
#'   (see 'Value' in \code{\link{run_nodesplit}}).}
#'   \code{diff} {The inconsistency factor of the split node (see 'Value' in
#'   \code{\link{run_nodesplit}}).}
#'   \code{phi} {The informative missingness parameter.}
#'   \code{beta} {The regression coefficient.}
#'
#'   \code{mcmc_diagnostics} also uses the
#'   \code{\link[mcmcplots:mcmcplot]{mcmcplot}} function of the R-package
#'   \href{https://CRAN.R-project.org/package=mcmcplots}{mcmcplots} to create an
#'   HTML file with a panel of diagnostic plots (trace, density, and
#'   autocorrelation) for each monitored parameter.
#'
#' @details For each monitored parameter, \code{mcmc_diagnostics} considers the
#'   maximum R-hat and MCMC error and compares them with the thresholds 1.1
#'   and 5\% of the posterior standard deviation (the rule of thumb),
#'   respectively.
#'   Convergence is achieved for the monitored parameter,
#'   when the maximum R-hat and MCMC error are below the corresponding
#'   thresholds; otherwise, the MCMC algorithm has not converged for that
#'   parameter. If the monitored parameter is a vector with the posterior
#'   results, there is only one R-hat and one MCMC error. If the monitored
#'   parameter is a matrix of the posterior results, there are as many R-hats
#'   and MCMC errors as the number of rows for that parameter. In that case, the
#'   maximum R-hat and MCMC error are considered, and the (MCMC) rule of thumb
#'   corresponds to the maximum of the 5\% of the posterior standard deviation
#'   of the investigated parameter.
#'
#' @author {Loukia M. Spineli}
#'
#' @seealso \code{\link[mcmcplots:mcmcplot]{mcmcplot}},
#'   \code{\link{run_metareg}}, \code{\link{run_model}},
#'   \code{\link{run_nodesplit}}, \code{\link{run_sensitivity}},
#'   \code{\link{run_series_meta}}, \code{\link{run_ume}}
#'
#' @references
#' Gelman, A, Rubin, DB. Inference from iterative simulation using multiple
#' sequences. \emph{Stat Sci} 1992;\bold{7}(4):457--72.
#' doi: 10.1214/ss/1177011136
#'
#' @examples
#' data("nma.baker2009")
#'
#' # Read results from 'run_nodesplit' (using the default arguments)
#' res <- readRDS(system.file('extdata/node_baker.rds', package = 'rnmamod'))
#'
#' # Check convergence based on R-hat
#' mcmc_diagnostics(net = res,
#'                  par = c("tau", "EM[2,1]", "EM.pred[2,1]"))
#'
#' @export
mcmc_diagnostics <- function(net, par = NULL) {

  if (!is.element(net$type,
                  c("nma", "nmr", "series", "node", "ume", "sens")) ||
      is.null(net$type)) {
    aa <- "'run_model', 'run_metareg', 'run_series_meta',"
    bb <- "'run_nodesplit', 'run_ume', or 'run_sensitivity'."
    stop(paste("'net' must be an object of S3 class", aa, bb), call. = FALSE)
  }
  cc <- "< 5% of the posterior standard deviation."
  message(paste("A parameter converges when R-hat < 1.10 *and* MCMC error", cc))

  par <- if (!is.null(net$jagsfit) & missing(par)) {
    stop("The argument 'par' needs to be defined.", call. = FALSE)
  } else if (!is.null(net$jagsfit) & !is.null(par)) {
    par
  } else if (is.null(net$jagsfit) & !is.null(par)) {
    aa <- "Note: The argument 'par' is ignored. It is used only"
    bb <- "with 'run_model', 'run_ume' and 'run_metareg'."
    message(paste(aa, bb))
    NULL
  }

  #save_res <- ((net$n_iter - net$n_burnin) * net$n_chains) / net$n_thin

  if (!is.null(net$jagsfit)) {
    jagsfit <- net$jagsfit
    # Turn results into a data-frame to select model parameters (using 'dplyr')
    get_results <- as.data.frame(t(jagsfit$BUGSoutput$summary))

    # Effect size of all unique pairwise comparisons
    EM0 <- t(get_results %>% select(starts_with("EM[")))
    EM <- max(EM0[, 8])
    EM_mcmc_error <- max(EM0[, 2]/sqrt(EM0[, 9]))
    EM_mcmc_rule <- max(0.05 * EM0[, 2])

    # Predictive effects of all unique pairwise comparisons
    if (net$model == "RE" & !is.null(net$EM_pred)) {
      EM_pred0 <- t(get_results %>% select(starts_with("EM.pred[")))
      EM_pred <- max(EM_pred0[, 8])
      EM_pred_mcmc_error <- max(EM_pred0[, 2]/sqrt(EM_pred0[, 9]))
      EM_pred_mcmc_rule <- max(0.05 * EM_pred0[, 2])
    } else if (net$model == "FE" || is.null(net$EM_pred)) {
      EM_pred <- EM_pred_mcmc_error <- EM_pred_mcmc_rule <- NA
    }

    # Within-trial effects size
    if (net$model == "RE" & !is.null(net$delta)) {
      delta0 <- t(get_results %>% select(starts_with("delta") &
                                           !ends_with(",1]")))
      delta <- max(delta0[, 8])
      delta_mcmc_error <- max(delta0[, 2]/sqrt(delta0[, 9]))
      delta_mcmc_rule <- max(0.05 * delta0[, 2])
    } else if (net$model == "FE" || is.null(net$delta)) {
      delta <- delta_mcmc_error <- delta_mcmc_rule <- NA
    }

    # Between-trial standard deviation
    if (net$model == "RE") {
      tau0 <- t(get_results %>% select(starts_with("tau")))
      tau <- tau0[8]
      tau_mcmc_error <- tau0[2]/sqrt(tau0[9])
      tau_mcmc_rule <- max(0.05 * tau0[2])
    } else {
      tau <- tau_mcmc_error <- tau_mcmc_rule <- NA
    }

    # Direct estimate from split nodes
    direct <- direct_mcmc_error <- direct_mcmc_rule <- NA

    # Indirect estimate from split nodes
    indirect <- indirect_mcmc_error <- indirect_mcmc_rule <- NA

    # Inconsistency factor estimate from split nodes
    diff <- diff_mcmc_error <- diff_mcmc_rule <- NA

    item <- data_preparation(net$data, net$measure)

    # Estimated missingness parameter
    if (!is.null(net$phi) & is.element(net$assumption,
                                       c("IDE-COMMON", "HIE-COMMON"))) {
      phi0 <- t(get_results %>% select(starts_with("phi") |
                                            starts_with("mean.phi")))
      phi <- phi0[8]
      phi_mcmc_error <- phi0[2]/sqrt(phi0[9])
      phi_mcmc_rule <- max(0.05 * phi0[2])
    } else if (!is.null(net$phi) & !is.element(net$assumption,
                                               c("IDE-COMMON", "HIE-COMMON"))) {
      phi0 <- t(get_results %>% select(starts_with("mean.phi[") |
                                         starts_with("phi[")))
      phi <- max(phi0[, 8])
      phi_mcmc_error <- max(phi0[, 2]/sqrt(phi0[, 9]))
      phi_mcmc_rule <- max(0.05 * phi0[, 2])
    } else if (is.null(net$phi)) {
      phi <- phi_mcmc_error <- phi_mcmc_rule <- NA
    }

    # Regression coefficient
    if (!is.null(net$beta_all)) {
      beta0 <- t(get_results %>% select(starts_with("beta.all[")))
      beta <- max(beta0[, 8])
      beta_mcmc_error <- max(beta0[, 2]/sqrt(beta0[, 9]))
      beta_mcmc_rule <- max(0.05 * beta0[, 2])
    } else if (is.null(net$beta)) {
      beta <- beta_mcmc_error <- beta_mcmc_rule <- NA
    }

  } else {
    if (!is.null(net$EM) & length(net$EM[1, ]) == 11) {
      # From 'run_model' function
      EM_pred <- EM_pred_mcmc_error <- EM_pred_mcmc_rule <- NA
      delta <- delta_mcmc_error <- delta_mcmc_rule <- NA
      phi <- phi_mcmc_error <- phi_mcmc_rule <- NA

      # From 'run_metareg' function
      beta <- beta_mcmc_error <- beta_mcmc_rule <- NA

      # From 'run_series_meta' function
      EM <- max(net$EM[, 10])
      EM_mcmc_error <- max(net$EM[, 4]/sqrt(net$EM[, 11]))
      EM_mcmc_rule <- max(0.05 * net$EM[, 4])
      if (!is.null(net$tau)) {
        tau <- max(net$tau[, 10])
        tau_mcmc_error <- max(net$tau[, 4]/sqrt(net$tau[, 11]))
        tau_mcmc_rule <- max(0.05 * net$tau[, 4])
      } else {
        tau <- tau_mcmc_error <- tau_mcmc_rule <- NA
      }

      # From 'run_nodesplit' function
      direct <- direct_mcmc_error <- direct_mcmc_rule <- NA
      indirect <- indirect_mcmc_error <- indirect_mcmc_rule <- NA
      diff <- diff_mcmc_error <- diff_mcmc_rule <- NA
    } else if (is.null(net$EM)) {
      # From 'run_model' function
      EM <- EM_mcmc_error <- EM_mcmc_rule <- NA
      EM_pred <- EM_pred_mcmc_error <- EM_pred_mcmc_rule <- NA
      delta <- delta_mcmc_error <- delta_mcmc_rule <- NA
      phi <- phi_mcmc_error <- phi_mcmc_rule <- NA

      # From 'run_metareg' function
      beta <- beta_mcmc_error <- beta_mcmc_rule <- NA

      # From 'run_nodesplit' function
      if (!is.null(net$tau)) {
        tau <- max(net$tau[, 7])
        tau_mcmc_error <- max(net$tau[, 4]/sqrt(net$tau[, 8]))
        tau_mcmc_rule <- max(0.05 * net$tau[, 4])
      } else {
        tau <- tau_mcmc_error <- tau_mcmc_rule <- NA
      }
      direct <- max(net$direct[, 7])
      direct_mcmc_error <- max(net$direct[, 4]/sqrt(net$direct[, 8]))
      direct_mcmc_rule <- max(0.05 * net$direct[, 4])
      indirect <- max(net$indirect[, 7])
      indirect_mcmc_error <- max(net$indirect[, 4]/sqrt(net$indirect[, 8]))
      indirect_mcmc_rule <- max(0.05 * net$indirect[, 4])
      diff <- max(net$diff[, 7])
      diff_mcmc_error <- max(net$diff[, 4]/sqrt(net$diff[, 8]))
      diff_mcmc_rule <- max(0.05 * net$diff[, 4])
    } else if (!is.null(net$EM) & length(net$EM[1, ]) == 9) {
      # From 'run_sensitivity' function
      EM <- max(net$EM[, 8])
      EM_mcmc_error <- max(net$EM[, 2]/sqrt(net$EM[, 9]))
      EM_mcmc_rule <- max(0.05 * net$EM[, 2])
      if (!is.null(net$tau)) {
        tau <- max(net$tau[, 8])
        tau_mcmc_error <- max(net$tau[, 2]/sqrt(net$tau[, 9]))
        tau_mcmc_rule <- max(0.05 * net$tau[, 2])
      } else {
        tau <- tau_mcmc_error <- tau_mcmc_rule <- NA
      }

      # From 'run_model' function
      delta <- delta_mcmc_error <- delta_mcmc_rule <- NA
      EM_pred <- EM_pred_mcmc_error <- EM_pred_mcmc_rule <- NA
      phi <- phi_mcmc_error <- phi_mcmc_rule <- NA

      # From 'run_metareg' function
      beta <- beta_mcmc_error <- beta_mcmc_rule <- NA

      # From 'run_nodesplit' function
      direct <- direct_mcmc_error <- direct_mcmc_rule <- NA
      indirect <- indirect_mcmc_error <- indirect_mcmc_rule <- NA
      diff <- diff_mcmc_error <- diff_mcmc_rule <- NA
    }
  }

  if (!is.null(net$jagsfit)) {
    # Turn 'R2jags' object into 'mcmc' object
    jagsfit_mcmc <- mcmcplots::as.mcmc.rjags(jagsfit)

    # An HTML file with a panel of diagnostic plots per monitored parameter
    mcmc_plot <- mcmcplots::mcmcplot(jagsfit_mcmc, parms = par)
  }

  r_hat_max <- c(EM,
                 EM_pred,
                 delta,
                 tau,
                 direct,
                 indirect,
                 diff,
                 phi,
                 beta)

  mcmc_error_max <- c(EM_mcmc_error,
                      EM_pred_mcmc_error,
                      delta_mcmc_error,
                      tau_mcmc_error,
                      direct_mcmc_error,
                      indirect_mcmc_error,
                      diff_mcmc_error,
                      phi_mcmc_error,
                      beta_mcmc_error)

  mcmc_rule_max <- c(EM_mcmc_rule,
                     EM_pred_mcmc_rule,
                     delta_mcmc_rule,
                     tau_mcmc_rule,
                     direct_mcmc_rule,
                     indirect_mcmc_rule,
                     diff_mcmc_rule,
                     phi_mcmc_rule,
                     beta_mcmc_rule)

  # Indicate whether each model parameter achieved or failed to converge
  conv <- rep(NA, length(r_hat_max))
  for (i in seq_len(length(r_hat_max))) {
    conv[i] <- ifelse(is.na(r_hat_max[i]), "Not applicable",
                      ifelse(!is.na(r_hat_max[i]) & r_hat_max[i] < 1.1 &
                               mcmc_error_max[i] < mcmc_rule_max[i],
                             "achieved", "failed"))
  }

  # A data-frame on convergence for all monitored parameters using the Rhat
  monitored_parameters <- c("Effect estimates (EM)",
                            "Predictions (EM_pred)",
                            "Within-trial estimates (delta)",
                            "Between-trial standard deviation (tau)",
                            "Direct effects (node-splitting; direct)",
                            "Indirect effect(s) (node-splitting; indirect)",
                            "Inconsistency factor(s) (node-splitting; diff)",
                            "Informative missingness parameter(s) (phi)",
                            "Regression coefficient(s) (beta)")
  convergence <- data.frame(monitored_parameters,
                            round(r_hat_max, 4),
                            round(mcmc_error_max, 4),
                            round(mcmc_rule_max, 4),
                            conv)
  rownames(convergence) <- NULL
  colnames(convergence) <- c("parameters",
                             "max R-hat",
                             "max MCMC error",
                             "MCMC rule",
                             "convergence")

  return(knitr::kable(convergence,
                      align = "lccl",
                      caption = "Markov Chain Monte Carlo diagnostic measures"))
}

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
#' @param par A vector of at least one character string that refers to the
#'   monitored parameters in \code{jagsfit} which is an object of S3 class
#'   \code{\link{run_metareg}}, \code{\link{run_model}}, and
#'   \code{\link{run_ume}}. The selected parameters will be considered in the
#'   diagnostic plots (see 'Value'). This argument will be ignored for objects
#'   of S3 class \code{\link{run_nodesplit}}, \code{\link{run_sensitivity}},
#'   and \code{\link{run_series_meta}}.
#'
#' @return \code{mcmc_diagnostics} considers the following monitored parameters:
#'   \item{EM}{The estimated summary effect measure.}
#'   \item{EM_pred}{The predicted summary effect measure.}
#'   \item{delta}{The estimated trial-specific effect measure.}
#'   \item{tau}{The between-trial standard deviation.}
#'   \item{direct}{The direct estimate of the split node (see 'Value' in
#'   \code{\link{run_nodesplit}}).}
#'   \item{indirect}{The indirect estimate of the split node
#'   (see 'Value' in \code{\link{run_nodesplit}}).}
#'   \item{diff}{The inconsistency factor of the split node (see 'Value' in
#'   \code{\link{run_nodesplit}}).}
#'   \item{phi}{The informative missingness parameter.}
#'   \item{beta}{The regression coefficient.}
#'
#'   For each monitored parameter mentioned above, \code{mcmc_diagnostics} also
#'   returns a barplot on the ratio of MCMC error to the posterior standard
#'   deviation and a barplot on the Gelman-Rubin R diagnostic. Bars that
#'   correspond to a ratio less than 5\% are indicated in green (the
#'   corresponding parameters have been estimated accurately); otherwise, the
#'   bars are indicated in red (inaccurate estimation). Furthermore, bars that
#'   correspond to an R value less than 1.10 are indicated in green (the
#'   corresponding parameters have been converged); otherwise, the bars are
#'   indicated in red (convergence is not achieved).
#'   \code{mcmc_diagnostics} returns histograms than barplots for \code{EM} when
#'   \code{\link{run_sensitivity}} is considered.
#'
#'   \code{mcmc_diagnostics} also uses the
#'   \code{\link[mcmcplots:mcmcplot]{mcmcplot}} function of the R-package
#'   \href{https://CRAN.R-project.org/package=mcmcplots}{mcmcplots} to create an
#'   HTML file with a panel of diagnostic plots (trace, density, and
#'   autocorrelation) for each monitored parameter.
#'
#' @details For each monitored parameter, \code{mcmc_diagnostics} considers the
#'   R-hat and MCMC error and compares them with the thresholds 1.1 and 5\% of
#'   the posterior standard deviation (the rule of thumb), respectively.
#'   Convergence is achieved for the monitored parameter, when the R-hat is
#'   below the corresponding threshold. Visual inspection of the trace plots
#'   and posterior density of the monitored parameters should also be considered
#'   when drawing conclusions about convergence.
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

  if (!is.element(class(net),
                  c("run_model", "run_metareg", "run_series_meta",
                    "run_nodesplit", "run_ume", "run_sensitivity")) ||
      is.null(net)) {
    aa <- "'run_model', 'run_metareg', 'run_series_meta',"
    bb <- "'run_nodesplit', 'run_ume', or 'run_sensitivity'."
    stop(paste("'net' must be an object of S3 class", aa, bb), call. = FALSE)
  }
  a <- "Visual inspection of the trace plot and posterior density of the"
  b <- "monitored parameters should *also* be considered when concluding"
  c <- "about convergence."
  message(paste("R-hat < 1.10 is an indication of convergence.", a, b, c))

  par <- if (!is.null(net$jagsfit) & missing(par)) {
    stop("The argument 'par' needs to be defined.", call. = FALSE)
  } else if (!is.null(net$jagsfit) & !is.null(par)) {
    par
  } else if (is.null(net$jagsfit) & !is.null(par)) {
    aa <- "Note: The argument 'par' is ignored. It is used only"
    bb <- "with 'run_metareg', 'run_model', and 'run_ume'."
    message(paste(aa, bb))
    NULL
  }

  if (!is.null(net$jagsfit)) {
    jagsfit <- net$jagsfit
    # Turn 'R2jags' object into 'mcmc' object
    jagsfit_mcmc <- mcmcplots::as.mcmc.rjags(jagsfit)

    # An HTML file with a panel of diagnostic plots per monitored parameter
    mcmc_plot <- mcmcplots::mcmcplot(jagsfit_mcmc, parms = par)
  }

  name <- ratio <- rhat <- NA

  ## Plot R diagnostic and MCMC rule per monitored parameter
  # Effect estimates (EM)
  EM_plot_ratio <- if (is.element(class(net),
                                  c("run_model", "run_metareg", "run_ume"))) {
    get_results <- as.data.frame(t(net$jagsfit$BUGSoutput$summary))
    #EM <- t(get_results %>% select(starts_with("EM[")))
    EM <- t(get_results)[startsWith(rownames(t(get_results)), "EM["), ]
    EM_mcmc_rule <- 1 / sqrt(EM[, 9])
    ggplot(data.frame(name = names(EM_mcmc_rule),
                      ratio = EM_mcmc_rule),
           aes(x = name,
               y = ratio,
               fill = factor(ifelse(ratio < 0.05, "Yes", "No"),
                             levels = c("Yes", "No")))) +
      geom_bar(stat = "identity", width = 0.5) +
      geom_hline(yintercept = 0.05, linetype = 2) +
      scale_fill_manual(values = c("Yes" = "#009E73", "No" = "#D55E00")) +
      labs(x = "",
           y = "MCSE to posterior standard deviation",
           fill = "Ratio < 0.05") +
      theme_classic() +
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_line(colour = "transparent"),
            axis.text = element_text(size = 12),
            axis.title = element_text(size = 12, face = "bold"),
            legend.position = "bottom",
            legend.text = element_text(size = 12),
            legend.title = element_text(size = 12, face = "bold"))
  } else if (inherits(net, "run_sensitivity")) {
    EM_mcmc_rule <- 1 / sqrt(net$EM[, 9])
    ggplot(data.frame(name = names(EM_mcmc_rule),
                      ratio = EM_mcmc_rule),
           aes(x = EM_mcmc_rule,
               fill = factor(ifelse(ratio < 0.05, "Yes", "No"),
                             levels = c("Yes", "No")))) +
      geom_histogram() +
      geom_vline(xintercept = 0.05, linetype = 2) +
      scale_fill_manual(values = c("Yes" = "#009E73", "No" = "#D55E00")) +
      labs(x = "MCSE to posterior standard deviation",
           y = "",
           fill = "Ratio < 0.05") +
      coord_cartesian(xlim = c(min(EM_mcmc_rule), max(EM_mcmc_rule))) +
      theme_classic() +
      theme(axis.text = element_text(size = 12),
            axis.title = element_text(size = 12, face = "bold"),
            legend.position = "bottom",
            legend.text = element_text(size = 12),
            legend.title = element_text(size = 12, face = "bold"))
  } else if (inherits(net, "run_series_meta")) {
    EM_mcmc_rule <- 1 / sqrt(net$EM[, 11])
    ggplot(data.frame(name = paste("EM:", net$EM[, 2], "vs", net$EM[, 1]),
                      ratio = EM_mcmc_rule),
           aes(x = name,
               y = ratio,
               fill = factor(ifelse(ratio < 0.05, "Yes", "No"),
                             levels = c("Yes", "No")))) +
      geom_bar(stat = "identity", width = 0.5) +
      geom_hline(yintercept = 0.05, linetype = 2) +
      scale_fill_manual(values = c("Yes" = "#009E73", "No" = "#D55E00")) +
      labs(x = "",
           y = "MCSE to posterior standard deviation",
           fill = "Ratio < 0.05") +
      theme_classic() +
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_line(colour = "transparent"),
            axis.text = element_text(size = 12),
            axis.title = element_text(size = 12, face = "bold"),
            legend.position = "bottom",
            legend.text = element_text(size = 12),
            legend.title = element_text(size = 12, face = "bold"))
  } else if (inherits(net, "run_nodesplit")) {
    NULL
  }

  EM_plot_rhat <- if (is.element(class(net),
                                 c("run_model", "run_metareg", "run_ume"))) {
    get_results <- as.data.frame(t(net$jagsfit$BUGSoutput$summary))
    #EM <- t(get_results %>% select(starts_with("EM[")))
    EM <- t(get_results)[startsWith(rownames(t(get_results)), "EM["), ]
    ggplot(data.frame(name = names(EM[, 8]), rhat = EM[, 8]),
           aes(x = name,
               y = rhat,
               fill = factor(ifelse(rhat < 1.10, "Yes", "No"),
                             levels = c("Yes", "No")))) +
      geom_bar(stat = "identity", width = 0.5) +
      geom_hline(yintercept = 1.10, linetype = 2) +
      scale_fill_manual(values = c("Yes" = "#009E73", "No" = "#D55E00")) +
      labs(x = "",
           y = bquote(Gelman-Rubin~hat(R)~diagnostic),
           fill = "R < 1.10") +
      scale_y_continuous(breaks = c(seq(0, max(EM[, 8]), by = 0.2), 1.10)) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
            axis.text = element_text(size = 12),
            axis.title = element_text(size = 12, face = "bold"),
            legend.position = "bottom",
            legend.text = element_text(size = 12),
            legend.title = element_text(size = 12, face = "bold"))
  } else if (inherits(net, "run_sensitivity")) {
    ggplot(data.frame(name = names(net$EM[, 8]), rhat = net$EM[, 8]),
           aes(x = rhat,
               fill = factor(ifelse(rhat < 1.10, "Yes", "No"),
                             levels = c("Yes", "No")))) +
      geom_histogram() +
      geom_vline(xintercept = 1.10, linetype = 2) +
      scale_fill_manual(values = c("Yes" = "#009E73", "No" = "#D55E00")) +
      labs(x = bquote(Gelman-Rubin~hat(R)~diagnostic),
           y = "",
           fill = "R < 1.10") +
      coord_cartesian(xlim = c(min(net$EM[, 8]), max(net$EM[, 8]))) +
      theme_classic() +
      theme(axis.text = element_text(size = 12),
            axis.title = element_text(size = 12, face = "bold"),
            legend.position = "bottom",
            legend.text = element_text(size = 12),
            legend.title = element_text(size = 12, face = "bold"))
  } else if (inherits(net, "run_series_meta")) {
    ggplot(data.frame(name = paste("EM:", net$EM[, 2], "vs", net$EM[, 1]),
                      rhat = net$EM[, 10]),
           aes(x = name,
               y = rhat,
               fill = factor(ifelse(rhat < 1.10, "Yes", "No"),
                             levels = c("Yes", "No")))) +
      geom_bar(stat = "identity", width = 0.5) +
      geom_hline(yintercept = 1.10, linetype = 2) +
      scale_fill_manual(values = c("Yes" = "#009E73", "No" = "#D55E00")) +
      labs(x = "",
           y = bquote(Gelman-Rubin~hat(R)~diagnostic),
           fill = "R < 1.10") +
      scale_y_continuous(breaks = c(seq(0, max(net$EM[, 10]), 0.2), 1.10)) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
            axis.text = element_text(size = 12),
            axis.title = element_text(size = 12, face = "bold"),
            legend.position = "bottom",
            legend.text = element_text(size = 12),
            legend.title = element_text(size = 12, face = "bold"))
  } else if (inherits(net, "run_nodesplit")) {
    NULL
  }

  EM_plot <- if(inherits(net, "run_nodesplit")) {
    NULL
  } else {
    suppressMessages({ggarrange(EM_plot_ratio, EM_plot_rhat, nrow = 2)})
  }

  # Predictions (EM_pred)
  EM_pred_plot_ratio <- if (any(is.element(class(net),
                                           c("run_model", "run_metareg")) &
                                net$model == "RE")) {
    get_results <- as.data.frame(t(net$jagsfit$BUGSoutput$summary))
    #EM_pred <- t(get_results %>% select(starts_with("EM.pred[")))
    EM_pred <-
      t(get_results)[startsWith(rownames(t(get_results)), "EM.pred["), ]
    EM_pred_mcmc_rule <- 1 / sqrt(EM_pred[, 9])
    ggplot(data.frame(name = names(EM_pred_mcmc_rule),
                      ratio = EM_pred_mcmc_rule),
           aes(x = name,
               y = ratio,
               fill = factor(ifelse(ratio < 0.05, "Yes", "No"),
                             levels = c("Yes", "No")))) +
      geom_bar(stat = "identity", width = 0.5) +
      geom_hline(yintercept = 0.05, linetype = 2) +
      scale_fill_manual(values = c("Yes" = "#009E73", "No" = "#D55E00")) +
      labs(x = "",
           y = "MCSE to posterior standard deviation",
           fill = "Ratio < 0.05") +
      theme_classic() +
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_line(colour = "transparent"),
            axis.text = element_text(size = 12),
            axis.title = element_text(size = 12, face = "bold"),
            legend.position = "bottom",
            legend.text = element_text(size = 12),
            legend.title = element_text(size = 12, face = "bold"))
  } else {
    NULL
  }

  EM_pred_plot_rhat <- if (any(is.element(class(net),
                                          c("run_model", "run_metareg")) &
                               net$model == "RE")) {
    get_results <- as.data.frame(t(net$jagsfit$BUGSoutput$summary))
    #EM_pred <- t(get_results %>% select(starts_with("EM.pred[")))
    EM_pred <-
      t(get_results)[startsWith(rownames(t(get_results)), "EM.pred["), ]
    ggplot(data.frame(name = names(EM_pred[, 8]), rhat = EM_pred[, 8]),
           aes(x = name,
               y = rhat,
               fill = factor(ifelse(rhat < 1.10, "Yes", "No"),
                             levels = c("Yes", "No")))) +
      geom_bar(stat = "identity", width = 0.5) +
      geom_hline(yintercept = 1.10, linetype = 2) +
      scale_fill_manual(values = c("Yes" = "#009E73", "No" = "#D55E00")) +
      labs(x = "",
           y = bquote(Gelman-Rubin~hat(R)~diagnostic),
           fill = "R < 1.10") +
      scale_y_continuous(breaks = c(seq(0, max(EM_pred[, 8]), 0.2), 1.10)) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
            axis.text = element_text(size = 12),
            axis.title = element_text(size = 12, face = "bold"),
            legend.position = "bottom",
            legend.text = element_text(size = 12),
            legend.title = element_text(size = 12, face = "bold"))
  } else {
    NULL
  }

  EM_pred_plot <-
    if(all(!is.element(class(net), c("run_model", "run_metareg")) |
           net$model == "FE")) {
      NULL
  } else {
    ggarrange(EM_pred_plot_ratio, EM_pred_plot_rhat, nrow = 2)
  }

  # Within-trial estimates (delta)
  delta_plot_ratio <- if (any(is.element(class(net),
                                         c("run_model", "run_metareg")) &
                          net$model == "RE")) {
    get_results <- as.data.frame(t(net$jagsfit$BUGSoutput$summary))
    #delta <- t(get_results %>% select(starts_with("delta") & !ends_with(",1]")))
    delta <-
      t(get_results)[startsWith(rownames(t(get_results)), "delta") &
                       !endsWith(rownames(t(get_results)), ",1]"), ]
    delta_mcmc_rule <- 1 / sqrt(delta[, 9])
    ggplot(data.frame(name = names(delta_mcmc_rule), ratio = delta_mcmc_rule),
           aes(x = name,
               y = ratio,
               fill = factor(ifelse(ratio < 0.05, "Yes", "No"),
                             levels = c("Yes", "No")))) +
      geom_bar(stat = "identity", width = 0.5) +
      geom_hline(yintercept = 0.05, linetype = 2) +
      scale_fill_manual(values = c("Yes" = "#009E73", "No" = "#D55E00")) +
      labs(x = "",
           y = "MCSE to posterior standard deviation",
           fill = "Ratio < 0.05") +
      theme_classic() +
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_line(colour = "transparent"),
            axis.text = element_text(size = 12),
            axis.title = element_text(size = 12, face = "bold"),
            legend.position = "bottom",
            legend.text = element_text(size = 12),
            legend.title = element_text(size = 12, face = "bold"))
  } else {
    NULL
  }

  delta_plot_rhat <-
    if (any(is.element(class(net), c("run_model", "run_metareg")) &
            net$model == "RE")) {
      get_results <- as.data.frame(t(net$jagsfit$BUGSoutput$summary))
      #delta <- t(get_results %>% select(starts_with("delta") &
      #                                    !ends_with(",1]")))
      delta <-
        t(get_results)[startsWith(rownames(t(get_results)), "delta") &
                         !endsWith(rownames(t(get_results)), ",1]"), ]
      ggplot(data.frame(name = names(delta[, 8]), rhat = delta[, 8]),
             aes(x = name,
                 y = rhat,
                 fill = factor(ifelse(rhat < 1.10, "Yes", "No"),
                               levels = c("Yes", "No")))) +
        geom_bar(stat = "identity", width = 0.5) +
        geom_hline(yintercept = 1.10, linetype = 2) +
        scale_fill_manual(values = c("Yes" = "#009E73", "No" = "#D55E00")) +
        labs(x = "",
             y = bquote(Gelman-Rubin~hat(R)~diagnostic),
             fill = "R < 1.10") +
        scale_y_continuous(breaks = c(seq(0, max(delta[, 8]), by = 0.2), 1.1)) +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
              axis.text = element_text(size = 12),
              axis.title = element_text(size = 12, face = "bold"),
              legend.position = "bottom",
              legend.text = element_text(size = 12),
              legend.title = element_text(size = 12, face = "bold"))
  } else {
    NULL
  }

  delta_plot <-
    if(all(!is.element(class(net), c("run_model", "run_metareg")) |
                   net$model == "FE")) {
      NULL
  } else {
    ggarrange(delta_plot_ratio, delta_plot_rhat, nrow = 2)
  }

  # Direct effects (node-splitting; direct)
  direct_plot_ratio <- if(!inherits(net, "run_nodesplit")) {
    NULL
  } else {
    direct_mcmc_rule <- 1 / sqrt(net$direct[, 8])
    ggplot(data.frame(name = paste("direct:",
                                   net$direct[, 1], "vs", net$direct[, 2]),
                      ratio = direct_mcmc_rule),
           aes(x = name,
               y = ratio,
               fill = factor(ifelse(ratio < 0.05, "Yes", "No"),
                             levels = c("Yes", "No")))) +
      geom_bar(stat = "identity", width = 0.5) +
      geom_hline(yintercept = 0.05, linetype = 2) +
      scale_fill_manual(values = c("Yes" = "#009E73", "No" = "#D55E00")) +
      labs(x = "",
           y = "MCSE to posterior standard deviation",
           fill = "Ratio < 0.05") +
      theme_classic() +
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_line(colour = "transparent"),
            axis.text = element_text(size = 12),
            axis.title = element_text(size = 12, face = "bold"),
            legend.position = "bottom",
            legend.text = element_text(size = 12),
            legend.title = element_text(size = 12, face = "bold"))
  }

  direct_plot_rhat <- if(!inherits(net, "run_nodesplit")) {
    NULL
  } else {
    ggplot(data.frame(name = paste("direct:",
                                   net$direct[, 1], "vs", net$direct[, 2]),
                      rhat = net$direct[, 7]),
           aes(x = name,
               y = rhat,
               fill = factor(ifelse(rhat < 1.10, "Yes", "No"),
                             levels = c("Yes", "No")))) +
      geom_bar(stat = "identity", width = 0.5) +
      geom_hline(yintercept = 1.10, linetype = 2) +
      scale_fill_manual(values = c("Yes" = "#009E73", "No" = "#D55E00")) +
      labs(x = "",
           y = bquote(Gelman-Rubin~hat(R)~diagnostic),
           fill = "R < 1.10") +
      scale_y_continuous(breaks = c(seq(0, max(net$direct[, 7]), 0.2), 1.10)) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
            axis.text = element_text(size = 12),
            axis.title = element_text(size = 12, face = "bold"),
            legend.position = "bottom",
            legend.text = element_text(size = 12),
            legend.title = element_text(size = 12, face = "bold"))
  }

  direct_plot <- if(!inherits(net, "run_nodesplit")) {
    NULL
  } else {
    ggarrange(direct_plot_ratio, direct_plot_rhat, nrow = 2)
  }

  # Indirect effects (node-splitting; indirect)
  indirect_plot_ratio <- if(!inherits(net, "run_nodesplit")) {
    NULL
  } else {
    indirect_mcmc_rule <- 1 / sqrt(net$indirect[, 8])
    ggplot(data.frame(name = paste("indirect:",
                                   net$indirect[, 1], "vs", net$indirect[, 2]),
                      ratio = indirect_mcmc_rule),
           aes(x = name,
               y = ratio,
               fill = factor(ifelse(ratio < 0.05, "Yes", "No"),
                             levels = c("Yes", "No")))) +
      geom_bar(stat = "identity", width = 0.5) +
      geom_hline(yintercept = 0.05, linetype = 2) +
      scale_fill_manual(values = c("Yes" = "#009E73", "No" = "#D55E00")) +
      labs(x = "",
           y = "MCSE to posterior standard deviation",
           fill = "Ratio < 0.05") +
      theme_classic() +
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_line(colour = "transparent"),
            axis.text = element_text(size = 12),
            axis.title = element_text(size = 12, face = "bold"),
            legend.position = "bottom",
            legend.text = element_text(size = 12),
            legend.title = element_text(size = 12, face = "bold"))
  }

  indirect_plot_rhat <- if(!inherits(net, "run_nodesplit")) {
    NULL
  } else {
    ggplot(data.frame(name = paste("indirect:",
                                   net$indirect[, 1], "vs", net$indirect[, 2]),
                      rhat = net$indirect[, 7]),
           aes(x = name,
               y = rhat,
               fill = factor(ifelse(rhat < 1.10, "Yes", "No"),
                             levels = c("Yes", "No")))) +
      geom_bar(stat = "identity", width = 0.5) +
      geom_hline(yintercept = 1.10, linetype = 2) +
      scale_fill_manual(values = c("Yes" = "#009E73", "No" = "#D55E00")) +
      labs(x = "",
           y = bquote(Gelman-Rubin~hat(R)~diagnostic),
           fill = "R < 1.10") +
      scale_y_continuous(breaks = c(seq(0, max(net$indirect[, 7]), 0.2), 1.1)) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
            axis.text = element_text(size = 12),
            axis.title = element_text(size = 12, face = "bold"),
            legend.position = "bottom",
            legend.text = element_text(size = 12),
            legend.title = element_text(size = 12, face = "bold"))
  }

  indirect_plot <- if(!inherits(net, "run_nodesplit")) {
    NULL
  } else {
    ggarrange(indirect_plot_ratio, indirect_plot_rhat, nrow = 2)
  }

  # Inconsistency factor(s) (node-splitting; diff)
  diff_plot_ratio <- if(!inherits(net, "run_nodesplit")) {
    NULL
  } else {
    diff_mcmc_rule <- 1 / sqrt(net$diff[, 8])
    ggplot(data.frame(name = paste("diff:", net$diff[, 1], "vs", net$diff[, 2]),
                      ratio = diff_mcmc_rule),
           aes(x = name,
               y = ratio,
               fill = factor(ifelse(ratio < 0.05, "Yes", "No"),
                             levels = c("Yes", "No")))) +
      geom_bar(stat = "identity", width = 0.5) +
      geom_hline(yintercept = 0.05, linetype = 2) +
      scale_fill_manual(values = c("Yes" = "#009E73", "No" = "#D55E00")) +
      labs(x = "",
           y = "MCSE to posterior standard deviation",
           fill = "Ratio < 0.05") +
      theme_classic() +
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_line(colour = "transparent"),
            axis.text = element_text(size = 12),
            axis.title = element_text(size = 12, face = "bold"),
            legend.position = "bottom",
            legend.text = element_text(size = 12),
            legend.title = element_text(size = 12, face = "bold"))
  }

  diff_plot_rhat <- if(!inherits(net, "run_nodesplit")) {
    NULL
  } else {
    ggplot(data.frame(name = paste("diff:", net$diff[, 1], "vs", net$diff[, 2]),
                      rhat = net$diff[, 7]),
           aes(x = name,
               y = rhat,
               fill = factor(ifelse(rhat < 1.10, "Yes", "No"),
                             levels = c("Yes", "No")))) +
      geom_bar(stat = "identity", width = 0.5) +
      geom_hline(yintercept = 1.10, linetype = 2) +
      scale_fill_manual(values = c("Yes" = "#009E73", "No" = "#D55E00")) +
      labs(x = "",
           y = bquote(Gelman-Rubin~hat(R)~diagnostic),
           fill = "R < 1.10") +
      scale_y_continuous(breaks = c(seq(0, max(net$diff[, 7]), 0.2), 1.10)) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
            axis.text = element_text(size = 12),
            axis.title = element_text(size = 12, face = "bold"),
            legend.position = "bottom",
            legend.text = element_text(size = 12),
            legend.title = element_text(size = 12, face = "bold"))
  }

  diff_plot <- if(!inherits(net, "run_nodesplit")) {
    NULL
  } else {
    ggarrange(diff_plot_ratio, diff_plot_rhat, nrow = 2)
  }

  # Informative missingness parameter(s) (phi)
  phi_plot_ratio <- if(any(!is.null(net$phi) &
                           !is.element(net$assumption,
                                       c("IDE-COMMON", "HIE-COMMON")))) {
    get_results <- as.data.frame(t(net$jagsfit$BUGSoutput$summary))
    #phi <- t(get_results %>% select(starts_with("mean.phi[") |
    #                                  starts_with("phi[")))
    phi <- t(get_results)[startsWith(rownames(t(get_results)), "mean.phi[") |
                            startsWith(rownames(t(get_results)), "phi["), ]
    phi_mcmc_rule <- 1 / sqrt(phi[, 9])
    ggplot(data.frame(name = names(phi_mcmc_rule), ratio = phi_mcmc_rule),
           aes(x = name,
               y = ratio,
               fill = factor(ifelse(ratio < 0.05, "Yes", "No"),
                             levels = c("Yes", "No")))) +
      geom_bar(stat = "identity", width = 0.5) +
      geom_hline(yintercept = 0.05, linetype = 2) +
      scale_fill_manual(values = c("Yes" = "#009E73", "No" = "#D55E00")) +
      labs(x = "",
           y = "MCSE to posterior standard deviation",
           fill = "Ratio < 0.05") +
      theme_classic() +
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_line(colour = "transparent"),
            axis.text = element_text(size = 12),
            axis.title = element_text(size = 12, face = "bold"),
            legend.position = "bottom",
            legend.text = element_text(size = 12),
            legend.title = element_text(size = 12, face = "bold"))
  } else {
    NULL
  }

  phi_plot_rhat <- if(any(!is.null(net$phi) &
                          !is.element(net$assumption,
                                      c("IDE-COMMON", "HIE-COMMON")))) {
    get_results <- as.data.frame(t(net$jagsfit$BUGSoutput$summary))
    #phi <- t(get_results %>% select(starts_with("mean.phi[") |
    #                                  starts_with("phi[")))
    phi <- t(get_results)[startsWith(rownames(t(get_results)), "mean.phi[") |
                            startsWith(rownames(t(get_results)), "phi["), ]
    ggplot(data.frame(name = names(phi[, 8]), rhat = phi[, 8]),
           aes(x = name,
               y = rhat,
               fill = factor(ifelse(rhat < 1.10, "Yes", "No"),
                             levels = c("Yes", "No")))) +
      geom_bar(stat = "identity", width = 0.5) +
      geom_hline(yintercept = 1.10, linetype = 2) +
      scale_fill_manual(values = c("Yes" = "#009E73", "No" = "#D55E00")) +
      labs(x = "",
           y = bquote(Gelman-Rubin~hat(R)~diagnostic),
           fill = "R < 1.10") +
      scale_y_continuous(breaks = c(seq(0, max(phi[, 8]), by = 0.2), 1.10)) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
            axis.text = element_text(size = 12),
            axis.title = element_text(size = 12, face = "bold"),
            legend.position = "bottom",
            legend.text = element_text(size = 12),
            legend.title = element_text(size = 12, face = "bold"))

  } else {
    NULL
  }

  phi_plot <-
    if(any(!is.null(net$phi) &
           !is.element(net$assumption, c("IDE-COMMON", "HIE-COMMON")))) {
      ggarrange(phi_plot_ratio, phi_plot_rhat, nrow = 2)
  } else {
    NULL
  }

  tabulate_phi <-
    if(any(!is.null(net$phi) & is.element(net$assumption,
                                          c("IDE-COMMON", "HIE-COMMON")))) {
      get_results <- as.data.frame(t(net$jagsfit$BUGSoutput$summary))
      #phi <- t(get_results %>% select(starts_with("phi") |
      #                                  starts_with("mean.phi")))
      phi <-
        t(get_results)[startsWith(rownames(t(get_results)), "phi") |
                         startsWith(rownames(t(get_results)), "mean.phi"), ]
      data.frame(R.hat = phi[8], MCMC.rule = 1 / sqrt(phi[9]))
  } else {
    data.frame(R.hat = "Not applicable", MCMC.rule = "Not applicable")
  }

  # Regression coefficient(s) (beta)
  beta_plot_ratio <- if (any(inherits(net, "run_metareg") &
                             net$covar_assumption != "common")) {
    get_results <- as.data.frame(t(net$jagsfit$BUGSoutput$summary))
    #beta <- t(get_results %>% select(starts_with("beta.all[")))
    beta <- t(get_results)[startsWith(rownames(t(get_results)), "beta.all["), ]
    beta_mcmc_rule <- 1 / sqrt(beta[, 9])
    ggplot(data.frame(name = names(beta_mcmc_rule),
                      ratio = beta_mcmc_rule),
           aes(x = name,
               y = ratio,
               fill = factor(ifelse(ratio < 0.05, "Yes", "No"),
                             levels = c("Yes", "No")))) +
      geom_bar(stat = "identity", width = 0.5) +
      geom_hline(yintercept = 0.05, linetype = 2) +
      scale_fill_manual(values = c("Yes" = "#009E73", "No" = "#D55E00")) +
      labs(x = "",
           y = "MCSE to posterior standard deviation",
           fill = "Ratio < 0.05") +
      theme_classic() +
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_line(colour = "transparent"),
            axis.text = element_text(size = 12),
            axis.title = element_text(size = 12, face = "bold"),
            legend.position = "bottom",
            legend.text = element_text(size = 12),
            legend.title = element_text(size = 12, face = "bold"))
  } else if (!inherits(net, "run_metareg")) {
    NULL
  }

  beta_plot_rhat <- if (any(inherits(net, "run_metareg") &
                            net$covar_assumption != "common")) {
    get_results <- as.data.frame(t(net$jagsfit$BUGSoutput$summary))
    #beta <- t(get_results %>% select(starts_with("beta.all[")))
    beta <- t(get_results)[startsWith(rownames(t(get_results)), "beta.all["), ]
    ggplot(data.frame(name = names(beta[, 8]), rhat = beta[, 8]),
           aes(x = name,
               y = rhat,
               fill = factor(ifelse(rhat < 1.10, "Yes", "No"),
                             levels = c("Yes", "No")))) +
      geom_bar(stat = "identity", width = 0.5) +
      geom_hline(yintercept = 1.10, linetype = 2) +
      scale_fill_manual(values = c("Yes" = "#009E73", "No" = "#D55E00")) +
      labs(x = "",
           y = bquote(Gelman-Rubin~hat(R)~diagnostic),
           fill = "R < 1.10") +
      scale_y_continuous(breaks = c(seq(0, max(beta[, 8]), by = 0.2), 1.10)) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
            axis.text = element_text(size = 12),
            axis.title = element_text(size = 12, face = "bold"),
            legend.position = "bottom",
            legend.text = element_text(size = 12),
            legend.title = element_text(size = 12, face = "bold"))
  } else if (!inherits(net, "run_metareg")) {
    NULL
  }

  beta_plot <- if (any(inherits(net, "run_metareg")  &
                       net$covar_assumption != "common")) {
    ggarrange(beta_plot_ratio, beta_plot_rhat, nrow = 2)
  } else {
    NULL
  }

  tabulate_beta <-
    if (any(inherits(net, "run_metareg") & net$covar_assumption == "common")) {
      get_results <- as.data.frame(t(net$jagsfit$BUGSoutput$summary))
      #beta <-
      #  t(get_results %>% select(starts_with("beta") & !starts_with("beta.all[")))
      beta <-
        t(get_results)[startsWith(rownames(t(get_results)), "beta") &
                         !startsWith(rownames(t(get_results)), "beta.all["), ]
      data.frame(R.hat = beta[8], MCMC.rule = 1 / sqrt(beta[9]))
    } else {
      data.frame(R.hat = "Not applicable", MCMC.rule = "Not applicable")
    }

  # Between-trial standard deviation (tau)
  tabulate_tau <-
    if(any(is.element(class(net), c("run_model", "run_metareg", "run_ume")) &
           net$model == "RE")) {
      get_results <- as.data.frame(t(net$jagsfit$BUGSoutput$summary))
      #tau <- t(get_results %>% select(starts_with("tau")))
      tau <- t(get_results)[startsWith(rownames(t(get_results)), "tau"), ]
      data.frame(R.hat = tau[8], MCMC.rule = 1 / sqrt(tau[9]))
  } else if (all(!is.element(class(net),
                             c("run_model", "run_metareg", "run_ume")) |
                 net$model == "FE")) {
    data.frame(R.hat = "Not applicable", MCMC.rule = "Not applicable")
  }

  tau_plot_ratio <- if (any(inherits(net, "run_nodesplit") &
                            net$model == "RE")) {
    tau_mcmc_rule <- 1 / sqrt(net$tau[, 8])
    ggplot(data.frame(name = paste("tau:", net$tau[, 1], "vs", net$tau[, 2]),
                      ratio = tau_mcmc_rule),
           aes(x = name,
               y = ratio,
               fill = factor(ifelse(ratio < 0.05, "Yes", "No"),
                             levels = c("Yes", "No")))) +
      geom_bar(stat = "identity", width = 0.5) +
      geom_hline(yintercept = 0.05, linetype = 2) +
      scale_fill_manual(values = c("Yes" = "#009E73", "No" = "#D55E00")) +
      labs(x = "",
           y = "MCSE to posterior standard deviation",
           fill = "Ratio < 0.05") +
      theme_classic() +
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_line(colour = "transparent"),
            axis.text = element_text(size = 12),
            axis.title = element_text(size = 12, face = "bold"),
            legend.position = "bottom",
            legend.text = element_text(size = 12),
            legend.title = element_text(size = 12, face = "bold"))
  } else if (any(inherits(net, "run_sensitivity") & net$model == "RE")) {
    tau_mcmc_rule <- 1 / sqrt(net$tau[, 9])
    ggplot(data.frame(name = paste("tau", 1:length(net$tau[, 8])),
                      ratio = tau_mcmc_rule),
           aes(x = name,
               y = ratio,
               fill = factor(ifelse(ratio < 0.05, "Yes", "No"),
                             levels = c("Yes", "No")))) +
      geom_bar(stat = "identity", width = 0.5) +
      geom_hline(yintercept = 0.05, linetype = 2) +
      scale_fill_manual(values = c("Yes" = "#009E73", "No" = "#D55E00")) +
      labs(x = "",
           y = "MCSE to posterior standard deviation",
           fill = "Ratio < 0.05") +
      theme_classic() +
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_line(colour = "transparent"),
            axis.text = element_text(size = 12),
            axis.title = element_text(size = 12, face = "bold"),
            legend.position = "bottom",
            legend.text = element_text(size = 12),
            legend.title = element_text(size = 12, face = "bold"))
  } else if (any(inherits(net, "run_series_meta") & net$model == "RE")) {
    tau_mcmc_rule <- 1 / sqrt(net$tau[, 11])
    ggplot(data.frame(name = paste("tau:", net$tau[, 2], "vs", net$tau[, 1]),
                      ratio = tau_mcmc_rule),
           aes(x = name,
               y = ratio,
               fill = factor(ifelse(ratio < 0.05, "Yes", "No"),
                             levels = c("Yes", "No")))) +
      geom_bar(stat = "identity", width = 0.5) +
      geom_hline(yintercept = 0.05, linetype = 2) +
      scale_fill_manual(values = c("Yes" = "#009E73", "No" = "#D55E00")) +
      labs(x = "",
           y = "MCSE to posterior standard deviation",
           fill = "Ratio < 0.05") +
      theme_classic() +
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_line(colour = "transparent"),
            axis.text = element_text(size = 12),
            axis.title = element_text(size = 12, face = "bold"),
            legend.position = "bottom",
            legend.text = element_text(size = 12),
            legend.title = element_text(size = 12, face = "bold"))
  } else {
    NULL
  }

  tau_plot_rhat <- if (any(inherits(net, "run_nodesplit") &
                           net$model == "RE")) {
    ggplot(data.frame(name = paste("tau:", net$tau[, 1], "vs", net$tau[, 2]),
                      rhat = net$tau[, 7]),
           aes(x = name,
               y = rhat,
               fill = factor(ifelse(rhat < 1.10, "Yes", "No"),
                             levels = c("Yes", "No")))) +
      geom_bar(stat = "identity", width = 0.5) +
      geom_hline(yintercept = 1.10, linetype = 2) +
      scale_fill_manual(values = c("Yes" = "#009E73", "No" = "#D55E00")) +
      labs(x = "",
           y = bquote(Gelman-Rubin~hat(R)~diagnostic),
           fill = "R < 1.10") +
      scale_y_continuous(breaks = c(seq(0, max(net$tau[, 7]), 0.2), 1.10)) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
            axis.text = element_text(size = 12),
            axis.title = element_text(size = 12, face = "bold"),
            legend.position = "bottom",
            legend.text = element_text(size = 12),
            legend.title = element_text(size = 12, face = "bold"))
  } else if (any(inherits(net, "run_sensitivity") & net$model == "RE")) {
    ggplot(data.frame(name = paste("tau", 1:length(net$tau[, 8])),
                      rhat = net$tau[, 8]),
           aes(x = name,
               y = rhat,
               fill = factor(ifelse(rhat < 1.10, "Yes", "No"),
                             levels = c("Yes", "No")))) +
      geom_bar(stat = "identity", width = 0.5) +
      geom_hline(yintercept = 1.10, linetype = 2) +
      scale_fill_manual(values = c("Yes" = "#009E73", "No" = "#D55E00")) +
      labs(x = "",
           y = bquote(Gelman-Rubin~hat(R)~diagnostic),
           fill = "R < 1.10") +
      scale_y_continuous(breaks = c(seq(0, max(net$tau[, 8]), 0.2), 1.10)) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
            axis.text = element_text(size = 12),
            axis.title = element_text(size = 12, face = "bold"),
            legend.position = "bottom",
            legend.text = element_text(size = 12),
            legend.title = element_text(size = 12, face = "bold"))
  } else if (any(inherits(net, "run_series_meta") & net$model == "RE")) {
    ggplot(data.frame(name = paste("tau:", net$tau[, 2], "vs", net$tau[, 1]),
                      rhat = net$tau[, 10]),
           aes(x = name,
               y = rhat,
               fill = factor(ifelse(rhat < 1.10, "Yes", "No"),
                             levels = c("Yes", "No")))) +
      geom_bar(stat = "identity", width = 0.5) +
      geom_hline(yintercept = 1.10, linetype = 2) +
      scale_fill_manual(values = c("Yes" = "#009E73", "No" = "#D55E00")) +
      labs(x = "",
           y = bquote(Gelman-Rubin~hat(R)~diagnostic),
           fill = "R < 1.10") +
      scale_y_continuous(breaks = c(seq(0, max(net$tau[, 10]), 0.2), 1.10)) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
            axis.text = element_text(size = 12),
            axis.title = element_text(size = 12, face = "bold"),
            legend.position = "bottom",
            legend.text = element_text(size = 12),
            legend.title = element_text(size = 12, face = "bold"))
  } else {
    NULL
  }

  tau_plot <-
    if (
      all(
        !is.element(class(net),
                    c("run_nodesplit", "run_sensitivity", "run_series_meta")) |
        net$model == "FE"
        )
      ) {
      NULL
  } else {
    ggarrange(tau_plot_ratio, tau_plot_rhat, nrow = 2)
  }

  results <- na.omit(list(Effect_estimates = EM_plot,
                          Predictions = EM_pred_plot,
                          Within_trial_estimates = delta_plot,
                          Between_trial_SD = tau_plot,
                          Direct_estimates = direct_plot,
                          Indirect_estimates = indirect_plot,
                          Inconsistency_factors = diff_plot,
                          Informative_missingness_parameters = phi_plot,
                          Regression_coefficients = beta_plot,
                          table_tau =
                          knitr::kable(
                            tabulate_tau,
                            caption =
                            "The (common) between-trial standard deviation"),
                          table_phi =
                          knitr::kable(
                            tabulate_phi,
                            caption =
                            "The common informative missingness parameter"),
                          table_beta =
                            knitr::kable(
                              tabulate_beta,
                              caption =
                                "The common regression coefficient")))

  return(Filter(Negate(is.null), results))
}

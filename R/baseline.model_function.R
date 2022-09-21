#' The baseline model for binary outcome
#'
#' @description
#'   To process the elements in the argument \code{base_risk} of the
#'   \code{\link{run_model}} function. It also runs the hierarchical baseline
#'   model, separately from the relative effects model as described in
#'   Dias et al. (2018) and Dias et al. (2013b). The output is to be passed to
#'   \code{\link{run_model}} and \code{\link{run_metareg}} to obtain the
#'   (unadjusted and adjusted, respectively) absolute risks for each
#'   intervention in the dataset.
#'
#' @param base_risk A scalar, a vector of length three with elements sorted in
#'   ascending order, or a matrix with two columns and number of rows equal to
#'   the number of relevant trials. In either choice the elements should be in
#'   the interval (0, 1). See 'Details' in \code{\link{run_model}}.
#'   This argument is only relevant for a binary outcome.
#' @param n_chains Positive integer specifying the number of chains for the
#'   MCMC sampling; an argument of the \code{\link[R2jags:jags]{jags}} function
#'   of the R-package \href{https://CRAN.R-project.org/package=R2jags}{R2jags}.
#'   The default argument is 2.
#' @param n_iter Positive integer specifying the number of Markov chains for the
#'   MCMC sampling; an argument of the \code{\link[R2jags:jags]{jags}} function
#'   of the R-package \href{https://CRAN.R-project.org/package=R2jags}{R2jags}.
#'   The default argument is 10000.
#' @param n_burnin Positive integer specifying the number of iterations to
#'   discard at the beginning of the MCMC sampling; an argument of the
#'   \code{\link[R2jags:jags]{jags}} function of the R-package
#'   \href{https://CRAN.R-project.org/package=R2jags}{R2jags}.
#'   The default argument is 1000.
#' @param n_thin Positive integer specifying the thinning rate for the
#'   MCMC sampling; an argument of the \code{\link[R2jags:jags]{jags}} function
#'   of the R-package \href{https://CRAN.R-project.org/package=R2jags}{R2jags}.
#'   The default argument is 1.
#'
#' @return When \code{base_risk} is scalar (fixed baseline), the function
#'   returns the user-defined baseline for the selected reference intervention.
#'   When \code{base_risk} is a vector (random baseline), the function returns
#'   a vector with the calculated logit of an event for the selected reference
#'   intervention and its precision. Finally, when \code{base_risk} is a matrix
#'   (predicted baseline), the function returns the following elements:
#'   \item{ref_base}{A vector with the posterior mean and precision of the
#'   predicted logit of an event for the selected reference intervention.
#'   This vector is be passed to \code{\link{run_model}} and
#'   \code{\link{run_metareg}}.}
#'   \item{mean_base_logit}{The posterior distribution of the summary mean of
#'   the random effects in the logit scale.}
#'   \item{tau_base_logit}{The posterior distribution of the between-trial
#'   standard deviation in the logit scale.}
#'
#'   When \code{base_risk} is a matrix, the function also returns a forest plot
#'   with the estimated trial-specific logit of an event and 95\% credible
#'   intervals (the random effects) alongside the corresponding observed logit
#'   of an event. A grey rectangular illustrates the summary mean and 95\%
#'   credible interval of the random effects.
#'
#' @details If \code{base_risk} is a matrix, \code{baseline_model} creates the
#'   hierarchical baseline model in the JAGS dialect of the BUGS language.
#'   The output of this function (see 'Value') constitutes the
#'   posterior mean and precision of the predicted logit of an event for the
#'   selected reference intervention and it is plugged in the WinBUGS code for
#'   the relative effects model (Dias et al., 2013a) via the
#'   \code{\link{prepare_model}} function. Following (Dias et al., 2013a), a
#'   uniform prior distribution is assigned on the between-trial standard
#'   deviation with upper and lower limit equal to 0 and 5, respectively.
#'
#' @author {Loukia M. Spineli}
#'
#' @seealso \code{\link{prepare_model}}, \code{\link[R2jags:jags]{jags}},
#'   \code{\link{run_metareg}}, \code{\link{run_model}}
#'
#' @references
#' Dias S, Ades AE, Welton NJ, Jansen JP, Sutton AJ. Network Meta-Analysis for
#' Decision Making. Chichester (UK): Wiley; 2018.
#'
#' Dias S, Sutton AJ, Ades AE, Welton NJ. Evidence synthesis for decision
#' making 2: a generalized linear modeling framework for pairwise and network
#' meta-analysis of randomized controlled trials. \emph{Med Decis Making}
#' 2013a;\bold{33}(5):607--17. doi: 10.1177/0272989X12458724
#'
#' Dias S, Welton NJ, Sutton AJ, Ades AE. Evidence synthesis for decision
#' making 5: the baseline natural history model. \emph{Med Decis Making}
#' 2013b;\bold{33}(5):657--70. doi: 10.1177/0272989X13485155
#'
#' @export
baseline_model <- function(base_risk,
                           n_chains,
                           n_iter,
                           n_burnin,
                           n_thin) {

  base_risk0 <- if (is.vector(base_risk) &
                    !is.element(length(base_risk), c(1, 3))) {
    stop("The argument 'base_risk' must be scalar or vector of three elements.",
         call. = FALSE)
  } else if (!is.vector(base_risk) & length(base_risk) != 2) {
    stop("The argument 'base_risk' must have two columns.", call. = FALSE)
  } else if (is.element(length(base_risk), c(1, 2, 3))) {
    base_risk
  }
  base_risk1 <- if (is.element(length(base_risk0), c(1, 3)) &
                    (min(base_risk0) <= 0 || max(base_risk0) >= 1)) {
    stop("The argument 'base_risk' must be defined in (0, 1).", call. = FALSE)
  } else if (length(base_risk0) == 3 & is.unsorted(base_risk0)) {
    aa <- "must be sort in ascending order."
    stop(paste("The elements in argument 'base_risk'" , aa), call. = FALSE)
  } else if (length(base_risk0) == 1) {
    rep(log(base_risk0/(1 - base_risk0)), 2)
  } else if (length(base_risk0) == 3) {
    # First element is mean & second element is precision, both in logit scale
    c(log(base_risk0[2]/(1 - base_risk0[2])),
      (3.92/((log(base_risk0[3])/(1 - log(base_risk0[3]))) -
         (log(base_risk0[1])/(1 - log(base_risk0[1])))))^2)
  } else if (length(base_risk0) == 2) {
    base_risk0
  }
  base_type <- if (length(base_risk0) == 1) {
    "fixed"
  } else if (length(base_risk0) == 3) {
    "random"
  } else if (!is.element(length(base_risk0), c(1, 3))) {
    "predicted"
  }
  n_chains <- if (missing(n_chains)) {
    2
  } else if (n_chains < 1) {
    stop("The argument 'n_chains' must be a positive integer.", call. = FALSE)
  } else {
    n_chains
  }
  n_iter <- if (missing(n_iter)) {
    10000
  } else if (n_iter < 1) {
    stop("The argument 'n_iter' must be a positive integer.", call. = FALSE)
  } else {
    n_iter
  }
  n_burnin <- if (missing(n_burnin)) {
    1000
  } else if (n_burnin < 1) {
    stop("The argument 'n_burnin' must be a positive integer.", call. = FALSE)
  } else {
    n_burnin
  }
  n_thin <- if (missing(n_thin)) {
    1
  } else if (n_thin < 1) {
    stop("The argument 'n_thin' must be a positive integer.", call. = FALSE)
  } else {
    n_thin
  }

  if (base_type == "predicted") {
    message("**Baseline model (predictions)**")
    # Data for the baseline model
    data_jag_base <- list("r.base" = base_risk1[, 1],
                          "n.base" = base_risk1[, 2],
                          "ns.base" = length(base_risk1[, 1]))

    # Parameters to monitor (baseline model)
    param_jags_base <- c("base.risk.logit",
                         "u.base",
                         "m.base",
                         "tau.base")

    # Run the baseline model
    jagsfit_base <- jags(data = data_jag_base,
                         parameters.to.save = param_jags_base,
                         model.file = textConnection('
                         model {
                            for (i in 1:ns.base) {
                              r.base[i] ~ dbin(p.base[i], n.base[i])
                              logit(p.base[i]) <- u.base[i]
                              u.base[i] ~ dnorm(m.base, prec.base)
                            }
                            # predicted baseline risk (logit scale)
                            base.risk.logit ~ dnorm(m.base, prec.base)
                            m.base ~ dnorm(0, .0001)
                            prec.base <- pow(tau.base, -2)
                            tau.base ~ dunif(0, 5)
                         }
                                                     '),
                         n.chains = n_chains,
                         n.iter = n_iter,
                         n.burnin = n_burnin,
                         n.thin = n_thin)

    # Turn R2jags object into a data-frame
    get_results_base <- as.data.frame(t(jagsfit_base$BUGSoutput$summary))

    # Effect size of all unique pairwise comparisons
    pred_base_logit <- t(get_results_base %>%
                           dplyr::select(starts_with("base.risk.logit")))
    mean_base_logit <- t(get_results_base %>%
                           dplyr::select(starts_with("m.base")))
    tau_base_logit <- t(get_results_base %>%
                           dplyr::select(starts_with("tau.base")))
    trial_base_logit <- t(get_results_base %>%
                            dplyr::select(starts_with("u.base[")))
  } else if (base_type == "fixed") {
    message("**Fixed baseline risk assigned**")
  } else if (base_type == "random") {
    message("**Random baseline risk assigned**")
  }

  ref_base <- if (is.element(base_type, c("fixed", "random"))) {
    base_risk1
  } else if (is.element(base_type, "predicted")) {
    c(pred_base_logit[1], 1 / (pred_base_logit[2])^2)
  }

  #Draw forest-plot of observed & estimated probabilities for "predicted"
  fig <- if (is.element(base_type, "predicted")) {
    # Back-transform to probability (trial-specific estimate)
    estim_prob <- exp(trial_base_logit[, c(1, 3, 7)]) /
      (1 + exp(trial_base_logit[, c(1, 3, 7)]))
    # Back-transform to probability (summary estimate)
    summary_prob <- round(exp(mean_base_logit[c(1, 3, 7)]) /
      (1 + exp(mean_base_logit[c(1, 3, 7)])) * 100, 0)


    # Create dataset for the forest-plot
    dataplot <- data.frame(rbind(matrix(rep(base_risk1[, 1] / base_risk1[, 2], 3),
                                        ncol = 3),
                                 estim_prob),
                           rep(c("Observed", "Estimated"),
                               each = data_jag_base$ns.base),
                           rep(as.factor(seq_len(data_jag_base$ns.base)), 2))
    colnames(dataplot) <- c("point", "lower", "upper", "type", "order")
    dataplot[, 1:3] <- round(dataplot[, 1:3] * 100, 0)

    # Rule to define x-axis label tick marks
    min_x <- ifelse(min(dataplot$lower) < 50, 0, min(dataplot$lower))

    # Present summary mean and between-trial StD
    caption <- paste0("Summary mean (%):", " ", summary_prob[1], " ",
                      "(", summary_prob[2], ",", " ", summary_prob[3],
                      ") \nBetween-trial SD:", " ", round(tau_base_logit[5], 2),
                      " ", "(", round(tau_base_logit[3], 2), ",", " ",
                      round(tau_base_logit[7], 2), ")")

    # Crate forest-plot
    ggplot(data = dataplot,
           aes(x = order,
               y = point,
               ymin = lower,
               ymax = upper)) +
      geom_hline(yintercept = summary_prob,
                 col = "grey",
                 size = 1,
                 lty = 2) +
      geom_rect(aes(xmin = -Inf,
                    xmax = Inf,
                    ymin = summary_prob[2],
                    ymax = summary_prob[3]),
                alpha = 0.01,
                fill = "grey") +
      geom_linerange(size = 1.5,
                     position = position_dodge(width = 0.5)) +
      geom_point(aes(colour = type),
                 stroke = 0.3,
                 size = 2.5) +
      geom_text(aes(x = order,
                    y = point,
                    label = point),
                color = "blue",
                hjust = 0,
                vjust = -0.5,
                size = 4.0,
                check_overlap = FALSE,
                position = position_dodge(width = 0.5)) +
      geom_text(data = subset(dataplot, type == "Estimated"),
                aes(x = order,
                    y = lower,
                    label = lower),
                color = "blue",
                hjust = 0,
                vjust = -0.5,
                size = 4.0,
                check_overlap = FALSE,
                position = position_dodge(width = 0.5)) +
      geom_text(data = subset(dataplot, type == "Estimated"),
                aes(x = order,
                    y = upper,
                    label = upper),
                color = "blue",
                hjust = 0,
                vjust = -0.5,
                size = 4.0,
                check_overlap = FALSE,
                position = position_dodge(width = 0.5)) +
      scale_colour_manual(values = c("Estimated" = "#D55E00",
                                     "Observed" = "#009E73")) +
      scale_y_continuous(breaks = round(seq(min_x, max(dataplot$upper,
                                                       summary_prob[3]),
                                            length.out = 11), 0),
                         limits = c(min_x, max(dataplot$upper,
                                               summary_prob[3])),
                         expand = c(0.01, 0.01)) +
      labs(x = "Trials in the baseline model",
           y = "Probability of an event in reference intervention (%)",
           colour = "",
           fill = "",
           caption = caption) +
      coord_flip() +
      theme_classic() +
      theme(axis.text.x = element_text(color = "black", size = 12),
            axis.text.y = element_text(color = "black", size = 12),
            axis.title.x = element_text(color = "black", face = "bold",
                                        size = 12),
            legend.position = "bottom",
            legend.text = element_text(color = "black", size = 12),
            legend.title = element_text(color = "black", face = "bold",
                                        size = 12),
            plot.caption = element_text(hjust = 0, face = "bold", size = 11,
                                        color = "grey47", lineheight = 1.1))
  }

  results <- list(ref_base = ref_base,
                  figure = fig)
  #if (!is.element(base_type, c("fixed", "random"))) {
  #  results <- append(results, list(mean_base_logit = mean_base_logit,
  #                                  tau_base_logit = tau_base_logit))
  #}

  return(results)
}

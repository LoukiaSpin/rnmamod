#' The baseline model for binary outcome
#'
#' @description
#'   XXXX
#'
#' @param base_risk Character string indicating the type of baseline model.
#'   Set \code{base_type} equal to one of the following: \code{"fixed"},
#'   \code{"random"}, or \code{"predicted"}.
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
#' @return An R character vector object to be passed to \code{\link{run_model}}
#'   and \code{\link{run_metareg}} through the
#'   \code{\link[base:textConnection]{textConnection}} function as the argument
#'   \code{object}.
#'
#' @details \code{prepare_model} creates the model in the JAGS dialect
#'   of the BUGS language. The output of this function constitutes the argument
#'   \code{model.file} of the \code{\link[R2jags:jags]{jags}} function (in the
#'   R-package \href{https://CRAN.R-project.org/package=R2jags}{R2jags}) via the
#'   \code{\link[base:textConnection]{textConnection}} function.
#'
#' @author {Loukia M. Spineli}
#'
#' @seealso \code{\link{run_metareg}}, \code{\link{run_model}},
#'    \code{\link[R2jags:jags]{jags}},
#'    \code{\link[base:textConnection]{textConnection}}
#'
#' @references
#' Cooper NJ, Sutton AJ, Morris D, Ades AE, Welton NJ. Addressing between-study
#' heterogeneity and inconsistency in mixed treatment comparisons: Application
#' to stroke prevention treatments in individuals with non-rheumatic atrial
#' fibrillation. \emph{Stat Med} 2009;\bold{28}(14):1861--81.
#' doi: 10.1002/sim.3594
#'
#' Dias S, Sutton AJ, Ades AE, Welton NJ. Evidence synthesis for decision
#' making 2: a generalized linear modeling framework for pairwise and network
#' meta-analysis of randomized controlled trials. \emph{Med Decis Making}
#' 2013;\bold{33}(5):607--17. doi: 10.1177/0272989X12458724
#'
#' Spineli LM, Kalyvas C, Papadimitropoulou K. Continuous(ly) missing outcome
#' data in network meta-analysis: a one-stage pattern-mixture model approach.
#' \emph{Stat Methods Med Res} 2021;\bold{30}(4):958--75.
#' doi: 10.1177/0962280220983544
#'
#' Spineli LM. An empirical comparison of Bayesian modelling strategies for
#' missing binary outcome data in network meta-analysis.
#' \emph{BMC Med Res Methodol} 2019;\bold{19}(1):86.
#' doi: 10.1186/s12874-019-0731-y
#'
#' Turner NL, Dias S, Ades AE, Welton NJ. A Bayesian framework to account
#' for uncertainty due to missing binary outcome data in pairwise meta-analysis.
#' \emph{Stat Med} 2015;\bold{34}(12):2062--80. doi: 10.1002/sim.6475
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
    message("*The baseline model also runs (using the predictive distribution)")
    # Data for the baseline model
    data_jag_base <- list("r.base" = base_risk1[, 1],
                          "n.base" = base_risk1[, 2],
                          "ns.base" = length(base_risk1[, 1]))

    # Parameters to monitor (baseline model)
    param_jags_base <- c("base.risk.logit", "tau.base")

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
  }

  ref_base <- if (is.element(base_type, c("fixed", "random"))) {
    base_risk1
  } else if (is.element(base_type, "predicted")) {
    c(pred_base_logit[1], 1/(pred_base_logit[2])^2)
  }

  return(ref_base)
}

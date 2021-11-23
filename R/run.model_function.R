#' Perform Bayesian pairwise or network meta-analysis
#'
#' @description
#'   Performs a one-stage pairwise or network meta-analysis while addressing
#'   aggregate binary or continuous missing participant outcome data via the
#'   pattern-mixture model.
#'
#' @param data A data-frame of the one-trial-per-row format with arm-level data.
#'   See 'Format' for the specification of the columns.
#' @param measure Character string indicating the effect measure with values
#'   \code{"OR"}, \code{"MD"}, \code{"SMD"}, or \code{"ROM"} for the odds ratio,
#'   mean difference, standardised mean difference and ratio of means,
#'   respectively.
#' @param model Character string indicating the analysis model with values
#'   \code{"RE"}, or \code{"FE"} for the random-effects and fixed-effect model,
#'   respectively. The default argument is \code{"RE"}.
#' @param assumption Character string indicating the structure of the
#'   informative missingness parameter.
#'   Set \code{assumption} equal to one of the following: \code{"HIE-COMMON"},
#'   \code{"HIE-TRIAL"}, \code{"HIE-ARM"}, \code{"IDE-COMMON"},
#'   \code{"IDE-TRIAL"}, \code{"IDE-ARM"}, \code{"IND-CORR"},
#'   or \code{"IND-UNCORR"}.
#'   The default argument is \code{"IDE-ARM"}. The abbreviations \code{"IDE"},
#'   \code{"HIE"}, and \code{"IND"} stand for identical, hierarchical and
#'   independent, respectively. \code{"CORR"} and \code{"UNCORR"} stand for
#'   correlated and uncorrelated, respectively.
#' @param heter_prior A list of three elements with the following order:
#'   1) a character string indicating the distribution with
#'   (currently available) values \code{"halfnormal"}, \code{"uniform"},
#'   \code{"lognormal"}, or \code{"logt"}; 2) two numeric values that refer to
#'   the parameters of the selected distribution.  For \code{"lognormal"}, and
#'   \code{"logt"} these numbers refer to the mean and precision, respectively.
#'   For \code{"halfnorm"}, these numbers refer to zero and the scale parameter
#'   (equal to 4 or 1 being the corresponding precision of the scale parameter
#'   0.5 or 1). For \code{"uniform"}, these numbers refer to the
#'   minimum and maximum value of the distribution.
#'   See 'Details' in \code{\link{heterogeneity_param_prior}}.
#' @param mean_misspar A numeric value or a vector of two numeric values for the
#'   mean of the normal distribution of the informative missingness parameter
#'   (see 'Details'). The default argument is 0 and corresponds to the
#'   missing-at-random assumption.
#'   See also 'Details' in \code{\link{missingness_param_prior}}.
#' @param var_misspar A positive non-zero number for the variance of the
#'   normal distribution of the informative missingness parameter.
#'   When the \code{measure} is \code{"OR"}, \code{"MD"}, or \code{"SMD"}
#'   the default argument is 1. When the \code{measure} is \code{"ROM"}
#'   the default argument is 0.04
#' @param D A binary number for the direction of the outcome.
#'   Set \code{D = 1} for beneficial outcome and \code{D = 0} for harmful
#'   outcome.
#' @param n_chains Positive integer specifying the number of chains for the
#'   MCMC sampling; an argument of the \code{\link[R2jags]{jags}} function
#'   of the R-package \href{https://CRAN.R-project.org/package=R2jags}{R2jags}.
#'   The default argument is 2.
#' @param n_iter Positive integer specifying the number of Markov chains for the
#'   MCMC sampling; an argument of the \code{\link[R2jags]{jags}} function
#'   of the R-package \href{https://CRAN.R-project.org/package=R2jags}{R2jags}.
#'   The default argument is 10000.
#' @param n_burnin Positive integer specifying the number of iterations to
#'   discard at the beginning of the MCMC sampling; an argument of the
#'   \code{\link[R2jags]{jags}} function of the R-package
#'   \href{https://CRAN.R-project.org/package=R2jags}{R2jags}.
#'   The default argument is 1000.
#' @param n_thin Positive integer specifying the thinning rate for the
#'   MCMC sampling; an argument of the \code{\link[R2jags]{jags}} function
#'   of the R-package \href{https://CRAN.R-project.org/package=R2jags}{R2jags}.
#'   The default argument is 1.
#'
#' @format The columns of the data-frame in the argument \code{data} refer
#'   to the following elements for a continuous outcome:
#'   \tabular{ll}{
#'    \strong{t} \tab An intervention identifier in each arm.\cr
#'    \tab \cr
#'    \strong{y} \tab The observed mean value of the outcome in each arm.\cr
#'    \tab \cr
#'    \strong{sd} \tab The observed standard deviation of the outcome in
#'    each arm.\cr
#'    \tab \cr
#'    \strong{m} \tab The number of missing participant outcome data in
#'    each arm.\cr
#'    \tab \cr
#'    \strong{n} \tab The number of randomised participants in each arm.\cr
#'   }
#'
#'   For a binary outcome, the columns of the data-frame in the argument
#'   \code{data} refer to the following elements:
#'   \tabular{ll}{
#'    \strong{t} \tab An intervention identifier in each arm.\cr
#'    \tab \cr
#'    \strong{r} \tab The observed number of events of the outcome in
#'    each arm.\cr
#'    \tab \cr
#'    \strong{m} \tab The number of missing participant outcome data in
#'    each arm.\cr
#'    \tab \cr
#'    \strong{n} \tab The number of randomised participants in each arm.\cr
#'   }
#'   The number of rows in \code{data} equals the number of collected trials.
#'   Each element appears in \code{data} as many times as the maximum number of
#'   interventions compared in a trial of the dataset.
#'   In pairwise meta-analysis (PMA), the maximum number of arms is inherently
#'   two. The same holds for a network meta-analysis (NMA) without multi-arm
#'   trials.
#'   In the case of NMA with multi-arm trials, the maximum number of arms
#'   exceeds two. See 'Examples' that illustrates the structure of \code{data}
#'   for a network with a maximum number of four arms.
#'   It is not a prerequisite of \code{run_model} that the multi-arm trials
#'   appear at the bottom of the dataset.
#'
#' @return A list of R2jags output on the summaries of the posterior
#'   distribution, and the Gelman-Rubin convergence diagnostic
#'   (Gelman et al., 1992) of the following monitored parameters for a
#'   fixed-effect PMA:
#'   \tabular{ll}{
#'    \code{EM} \tab The estimated summary effect measure (according to the
#'    argument \code{measure}).\cr
#'    \tab \cr
#'    \code{dev_o} \tab The deviance contribution of each trial-arm based
#'    on the observed outcome.\cr
#'    \tab \cr
#'    \code{hat_par} \tab The fitted outcome at each trial-arm.\cr
#'    \tab \cr
#'    \code{phi} \tab The informative missingness parameter.\cr
#'   }
#'
#'     For a fixed-effect NMA, the output additionally includes:
#'   \tabular{ll}{
#'    \code{EM_ref} \tab The estimated summary effect measure
#'    (according to the argument \code{measure}) of all comparisons
#'    with the reference intervention.\cr
#'    \tab \cr
#'    \code{SUCRA} \tab The surface under the cumulative ranking curve
#'    for each intervention.\cr
#'    \tab \cr
#'    \code{effectiveneness} \tab The ranking probability of each intervention
#'     for every rank.\cr
#'   }
#'
#'   For a random-effects PMA, the output additionally includes the
#'   following elements:
#'   \tabular{ll}{
#'    \code{EM_pred} \tab The predicted summary effect measure
#'    (according to the argument \code{measure}).\cr
#'    \tab \cr
#'    \code{delta} \tab The estimated trial-specific effect measure
#'    (according to the argument \code{measure}).\cr
#'    \tab \cr
#'    \code{tau} \tab The between-trial standard deviation.\cr
#'   }
#'
#'   For a random-effects NMA, the output additionally includes:
#'   \tabular{ll}{
#'    \code{pred_ref} \tab The predicted summary effect measure
#'    (according to the argument \code{measure}) of all comparisons
#'    with the reference intervention.\cr
#'   }
#'   In NMA, \code{EM} and \code{EM_pred} refer to all possible pairwise
#'   comparisons of interventions in the network. Furthermore, \code{tau} is
#'   typically assumed to be common for all observed comparisons in the network.
#'   For a multi-arm trial, we estimate a total of \emph{T-1} \code{delta} for
#'   comparisons with the baseline intervention of the trial (found in the first
#'   column of the element \bold{t}), with \emph{T} being the number of
#'   interventions in the trial.
#'
#'   Furthermore, the output includes the following elements:
#'   \tabular{ll}{
#'    \code{leverage_o} \tab The leverage for the observed outcome
#'    at each trial-arm.\cr
#'    \tab \cr
#'    \code{sign_dev_o} \tab The sign of the difference between
#'    observed and fitted outcome at each trial-arm.\cr
#'    \tab \cr
#'    \code{model_assessment} \tab A data-frame on the measures of
#'    model assessment: deviance information criterion,
#'    number of effective parameters, and total residual deviance.\cr
#'    \tab \cr
#'    \code{jagsfit} \tab An object of S3 class \code{\link[R2jags]{jags}}
#'    with the posterior results on all monitored parameters to be used
#'    in the \code{\link{mcmc_diagnostics}} function.\cr
#'   }
#'   The \code{run_model} function also returns the arguments \code{data},
#'   \code{measure}, \code{model}, \code{assumption}, \code{heter_prior},
#'   \code{mean_misspar}, \code{var_misspar}, and \code{D} as specified by the
#'   user to be considered in other functions of the package.
#'
#' @details The model runs in \code{JAGS} and the progress of the simulation
#'   appears on the R console. The output of \code{run_model} is used as an S3
#'   object by other functions of the package to be processed further and
#'   provide an end-user-ready output.
#'
#'   The \code{\link{data_preparation}} function is called to prepare the data
#'   for the Bayesian analysis. \code{\link{data_preparation}} checks whether
#'   the element \strong{m} exists in the \code{data}. If this element is
#'   missing, \code{\link{data_preparation}} creates a pseudo-data-frame for
#'   \strong{m} that has the zero value for the observed trial-arms, and
#'   \code{NA} for the unobserved trial-arms, and the pseudo-data-frame \code{I}
#'   that is identical with the pseudo-data-frame for \code{m}. If the element
#'   \strong{m} exists in the \code{data} and has values only for some
#'   trial-arms, the pseudo-data-frame for \strong{m} is identical to \strong{m}
#'   for the corresponding trial-arms, and the pseudo-data-frame \code{I} has
#'   the value one for these trial-arms. Both pseudo-data-frames aim to retain
#'   the trials without information on missing participant outcome data.
#'
#'   Furthermore, \code{\link{data_preparation}} sorts the interventions across
#'   the arms of each trial in an ascending order and correspondingly the
#'   remaining elements in \code{data} (see 'Format').
#'   \code{\link{data_preparation}} considers the
#'   first column in \strong{t} as being the control arm for every trial. Thus,
#'   this sorting ensures that interventions with a lower identifier are
#'   consistently treated as the control arm in each trial. This case is
#'   relevant in non-star-shaped networks. By default, \code{run_model} treats
#'   the intervention with identifier equal to one as the reference intervention
#'   of the network.
#'
#'   To perform a Bayesian PMA or NMA, the \code{\link{prepare_model}} function
#'   is called which contains the WinBUGS code as written by Dias et al., (2013)
#'   for binomial and normal likelihood to analyse binary and continuous data,
#'   respectively. \code{\link{prepare_model}} uses the consistency model (as
#'   described in Lu and Ades (2006)) to estimate all possible comparisons in
#'   the network.
#'   It also accounts for the multi-arm trials by assigning conditional
#'   univariate normal distributions on the basic parameters of these trials,
#'   namely, effect parameters between the non-baseline arms and the baseline
#'   arm of the multi-arm trial (Dias et al., 2013).
#'
#'   The code of Dias et al., (2013) has been extended to incorporate the
#'   pattern-mixture model to adjust the underlying outcome in each arm of
#'   every trial for missing participant outcome data (Turner et al., 2015;
#'   Spineli, 2019a; Spineli et al., 2021). The assumptions about the
#'   missingness parameter are specified using the arguments \code{mean_misspar}
#'   and \code{var_misspar}. Specifically, \code{run_model} considers the
#'   informative missingness odds ratio in the logarithmic scale for binary
#'   outcome data (White et al., 2008; Turner et al., 2015; Spineli, 2019a), the
#'   informative missingness difference of means when \code{measure} is
#'   \code{"MD"} or \code{"SMD"}, and the informative missingness ratio of means
#'   in the logarithmic scale when \code{measure} is \code{"ROM"}
#'   (Mavridis et al., 2015; Spineli et al., 2021).
#'
#'   When \code{assumption} is trial-specific (i.e., \code{"IDE-TRIAL"} or
#'   \code{"HIE-TRIAL"}), or independent (i.e., \code{"IND-CORR"} or
#'   \code{"IND-UNCORR"}), only one numeric value can be assigned to
#'   \code{mean_misspar} because the same missingness scenario is applied to all
#'   trials and trial-arms of the dataset, respectively. When \code{assumption}
#'   is \code{"IDE-ARM"} or \code{"HIE-ARM"}, a maximum of two
#'   \emph{different} or \emph{identical} numeric values can be assigned as a
#'   vector to \code{mean_misspars}: the first value refers to the experimental
#'   arm, and the second value refers to the control arm of a trial.
#'   In the case of a network, the first value is considered for all
#'   non-reference interventions and the second value is considered for the
#'   reference intervention of the network (i.e., the intervention with
#'   identifier equal to one). This is necessary to ensure transitivity in the
#'   assumptions for the missingness parameter across the network (Spineli,
#'   2019b).
#'
#'   Currently, there are no empirically-based prior distributions for the
#'   informative missingness parameters. The user may refer to
#'   White et al., (2008); Mavridis et al., (2015); Turner et al., (2015) and
#'   Spineli (2019) to determine \code{mean_misspar} and select a proper value
#'   for \code{var_misspar}.
#'
#' @author {Loukia M. Spineli}
#'
#' @seealso \code{\link{data_preparation}},
#'   \code{\link{heterogeneity_param_prior}}, \code{\link[R2jags]{jags}}
#'   \code{\link{missingness_param_prior}}, \code{\link{prepare_model}}
#'
#' @references
#' Spineli LM, Kalyvas C, Papadimitropoulou K. Continuous(ly) missing outcome
#' data in network meta-analysis: a one-stage pattern-mixture model approach.
#' \emph{Stat Methods Med Res} 2021. \doi{10.1177/0962280220983544}
#'
#' Spineli LM. An empirical comparison of Bayesian modelling strategies for
#' missing binary outcome data in network meta-analysis.
#' \emph{BMC Med Res Methodol} 2019a;\bold{19}(1):86.
#' \doi{10.1186/s12874-019-0731-y}
#'
#' Spineli LM. Modeling missing binary outcome data while preserving
#' transitivity assumption yielded more credible network meta-analysis
#' results. \emph{J Clin Epidemiol} 2019b;\bold{105}:19--26.
#'
#' Mavridis D, White IR, Higgins JP, Cipriani A, Salanti G. Allowing for
#' uncertainty due to missing continuous outcome data in pairwise and
#' network meta-analysis. \emph{Stat Med} 2015;\bold{34}(5):721--741.
#' \doi{10.1002/sim.6365}
#'
#' Turner NL, Dias S, Ades AE, Welton NJ. A Bayesian framework to account
#' for uncertainty due to missing binary outcome data in pairwise
#' meta-analysis. \emph{Stat Med} 2015;\bold{34}(12):2062--2080.
#' \doi{10.1002/sim.6475}
#'
#' Dias S, Sutton AJ, Ades AE, Welton NJ. Evidence synthesis for decision
#' making 2: a generalized linear modeling framework for pairwise and network
#' meta-analysis of randomized controlled trials. \emph{Med Decis Making}
#' 2013;\bold{33}(5):607--617. \doi{10.1177/0272989X12458724}
#'
#' Cooper NJ, Sutton AJ, Morris D, Ades AE, Welton NJ. Addressing
#' between-study heterogeneity and inconsistency in mixed treatment
#' comparisons: Application to stroke prevention treatments in individuals
#' with non-rheumatic atrial fibrillation.
#' \emph{Stat Med} 2009;\bold{28}(14):1861--81. \doi{10.1002/sim.3594}
#'
#' White IR, Higgins JP, Wood AM. Allowing for uncertainty due to missing data
#' in meta-analysis--part 1: two-stage methods. \emph{Stat Med}
#' 2008;\bold{27}(5):711--727. \doi{10.1002/sim.3008}
#'
#' Lu G, Ades AE. Assessing evidence inconsistency in mixed treatment
#' comparisons. \emph{J Am Stat Assoc} 2006;\bold{101}:447--459.
#' \doi{10.1198/016214505000001302}
#'
#' Gelman A, Rubin DB. Inference from iterative simulation using multiple
#' sequences. \emph{Stat Sci} 1992;\bold{7}:457--472.
#'
#' @examples
#' data("nma.baker2009")
#'
#' # Show the first six trials of the dataset
#' head(nma.baker2009)
#'
#' \donttest{
#' # Perform a random-effects network meta-analysis
#' # Note: Ideally, set 'n_iter' to 10000 and 'n_burnin' to 1000
#' run_model(data = nma.baker2009,
#'           measure = "OR",
#'           model = "RE",
#'           assumption = "IDE-ARM",
#'           heter_prior = list("halfnormal", 0, 1),
#'           mean_misspar = c(0, 0),
#'           var_misspar = 1,
#'           D = 0,
#'           n_chains = 3,
#'           n_iter = 1000,
#'           n_burnin = 100,
#'           n_thin = 1)
#' }
#'
#' @export
run_model <- function(data,
                      measure,
                      model,
                      assumption,
                      heter_prior,
                      mean_misspar,
                      var_misspar,
                      D,
                      n_chains,
                      n_iter,
                      n_burnin,
                      n_thin) {


  # Prepare the dataset for the R2jags
  item <- data_preparation(data, measure)

  # Missing and default arguments
  model <- if (missing(model)) {
    "RE"
  } else if (!is.element(model, c("RE", "FE"))) {
    stop("Insert 'RE', or 'FE'", call. = FALSE)
  } else {
    model
  }
  assumption <- if (missing(assumption)) {
    "IDE-ARM"
  } else {
    assumption
  }
  D <- if (missing(D)) {
    stop("The argument 'D' needs to be defined", call. = FALSE)
  } else {
    D
  }
  mean_misspar <- missingness_param_prior(assumption, mean_misspar)
  heterog_prior <- heterogeneity_param_prior(measure, model, heter_prior)
  var_misspar <- if (missing(var_misspar) &
                     is.element(measure, c("OR", "MD", "SMD"))) {
    1
  } else if (missing(var_misspar) & measure == "ROM") {
    0.2^2
  } else {
    var_misspar
  }
  n_chains <- ifelse(missing(n_chains), 2, n_chains)
  n_iter <- ifelse(missing(n_iter), 10000, n_iter)
  n_burnin <- ifelse(missing(n_burnin), 1000, n_burnin)
  n_thin <- ifelse(missing(n_thin), 1, n_thin)

  # Data in list format for R2jags
  data_jag <- list("m" = item$m,
                   "N" = item$N,
                   "t" = item$t,
                   "na" = item$na,
                   "nt" = item$nt,
                   "ns" = item$ns,
                   "ref" = item$ref,
                   "I" = item$I,
                   "meand.phi" = mean_misspar,
                   "precd.phi" = 1 / var_misspar,
                   "D" = D,
                   "heter.prior" = heterog_prior)

  data_jag <- if (is.element(measure, c("MD", "SMD", "ROM"))) {
    append(data_jag, list("y.o" = item$y0, "se.o" = item$se0))
  } else if (measure == "OR") {
    append(data_jag, list("r" = item$r))
  }

  data_jag <- if (is.element(assumption, c("IND-CORR", "IND-UNCORR"))) {
    append(data_jag, list("M" = ifelse(!is.na(item$m), mean_misspar, NA),
                          "cov.phi" = 0.5 * var_misspar,
                          "var.phi" = var_misspar,))
  } else {
    data_jag
  }

  param_jags <- c("delta",
                  "EM",
                  "EM.ref",
                  "EM.pred",
                  "pred.ref",
                  "tau",
                  "SUCRA",
                  "effectiveness",
                  "dev.o",
                  "totresdev.o",
                  "hat.par")

  param_jags <- if (is.element(assumption,
                               c("HIE-COMMON", "HIE-TRIAL", "HIE-ARM"))) {
    append(param_jags, "mean.phi")
  } else {
    append(param_jags, "phi")
  }

  param_jags <- if (model == "RE") {
    param_jags
  } else {
    param_jags[!is.element(param_jags,
                           c("EM.pred", "pred.ref", "tau", "delta"))]
  }

  # Run the Bayesian analysis
  jagsfit <- jags(data = data_jag,
                  parameters.to.save = param_jags,
                  model.file = textConnection(
                    prepare_model(measure,
                                  model,
                                  covar_assumption = "NO",
                                  assumption)
                    ),
                  n.chains = n_chains,
                  n.iter = n_iter,
                  n.burnin = n_burnin,
                  n.thin = n_thin)

  # Turn R2jags object into a data-frame
  get_results <- as.data.frame(t(jagsfit$BUGSoutput$summary))

  # Effect size of all unique pairwise comparisons
  EM <- t(get_results %>% dplyr::select(starts_with("EM[")))

  # Effect size of all comparisons with the reference intervention
  EM_ref <- t(get_results %>% dplyr::select(starts_with("EM.ref[")))

  # Predictive effects of all unique pairwise comparisons
  EM_pred <- t(get_results %>% dplyr::select(starts_with("EM.pred[")))

  # Predictive effects of all comparisons with the reference intervention
  pred_ref <- t(get_results %>% dplyr::select(starts_with("pred.ref[")))

  # Between-trial standard deviation
  tau <- t(get_results %>% dplyr::select(starts_with("tau")))

  # SUrface under the Cumulative RAnking curve values
  SUCRA <- t(get_results %>% dplyr::select(starts_with("SUCRA")))

  # Within-trial effects size
  delta <- t(get_results %>% dplyr::select(starts_with("delta") &
                                             !ends_with(",1]")))

  # Ranking probability of each intervention for every rank
  effectiveness <- t(get_results %>% dplyr::select(
    starts_with("effectiveness")))

  # Estimated missingness parameter
  phi <- if (length(unique(unlist(item$m))) > 2) {
    t(get_results %>% dplyr::select(starts_with("phi") |
                                     starts_with("mean.phi") |
                                     starts_with("mean.phi[") |
                                     starts_with("phi[")))
  } else {
    NA
  }

  # Trial-arm deviance contribution for observed outcome
  dev_o <- t(get_results %>% dplyr::select(starts_with("dev.o")))

  # Fitted/predicted outcome
  hat_par <- t(get_results %>% dplyr::select(starts_with("hat.par")))

  # Total residual deviance
  dev <- jagsfit$BUGSoutput$summary["totresdev.o", "mean"]

  # Calculate the deviance at posterior mean of fitted values
  # Turn 'N' and 'm' into a vector (first column, followed by second, etc)
  m_new <- suppressMessages({
    as.vector(na.omit(melt(item$m)[, 2]))
    })
  N_new <- suppressMessages({
    as.vector(na.omit(melt(item$N)[, 2]))
    })
  obs <- N_new - m_new

  if (is.element(measure, c("MD", "SMD", "ROM"))) {
    # Turn 'y0', 'se0' into a vector as above
    y0_new <- suppressMessages({
      as.vector(na.omit(melt(item$y0)[, 2]))
      })
    se0_new <- suppressMessages({
      as.vector(na.omit(melt(item$se0)[, 2]))
      })

    # Deviance at the posterior mean of the fitted mean outcome
    dev_post_o <- (y0_new -
                     as.vector(hat_par[, 1])) *
      (y0_new - as.vector(hat_par[, 1])) * (1 / se0_new^2)

    # Sign of the difference between observed and fitted mean outcome
    sign_dev_o <- sign(y0_new - as.vector(hat_par[, 1]))
  } else {
    # Turn 'r' and number of observed into a vector as above
    r_new <- suppressMessages({
      as.vector(na.omit(melt(item$r)[, 2]))
      })

    # Correction for zero events in trial-arm
    r0 <- ifelse(r_new == 0, r_new + 0.01,
                 ifelse(r_new == obs, r_new - 0.01, r_new))

    # Deviance at the posterior mean of the fitted response
    dev_post_o <- 2 * (r0 * (log(r0) -
                               log(as.vector(hat_par[, 1]))) +
                         (obs - r0) * (log(obs - r0) -
                                         log(obs - as.vector(hat_par[, 1]))))

    # Sign of the difference between observed and fitted response
    sign_dev_o <- sign(r0 - as.vector(hat_par[, 1]))
  }

  # Obtain the leverage for observed outcomes
  leverage_o <- as.vector(dev_o[, 1]) - dev_post_o

  # Number of effective parameters
  pD <- dev - sum(dev_post_o)

  # Deviance information criterion
  DIC <- pD + dev

  # A data-frame on the measures of model assessment:
  # DIC, pD, and total residual deviance
  model_assessment <- data.frame(DIC, pD, dev)

  # Return a list of results
  if (model == "RE") {
    ma_results <- list(EM = EM,
                       EM_pred = EM_pred,
                       tau = tau,
                       delta = delta,
                       dev_o = dev_o,
                       hat_par = hat_par,
                       leverage_o = leverage_o,
                       sign_dev_o = sign_dev_o,
                       phi = phi,
                       model_assessment = model_assessment,
                       data = data,
                       measure = measure,
                       model = model,
                       assumption = assumption,
                       heter_prior = heterog_prior,
                       mean_misspar = mean_misspar,
                       var_misspar = var_misspar,
                       D = D,
                       jagsfit = jagsfit)
    nma_results <- append(ma_results, list(EM_ref = EM_ref,
                                           pred_ref = pred_ref,
                                           SUCRA = SUCRA,
                                           effectiveness = effectiveness))
  } else {
    ma_results <- list(EM = EM,
                       dev_o = dev_o,
                       hat_par = hat_par,
                       leverage_o = leverage_o,
                       sign_dev_o = sign_dev_o,
                       phi = phi,
                       model_assessment = model_assessment,
                       data = data,
                       measure = measure,
                       model = model,
                       assumption = assumption,
                       mean_misspar = mean_misspar,
                       var_misspar = var_misspar,
                       D = D,
                       jagsfit = jagsfit)
    nma_results <- append(ma_results, list(EM_ref = EM_ref,
                                           SUCRA = SUCRA,
                                           effectiveness = effectiveness))
  }

  ifelse(item$nt > 2, return(nma_results), return(ma_results))
}

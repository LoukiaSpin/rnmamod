#' Pairwise or network meta-analysis with missing participant outcome data
#'
#' @description This function performs a one-stage pairwise or network meta-analysis while addressing aggregate binary or continuous missing participant outcome data via the pattern-mixture model.
#'
#' @param data A data-frame of the one-trial-per-row format containing arm-level data of each trial. This format is widely used for BUGS models. See 'Format' for the specification of the columns.
#' @param measure Character string indicating the effect measure with values \code{"OR"}, \code{"MD"}, \code{"SMD"}, or \code{"ROM"} for the odds ratio, mean difference,
#'   standardised mean difference and ratio of means, respectively.
#' @param model Character string indicating the analysis model with values \code{"RE"}, or \code{"FE"} for the random-effects and fixed-effect model, respectively. The default argument is \code{"RE"}.
#' @param covar.assumption Character string indicating the structure of the slope for the intervention by covariate interaction, as described in Cooper et al., (2009).
#'  Set \code{covar.assumption} equal to one of the following: \code{"NO"}, when no meta-regression is performed; otherwise, \code{"exchangeable"} \code{"independent"}, and \code{"common"}.
#'  See the \code{ru.metareg} function.
#' @param assumption Character string indicating the structure of the informative missingness parameter.
#'   Set \code{assumption} equal to one of the following: \code{"HIE-COMMON"}, \code{"HIE-TRIAL"}, \code{"HIE-ARM"}, \code{"IDE-COMMON"}, \code{"IDE-TRIAL"}, \code{"IDE-ARM"}, \code{"IND-CORR"}, or \code{"IND-UNCORR"}.
#'   The default argument is \code{"IDE-ARM"}. The abbreviations \code{"IDE"}, \code{"HIE"}, and \code{"IND"} stand for identical, hierarchical and independent, respectively. \code{"CORR"} and \code{"UNCORR"} stand for correlated and uncorrelated, respectively.
#' @param heter.prior A list of three elements with the following order: 1) a character string indicating the distribution with (currently available) values \code{"halfnormal"},
#'   \code{"uniform"}, \code{"lognormal"}, or \code{"logt"}; 2) two numeric values that refer to the parameters of the selected distribution. For \code{"halfnormal"}, \code{"lognormal"}, and \code{"logt"}
#'   these numbers refer to the mean and precision, respectively. For \code{"uniform"}, these numbers refer to the minimum and maximum value of the distribution.
#'   See the function \code{\link{heterogeneity.param.prior}}.
#' @param mean.misspar A numeric value or a vector of two numeric values for the mean of the normal distribution of the informative missingness parameter (see 'Details'). The default argument is 0 and corresponds to the missing-at-random assumption.
#'   See the function \code{\link{missingness.param.prior}}.
#' @param var.misspar A positive non-zero number for the variance of the normal distribution of the informative missingness parameter. When the \code{measure} is \code{"OR"}, \code{"MD"}, or \code{"SMD"}
#'   the default argument is 1; When the \code{measure} is \code{"ROM"} the default argument is 0.04
#' @param D A binary number for the direction of the outcome. Set \code{D = 1} for a beneficial outcome and \code{D = 0} for a harmful outcome.
#' @param n.chains Positive integer specifying the number of chains for the MCMC sampling; an argument of the \code{\link[R2jags]{jags}} function of the R-package \href{https://CRAN.R-project.org/package=R2jags}{R2jags}.
#'   The default argument is 2.
#' @param n.iter Positive integer specifying the number of Markov chains for the MCMC sampling; an argument of the \code{\link[R2jags]{jags}} function of the R-package \href{https://CRAN.R-project.org/package=R2jags}{R2jags}.
#'   The default argument is 10000.
#' @param n.burnin Positive integer specifying the number of iterations to discard at the beginning of the MCMC sampling; an argument of the \code{\link[R2jags]{jags}} function of the R-package \href{https://CRAN.R-project.org/package=R2jags}{R2jags}.
#'   The default argument is 1000.
#' @param n.thin Positive integer specifying the thinning rate for the MCMC sampling; an argument of the \code{\link[R2jags]{jags}} function of the R-package \href{https://CRAN.R-project.org/package=R2jags}{R2jags}.
#'   The default argument is 1.
#'
#' @format The columns of the data frame in the argument \code{data} refer to the following ordered elements for a continuous outcome:
#' \tabular{ll}{
#'  \strong{t} \tab An intervention identifier in each arm.\cr
#'  \tab \cr
#'  \strong{y} \tab The observed mean value of the outcome in each arm.\cr
#'  \tab \cr
#'  \strong{sd} \tab The observed standard deviation of the outcome in each arm.\cr
#'  \tab \cr
#'  \strong{m} \tab The number of missing participant outcome data (MOD) in each arm.\cr
#'  \tab \cr
#'  \strong{n} \tab The number of participants randomised on the assigned intervention in each arm.\cr
#' }
#'
#' For a binary outcome, the columns of the data-frame in the argument \code{data} refer to the following ordered elements:
#' \tabular{ll}{
#'  \strong{t} \tab An intervention identifier in each arm.\cr
#'  \tab \cr
#'  \strong{r} \tab The observed number of events of the outcome in each arm.\cr
#'  \tab \cr
#'  \strong{m} \tab The number of MOD in each arm.\cr
#'  \tab \cr
#'  \strong{n} \tab The number of participants randomised on the assigned intervention in each arm.\cr
#' }
#' The number of rows in \code{data} equals the number of collected trials. Each element appears in \code{data} as many times as the maximum number of interventions compared in a trial of the dataset. In pairwise meta-analysis (PMA), the maximum number of arms is inherently two.
#' The same holds for NMA without multi-arm trials. In the case of NMA with multi-arm trials, the maximum number of arms exceeds two. See 'Examples' that illustrates the structure of \code{data} for a network with a maximum number of four arms.
#' In this example, the multi-arm trials appear at the bottom of the dataset; however, this structure is not necessary prerequisite \code{run.model}.
#'
#' @return A list of R2jags outputs on the summaries of the posterior distribution, and the Gelman-Rubin convergence diagnostic (Gelman et al., 1992) of the following monitored parameters for a fixed-effect PMA:
#' \tabular{ll}{
#'  \code{EM} \tab The estimated summary effect measure (according to the argument \code{measure}).\cr
#'  \tab \cr
#'  \code{dev.o} \tab The deviance contribution of each trial-arm based on the observed outcome.\cr
#'  \tab \cr
#'  \code{hat.par} \tab The fitted outcome at each trial-arm.\cr
#'  \tab \cr
#'  \code{phi} \tab The informative missingness parameter.\cr
#' }
#'
#' For a random-effects PMA, the output additionally includes the following elements:
#' \tabular{ll}{
#'  \code{EM.pred} \tab The predicted summary effect measure (according to the argument \code{measure}).\cr
#'  \tab \cr
#'  \code{delta} \tab The estimated trial-specific effect measure (according to the argument \code{measure}).
#'  For a multi-arm trial, we estimate \emph{T-1} effects, where \emph{T} is the number of interventions in the trial.\cr
#'  \tab \cr
#'  \code{tau} \tab The between-trial standard deviation.\cr
#' }
#'
#' For a random-effects NMA, the output additionally includes:
#' \tabular{ll}{
#'  \code{EM.ref} \tab The estimated summary effect measure (according to the argument \code{measure}) of all comparisons with the reference intervention.\cr
#'  \tab \cr
#'  \code{pred.ref} \tab The predicted summary effect measure (according to the argument \code{measure}) of all comparisons with the reference intervention.\cr
#'  \tab \cr
#'  \code{SUCRA} \tab The surface under the cumulative ranking curve for each intervention.\cr
#'  \tab \cr
#'  \code{effectiveneness} \tab The ranking probability of each intervention for every rank.\cr
#' }
#' In NMA, \code{EM} and \code{EM.pred} refer to all possible pairwise comparisons of interventions in the network. Furthermore, \code{tau} is typically assumed to be common for all observed comparisons in the network.
#'
#' Furthermore, the output includes the following elements - the first three resulting from relevant monitored parameters:
#' \tabular{ll}{
#'  \code{leverage.o} \tab The leverage for the observed outcome at each trial-arm.\cr
#'  \tab \cr
#'  \code{sign.dev.o} \tab The sign of the difference between observed and fitted outcome at each trial-arm.\cr
#'  \tab \cr
#'  \code{model.assessment} \tab A data-frame on the measures of model assessment: deviance information criterion, number of effective parameters, and total residual deviance.\cr
#'  \tab \cr
#'  \code{jagsfit} \tab An object of S3 class \code{\link[R2jags]{jags}} with the posterior results on all monitored parameters to be used in the \code{mcmc.diagnostics} function.\cr
#' }
#' The \code{run.model} function also returns the arguments \code{data}, \code{measure}, \code{model}, \code{assumption}, \code{heter.prior}, \code{mean.misspar}, \code{var.misspar}, and \code{D}
#' as specified by the user to be considered in other functions of the package.
#'
#' @details The model as specified by the arguments of \code{run.model} runs in \code{JAGS} and the progress of the simulation appears in the R console.
#'   The output of \code{run.model} is used as an S3 object by other functions of the package function to be processed further and provide an end-user-ready output.
#'
#'   The \code{data.preparation} function is called to prepare the data for the Bayesian analysis. The data preparation includes the following actions. First, \code{data.preparation} checks
#'   whether the element \strong{m} exists in \code{data}. If this element is missing, the \code{data.preparation} function creates a pseudo-data-frame for \code{m} that has the zero value for the observed trial-arms, and \code{NA} for the unobserved trial-arms,
#'   and the pseudo-data-frame \code{I} that has the same values with the pseudo-data-frame for \code{m}. If the element \strong{m} exists in \code{data} and has values only for some trials, the pseudo-data-frame for \code{m} has the same values with
#'   \code{m} for the corresponding trial-arms, and the pseudo-data-frame \code{I} has the value one for these trial-arms. Both pseudo-data-frames aim to retain the trials without information on MOD in \code{data} and allow running the Bayesian analysis via the \code{prepare.model}
#'   function when MOD have not been extracted for any trial of the dataset. Second, the \code{data.preparation} function sorts the interventions across the arms of each trial in an ascending order and correspondingly sorts within each trial the remaining elements in \code{data} (see 'Format').
#'   Since the \code{prepare.model} function considers the first column in \strong{t} to include the control arm for every trial, this sorting ensures that interventions with a lower identifier are consistently treated as the control arm in each trial. This case is relevant in
#'   non-star-shaped networks. By default, \code{run.model} function treats the intervention with identifier equal to one as the reference intervention of the network.
#'
#'   To perform the Bayesian PMA or NMA, the \code{prepare.model} function is called which contains the WinBUGS code as written by Dias et al. (2013) for binomial and normal likelihood to analyse binary and continuous outcome data, respectively. \code{prepare.model} uses consistency models (as
#'   described in Lu and Ades (2006)) to estimate all possible comparisons in the network. It also accounts for the multi-arm trials by assigning conditional univariate normal distributions on the basic parameters of these trials, that is, effect parameters between non-baseline and baseline arms (Dias et al., 2013).
#'
#'   The code of Dias et al. (2013) has been extended to incorporate the pattern-mixture model to adjust the underlying outcome in each arm of every trial for MOD (Turner et al., 2015; Spineli et al., 2021). Assumptions about the missingness parameter is specified using the arguments \code{mean.misspar} and \code{var.misspar}.
#'   Specifically, \code{run.model} considers the informative missingness odds ratio in the logarithmic scale for binary outcome data (White et al., 2008; Turner et al., 2015; Spineli, 2019), the informative missingness difference of means when \code{measure} is \code{"MD"} or \code{"SMD"},
#'   and the informative missingness ratio of means in the logarithmic scale when \code{measure} is \code{"ROM"} (Mavridis et al., 2015; Spineli et al., 2021).
#'
#'   When \code{assumption} is trial-specific (i.e., \code{"IDE-TRIAL"} or \code{"HIE-TRIAL"}), or independent (i.e., \code{"IND-CORR"} or \code{"IND-UNCORR"}),
#'   only one numeric value can be assigned to \code{mean.misspar} as the same missingness scenario is applied to all trials and trial-arms of the dataset, respectively. When \code{assumption} is \code{"IDE-ARM"} or \code{"HIE-ARM"}, a maximum of two
#'   \emph{different} numeric values can be assigned as a vector to \code{mean.misspars}: the first value refers to the experimental arm, and the second value refers to the control arm of each trial. In the case of a network, the first value is considered for all non-reference interventions
#'   and the second value is considered for the reference intervention of the network (i.e., the intervention with identifier equal to one). This is necessary to ensure transitivity in the assumed values for the missingness parameter across the comparisons in the network (Spineli, 2019).
#'   When one numeric value is considered, the same missingness scenario is applied to all interventions in the dataset.
#'
#'   Currently, there are no empirically-based prior distributions for the informative missingness parameters. The users may refer to White et al. (2008), Mavridis et al. (2015), Turner et al. (2015) and Spineli (2019) to determine \code{mean.misspar} for an informative missingness mechanism and select a proper value for \code{var.misspar}.
#'
#' @author {Loukia M. Spineli}
#'
#' @seealso \code{\link{data.preparation}}, \code{\link{prepare.model}}, \code{\link{missingness.param.prior}}, \code{\link{heterogeneity.param.prior}}, \code{\link[R2jags]{jags}}
#'
#' @references
#' Spineli LM, Kalyvas C, Papadimitropoulou K. Continuous(ly) missing outcome data in network meta-analysis: a one-stage pattern-mixture model approach. \emph{Stat Methods Med Res} 2021. [\doi{10.1177/0962280220983544}]
#'
#' Spineli LM. An empirical comparison of Bayesian modelling strategies for missing binary outcome data in network meta-analysis. \emph{BMC Med Res Methodol} 2019;\bold{19}(1):86. [\doi{10.1186/s12874-019-0731-y}]
#'
#' Spineli LM. Modeling missing binary outcome data while preserving transitivity assumption yielded more credible network meta-analysis results. \emph{J Clin Epidemiol} 2019;\bold{105}:19--26. [\doi{10.1016/j.jclinepi.2018.09.002}]
#'
#' Mavridis D, White IR, Higgins JP, Cipriani A, Salanti G. Allowing for uncertainty due to missing continuous outcome data in pairwise and network meta-analysis. \emph{Stat Med} 2015;\bold{34}(5):721--741. [\doi{10.1002/sim.6365}]
#'
#' Turner NL, Dias S, Ades AE, Welton NJ. A Bayesian framework to account for uncertainty due to missing binary outcome data in pairwise meta-analysis. \emph{Stat Med} 2015;\bold{34}(12):2062--2080. [\doi{10.1002/sim.6475}]
#'
#' Dias S, Sutton AJ, Ades AE, Welton NJ. Evidence synthesis for decision making 2: a generalized linear modeling framework for pairwise and network meta-analysis of randomized controlled trials. \emph{Med Decis Making} 2013;\bold{33}(5):607--617. [\doi{10.1177/0272989X12458724}]
#'
#' Cooper NJ, Sutton AJ, Morris D, Ades AE, Welton NJ. Addressing between-study heterogeneity and inconsistency in mixed treatment comparisons: Application to stroke prevention treatments in individuals with non-rheumatic atrial fibrillation. \emph{Stat Med} 2009;\bold{28}(14):1861--81. [\doi{10.1002/sim.3594}]
#'
#' White IR, Higgins JP, Wood AM. Allowing for uncertainty due to missing data in meta-analysis--part 1: two-stage methods. \emph{Stat Med} 2008;\bold{27}(5):711--727. [\doi{10.1002/sim.3008}]
#'
#' Lu G, Ades AE. Assessing evidence inconsistency in mixed treatment comparisons. \emph{J Am Stat Assoc} 2006;\bold{101}:447--459. [\doi{10.1198/016214505000001302}]
#'
#' Gelman A, Rubin DB. Inference from iterative simulation using multiple sequences. \emph{Stat Sci} 1992;\bold{7}:457--472. [\doi{10.1214/ss/1177011136}]
#'
#' @examples
#' data("nma.baker2009")
#'
#' # Show the first six trials of the dataset (one-trial-per-row format)
#' head(nma.baker2009)
#'
#' \dontrun{
#' # Perform a random-effects network meta-analysis
#' run.model(data = nma.baker2009,
#'           measure = "OR",
#'           model = "RE",
#'           assumption = "IDE-ARM",
#'           heter.prior = list("halfnormal", 0, 1),
#'           mean.misspar = 0,
#'           var.misspar = 1,
#'           D = 1,
#'           n.chains = 2,
#'           n.iter = 1000,
#'           n.burnin = 100,
#'           n.thin = 1)
#' }
#'
#' \dontshow{
#' closeAllConnections()
#' }
#'
#' @export
run.model <- function(data, measure, model, covar.assumption, assumption, heter.prior, mean.misspar, var.misspar, D, n.chains, n.iter, n.burnin, n.thin) {


  ## Turn off warning when variables in the 'data.jag' are not used
  options(warn = -1)


  ## Prepare the dataset for the R2jags
  item <- data.preparation(data, measure)


  ## Missing and default arguments
  model <- if (missing(model)) {
    "RE"
  } else if (!is.element(model, c("RE", "FE"))) {
    stop("Insert 'RE', or 'FE'", call. = F)
  } else {
    model
  }
  covar.assumption <- if (missing(covar.assumption)) {
    "NO"
  } else if (!is.element(covar.assumption,  c("NO", "exchangeable", "independent", "common"))) {
    stop("Insert 'NO', 'exchangeable', 'independent', or 'common'", call. = F)
  } else {
    covar.assumption
  }
  assumption <- if (missing(assumption)) {
    "IDE-ARM"
  } else if (!is.element(assumption,  c("IDE-ARM", "IDE-TRIAL", "IDE-COMMON", "HIE-ARM", "HIE-TRIAL", "HIE-COMMON", "IND-CORR", "IND-UNCORR"))) {
    stop("Insert 'IDE-ARM', 'IDE-TRIAL', 'IDE-COMMON', 'HIE-ARM', 'HIE-TRIAL', 'HIE-COMMON', 'IND-CORR', or 'IND-UNCORR'", call. = F)
  } else {
    assumption
  }
  D <- if (missing(D)) {
    stop("The argument 'D' needs to be defined", call. = F)
  } else {
    D
  }
  mean.misspar <- missingness.param.prior(assumption, mean.misspar)
  heterog.prior <- heterogeneity.param.prior(measure, model, heter.prior)
  var.misspar <- ifelse(missing(var.misspar) & (is.element(measure, c("OR", "MD", "SMD"))), 1, ifelse(missing(var.misspar) & measure == "ROM", 0.2^2, var.misspar))
  n.chains <- ifelse(missing(n.chains), 2, n.chains)
  n.iter <- ifelse(missing(n.iter), 10000, n.iter)
  n.burnin <- ifelse(missing(n.burnin), 1000, n.burnin)
  n.thin <- ifelse(missing(n.thin), 1, n.thin)


  ## Data in list format for R2jags
  data.jag <- list("m" = item$m,
                   "N" = item$N,
                   "t" = item$t,
                   "na" = item$na,
                   "nt" = item$nt,
                   "ns" = item$ns,
                   "ref" = item$ref,
                   "I" = item$I,
                   "M" = ifelse(!is.na(item$m), mean.misspar, NA),
                   "cov.phi" = 0.5*var.misspar,
                   "var.phi" = var.misspar,
                   "meand.phi" = mean.misspar,
                   "precd.phi" = 1/var.misspar,
                   "D" = D,
                   "heter.prior" = heterog.prior)


  if (is.element(measure, c("MD", "SMD", "ROM"))) {
    data.jag <- append(data.jag, list("y.o" = item$y0, "se.o" = item$se0))
  } else if (measure == "OR") {
    data.jag <- append(data.jag, list("r" = item$r))
  }


  param.jags <- c("delta", "EM", "EM.ref", "EM.pred", "pred.ref", "tau", "SUCRA",  "effectiveness", "dev.o", "totresdev.o", "hat.par")
  if (is.element(assumption, c("HIE-COMMON", "HIE-TRIAL", "HIE-ARM"))) {
    param.jags <- append(param.jags, "mean.phi")
  } else {
    param.jags <- append(param.jags, "phi")
  }


  param.jags <- if (model == "RE") {
    param.jags
  } else {
    param.jags[!is.element(param.jags, c("EM.pred", "pred.ref", "tau","delta"))]
  }


  ## Run the Bayesian analysis
  jagsfit <- jags(data = data.jag,
                  parameters.to.save = param.jags,
                  model.file = textConnection(prepare.model(measure, model, covar.assumption, assumption)),
                  n.chains = n.chains,
                  n.iter = n.iter,
                  n.burnin = n.burnin,
                  n.thin = n.thin)


  ## Turn summary of posterior results (R2jags object) into a data-frame to select model parameters (using 'dplyr')
  getResults <- as.data.frame(t(jagsfit$BUGSoutput$summary))

  # Effect size of all unique pairwise comparisons
  EM <- t(getResults %>% dplyr::select(starts_with("EM[")))

  # Effect size of all comparisons with the reference intervention
  EM.ref <- t(getResults %>% dplyr::select(starts_with("EM.ref[")))

  # Predictive effects of all unique pairwise comparisons
  EM.pred <- t(getResults %>% dplyr::select(starts_with("EM.pred[")))

  # Predictive effects of all comparisons with the reference intervention
  pred.ref <- t(getResults %>% dplyr::select(starts_with("pred.ref[")))

  # Between-trial standard deviation
  tau <- t(getResults %>% dplyr::select(starts_with("tau")))

  # SUrface under the Cumulative RAnking curve values
  SUCRA <- t(getResults %>% dplyr::select(starts_with("SUCRA")))

  # Within-trial effects size (multi-arm trials with T interventions provide T-1 such effect sizes)
  delta <- t(getResults %>% dplyr::select(starts_with("delta") & !ends_with(",1]")))

  # Ranking probability of each intervention for every rank
  effectiveness <- t(getResults %>% dplyr::select(starts_with("effectiveness")))

  # Estimated missingness parameter
  phi <- if (length(unique(unlist(item$m))) > 2) {
    t(getResults %>% dplyr::select(starts_with("phi") | starts_with("mean.phi") | starts_with("mean.phi[") | starts_with("phi[")))
  } else {
    NA
  }

  # Trial-arm deviance contribution for observed outcome
  dev.o <- t(getResults %>% dplyr::select(starts_with("dev.o")))

  # Fitted/predicted outcome
  hat.par <- t(getResults %>% dplyr::select(starts_with("hat.par")))

  # Total residual deviance
  dev <- jagsfit$BUGSoutput$summary["totresdev.o", "mean"]


  ## Calculate the deviance at posterior mean of fitted values
  # Turn 'number of observed' and 'm' into a vector (first column, followed by second column, and so on)
  m.new <- suppressMessages({as.vector(na.omit(melt(item$m)[, 2]))})
  N.new <- suppressMessages({as.vector(na.omit(melt(item$N)[, 2]))})
  obs <- N.new - m.new


  if (is.element(measure, c("MD", "SMD", "ROM"))) {

    # Turn 'y0', 'se0'into a vector (first column, followed by second column, and so on)
    y0.new <- suppressMessages({as.vector(na.omit(melt(item$y0)[, 2]))})
    se0.new <- suppressMessages({as.vector(na.omit(melt(item$se0)[, 2]))})

    # Deviance at the posterior mean of the fitted mean outcome
    dev.post.o <- (y0.new - as.vector(hat.par[, 1]))*(y0.new - as.vector(hat.par[, 1]))*(1/se0.new^2)

    # Sign of the difference between observed and fitted mean outcome
    sign.dev.o <- sign(y0.new - as.vector(hat.par[, 1]))

  } else {

    # Turn 'r' and number of observed into a vector (first column, followed by second column, and so on)
    r.new <- suppressMessages({as.vector(na.omit(melt(item$r)[, 2]))})

    # Correction for zero events in trial-arm
    r0 <- ifelse(r.new == 0, r.new + 0.01, ifelse(r.new == obs, r.new - 0.01, r.new))

    # Deviance at the posterior mean of the fitted response
    dev.post.o <- 2*(r0*(log(r0) - log(as.vector(hat.par[, 1]))) + (obs - r0)*(log(obs - r0) - log(obs - as.vector(hat.par[, 1]))))

    # Sign of the difference between observed and fitted response
    sign.dev.o <- sign(r0 - as.vector(hat.par[, 1]))

  }


  ## Obtain the leverage for observed outcomes
  leverage.o <- as.vector(dev.o[, 1]) - dev.post.o

  # Number of effective parameters
  pD <- dev - sum(dev.post.o)

  # Deviance information criterion
  DIC <- pD + dev

  # A data-frame on the measures of model assessment: DIC, pD, and total residual deviance
  model.assessment <- data.frame(DIC, pD, dev)


  ## Return a list of results
  if (model == "RE") {
    ma.results <- list(EM = EM,
                       EM.pred = EM.pred,
                       tau = tau,
                       delta = delta,
                       dev.o = dev.o,
                       hat.par = hat.par,
                       leverage.o = leverage.o,
                       sign.dev.o = sign.dev.o,
                       phi = phi,
                       model.assessment = model.assessment,
                       data = data,
                       measure = measure,
                       model = model,
                       assumption = assumption,
                       heter.prior = heterog.prior,
                       mean.misspar = mean.misspar,
                       var.misspar = var.misspar,
                       D = D,
                       jagsfit = jagsfit)
    nma.results <- append(ma.results, list(EM.ref = EM.ref, pred.ref = pred.ref, SUCRA = SUCRA, effectiveness = effectiveness))
  } else {
    ma.results <- list(EM = EM,
                       dev.o = dev.o,
                       hat.par = hat.par,
                       leverage.o = leverage.o,
                       sign.dev.o = sign.dev.o,
                       phi = phi,
                       model.assessment = model.assessment,
                       data = data,
                       measure = measure,
                       model = model,
                       mean.misspar = mean.misspar,
                       var.misspar = var.misspar,
                       D = D,
                       jagsfit = jagsfit)
    nma.results <- append(ma.results, list(EM.ref = EM.ref, SUCRA = SUCRA, effectiveness = effectiveness))
  }

  ifelse(item$nt > 2, return(nma.results), return(ma.results))
}


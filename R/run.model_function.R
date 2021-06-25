#' Perform pairwise or network meta-analysis with missing participant outcome data
#'
#' @description Perform one-stage pairwise or network meta-analysis with incorporation of pattern-mixture model to address aggregate binary or continuous missing participant outcome data
#'
#' @param data A data-frame of a one-trial-per-row format containing arm-level data of each trial. This format is widely used for BUGS models. See 'Format' for the specification of the columns.
#' @param measure Character string indicating the effect measure with values \code{"OR"}, \code{"MD"}, \code{"SMD"}, or \code{"ROM"} for the odds ratio, mean difference,
#'   standardised mean difference and ratio of means, respectively.
#' @param assumption Character string indicating the structure of the informative missingness parameter.
#'   Set \code{assumption} equal to one of the following: \code{"HIE-COMMON"}, \code{"HIE-TRIAL"}, \code{"HIE-ARM"}, \code{"IDE-COMMON"}, \code{"IDE-TRIAL"}, \code{"IDE-ARM"}, \code{"IND-CORR"}, or \code{"IND-UNCORR"}.
#' @param model Character string indicating the analysis model with values \code{"RE"}, or \code{"FE"} for the random-effects and fixed-effect model, respectively.
#' @param assumption Character string indicating the structure of the informative missingness parameter. Set \code{assumption} equal to one of the following:  \code{"IDE-ARM"}, \code{"IDE-TRIAL"},
#' \code{"IDE-COMMON"}, \code{"HIE-ARM"}, \code{"HIE-TRIAL"}, \code{"HIE-COMMON"}, \code{"IND-CORR"}, or \code{"IND-UNCORR"}.
#' @param heter.prior A list of three elements with the following order: 1) a character string indicating the distribution with (currently available) values \code{"halfnormal"},
#'   \code{"uniform"}, \code{"lognormal"}, or \code{"logt"}; 2) two numeric values that refer to the parameters of the selected distribution. For \code{"halfnormal"}, \code{"lognormal"}, and \code{"logt"}
#'   these numbers refer to the mean and precision, respectively. For \code{"uniform"}, these numbers refer to the minimum and maximum value of the distribution.
#' @param mean.misspar A real number for the mean of the normal distribution of the selected informative missingness parameter (see \code{assumption}).
#' @param var.misspar A positive non-zero number for the variance of the normal distribution of the selected informative missingness parameter (see \code{assumption}).
#' @param D A binary number for the direction of the outcome. Set \code{D = 1} for a positive outcome and \code{D = 0} for a negative outcome.
#' @param n.chains Integer specifying the number of chains for the MCMC sampling; an argument of the \code{\link[R2jags]{jags}} function.
#' @param n.iter Integer specifying the number of Markov chains for the MCMC sampling; an argument of the \code{\link[R2jags]{jags}} function.
#' @param n.burnin Integer specifying the number of iterations to discard at the beginning of the MCMC sampling; an argument of the \code{\link[R2jags]{jags}} function.
#' @param n.thin Integer specifying the thinning rate for the MCMC sampling; an argument of the \code{\link[R2jags]{jags}} function.
#'
#' @return An R2jags output on the summaries of the posterior distribution, and the GelmanRubin convergence diagnostic of the following parameters for a fixed-effect pairwise meta-analysis model:
#' \tabular{ll}{
#'  \code{EM} \tab The estimated effect measure (according to \code{measure}).\cr
#'  \tab \cr
#'  \code{dev.o} \tab The deviance contribution of each trial-arm based on the observed outcome.\cr
#'  \tab \cr
#'  \code{hat.par} \tab The fitted outcome at each trial-arm.\cr
#'  \tab \cr
#'  \code{phi} \tab The informative missingness parameter.\cr
#' }
#'
#' For a random-effects pairwise meta-analysis model, the output additionally includes the following elements:
#' \tabular{ll}{
#'  \code{EM.pred} \tab The predicted effect measure (according to \code{measure}).\cr
#'  \tab \cr
#'  \code{delta} \tab The underlying trial-specific estimated effect measure (according to \code{measure}).
#'  For a multi-arm trial, we estimate \emph{T-1} trial-specific effect estimates, where \emph{T} is the number of interventions in the trial.\cr
#'  \tab \cr
#'  \code{tau} \tab The between-trial standard deviation.\cr
#' }
#'
#' In the case of network meta-analysis (NMA), the output additionally includes:
#' \tabular{ll}{
#'  \code{EM.ref} \tab The estimated effect measure (according to \code{measure}) of all comparisons with the reference intervention.\cr
#'  \tab \cr
#'  \code{pred.ref} \tab The predicted effect measure (according to \code{measure}) of all comparisons with the reference intervention.\cr
#'  \tab \cr
#'  \code{SUCRA} \tab The surface under the cumulative ranking curve for each intervention.\cr
#'  \tab \cr
#'  \code{effectiveneness} \tab The ranking probability of each intervention for every rank.\cr
#' }
#' In NMA, \code{EM} and \code{EM.pred} refer to all possible comparisons in the network. Furthermore, \code{tau} is typically assumed to be common for all observed comparisons in the network.
#'
#' Furthermore, the output includes the following elements - the first three resulting from monitored parameters:
#' \tabular{ll}{
#'  \code{leverage.o} \tab The leverage for the observed outcome at each trial-arm.\cr
#'  \tab \cr
#'  \code{sign.dev.o} \tab The sign of the difference between observed and fitted outcome at each trial-arm.\cr
#'  \tab \cr
#'  \code{model.assessment} \tab A data-frame on the measures of model assessment: deviance information criterion, number of effective parameters, and total residual deviance.\cr
#'  \tab \cr
#'  \code{measure} \tab The effect measured as defined in \code{measure} to be used in other functions of the package.\cr
#'  \tab \cr
#'  \code{model} \tab The analysis model as defined in \code{model} to be used in other functions of the package.\cr
#'  \tab \cr
#'  \code{jagsfit} \tab An object of S3 class \code{\link[R2jags]{jags}} with the posterior results on all monitored parameters to be used in \code{mcmc.diagnostics} function.\cr
#' }
#'
#' @format The columns of the data frame \code{data} refer to the following ordered elements for a continuous outcome:
#' \tabular{ll}{
#'  \strong{t} \tab An intervention identifier.\cr
#'  \tab \cr
#'  \strong{y} \tab The observed mean value of the outcome.\cr
#'  \tab \cr
#'  \strong{sd} \tab The observed standard deviation of the outcome.\cr
#'  \tab \cr
#'  \strong{m} \tab The number of missing outcome data.\cr
#'  \tab \cr
#'  \strong{n} \tab The number of participants randomised on the assigned intervention.\cr
#' }
#'
#' For a binary outcome, the columns of the data-frame \code{data} refer to the following ordered elements:
#' \tabular{ll}{
#'  \strong{t} \tab An intervention identifier.\cr
#'  \tab \cr
#'  \strong{r} \tab The observed number of events of the outcome.\cr
#'  \tab \cr
#'  \strong{m} \tab The number of missing outcome data.\cr
#'  \tab \cr
#'  \strong{n} \tab The number of participants randomised on the assigned intervention.\cr
#' }
#' All elements appear in \code{data} as many times as the maximum number of interventions compared in a trial.
#'
#' @author {Loukia M. Spineli}
#'
#' @examples
#' \dontshow{load("./data/cipriani2011.bin.RData")}
#' ### Show the data (one-trial-per-row format)
#' (data <- as.data.frame(cipriani.bin.new))
#'
#' ### Run a random-effects network meta-analysis with consistency equations for the odds ratio (in the logarithmic scale)
#' ### assuming missing at random for identical, intervention-specific informative missingness odds ratio.
#' run.model(data = data, measure = "OR", model = "RE", assumption = "IDE-ARM", heter.prior = list("halfnormal", 0, 1), mean.misspar = 0, var.misspar = 1, D = 0, n.chains = 3, n.iter = 10000, n.burnin = 1000, n.thin = 1)
#'
#' @export
run.model <- function(data, measure, model, assumption, heter.prior, mean.misspar, var.misspar, D, n.chains, n.iter, n.burnin, n.thin) {


  ## Turn off warning when variables in the 'data.jag' are not used
  options(warn = -1)


  ## Prepare the dataset for the R2jags
  item <- data.preparation(data, measure)


  ## Default arguments
  model <- if (missing(model)) {
    "RE"
  } else if (!is.element(model, c("RE", "FE"))) {
    stop("Insert 'RE', or 'FE'")
  } else {
    model
  }
  assumption <- if (missing(assumption)) {
    "IDE-ARM"
  } else if (!is.element(assumption,  c("IDE-ARM", "IDE-TRIAL", "IDE-COMMON", "HIE-ARM", "HIE-TRIAL", "HIE-COMMON", "IND-CORR", "IND-UNCORR"))) {
    stop("Insert 'IDE-ARM', 'IDE-TRIAL', 'IDE-COMMON', 'HIE-ARM', 'HIE-TRIAL', 'HIE-COMMON', 'IND-CORR', or 'IND-UNCORR'")
  } else {
    assumption
  }
  D <- if (missing(D)) {
    stop("The 'D' needs to be defined")
  } else {
    D
  }
  mean.misspar <- missingness.param.prior(assumption, mean.misspar)
  heter.prior <- heterogeneity.param.prior(measure, model, heter.prior)
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
                   "heter.prior" = heter.prior,
                   "eff.mod2" = matrix(0, nrow = item$ns, ncol = max(item$na)),
                   "eff.mod" = rep(0, item$ns))


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
                  model.file = textConnection(prepare.model(measure, model, assumption)),
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
  phi <- t(getResults %>% dplyr::select(starts_with("phi") | starts_with("mean.phi") | starts_with("mean.phi[") | starts_with("phi[")))

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
                       measure = measure,
                       model = model,
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
                       measure = measure,
                       model = model,
                       jagsfit = jagsfit)
    nma.results <- append(ma.results, list(EM.ref = EM.ref, SUCRA = SUCRA, effectiveness = effectiveness))
  }

  ifelse(item$nt > 2, return(nma.results), return(ma.results))
}




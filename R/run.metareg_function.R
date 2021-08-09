#' Pairwise or network meta-regression with missing participant outcome data
#'
#' @description This function performs a one-stage pairwise or network meta-regression while addressing aggregate binary or continuous missing participant outcome data via the pattern-mixture model.
#'
#' @param full An object of S3 class \code{\link{run.model}}. See 'Value' in \code{\link{run.model}}.
#' @param covariate A numeric vector or a matrix for a trial-specific covariate that is a potential effect modifier. See 'Details'.
#' @param covar.assumption Character string indicating the structure of the slope for the intervention by covariate interaction, as described in Cooper et al., (2009).
#'  Set \code{covar.assumption} equal to one of the following: \code{"NO"}, when no meta-regression is performed; otherwise, \code{"exchangeable"} \code{"independent"}, and \code{"common"}.
#'  See the \code{ru.metareg} function.
#' @param n.chains Positive integer specifying the number of chains for the MCMC sampling; an argument of the \code{\link[R2jags]{jags}} function of the R-package \href{https://CRAN.R-project.org/package=R2jags}{R2jags}.
#'   The default argument is 2.
#' @param n.iter Positive integer specifying the number of Markov chains for the MCMC sampling; an argument of the \code{\link[R2jags]{jags}} function of the R-package \href{https://CRAN.R-project.org/package=R2jags}{R2jags}.
#'   The default argument is 10000.
#' @param n.burnin Positive integer specifying the number of iterations to discard at the beginning of the MCMC sampling; an argument of the \code{\link[R2jags]{jags}} function of the R-package \href{https://CRAN.R-project.org/package=R2jags}{R2jags}.
#'   The default argument is 1000.
#' @param n.thin Positive integer specifying the thinning rate for the MCMC sampling; an argument of the \code{\link[R2jags]{jags}} function of the R-package \href{https://CRAN.R-project.org/package=R2jags}{R2jags}.
#'   The default argument is 1.
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
#' @details \code{run.metareg} does not contain the arguments \code{data}, \code{measure}, \code{model}, \code{assumption}, \code{heter.prior}, \code{mean.misspar}, and \code{var.misspar} that are found in \code{run.model}.
#'   This is to prevent misspecifying the Bayesian model as it would make the comparison of the consistency model with and a covariate meaningless.
#'   Instead, these arguments are contained in the argument \code{full} of the function. Therefore, the user needs first to apply \code{run.model}, and then use \code{run.UME} (see, 'Examples').
#'
#'   The model as specified by the arguments of \code{run.model} runs in \code{JAGS} and the progress of the simulation appears in the R console.
#'   The output of \code{run.model} is used as an S3 object by other functions of the package function to be processed further and provide an end-user-ready output.
#'
#'   The code of Spineli, (2019) and Spineli et al. (2021) has been extended to incorporate one \emph{study-level covariate} variable following the assumptions of Cooper et al. (2009) for the structure of the intervention by covariate interaction.
#'   The covariate can be either a numeric vector or a matrix (the number of columns equals the maximum number of arms in the dataset).
#'
#' @author {Loukia M. Spineli}
#'
#' @seealso  \code{\link{run.model}}, \code{\link[R2jags]{jags}}
#'
#' @references
#' Spineli LM, Kalyvas C, Papadimitropoulou K. Continuous(ly) missing outcome data in network meta-analysis: a one-stage pattern-mixture model approach. \emph{Stat Methods Med Res} 2021. [\doi{10.1177/0962280220983544}]
#'
#' Spineli LM. An empirical comparison of Bayesian modelling strategies for missing binary outcome data in network meta-analysis. \emph{BMC Med Res Methodol} 2019;\bold{19}(1):86. [\doi{10.1186/s12874-019-0731-y}]
#'
#' Cooper NJ, Sutton AJ, Morris D, Ades AE, Welton NJ. Addressing between-study heterogeneity and inconsistency in mixed treatment comparisons: Application to stroke prevention treatments in individuals with non-rheumatic atrial fibrillation. \emph{Stat Med} 2009;\bold{28}(14):1861--81. [\doi{10.1002/sim.3594}]
#'
#' Gelman A, Rubin DB. Inference from iterative simulation using multiple sequences. \emph{Stat Sci} 1992;\bold{7}:457--472. [\doi{10.1214/ss/1177011136}]
#'
#' @examples
#' data("nma.baker2009")
#'
#' # Show the first six trials of the dataset (one-trial-per-row format)
#' head(nma.baker2009)
#'
#' # Perform a random-effects network meta-analysis
#' res1 <- run.model(data = nma.baker2009,
#'                   measure = "OR",
#'                   model = "RE",
#'                   assumption = "IDE-ARM",
#'                   heter.prior = list("halfnormal", 0, 1),
#'                   mean.misspar = 0,
#'                   var.misspar = 1,
#'                   D = 1,
#'                   n.chains = 3,
#'                   n.iter = 10000,
#'                   n.burnin = 1000,
#'                   n.thin = 1)
#'
#' # Whether a trial is placebo-controlled.
#' covar.binary <- ifelse(nma.baker2009[, "t1"] == 1, 1, 0)
#'
#' # Perform a random-effects network meta-regression (exchangeable structure)
#' run.metareg(full = res1,
#'             covariate = covar.binary,
#'             covar.assumption = "exchangeable",
#'             n.chains = 3,
#'             n.iter = 10000,
#'             n.burnin = 1000,
#'             n.thin = 1)
#'
#' \dontshow{
#' closeAllConnections()
#' }
#'
#' @export
run.metareg <- function(full, covariate, covar.assumption, n.chains, n.iter, n.burnin, n.thin){


  ## Turn off warning when variables in the 'data.jag' are not used
  options(warn = -1)

  data <- full$data
  measure <- full$measure
  model <- full$model
  assumption <- full$assumption
  heter.prior <- full$heter.prior
  mean.misspar <- full$mean.misspar
  var.misspar <- full$var.misspar
  D <- full$D


  ## Prepare the dataset for the R2jags
  item <- data.preparation(data, measure)

  covariate <- if (missing(covariate)) {
    stop("The argument 'covariate' needs to be defined", call. = F)
  } else {
    covariate
  }

  ## Default arguments
  covar.assumption <- if (missing(covar.assumption)) {
    "NO"
  } else if (!is.element(covar.assumption,  c("NO", "exchangeable", "independent", "common"))) {
    stop("Insert 'NO', 'exchangeable', 'independent', or 'common'", call. = F)
  } else {
    covar.assumption
  }
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
                   "heter.prior" = heter.prior)


  if (is.element(measure, c("MD", "SMD", "ROM"))) {
    data.jag <- append(data.jag, list("y.o" = item$y0, "se.o" = item$se0))
  } else if (measure == "OR") {
    data.jag <- append(data.jag, list("r" = item$r))
  }


  ## Center covariate if metric and not arm-specific
  ## Whether covariate is a vector (trial-specific) or matrix (arm-specific)
  if (!is.factor(covariate) & is.vector(covariate)) {
    data.jag <- append(data.jag, list("cov.vector" = covariate - mean(covariate), "cov.matrix" = matrix(0, nrow = item$ns, ncol = max(item$na))))
  } else if (is.factor(covariate) & is.vector(covariate)) {
    data.jag <- append(data.jag, list("cov.vector" = covariate, "cov.matrix" = matrix(0, nrow = item$ns, ncol = max(item$na))))
  } else if (!is.vector(covariate) & !is.factor(covariate)) {
    data.jag <- append(data.jag, list("cov.vector" = rep(0, item$ns), "cov.matrix" = covariate))
  }


  param.jags <- c("delta", "EM", "EM.ref", "EM.pred", "pred.ref", "tau", "beta", "SUCRA",  "effectiveness", "dev.o", "totresdev.o", "hat.par")
  if (is.element(assumption, c("HIE-COMMON", "HIE-TRIAL", "HIE-ARM"))) {
    param.jags <- append(param.jags, "mean.phi")
  } else {
    param.jags <- append(param.jags, "phi")
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

  # Regression coefficient
  beta <- t(getResults %>% dplyr::select(starts_with("beta[") | starts_with("beta")))

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


  ## Obtain the leverage for observed and missing outcomes
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
                       beta = beta,
                       dev.o = dev.o,
                       hat.par = hat.par,
                       leverage.o = leverage.o,
                       sign.dev.o = sign.dev.o,
                       phi = phi,
                       model.assessment = model.assessment,
                       measure = measure,
                       model = model,
                       covariate = covariate,
                       covar.assumption = covar.assumption,
                       jagsfit = jagsfit)
    nma.results <- append(ma.results, list(EM.ref = EM.ref, pred.ref = pred.ref, SUCRA = SUCRA, effectiveness = effectiveness))
  } else {
    ma.results <- list(EM = EM,
                       beta = beta,
                       dev.o = dev.o,
                       hat.par = hat.par,
                       leverage.o = leverage.o,
                       sign.dev.o = sign.dev.o,
                       phi = phi,
                       model.assessment = model.assessment,
                       measure = measure,
                       model = model,
                       covariate = covariate,
                       covar.assumption = covar.assumption,
                       jagsfit = jagsfit)
    nma.results <- append(ma.results, list(EM.ref = EM.ref, SUCRA = SUCRA, effectiveness = effectiveness))
  }

  ifelse(item$nt > 2, return(nma.results), return(ma.results))

}


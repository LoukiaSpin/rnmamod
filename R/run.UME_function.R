#' A function to perform Bayesian unrelated mean effects model for aggregate binary or continuous outcomes
#'
#' @param data A data-frame of a one-trial-per-row format containing arm-level data of each trial. This format is widely used for BUGS models. See 'Format' for the specification of the columns.
#' @param measure Character string indicating the effect measure with values \code{"OR"}, \code{"MD"}, \code{"SMD"}, or \code{"ROM"}.
#' @param assumption Character string indicating the structure of the informative missingness parameter. Set \code{assumption} equal to one of the following: \code{"HIE-COMMON"}, \code{"HIE-TRIAL"}, \code{"HIE-ARM"}, \code{"IDE-COMMON"}, \code{"IDE-TRIAL"}, \code{"IDE-ARM"}, \code{"IND-CORR"}, or \code{"IND-UNCORR"}.
#' @param heter.prior A vector of length equal to two with the following values: \code{rep(1, 2)}, \code{rep(2, 2)}, and \code{rep(3, 2)} refers to half-normal distribution with variance 1 or 0.5, and uniform distribution with interval [0, 5], respectively,
#' for the between-trial standard deviation. To indicate an empirically-based prior distribution for the between-trial variance, the first and second values of the vector should be the mean and precision
#' of the selected prior distribution. The empirically-based prior distribution for the between-trial variance is applicable only when \code{"OR"} or \code{"SMD"} is considered.
#' @param mean.misspar A positive non-zero number for the mean of the normal distribution of the informative missingness parameter.
#' @param var.misspar A positive non-zero number for the variance of the normal distribution of the informative missingness parameter.
#' @param n.chains Integer specifying the number of chains for the MCMC sampling; an argument of the \code{\link[R2jags]{jags}} function.
#' @param n.iter Integer specifying the number of Markov chains for the MCMC sampling; an argument of the \code{\link[R2jags]{jags}} function.
#' @param n.burnin Integer specifying the number of iterations to discard at the beginning of the MCMC sampling; an argument of the \code{\link[R2jags]{jags}} function.
#' @param n.thin Integer specifying the thinning rate for the MCMC sampling; an argument of the \code{\link[R2jags]{jags}} function.
#'
#' @return An R2jags output on the summaries of the posterior distribution, and the Gelman–Rubin convergence diagnostic of the following parameters:
#' \describe{
#'  \item{\code{EM}}{The effect estimate of all observed comparisons of interventions in the network.}
#'  \item{\code{dev.o}}{The deviance contribution of each trial-arm based on the observed outcomes.}
#'  \item{\code{resdev.o}}{The total residual deviance contribution of every trial based on the observed outcomes.}
#'  \item{\code{totresdev.o}}{The total residual deviance based on the observed outcomes.}
#'  \item{\code{dev.m}}{The deviance contribution of each trial-arm based on the missing outcomes.}
#'  \item{\code{resdev.m}}{The total residual deviance contribution of every trial based on the missing outcomes.}
#'  \item{\code{resdev.m}}{The total residual deviance based on the missing outcomes.}
#'  \item{\code{tau}}{The between-trial standard deviation assumed to be common for all observed comparisons.}
#' }
#'
#' @format The columns of the data frame \code{data} refer to the following ordered elements for a continuous outcome:
#' \describe{
#'  \item{\strong{t}}{An intervention identifier.}
#'  \item{\strong{y}}{The observed mean value of the outcome.}
#'  \item{\strong{sd}}{The observed standard deviation of the outcome.}
#'  \item{\strong{m}}{The number of missing outcome data.}
#'  \item{\strong{c}}{The number of participants completing the assigned intervention.}
#'  \item{\strong{na}}{The number of compared interventions.}
#' }
#' Apart from \strong{na}, all other elements appear in \code{data} as many times as the maximum number of interventions compared in a trial. See, 'Example'.
#'
#' @seealso \code{\link{R2jags}}
#'
#' @references
#' Gelman, A, Rubin, DB. Inference from iterative simulation using multiple sequences. Stat Sci. 1992;7:457–472.
#'
#' \dontshow{load("./data/NMA Dataset Continuous.RData")}
#' @examples
#' ### Show the data (one-trial-per-row format)
#' (data <- as.data.frame(one.stage.dataset.NMA[[3]]))
#'
#' ### Run a random-effects network meta-analysis with consistency equations for the standardised mean difference
#' ### assuming missing at random for identical, common informative missingness difference of means.
#' run.UME(data = data, measure = "SMD", assumption = "IDE-COMMON", mean.misspar = 0, var.misspar = 1, n.chains = 3, n.iter = 10000, n.burnin = 1000, n.thin = 1)
#'
#' @export
run.UME <- function(data, measure, model, assumption, heter.prior, mean.misspar, var.misspar, n.chains, n.iter, n.burnin, n.thin) {


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
  mean.misspar <- missingness.param.prior(assumption, mean.misspar)
  heter.prior <- heterogeneity.param.prior(measure, model, heter.prior)
  var.misspar <- ifelse(missing(var.misspar) & (is.element(measure, c("OR", "MD", "SMD"))), 1, ifelse(missing(var.misspar) & measure == "ROM", 0.2^2, var.misspar))
  n.chains <- ifelse(missing(n.chains), 2, n.chains)
  n.iter <- ifelse(missing(n.iter), 10000, n.iter)
  n.burnin <- ifelse(missing(n.burnin), 1000, n.burnin)
  n.thin <- ifelse(missing(n.thin), 1, n.thin)


  ## Unique comparisons with the baseline intervention
  ## A function to extract numbers from a character. Source: http://stla.github.io/stlapblog/posts/Numextract.html
  Numextract <- function(string) {
    unlist(regmatches(string,gregexpr("[[:digit:]]+\\.*[[:digit:]]*", string)))
  }


  ## Sort by the number of arms in ascending order
  t <- item$t[order(item$na, na.last = T), ]
  m <- item$m[order(item$na, na.last = T), ]
  N <- item$N[order(item$na, na.last = T), ]
  na <- sort(item$na)

  ## Observed comparisons in the network
  observed.comp0 <- improved.UME(t, m, N, item$ns, na)$obs.comp
  observed.comp <- matrix(Numextract(observed.comp0[, 1]), nrow = length(observed.comp0[, 1]), ncol = 2, byrow = T)
  t1.obs.com <- as.numeric(as.character(observed.comp[, 1]))
  t2.obs.com <- as.numeric(as.character(observed.comp[, 2]))
  obs.comp <- paste0(t2.obs.com, "vs", t1.obs.com)


  ## Keep only comparisons with the baseline intervention
  indic0 <- list()
  for(i in 1:item$ns) {
    indic0[[i]] <- combn(t(na.omit(t(t[i, ]))), 2)[, 1:(item$na[i] - 1)]
  }
  (indic <- unique(t(do.call(cbind, indic0))))
  t1.indic <- indic[, 1]
  t2.indic <- indic[, 2]
  N.obs <- length(t1.indic)


  ## Keep only comparisons with the baseline intervention in multi-arm trials
  ns.multi <- length(item$na[na > 2])
  if (ns.multi < 1) {
    N.obs.multi <- 0
    t1.indic.multi <- 0
    t2.indic.multi <- 0

  } else {
    indic.multi0 <- list()
    for(i in (item$ns - ns.multi + 1):item$ns) {
      indic.multi0[[i]] <- combn(t(na.omit(t(t[i, ]))), 2)[, 1:(na[i] - 1)]
    }
    (indic.multi <- unique(t(do.call(cbind, indic.multi0))))
    t1.indic.multi <- indic.multi[, 1]
    t2.indic.multi <- indic.multi[, 2]
    N.obs.multi <- length(t1.indic.multi)
  }


  ## Data in list format for R2jags
  data.jag <- list("m" = m,
                   "N" = N,
                   "t" = t,
                   "na" = na,
                   "nt" = item$nt,
                   "ns" = item$ns,
                   "ref" = ifelse(is.element(assumption, c("HIE-ARM", "IDE-ARM")), item$ref, NA),
                   "I" = item$I[order(item$na, na.last = T), ],
                   "M" = ifelse(!is.na(m), mean.misspar, NA),
                   "cov.phi" = 0.5*var.misspar,
                   "var.phi" = var.misspar,
                   "meand.phi" = mean.misspar,
                   "precd.phi" = 1/var.misspar,
                   "heter.prior" = heter.prior,
                   "t1" = t1.indic,
                   "t2" = t2.indic,
                   "N.obs" = N.obs)


  if (is.element(measure, c("MD", "SMD", "ROM"))) {
    data.jag <- append(data.jag, list("y.o" = item$y0[order(item$na, na.last = T), ], "se.o" = item$se0[order(item$na, na.last = T), ]))
  } else if (measure == "OR") {
    data.jag <- append(data.jag, list("r" = item$r[order(item$na, na.last = T), ]))
  }


  if (max(na) > 2 & has_error(improved.UME(t, m, N, item$ns, na), silent = T) == F) {
    impr.UME <- improved.UME(t, m, N, item$ns, na)
    data.jag <- append(data.jag, list("t1.m" = t1.indic.multi,
                                      "t2.m" = t2.indic.multi,
                                      "N.obs.multi" = N.obs.multi,
                                      "ns.multi" = ns.multi,
                                      "t1.bn" = impr.UME$t1.bn,
                                      "t2.bn" = impr.UME$t2.bn,
                                      "base" = impr.UME$base,
                                      "nbase.multi" = impr.UME$nbase.multi))
  } else if (max(na) < 3 || has_error(improved.UME(t, m, N, item$ns, na), silent = T) == T) {
    data.jag <- append(data.jag, list("t1.m" = t1.indic.multi,
                                      "t2.m" = t2.indic.multi,
                                      "N.obs.multi" = N.obs.multi,
                                      "ns.multi" = ns.multi,
                                      "t1.bn" = 0,
                                      "t2.bn" = 0,
                                      "base" = 0,
                                      "nbase.multi" = 0))
  }


  ## Define the nodes to be monitored
  param.jags <- if (model == "RE") {
    c("EM", "EM.m", "dev.o", "resdev.o", "totresdev.o", "tau", "hat.par")
  } else {
    c("EM", "EM.m", "dev.o", "resdev.o", "totresdev.o", "hat.par")
  }


  ## Run the Bayesian analysis
  jagsfit <- jags(data = data.jag,
                  parameters.to.save = param.jags,
                  model.file = textConnection(prepare.UME(measure, model, assumption)),
                  n.chains = n.chains,
                  n.iter = n.iter,
                  n.burnin = n.burnin,
                  n.thin = n.thin)


  ## Turn summary of posterior results (R2jags object) into a data-frame to select model parameters (using 'dplyr')
  getResults <- as.data.frame(t(jagsfit$BUGSoutput$summary))

  # Effect size of all observed pairwise comparisons
  EM <- t(getResults %>% dplyr::select(starts_with("EM[")))

  # Effect size of observed pairwise comparisons with the baseline intervention in multi-arm trials
  EM.m <- t(getResults %>% dplyr::select(starts_with("EM.m[")))

  # Between-trial standard deviation
  tau <- t(getResults %>% dplyr::select(starts_with("tau")))

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

  # Collect the minimum results at common
  results <- if (model == "RE") {
    list(EM = EM,
         dev.o = dev.o,
         hat.par = hat.par,
         leverage.o = leverage.o,
         sign.dev.o = sign.dev.o,
         tau = tau,
         model.assessment = model.assessment,
         measure = measure,
         model = model,
         obs.comp = obs.comp,
         jagsfit = jagsfit)
  } else {
    list(EM = EM,
         dev.o = dev.o,
         hat.par = hat.par,
         leverage.o = leverage.o,
         sign.dev.o = sign.dev.o,
         model.assessment = model.assessment,
         measure = measure,
         model = model,
         obs.comp = obs.comp,
         jagsfit = jagsfit)
  }

  # Return different list of results according to a condition
  if (max(na) < 3 || has_error(improved.UME(t, m, N, ns, na), silent = T) == T) {
    return(results)
  } else {
    return(append(results, list(EM.m = EM.m, frail.comp = paste0(impr.UME$t2.bn, "vs", impr.UME$t1.bn))))
  }

}


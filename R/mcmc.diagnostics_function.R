#' Markov Chain Monte Carlo Diagnostics
#'
#' @param data A data-frame of a one-trial-per-row format containing arm-level data of each trial. This format is widely used for BUGS models. See 'Format' for the specification of the columns.
#' @param measure Character string indicating the effect measure with values \code{"OR"}, \code{"MD"}, \code{"SMD"}, or \code{"ROM"}.
#' @param assumption Character string indicating the structure of the informative missingness parameter. Set \code{assumption} equal to one of the following: \code{"HIE-COMMON"}, \code{"HIE-TRIAL"}, \code{"HIE-ARM"}, \code{"IDE-COMMON"}, \code{"IDE-TRIAL"}, \code{"IDE-ARM"}, \code{"IND-CORR"}, or \code{"IND-UNCORR"}.
#' @param heter.prior A vector of length equal to two with the following values: \code{rep(1, 2)}, \code{rep(2, 2)}, and \code{rep(3, 2)} refers to half-normal distribution with variance 1 or 0.5, and uniform distribution with interval [0, 5], respectively,
#' for the between-trial standard deviation. To indicate an empirically-based prior distribution for the between-trial variance, the first and second values of the vector should be the mean and precision
#' of the selected prior distribution. The empirically-based prior distribution for the between-trial variance is applicable only when \code{"OR"} or \code{"SMD"} is considered.
#' @param net.ref Integer specifying the reference intervention of the network. The default is the most frequently appeared intervention in the network.
#' @param mean.misspar A positive non-zero number for the mean of the normal distribution of the informative missingness parameter.
#' @param var.misspar A positive non-zero number for the variance of the normal distribution of the informative missingness parameter.
#' @param D A binary number for the direction of the outcome. Set \code{D = 1} for a positive outcome and \code{D = 0} for a negative outcome.
#' @param n.chains Integer specifying the number of chains for the MCMC sampling; an argument of the \code{\link[R2jags]{jags}} function.
#' @param n.iter Integer specifying the number of Markov chains for the MCMC sampling; an argument of the \code{\link[R2jags]{jags}} function.
#' @param n.burnin Integer specifying the number of iterations to discard at the beginning of the MCMC sampling; an argument of the \code{\link[R2jags]{jags}} function.
#' @param n.thin Integer specifying the thinning rate for the MCMC sampling; an argument of the \code{\link[R2jags]{jags}} function.
#'
#' @return A panel of autocorrelation plots where the rows correspond to the chains and the columns correspond to the monitor parameters (maximum three).
#' Additionally, it uses the \code{\link[mcmcplots]{mcmcplot}} function to create an HTML file with a panel of diagnostic plots (trace, density, and autocorrelation) for each monitored parameter.
#'
#' @format See, function \code{nma.continuous.full.model}.
#'
#' @seealso \code{\link{mcmcplots}}
#'
#' @references
#' Gelman, A, Rubin, DB. Inference from iterative simulation using multiple sequences. Stat Sci. 1992;7:457–472.
#'
#' \dontshow{load("netmodr/data/One-stage model_NMA Dataset.RData")}
#' @examples
#' ### Obtain the diagnostic plots and check convergence for all monitored parameters using the R.hat
#' mcmc.diagnostics(par = c("tau2", "EM[3,1]", "EM[3,2]"), net = res1)
#'
#' @export
mcmc.diagnostics <- function(par, data, measure, assumption, heter.prior, net.ref, mean.misspar, var.misspar, D, n.chains, n.iter, n.burnin, n.thin){



  ## Default arguments
  assumption <- ifelse(missing(assumption), "IDE-ARM", assumption)
  var.misspar <- ifelse(missing(var.misspar) & (measure == "OR" || measure == "MD"|| measure == "SMD"), 1, ifelse(missing(var.misspar) & measure == "ROM", 0.2^2, var.misspar))
  n.chains <- ifelse(missing(n.chains), 2, n.chains)
  n.iter <- ifelse(missing(n.iter), 10000, n.iter)
  n.burnin <- ifelse(missing(n.burnin), 1000, n.burnin)
  n.thin <- ifelse(missing(n.thin), 1, n.thin)


  ## For a continuous outcome
  if(measure == "MD" || measure == "SMD"|| measure == "ROM"){


    ## Continuous: arm-level, wide-format dataset
    (y.obs <- data %>% dplyr::select(starts_with("y")))                             # Observed mean value in each arm of every trial
    (sd.obs <- data %>% dplyr::select(starts_with("sd")))                           # Observed standard deviation in each arm of every trial
    (mod <- data %>% dplyr::select(starts_with("m")))                               # Number of missing participants in each arm of every trial
    (c <- data %>% dplyr::select(starts_with("c")))                                 # Number of completers in each arm of every trial
    (se.obs <- sd.obs/sqrt(c))                                                      # Observed standard error in each arm of every trial
    (rand <- mod + c)                                                               # Number of randomised participants in each arm of every trial
    (treat <- data %>% dplyr::select(starts_with("t")))                             # Intervention studied in each arm of every trial
    na <- apply(treat, 1, function(x) length(which(!is.na(x))))                     # Number of interventions investigated in every trial per network
    nt <- length(table(as.matrix(treat)))                                           # Total number of interventions per network
    ns <- length(y.obs[, 1])                                                        # Total number of included trials per network
    ref <- ifelse(missing(net.ref), which.max(table(as.matrix(treat))), net.ref)    # Reference intervention per network. If not specify, the most frequently appeared intervention in the network is selected
    # Trial-specific observed pooled standard deviation
    (sigma <- sqrt(apply((sd.obs^2)*(c - 1), 1, sum, na.rm = T)/(apply(c, 1, sum, na.rm = T) - na)))


    ## Order by 'id of t1' < 'id of t1'
    y0 <- se0 <- m <- N <- t <- t0 <- treat
    for(i in 1:ns){

      t0[i, ] <- order(treat[i, ], na.last = T)
      y0[i, ] <- y.obs[i, order(t0[i, ], na.last = T)]
      se0[i, ] <- se.obs[i, order(t0[i, ], na.last = T)]
      m[i, ] <- mod[i, order(t0[i, ], na.last = T)]
      N[i, ] <- rand[i, order(t0[i, ], na.last = T)]
      t[i, ] <- sort(treat[i, ], na.last = T)
    }



    ## Condition regarding the specification of the prior mean ('mean.misspar') for the missingness parameter
    if(missing(mean.misspar)) {

      mean.misspar <- rep(0, 2)

    } else if(!missing(mean.misspar) & (assumption == "HIE-ARM" || assumption == "IDE-ARM" ) & !is.null(dim(mean.misspar))) {

      mean.misspar <- as.vector(mean.misspar)

    } else if(!missing(mean.misspar) & (assumption == "HIE-ARM" || assumption == "IDE-ARM" ) & is.null(dim(mean.misspar))) {

      mean.misspar <- rep(mean.misspar, 2)

    } else {

      mean.misspar <- mean.misspar

    }



    ## Information for the prior distribution on the missingness parameter (IMDOM or logIMROM)
    M <- ifelse(!is.na(y0), mean.misspar, NA)  # Vector of the mean value of the normal distribution of the informative missingness parameter as the number of arms in trial i (independent structure)
    prec.misspar <- 1/var.misspar
    psi.misspar <- sqrt(var.misspar)           # the lower bound of the uniform prior distribution for the prior standard deviation of the missingness parameter (hierarchical structure)
    cov.misspar <- 0.5*var.misspar             # covariance of pair of missingness parameters in a trial (independent structure)



    ## Specification of the prior distribution for the between-trial parameter
    if (heter.prior[[1]] == "halfnormal") {

      heter.prior <- as.numeric(c(0, heter.prior[[3]], 1))

    } else if (heter.prior[[1]] == "uniform") {

      heter.prior <- as.numeric(c(0, heter.prior[[3]], 2))

    } else if (measure == "SMD" & heter.prior[[1]] == "logt") {

      heter.prior <- as.numeric(c(heter.prior[[2]], heter.prior[[3]], 3))

    } else if (measure != "SMD" & heter.prior[[1]] == "logt") {

      stop("There are currently no empirically-based prior distributions for MD and ROM. Choose a half-normal or a uniform prior distribution, instead")

    }



    # Under the Independent structure with or without SMD as effect measure
    if (measure == "SMD" & assumption != "IND-CORR") {

      data.jag <- list("y.o" = y0, "se.o" = se0, "m" = m, "N" = N, "t" = t, "na" = na, "nt" = nt, "ns" = ns, "ref" = ref, "sigma" = sigma, "meand.phi" = mean.misspar, "precd.phi" = prec.misspar, "D" = D, "heter.prior" = heter.prior)

    } else if (measure == "SMD" & assumption == "IND-CORR"){

      data.jag <- list("y.o" = y0, "se.o" = se0, "m" = m, "N" = N, "t" = t, "na" = na, "nt" = nt, "ns" = ns, "ref" = ref, "sigma" = sigma, "M" = M, "cov.phi" = cov.misspar, "var.phi" = var.misspar, "D" = D, "heter.prior" = heter.prior)

    } else if (measure != "SMD" & assumption == "IND-CORR") {

      data.jag <- list("y.o" = y0, "se.o" = se0, "m" = m, "N" = N, "t" = t, "na" = na, "nt" = nt, "ns" = ns, "ref" = ref, "M" = M, "cov.phi" = cov.misspar, "var.phi" = var.misspar, "D" = D, "heter.prior" = heter.prior)

    } else {

      data.jag <- list("y.o" = y0, "se.o" = se0, "m" = m, "N" = N, "t" = t, "na" = na, "nt" = nt, "ns" = ns, "ref" = ref, "meand.phi" = mean.misspar, "precd.phi" = prec.misspar, "D" = D, "heter.prior" = heter.prior)

    }


  } else {



    ## Binary: arm-level, wide-format dataset
    (event <- data %>% dplyr::select(starts_with("r")))                             # Number of observed events in each arm of every trial
    (mod <- data %>% dplyr::select(starts_with("m")))                               # Number of missing participants in each arm of every trial
    (rand <- data %>% dplyr::select(starts_with("n")))                              # Number randomised participants in each arm of every trial
    (treat <- data %>% dplyr::select(starts_with("t")))                             # Intervention studied in each arm of every trial
    na <- apply(treat, 1, function(x) length(which(!is.na(x))))                     # Number of interventions investigated in every trial per network
    nt <- length(table(as.matrix(treat)))                                           # Total number of interventions per network
    ns <- length(event[, 1])                                                        # Total number of included trials per network
    ref <- ifelse(missing(net.ref), which.max(table(as.matrix(treat))), net.ref)    # Reference intervention per network. If not specify, the most frequently appeared intervention in the network is selected


    ## Order by 'id of t1' < 'id of t1'
    r <- m <- N <- t <- t0 <- treat
    for(i in 1:ns){

      t0[i, ] <- order(treat[i, ], na.last = T)
      r[i, ] <- event[i, order(t0[i, ], na.last = T)]
      m[i, ] <- mod[i, order(t0[i, ], na.last = T)]
      N[i, ] <- rand[i, order(t0[i, ], na.last = T)]
      t[i, ] <- sort(treat[i, ], na.last = T)
    }


    ## Condition regarding the specification of the prior mean ('mean.misspar') for the missingness parameter
    if(missing(mean.misspar)) {

      mean.misspar <- rep(0.0001, 2)

    } else if(!missing(mean.misspar) & (assumption == "HIE-ARM" || assumption == "IDE-ARM" ) & !is.null(dim(mean.misspar))) {

      mean.misspar <- as.vector(mean.misspar)
      mean.misspar[1] <- ifelse(mean.misspar[1] == 0, 0.0001, mean.misspar[1])
      mean.misspar[2] <- ifelse(mean.misspar[2] == 0, 0.0001, mean.misspar[2])

    } else if(!missing(mean.misspar) & (assumption == "HIE-ARM" || assumption == "IDE-ARM" ) & is.null(dim(mean.misspar))) {

      mean.misspar <- rep(ifelse(mean.misspar == 0, 0.0001, mean.misspar), 2)

    } else if(!missing(mean.misspar) & (assumption != "HIE-ARM" || assumption != "IDE-ARM")) {

      mean.misspar <- ifelse(mean.misspar == 0, 0.0001, mean.misspar)

    }


    M <- ifelse(!is.na(r), mean.misspar, NA)   # Vector of the mean value of the normal distribution of the informative missingness parameter as the number of arms in trial i (independent structure)
    prec.misspar <- 1/var.misspar
    psi.misspar <- sqrt(var.misspar)           # the lower bound of the uniform prior distribution for the prior standard deviation of the missingness parameter (hierarchical structure)
    cov.misspar <- 0.5*var.misspar             # covariance of pair of missingness parameters in a trial (independent structure)



    ## Specification of the prior distribution for the between-trial parameter
    if (heter.prior[[1]] == "halfnormal") {

      heter.prior <- as.numeric(c(0, heter.prior[[3]], 1))

    } else if (heter.prior[[1]] == "uniform") {

      heter.prior <- as.numeric(c(0, heter.prior[[3]], 2))

    } else if (heter.prior[[1]] == "lognormal")  {

      heter.prior <- as.numeric(c(heter.prior[[2]], heter.prior[[3]], 3))

    }



    ## Condition for the Independent structure
    if (assumption != "IND-CORR") {

      data.jag <- list("r" = r, "m" = m, "N" = N, "t" = t, "na" = na, "nt" = nt, "ns" = ns, "ref" = ref, "meand.phi" = mean.misspar, "precd.phi" = prec.misspar, "D" = D, "heter.prior" = heter.prior)

    } else {

      data.jag <- list("r" = r, "m" = m, "N" = N, "t" = t, "na" = na, "nt" = nt, "ns" = ns, "ref" = ref, "M" = M, "cov.phi" = cov.misspar, "var.phi" = var.misspar, "D" = D, "heter.prior" = heter.prior)

    }

  }


  ## Condition for the hierarchical structure of the missingness parameter
  if (assumption == "HIE-COMMON" || assumption == "HIE-TRIAL" || assumption == "HIE-ARM") {

    param.jags <- c("delta", "EM", "EM.ref", "EM.pred", "pred.ref", "tau", "SUCRA", "mean.phi", "effectiveness", "dev.m", "dev.o", "hat.par", "hat.m")

  } else {

    param.jags <- c("delta", "EM", "EM.ref", "EM.pred", "pred.ref", "tau", "SUCRA", "phi", "effectiveness", "dev.m", "dev.o", "hat.par", "hat.m")

  }



  ## Run the Bayesian analysis
  jagsfit <- jags(data = data.jag, parameters.to.save = param.jags, model.file = paste0("./model/Full RE-NMA/Full RE-NMA_", measure, "_Pattern-mixture_", assumption, ".txt"),
                  n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin, DIC = T)



  ## Turn results into a data-frame to select model parameters (using 'dplyr')
  getResults <- as.data.frame(t(jagsfit$BUGSoutput$summary))

  # Effect size of all unique pairwise comparisons
  (EM <- t(getResults %>% dplyr::select(starts_with("EM["))))

  # Predictive effects of all unique pairwise comparisons
  (EM.pred <- t(getResults %>% dplyr::select(starts_with("EM.pred["))))

  # SUrface under the Cumulative RAnking curve values
  (SUCRA <- t(getResults %>% dplyr::select(starts_with("SUCRA"))))

  # Within-trial effects size (multi-arm trials with T interventions provide T-1 such effect sizes)
  (delta <- t(getResults %>% dplyr::select(starts_with("delta"))))

  # Ranking probability of each intervention for every rank
  (effectiveness <- t(getResults %>% dplyr::select(starts_with("effectiveness"))))

  # Estimated missingness parameter
  (phi <- t(getResults %>% dplyr::select(starts_with("mean.phi"))))

  # Between-trial standard deviation
  (tau <- t(getResults %>% dplyr::select(starts_with("tau"))))



  ## Turn 'R2jags' object into 'mcmc.plot' object
  jagsfit.mcmc <- as.mcmc(jagsfit)



  ## A panel of autocorrelation plots for each chain and every monitored parameter
  autocorrelation <- par(mfrow = c(3, n.chains))
  for(i in 1:n.chains){

    autplot1(jagsfit.mcmc[, par[1]], chain = i, main = paste(par[1], "-", "chain", i))
    autplot1(jagsfit.mcmc[, par[2]], chain = i, main = paste(par[2], "-","chain", i))
    autplot1(jagsfit.mcmc[, par[3]], chain = i, main = paste(par[3], "-","chain", i))

  }



  ## An HTML file with a panel of diagnostic plots per monitored paraemter
  mcmcplot <- mcmcplot(jagsfit.mcmc, parms = par)



  ## Keep results on the maximum Rhat for the selected monitored model parameters
  if(is.null(dim(phi))){

    R.hat.max <- c(max(EM[, 5]), max(EM.pred[, 5]), max(delta[, 5]), max(tau[5]), max(SUCRA[, 5]), max(effectiveness[, 5]), phi[5])

  } else {

    R.hat.max <- c(max(EM[, 5]), max(EM.pred[, 5]), max(delta[, 5]), max(tau[5]), max(SUCRA[, 5]), max(effectiveness[, 5]), max(phi[, 5]))

  }



  ## Indicate whether each model parameter achieved or failed to achieve convergence
  conv <- rep(NA, length(R.hat.max))
  for(i in 1:length(R.hat.max)) {

    conv[i] <- ifelse(R.hat.max[i] < 1.1, "achieved", "failed")

  }



  ## A data-frame with results on convergence for all monitored parameters using the Rhat
  convergence <- data.frame(R.hat.max, conv)
  rownames(convergence) <- c("EM", "Pred", "delta", "tau", "SUCRA", "effectiveness", "phi")
  colnames(convergence) <- c("R.hat max", "convergence status")


  return(list(convergence = convergence))
}


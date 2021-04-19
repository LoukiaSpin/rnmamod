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
run.UME <- function(data, measure, rho, assumption, heter.prior, mean.misspar, var.misspar, n.chains, n.iter, n.burnin, n.thin) {


  options(warn = -1)

  ## Default arguments
  rho <- ifelse(missing(rho) & (measure == "MD" || measure == "SMD"|| measure == "ROM"), 0.5, ifelse(missing(rho) & measure == "OR", NA, rho))
  assumption <- ifelse(missing(assumption), "IDE-ARM", assumption)
  var.misspar <- ifelse(missing(var.misspar) & (measure == "OR" || measure == "MD"|| measure == "SMD"), 1, ifelse(missing(var.misspar) & measure == "ROM", 0.2^2, var.misspar))
  n.chains <- ifelse(missing(n.chains), 2, n.chains)
  n.iter <- ifelse(missing(n.iter), 10000, n.iter)
  n.burnin <- ifelse(missing(n.burnin), 1000, n.burnin)
  n.thin <- ifelse(missing(n.thin), 1, n.thin)


  ## For a continuous outcome
  if(measure == "MD" || measure == "SMD"|| measure == "ROM"){


    ## Continuous: arm-level, wide-format dataset
    y.obs <- data %>% dplyr::select(starts_with("y"))                               # Observed mean value in each arm of every trial
    sd.obs <- data %>% dplyr::select(starts_with("sd"))                             # Observed standard deviation in each arm of every trial
    y.bas <- data %>% dplyr::select(starts_with("bas.y"))                           # Mean value at baseline in each arm of every trial
    sd.bas <- data %>% dplyr::select(starts_with("bas.sd"))                         # Standard deviation at baseline in each arm of every trial
    ind <- data %>% dplyr::select(starts_with("ind"))                               # Trial indicator (1: final measurement; 2: change from baseline and baseline measurement
    mod <- data %>% dplyr::select(starts_with("m"))                                 # Number of missing participants in each arm of every trial
    c <- data %>% dplyr::select(starts_with("c"))                                   # Number of completers in each arm of every trial
    se.obs <- sd.obs/sqrt(c)                                                        # Observed standard error in each arm of every trial
    rand <- mod + c                                                                 # Number of randomised participants in each arm of every trial
    if (dim(ind)[2] == 0) {
      se.bas <- NA                                                                  # Standard error at baseline in each arm of every trial
    } else {
      se.bas <- sd.bas/sqrt(rand)                                                   # Standard error at baseline in each arm of every trial
    }
    treat <- data %>% dplyr::select(starts_with("t"))                               # Intervention studied in each arm of every trial
    na <- apply(treat, 1, function(x) length(which(!is.na(x))))                     # Number of interventions investigated in every trial per network
    nt <- length(table(as.matrix(treat)))                                           # Total number of interventions per network
    ns <- length(y.obs[, 1])                                                        # Total number of included trials per network
    n1 <- ifelse(dim(ind)[2] == 0, ns, table(ind)[1])                               # Number of trials reporting data at final point
    n2 <- ifelse(dim(ind)[2] == 0, 0, table(ind)[2])                                # Number of trials reporting data at change from baseline
    ref <- 1                                                                        # The first intervention (t1 = 1) is the reference of the network


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


    ## Order by 'id of t1' < 'id of t1' - Trials reporting change from baseline and baseline per arm
    y.b <- se.b <- treat
    if (dim(ind)[2] == 0) {

      for(i in 1:ns){
        y.b[i, ] <- NA
        se.b[i, ] <- NA
      }

    } else {

      for(i in 1:ns){
        y.b[i, ] <- y.bas[i, order(t0[i, ], na.last = T)]
        se.b[i, ] <- se.bas[i, order(t0[i, ], na.last = T)]
      }

    }



    ## Unique comparisons with the baseline intervention
    ## A function to extract numbers from a character. Source: http://stla.github.io/stlapblog/posts/Numextract.html
    Numextract <- function(string){
      unlist(regmatches(string,gregexpr("[[:digit:]]+\\.*[[:digit:]]*",string)))
    }


    ## Observed comparisons in the network
    observed.comp0 <- improved.UME(t, m, N, ns, na)$obs.comp
    observed.comp <- matrix(Numextract(observed.comp0[, 1]), nrow = length(observed.comp0[, 1]), ncol = 2, byrow = T)
    t1.obs.com <- as.numeric(as.character(observed.comp[, 1]))
    t2.obs.com <- as.numeric(as.character(observed.comp[, 2]))
    obs.comp <- paste0(t2.obs.com, "vs", t1.obs.com)


    ## Keep only comparisons with the baseline intervention
    indic0 <- list()
    for(i in 1:ns) {
      indic0[[i]] <- combn(t(na.omit(t(t[i, ]))), 2)[, 1:(na[i] - 1)]
    }
    (indic <- unique(t(do.call(cbind, indic0))))
    t1.indic <- indic[, 1]
    t2.indic <- indic[, 2]
    N.obs <- length(t1.indic)



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



    if (max(na) > 2 & has_error(improved.UME(t, m, N, ns, na), silent = T) == F) {

      impr.UME <- improved.UME(t, m, N, ns, na)

      data.jag <- list("y.o" = y0, "se.o" = se0, "y.b" = y.b, "se.b" = se.b, "m" = m, "N" = N, "t" = t, "na" = na, "nt" = nt, "ns" = ns, "n1" = n1, "n2" = n2, rho = "rho", "ref" = ref, "M" = M, "cov.phi" = cov.misspar, "var.phi" = var.misspar, "meand.phi" = mean.misspar, "precd.phi" = prec.misspar, "heter.prior" = heter.prior, "t1" = t1.indic, "t2" = t2.indic, "N.obs" = N.obs, "t1.bn" = impr.UME$t1.bn, "t2.bn" = impr.UME$t2.bn, "base" = impr.UME$base, "nbase.multi" = impr.UME$nbase.multi)


    } else if (max(na) < 3 || has_error(improved.UME(t, m, N, ns, na), silent = T) == T) {

      ## Condition for the data specification based on the assumption about the structure of the missingness parameter
     if (assumption == "IND-CORR"){

        data.jag <- list("y.o" = y0, "se.o" = se0, "m" = m, "N" = N, "t" = t, "na" = na, "ns" = ns, "M" = M, "cov.phi" = cov.misspar, "var.phi" = var.misspar, "heter.prior" = heter.prior, "t1" = t1.indic, "t2" = t2.indic, "N.obs" = N.obs, "nbase.multi" = 0)

      } else if (assumption == "HIE-ARM" || assumption == "IDE-ARM") {

        data.jag <- list("y.o" = y0, "se.o" = se0, "y.b" = y.b, "se.b" = se.b, "m" = m, "N" = N, "t" = t, "na" = na, "nt" = nt, "ns" = ns, "n1" = n1, "n2" = n2, "rho" = rho, "ref" = ref, "meand.phi" = mean.misspar, "precd.phi" = prec.misspar, "heter.prior" = heter.prior, "t1" = t1.indic, "t2" = t2.indic, "N.obs" = N.obs, "nbase.multi" = 0)

      } else if (assumption != "IND-CORR" || assumption != "HIE-ARM" || assumption != "IDE-ARM") {

        data.jag <- list("y.o" = y0, "se.o" = se0, "m" = m, "N" = N, "t" = t, "na" = na, "ns" = ns, "meand.phi" = mean.misspar, "precd.phi" = prec.misspar, "heter.prior" = heter.prior, "t1" = t1.indic, "t2" = t2.indic, "N.obs" = N.obs, "nbase.multi" = 0)

      }

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
    ref <- 1                                                                        # The first intervention (t1 = 1) is the reference of the network


    ## Order by 'id of t1' < 'id of t1'
    r <- m <- N <- t <- t0 <- treat
    for(i in 1:ns){

      t0[i, ] <- order(treat[i, ], na.last = T)
      r[i, ] <- event[i, order(t0[i, ], na.last = T)]
      m[i, ] <- mod[i, order(t0[i, ], na.last = T)]
      N[i, ] <- rand[i, order(t0[i, ], na.last = T)]
      t[i, ] <- sort(treat[i, ], na.last = T)
    }



    ## Unique comparisons with the baseline intervention
    ## A function to extract numbers from a character. Source: http://stla.github.io/stlapblog/posts/Numextract.html
    Numextract <- function(string){
      unlist(regmatches(string,gregexpr("[[:digit:]]+\\.*[[:digit:]]*",string)))
    }


    ## Observed comparisons in the network
    observed.comp0 <- improved.UME(t, m, N, ns, na)$obs.comp
    observed.comp <- matrix(Numextract(observed.comp0[, 1]), nrow = length(observed.comp0[, 1]), ncol = 2, byrow = T)
    t1.obs.com <- as.numeric(as.character(observed.comp[, 1]))
    t2.obs.com <- as.numeric(as.character(observed.comp[, 2]))
    obs.comp <- paste0(t2.obs.com, "vs", t1.obs.com)


    ## Keep only comparisons with the baseline intervention
    indic0 <- list()
    for(i in 1:ns) {
      indic0[[i]] <- combn(t(na.omit(t(t[i, ]))), 2)[, 1:(na[i] - 1)]
    }
    (indic <- unique(t(do.call(cbind, indic0))))
    t1.indic <- indic[, 1]
    t2.indic <- indic[, 2]
    N.obs <- length(t1.indic)


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



    if (max(na) > 2 & has_error(improved.UME(t, m, N, ns, na), silent = T) == F) {

      impr.UME <- improved.UME(t, m, N, ns, na)

      ## Condition for the data specification based on the assumption about the structure of the missingness parameter
      if (assumption == "IND-CORR") {

        data.jag <- list("r" = r, "m" = m, "N" = N, "t" = t, "na" = na, "ns" = ns, "M" = M, "cov.phi" = cov.misspar, "var.phi" = var.misspar, "heter.prior" = heter.prior, "t1" = t1.indic, "t2" = t2.indic, "N.obs" = N.obs, "t1.bn" = impr.UME$t1.bn, "t2.bn" = impr.UME$t2.bn, "base" = impr.UME$base, "nbase.multi" = impr.UME$nbase.multi)

      } else if (assumption == "HIE-ARM" || assumption == "IDE-ARM") {

        data.jag <- list("r" = r, "m" = m, "N" = N, "t" = t, "na" = na, "nt" = nt, "ns" = ns, "ref" = ref, "meand.phi" = mean.misspar, "precd.phi" = prec.misspar, "heter.prior" = heter.prior, "t1" = t1.indic, "t2" = t2.indic, "N.obs" = N.obs, "t1.bn" = impr.UME$t1.bn, "t2.bn" = impr.UME$t2.bn, "base" = impr.UME$base, "nbase.multi" = impr.UME$nbase.multi)

      } else if (assumption != "HIE-ARM" || assumption != "IDE-ARM" || assumption != "IND-CORR") {

        data.jag <- list("r" = r, "m" = m, "N" = N, "t" = t, "na" = na, "ns" = ns, "meand.phi" = mean.misspar, "precd.phi" = prec.misspar, "heter.prior" = heter.prior, "t1" = t1.indic, "t2" = t2.indic, "N.obs" = N.obs, "t1.bn" = impr.UME$t1.bn, "t2.bn" = impr.UME$t2.bn, "base" = impr.UME$base, "nbase.multi" = impr.UME$nbase.multi)

      }

    } else if (max(na) < 3 || has_error(improved.UME(t, m, N, ns, na), silent = T) == T) {

      if (assumption == "IND-CORR") {

        data.jag <- list("r" = r, "m" = m, "N" = N, "t" = t, "na" = na, "ns" = ns, "M" = M, "cov.phi" = cov.misspar, "var.phi" = var.misspar, "heter.prior" = heter.prior, "t1" = t1.indic, "t2" = t2.indic, "N.obs" = N.obs, "nbase.multi" = 0)

      } else if (assumption == "HIE-ARM" || assumption == "IDE-ARM") {

        data.jag <- list("r" = r, "m" = m, "N" = N, "t" = t, "na" = na, "nt" = nt, "ns" = ns, "ref" = ref, "meand.phi" = mean.misspar, "precd.phi" = prec.misspar, "heter.prior" = heter.prior, "t1" = t1.indic, "t2" = t2.indic, "N.obs" = N.obs, "nbase.multi" = 0)

      } else if (assumption != "HIE-ARM" || assumption != "IDE-ARM" || assumption != "IND-CORR") {

        data.jag <- list("r" = r, "m" = m, "N" = N, "t" = t, "na" = na, "ns" = ns, "meand.phi" = mean.misspar, "precd.phi" = prec.misspar, "heter.prior" = heter.prior, "t1" = t1.indic, "t2" = t2.indic, "N.obs" = N.obs, "nbase.multi" = 0)

      }

    }


  }



  ## Define the nodes to be monitored
  param.jags <- c("EM", "dev.o", "resdev.o", "totresdev.o", "dev.m", "resdev.m", "totresdev.m", "tau", "hat.par", "hat.m")



  ## Run the Bayesian analysis
  jagsfit <- jags(data = data.jag, parameters.to.save = param.jags, model.file = paste0("./model/UME RE model/RE-UME_", measure, "_Pattern-mixture_", assumption, ".txt"),
                  n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin, DIC = T)




  ## Turn summary of posterior results (R2jags object) into a data-frame to select model parameters (using 'dplyr')
  getResults <- as.data.frame(t(jagsfit$BUGSoutput$summary))

  # Effect size of all unique pairwise comparisons
  EM <- t(getResults %>% dplyr::select(starts_with("EM[")))

  # Between-trial standard deviation
  tau <- t(getResults %>% dplyr::select(starts_with("tau")))

  # Trial-arm deviance contribution for observed outcome
  dev.o <- t(getResults %>% dplyr::select(starts_with("dev.o")))

  # Trial-arm deviance contribution for missing outcome data
  dev.m <- t(getResults %>% dplyr::select(starts_with("dev.m")))

  # Fitted/predicted number of missing outcome data
  hat.m <- t(getResults %>% dplyr::select(starts_with("hat.m")))

  # Fitted/predicted outcome
  hat.par <- t(getResults %>% dplyr::select(starts_with("hat.par")))

  # Deviance information criterion (as obtained from 'R2jags')
  DIC <- jagsfit$BUGSoutput$DIC

  # Number of effective parameters obtained as pD = var(deviance)/2
  pD <- jagsfit$BUGSoutput$pD

  # Total residual deviance (as obtained from 'R2jags')
  dev <- jagsfit$BUGSoutput$summary["deviance", "50%"]

  # A data-frame on the measures of model assessment: DIC, pD, and total residual deviance
  model.assessment <- data.frame(DIC, pD, dev)



  ## Calculate the deviance at posterior mean of fitted values
  # Turn 'number of observed' and 'm' into a vector (first column, followed by second column, and so on)
  m.new <- suppressMessages({as.vector(na.omit(melt(m)[, 2]))})
  N.new <- suppressMessages({as.vector(na.omit(melt(N)[, 2]))})
  obs <- N.new - m.new

  # Correction for zero MOD in trial-arm
  m0 <- ifelse(m.new == 0, m.new + 0.01, m.new)

  # Deviance at the posterior mean of the fitted MOD
  dev.post.m <- 2*(m0*(log(m0) - log(as.vector(hat.m[, 1]))) + (N.new - m0)*(log(N.new - m0) - log(N.new - as.vector(hat.m[, 1]))))

  # Sign of the difference between observed and fitted MOD
  sign.dev.m <- sign(m0 - as.vector(hat.m[, 1]))

  if(measure == "MD" || measure == "SMD"|| measure == "ROM") {

    # Turn 'y0', 'se0'into a vector (first column, followed by second column, and so on)
    y0.new <- suppressMessages({as.vector(na.omit(melt(y0)[, 2]))})
    se0.new <- suppressMessages({as.vector(na.omit(melt(se0)[, 2]))})

    # Deviance at the posterior mean of the fitted mean outcome
    dev.post.o <- (y0.new - as.vector(hat.par[, 1]))*(y0.new - as.vector(hat.par[, 1]))*(1/se0.new^2)

    # Sign of the difference between observed and fitted mean outcome
    sign.dev.o <- sign(y0.new - as.vector(hat.par[, 1]))

  } else {

    # Turn 'r' and number of observed into a vector (first column, followed by second column, and so on)
    r.new <- suppressMessages({as.vector(na.omit(melt(r)[, 2]))})

    # Correction for zero events in trial-arm
    r0 <- ifelse(r.new == 0, r.new + 0.01, ifelse(r.new == obs, r.new - 0.01, r.new))

    # Deviance at the posterior mean of the fitted response
    dev.post.o <- 2*(r0*(log(r0) - log(as.vector(hat.par[, 1]))) + (obs - r0)*(log(obs - r0) - log(obs - as.vector(hat.par[, 1]))))

    # Sign of the difference between observed and fitted response
    sign.dev.o <- sign(r0 - as.vector(hat.par[, 1]))

  }


  ## Obtain the leverage for observed and missing outcomes
  leverage.o <- as.vector(dev.o[, 1]) - dev.post.o
  leverage.m <- as.vector(dev.m[, 1]) - dev.post.m


  return(list(EM = EM, dev.m = dev.m, dev.o = dev.o, hat.m = hat.m, hat.par = hat.par, leverage.o = leverage.o, sign.dev.o = sign.dev.o, leverage.m = leverage.m, sign.dev.m = sign.dev.m, tau = tau, model.assessment = model.assessment, obs.comp = obs.comp))

}


#' Perform sensitivity meta-analysis for an aggregate binary or continuous outcome with missing participant data
#'
#' @param data A data-frame of a one-trial-per-row format containing arm-level data of each trial. This format is widely used for BUGS models. See 'Format' for the specification of the columns.
#' @param measure Character string indicating the effect measure with values \code{"OR"}, \code{"MD"}, \code{"SMD"}, or \code{"ROM"}.
#' @param assumption Character string indicating the structure of the informative missingness parameter. Set \code{assumption} equal to one of the following:  \code{"HIE-ARM"}, or \code{"IDE-ARM"}.
#' @param mean.misspar A positive non-zero number for the mean of the normal distribution of the informative missingness parameter.
#' @param var.misspar A positive non-zero number for the variance of the normal distribution of the informative missingness parameter.
#' @param D A binary number for the direction of the outcome. Set \code{D = 1} for a positive outcome and \code{D = 0} for a negative outcome.
#' @param n.chains Integer specifying the number of chains for the MCMC sampling; an argument of the \code{\link[R2jags]{jags}} function.
#' @param n.iter Integer specifying the number of Markov chains for the MCMC sampling; an argument of the \code{\link[R2jags]{jags}} function.
#' @param n.burnin Integer specifying the number of iterations to discard at the beginning of the MCMC sampling; an argument of the \code{\link[R2jags]{jags}} function.
#' @param n.thin Integer specifying the thinning rate for the MCMC sampling; an argument of the \code{\link[R2jags]{jags}} function.
#'
#' @return An R2jags output on the summaries of the posterior distribution, and the Gelman–Rubin convergence diagnostic of the following parameters:
#' \describe{
#'  \item{\code{EM}}{The effect estimate of all possible comparisons of interventions.}
#'  \item{\code{SUCRA}}{The surface under the cumulative ranking curve for each intervention.}
#'  \item{\code{phi}}{The informative missingness parameter.}
#'  \item{\code{delta}}{The underlying trial-specific effect estimate. For a multi-arm trial, we estimate \emph{T-1} trial-specific effect estimates, where \emph{T} is the number of interventions in the trial.}
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
#' run.sensitivity(data = data, measure = "SMD", assumption = "IDE-COMMON", mean.misspar = 0, var.misspar = 1, D = 0, n.chains = 3, n.iter = 10000, n.burnin = 1000, n.thin = 1)
#'
#' @export
run.sensitivity <- function(data, measure, assumption, heter.prior, var.misspar, D, n.chains, n.iter, n.burnin, n.thin){



  if(measure == "MD" || measure == "SMD"|| measure == "ROM"){

    ## Continuous: arm-level, wide-format dataset
    (y0 <- data %>% dplyr::select(starts_with("y")))            # Observed mean value in each arm of every trial
    (sd0 <- data %>% dplyr::select(starts_with("sd")))          # Observed standard deviation in each arm of every trial
    (m <- data %>% dplyr::select(starts_with("m")))             # Number of missing participants in each arm of every trial
    (c <- data %>% dplyr::select(starts_with("c")))             # Number of completers in each arm of every trial
    (se0 <- sd0/sqrt(c))                                        # Observed standard error in each arm of every trial
    (N <- m + c)                                                # Number of randomised participants in each arm of every trial
    (t <- data %>% dplyr::select(starts_with("t")))             # Intervention studied in each arm of every trial
    na <- apply(t, 1, function(x) length(which(!is.na(x))))     # Number of interventions investigated in every trial per network
    nt <- length(table(as.matrix(t)))                           # Total number of interventions per network
    ns <- length(y0[, 1])                                       # Total number of included trials per network
    ref <- 1                                                                        # The first intervention (t1 = 1) is the reference of the network
    # Trial-specific observed pooled standard deviation
    (sigma <- sqrt(apply((sd0^2)*(c - 1), 1, sum, na.rm = T)/(apply(c, 1, sum, na.rm = T) - na)))



    ## A 2x2 matrix of 25 reference-specific scenarios (PMID: 30223064)
    (scenarios <- c(-2, -1, 0, 1, 2))
    (mean.misspar <- as.matrix(cbind(rep(scenarios, each = 5), rep(scenarios, 5)))) # 2nd column refers to the reference intervention (control in MA)


    ## Information for the prior distribution on the missingness parameter (IMDOM or logIMROM)
    prec.misspar <- 1/var.misspar
    psi.misspar <- sqrt(var.misspar)           # the lower bound of the uniform prior distribution for the prior standard deviation of the missingness parameter (hierarchical structure)



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



    ## Prepare parameters for JAGS
    jagsfit <- data.jag <- list()


    ## Parameters to save
    param.jags <- c("EM", "tau")


    ## Calculate time needed for all models
    start.time <- Sys.time()

    set.seed(123)
    memory.limit(size = 40000)
    for(i in 1:length(mean.misspar[, 1])){

      if (measure == "SMD") {

        data.jag[[i]] <- list("y.o" = y0, "se.o" = se0, "m" = m, "N" = N, "t" = t, "na" = na, "nt" = nt, "ns" = ns, "ref" = ref, "sigma" = sigma, "meand.phi" = mean.misspar[i, ], "precd.phi" = prec.misspar,
                              "D" = D, "heter.prior" = heter.prior, "eff.mod" = rep(0, ns))


      } else {

        data.jag[[i]] <- list("y.o" = y0, "se.o" = se0, "m" = m, "N" = N, "t" = t, "na" = na, "nt" = nt, "ns" = ns, "ref" = ref, "meand.phi" = mean.misspar[i, ], "precd.phi" = prec.misspar,
                              "D" = D, "heter.prior" = heter.prior, "eff.mod" = rep(0, ns))

      }


      jagsfit[[i]] <- jags(data = data.jag[[i]], parameters.to.save = param.jags, model.file = paste0("./model/Full RE-NMA/Full RE-NMA_", measure, "_Pattern-mixture_", assumption, ".txt"),
                           n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin, DIC = F)
    }

    end.time <- Sys.time()
    time.taken <- end.time - start.time



  } else {

    ## Binary: arm-level, wide-format dataset
    (r <- data %>% dplyr::select(starts_with("r")))             # Number of observed events in each arm of every trial
    (m <- data %>% dplyr::select(starts_with("m")))             # Number of missing participants in each arm of every trial
    (N <- data %>% dplyr::select(starts_with("n")))             # Number randomised participants in each arm of every trial
    (t <- data %>% dplyr::select(starts_with("t")))             # Intervention studied in each arm of every trial
    na <- apply(t, 1, function(x) length(which(!is.na(x))))     # Number of interventions investigated in every trial per network
    nt <- length(table(as.matrix(t)))                           # Total number of interventions per network
    ns <- length(r[, 1])                                        # Total number of included trials per network
    ref <- 1                                                                        # The first intervention (t1 = 1) is the reference of the network



    ## A 2x2 matrix of 25 reference-specific scenarios (PMID: 30223064)
    (scenarios <- c(-log(3), -log(2), log(0.9999), log(2), log(3)))
    (mean.misspar <- as.matrix(cbind(rep(scenarios, each = 5), rep(scenarios, 5)))) # 2nd column refers to the reference intervention
    prec.misspar <- 1/var.misspar
    psi.misspar <- sqrt(var.misspar)           # the lower bound of the uniform prior distribution for the prior standard deviation of the missingness parameter (hierarchical structure)



    ## Specification of the prior distribution for the between-trial parameter
    if (heter.prior[[1]] == "halfnormal") {

      heter.prior <- as.numeric(c(0, heter.prior[[3]], 1))

    } else if (heter.prior[[1]] == "uniform") {

      heter.prior <- as.numeric(c(0, heter.prior[[3]], 2))

    } else if (heter.prior[[1]] == "lognormal")  {

      heter.prior <- as.numeric(c(heter.prior[[2]], heter.prior[[3]], 3))

    }



    ## Prepare parameters for JAGS
    jagsfit <- data.jag <- list()


    ## Parameters to save
    param.jags <- c("EM", "tau")


    ## Calculate time needed for all models
    start.time <- Sys.time()

    set.seed(123)
    memory.limit(size = 40000)
    for(i in 1:length(mean.misspar[, 1])){

      data.jag[[i]] <- list("r" = r, "m" = m, "N" = N, "t" = t, "na" = na, "nt" = nt, "ns" = ns, "ref" = ref, "meand.phi" = mean.misspar[i, ], "precd.phi" = prec.misspar, "D" = D,
                            "heter.prior" = heter.prior, "eff.mod" = rep(0, ns))


      jagsfit[[i]] <- jags(data = data.jag[[i]], parameters.to.save = param.jags, model.file = paste0("./model/Full RE-NMA/Full RE-NMA_", measure, "_Pattern-mixture_", assumption, ".txt"),
                           n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin, DIC = F)
    }

    end.time <- Sys.time()
    time.taken <- end.time - start.time

  }


  ## Obtain the posterior distribution of the necessary model paramters
  EM <- do.call(rbind,lapply(1:length(mean.misspar[, 1]), function(i) jagsfit[[i]]$BUGSoutput$summary[1:(nt*(nt - 1)*0.5), c("mean", "sd", "2.5%", "97.5%", "Rhat", "n.eff")]))
  tau <- do.call(rbind,lapply(1:length(mean.misspar[, 1]), function(i) jagsfit[[i]]$BUGSoutput$summary["tau", c("50%", "sd", "2.5%", "97.5%", "Rhat", "n.eff")]))


  return(list(EM = EM, tau = tau))

}


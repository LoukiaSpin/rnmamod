#' A function to perform separate Bayesian random-effects pairwise meta-analyses for aggregate binary or continuous outcomes based on a connected network of interventions
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
run.separate.meta <- function(data, measure, assumption, heter.prior, mean.misspar, var.misspar, n.chains, n.iter, n.burnin, n.thin){


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
    (rand <- mod + c)                                                               # Number of randomised participants in each arm of every trial
    (treat <- data %>% dplyr::select(starts_with("t")))                             # Intervention studied in each arm of every trial
    na <- apply(treat, 1, function(x) length(which(!is.na(x))))                     # Number of interventions investigated in every trial per network
    nt <- length(table(as.matrix(treat)))                                           # Total number of interventions per network
    ns <- length(y.obs[, 1])                                                        # Total number of included trials per network
    # Trial-specific observed pooled standard deviation
    (sigma <- sqrt(apply((sd.obs^2)*(c - 1), 1, sum, na.rm = T)/(apply(c, 1, sum, na.rm = T) - na)))


    ## Order by 'id of t1' < 'id of t1'
    y0 <- sd0 <- m <- N <- t <- t0 <- treat
    for(i in 1:ns){

      t0[i, ] <- order(treat[i, ], na.last = T)
      y0[i, ] <- y.obs[i, order(t0[i, ], na.last = T)]
      sd0[i, ] <- sd.obs[i, order(t0[i, ], na.last = T)]
      m[i, ] <- mod[i, order(t0[i, ], na.last = T)]
      N[i, ] <- rand[i, order(t0[i, ], na.last = T)]
      t[i, ] <- sort(treat[i, ], na.last = T)
    }


    ## Turn into contrast-level data: one row per possible comparison in each trial ('netmeta')
    # Maintain study-id, intervention, observed mean outcome, observed standard deviation, and number randomised
    (pairwise.observed <- pairwise(as.list(t), mean = as.list(y0), sd = as.list(sd0), n = as.list(N), data = data, studlab = 1:ns)[, c(3:5, 7, 10, 8, 11, 6, 9)])
    colnames(pairwise.observed) <- c("study", "arm1", "arm2", "y1", "y2", "sd1", "sd2", "n1", "n2")


    # Maintain MOD and merge with 'pairwise.observed'
    (pairwise.mod <- pairwise(as.list(t), mean = as.list(y0), sd = as.list(sd0), n = as.list(m), data = data, studlab = 1:ns)[, c(6, 9)])
    colnames(pairwise.mod) <- c("m1", "m2")


  } else {


    ## Binary: arm-level, wide-format dataset
    (event <- data %>% dplyr::select(starts_with("r")))                             # Number of observed events in each arm of every trial
    (mod <- data %>% dplyr::select(starts_with("m")))                               # Number of missing participants in each arm of every trial
    (rand <- data %>% dplyr::select(starts_with("n")))                              # Number randomised participants in each arm of every trial
    (treat <- data %>% dplyr::select(starts_with("t")))                             # Intervention studied in each arm of every trial
    na <- apply(treat, 1, function(x) length(which(!is.na(x))))                     # Number of interventions investigated in every trial per network
    nt <- length(table(as.matrix(treat)))                                           # Total number of interventions per network
    ns <- length(event[, 1])                                                        # Total number of included trials per network


    ## Order by 'id of t1' < 'id of t1'
    r <- m <- N <- t <- t0 <- treat
    for(i in 1:ns){

      t0[i, ] <- order(treat[i, ], na.last = T)
      r[i, ] <- event[i, order(t0[i, ], na.last = T)]
      m[i, ] <- mod[i, order(t0[i, ], na.last = T)]
      N[i, ] <- rand[i, order(t0[i, ], na.last = T)]
      t[i, ] <- sort(treat[i, ], na.last = T)
    }


    ## Turn into contrast-level data: one row per possible comparison in each trial ('netmeta')
    # Maintain study-id, intervention, observed events, and number randomised
    (pairwise.observed <- pairwise(as.list(t), event = as.list(r), n = as.list(N), data = data, studlab = 1:ns)[, c(3:6, 8, 7, 9)])
    colnames(pairwise.observed) <- c("study", "arm1", "arm2", "r1", "r2", "n1", "n2")


    # Maintain MOD and merge with 'pairwise.observed'
    (pairwise.mod <- pairwise(as.list(t), event = as.list(m), n = as.list(N), data = data, studlab = 1:ns)[, c(6, 8)])
    colnames(pairwise.mod) <- c("m1", "m2")

  }



  ## The dataset for the analysis
  pairwise <- data.frame(pairwise.observed, pairwise.mod)
  pairwise$t1 <- rep(1, dim(pairwise)[1])
  pairwise$t2 <- rep(2, dim(pairwise)[1])


  ## Unique comparisons with the baseline intervention
  ## A function to extract numbers from a character. Source: http://stla.github.io/stlapblog/posts/Numextract.html
  Numextract <- function(string){
    unlist(regmatches(string,gregexpr("[[:digit:]]+\\.*[[:digit:]]*",string)))
  }


  ## Observed comparisons in the network
  comp <- as.data.frame(table(paste0(pairwise$arm1, "vs", pairwise$arm2)))
  colnames(comp) <- c("comparison", "frequency")


  ## Keep comparisons withat least two trials
  keep.comp0 <- subset(comp, frequency > 1)
  keep.comp <- matrix(as.numeric(Numextract(keep.comp0[, 1])), nrow = dim(keep.comp0)[1], ncol = 2, byrow = T)
  N.comp <- dim(keep.comp)[1]


  ## Run each random-effects paiwise meta-analysis
  meta <- list()
  for(i in 1:N.comp) {

    meta[[i]] <- run.model(data = pairwise[pairwise$arm1 == keep.comp[i, 1] & pairwise$arm2 == keep.comp[i, 2], -c(1:3)], measure, assumption, heter.prior, net.ref = 1, mean.misspar, var.misspar, D = 1, n.chains, n.iter, n.burnin, n.thin) # 'D' does not matter in pairwise meta-analysis

  }

  EM <- data.frame(keep.comp, do.call(rbind, lapply(1:N.comp, function(i) meta[[i]]$EM)))
  colnames(EM) <- c("t1", "t2", "mean", "sd", "2.5%", "25%", "50%", "75%", "97.5%", "Rhat", "n.eff")
  rownames(EM) <- NULL

  tau <- data.frame(keep.comp, do.call(rbind, lapply(1:N.comp, function(i) meta[[i]]$tau)))
  colnames(tau) <- c("t1", "t2", "median", "sd", "2.5%", "25%", "50%", "75%", "97.5%", "Rhat", "n.eff")
  rownames(tau) <- NULL

  return(list(EM = EM, tau = tau))


}


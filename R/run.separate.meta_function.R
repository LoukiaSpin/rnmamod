#' A function to perform separate Bayesian random-effects pairwise meta-analyses for aggregate binary or continuous outcomes based on a connected network of interventions
#'
#' @param data A data-frame of a one-trial-per-row format containing arm-level data of each trial. This format is widely used for BUGS models. See 'Format' for the specification of the columns.
#' @param net An object of S3 class \code{\link{run.model}}. See 'Value' in \code{\link{run.model}}.
#' @param n.chains Integer specifying the number of chains for the MCMC sampling; an argument of the \code{\link[R2jags]{jags}} function of the R-package \href{https://CRAN.R-project.org/package=R2jags}{R2jags}.
#'   The default argument is 2.
#' @param n.iter Integer specifying the number of Markov chains for the MCMC sampling; an argument of the \code{\link[R2jags]{jags}} function of the R-package \href{https://CRAN.R-project.org/package=R2jags}{R2jags}.
#'   The default argument is 10000.
#' @param n.burnin Integer specifying the number of iterations to discard at the beginning of the MCMC sampling; an argument of the \code{\link[R2jags]{jags}} function of the R-package \href{https://CRAN.R-project.org/package=R2jags}{R2jags}.
#'   The default argument is 1000.
#' @param n.thin Integer specifying the thinning rate for the MCMC sampling; an argument of the \code{\link[R2jags]{jags}} function of the R-package \href{https://CRAN.R-project.org/package=R2jags}{R2jags}.
#'   The default argument is 1.
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
#' @author {Loukia M. Spineli}
#'
#' @seealso \code{\link{R2jags}}
#'
#' @examples
#' \dontshow{
#' load("./data/nma.baker2009.RData")
#' }
#'
#' # Perform a random-effects NMA with consistency equations for the odds ratio (in the logarithmic scale) assuming missing at random for identical, intervention-specific informative missingness odds ratio.
#' res1 <- run.model(data = nma.baker2009, measure = "OR", model = "RE", assumption = "IDE-ARM", heter.prior = list("halfnormal", 0, 1), mean.misspar = 0, var.misspar = 1, D = 1, n.chains = 3, n.iter = 10000, n.burnin = 1000, n.thin = 1)
#'
#' # Run separate random-effects pairwise meta-analyses using the same arguments with the network meta-analysis above.
#' run.separate.meta(data = data1, net = res1, n.chains = 3, n.iter = 10000, n.burnin = 1000, n.thin = 1)
#'
#' @export
run.separate.meta <- function(data, net, n.chains, n.iter, n.burnin, n.thin){


  measure <- net$measure
  model <- net$model
  assumption <- net$assumption
  heter.prior <- net$heter.prior
  mean.misspar <- net$mean.misspar
  var.misspar <- net$var.misspar

  ## Turn off warning when variables in the 'data.jag' are not used
  options(warn = -1)


  ## Prepare the dataset for the R2jags
  item <- data.preparation(data, measure)
  if(item$nt < 3) {
    stop("This function is *not* relevant for a pairwise meta-analysis", call. = F)
  }



  ## For a continuous outcome
  if (is.element(measure, c("MD", "SMD", "ROM"))) {

    ## Turn into contrast-level data: one row per possible comparison in each trial ('netmeta')
    # Maintain study-id, intervention, observed mean outcome, observed standard deviation, and number randomised
    (pairwise.observed <- pairwise(as.list(item$t), mean = as.list(item$y0), sd = as.list(item$sd0), n = as.list(item$N), data = data, studlab = 1:item$ns)[, c(3:5, 7, 10, 8, 11, 6, 9)])
    colnames(pairwise.observed) <- c("study", "arm1", "arm2", "y1", "y2", "sd1", "sd2", "n1", "n2")


    # Maintain MOD and merge with 'pairwise.observed'
    (pairwise.mod <- pairwise(as.list(item$t), mean = as.list(item$y0), sd = as.list(item$sd0), n = as.list(item$m), data = data, studlab = 1:item$ns)[, c(6, 9)])
    colnames(pairwise.mod) <- c("m1", "m2")


  } else {

    ## Turn into contrast-level data: one row per possible comparison in each trial ('netmeta')
    # Maintain study-id, intervention, observed events, and number randomised
    (pairwise.observed <- pairwise(as.list(item$t), event = as.list(item$r), n = as.list(item$N), data = data, studlab = 1:item$ns)[, c(3:6, 8, 7, 9)])
    colnames(pairwise.observed) <- c("study", "arm1", "arm2", "r1", "r2", "n1", "n2")


    # Maintain MOD and merge with 'pairwise.observed'
    (pairwise.mod <- pairwise(as.list(item$t), event = as.list(item$m), n = as.list(item$N), data = data, studlab = 1:item$ns)[, c(6, 8)])
    colnames(pairwise.mod) <- c("m1", "m2")

  }


  ## The dataset for the analysis
  pairwise.data <- data.frame(pairwise.observed, pairwise.mod)
  pairwise.data$t1 <- rep(1, dim(pairwise.data)[1])
  pairwise.data$t2 <- rep(2, dim(pairwise.data)[1])


  ## Unique comparisons with the baseline intervention
  ## A function to extract numbers from a character. Source: http://stla.github.io/stlapblog/posts/Numextract.html
  Numextract <- function(string){
    unlist(regmatches(string,gregexpr("[[:digit:]]+\\.*[[:digit:]]*",string)))
  }


  ## Observed comparisons in the network
  comp <- as.data.frame(table(paste0(pairwise.data$arm1, "vs", pairwise.data$arm2)))
  colnames(comp) <- c("comparison", "frequency")


  ## Keep comparisons withat least two trials
  keep.comp0 <- subset(comp, frequency > 1)
  keep.comp <- matrix(as.numeric(Numextract(keep.comp0[, 1])), nrow = dim(keep.comp0)[1], ncol = 2, byrow = T)
  N.comp <- dim(keep.comp)[1]


  ## Run each random-effects pairwise meta-analysis
  meta <- list()
  for (i in 1:N.comp) {
    message(paste(i, "out of", N.comp, "observed comparisons"))
    meta[[i]] <- run.model(data = pairwise.data[pairwise.data$arm1 == keep.comp[i, 1] & pairwise.data$arm2 == keep.comp[i, 2], ], measure, model, assumption, heter.prior, mean.misspar, var.misspar, D = 1, n.chains, n.iter, n.burnin, n.thin) # 'D' does not matter in pairwise meta-analysis
  }

  EM <- data.frame(keep.comp, do.call(rbind, lapply(1:N.comp, function(i) meta[[i]]$EM)))
  colnames(EM) <- c("t1", "t2", "mean", "sd", "2.5%", "25%", "50%", "75%", "97.5%", "Rhat", "n.eff")
  rownames(EM) <- NULL

  if (model == "RE") {
    tau <- data.frame(keep.comp, do.call(rbind, lapply(1:N.comp, function(i) meta[[i]]$tau)))
    colnames(tau) <- c("t1", "t2", "median", "sd", "2.5%", "25%", "50%", "75%", "97.5%", "Rhat", "n.eff")
    rownames(tau) <- NULL
  } else {
    tau <- NA
  }


  ## Return results based on the model
  return.results <- if (model == "RE") {
    list(EM = EM, tau = tau, measure = measure, model = model)
  } else {
    list(EM = EM, measure = measure, model = model)
  }

  return(return.results)
}


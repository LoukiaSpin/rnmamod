#' A series of Bayesian pairwise meta-analyses from a network of interventions
#'
#' @description This function performs a Bayesian pairwise meta-analysis separately for pairwise comparisons with at least two trials observed in the investigated network of interventions.
#'
#' @param data A data-frame of a one-trial-per-row format containing arm-level data of each trial. This format is widely used for BUGS models. See 'Format' in \code{\link{run.model}} for the specification of the columns.
#' @param full An object of S3 class \code{\link{run.model}}. See 'Value' in \code{\link{run.model}}.
#' @param n.chains Integer specifying the number of chains for the MCMC sampling; an argument of the \code{\link[R2jags]{jags}} function of the R-package \href{https://CRAN.R-project.org/package=R2jags}{R2jags}.
#'   The default argument is 2.
#' @param n.iter Integer specifying the number of Markov chains for the MCMC sampling; an argument of the \code{\link[R2jags]{jags}} function of the R-package \href{https://CRAN.R-project.org/package=R2jags}{R2jags}.
#'   The default argument is 10000.
#' @param n.burnin Integer specifying the number of iterations to discard at the beginning of the MCMC sampling; an argument of the \code{\link[R2jags]{jags}} function of the R-package \href{https://CRAN.R-project.org/package=R2jags}{R2jags}.
#'   The default argument is 1000.
#' @param n.thin Integer specifying the thinning rate for the MCMC sampling; an argument of the \code{\link[R2jags]{jags}} function of the R-package \href{https://CRAN.R-project.org/package=R2jags}{R2jags}.
#'   The default argument is 1.
#'
#' @return An R2jags output on the summaries of the posterior distribution, and the Gelman-Rubin convergence diagnostic of the following monitored parameters:
#' \tabular{ll}{
#'  \code{EM} \tab The summary effect estimate of each pairwise comparison with at least two trials observed in the network.\cr
#'  \tab \cr
#'  \code{tau} \tab The between-trial standard deviation for pairwise comparisons with at least two trials, when the random-effects model has been specified.\cr
#' }
#'
#' @details \code{run.series.meta} does not contain the arguments \code{measure}, \code{model}, \code{assumption}, \code{heter.prior}, \code{mean.misspar}, and \code{var.misspar} that are found in \code{run.model}.
#'   This is to prevent misspecifying the Bayesian model as it would make the comparison of the consistency model (via \code{run.model}) with the separate pairwise meta-analyses for the observed comparisons meaningless.
#'   Instead, these arguments are contained in the argument \code{full} of the function. Therefore, the user needs first to apply \code{run.model}, and then use \code{run.series.meta} (see, 'Examples').
#'
#'   \code{run.series.meta} runs Bayesian pairwise meta-analysis in \code{JAGS}. The progress of the simulation appears in the R console. The number of times \code{run.series.meta} is used appears in the R console as a text in red
#'   and it equals the number of pairwise comparisons observed in the network of interventions (see 'Examples').
#'
#'   The output of \code{run.series.meta} is not end-user-ready. The \code{series.meta.plot} function uses the output of \code{run.series.meta} as an S3 object and processes it further to provide an end-user-ready output.
#'
#'   \code{run.series.meta} can be used only for a network of interventions. In the case of two interventions, the execution of the function will be stopped and an error message will be printed in the R console.
#'
#' @author {Loukia M. Spineli}
#'
#' @seealso \code{\link{run.model}}, \code{\link{series.meta.plot}}, \code{\link[R2jags]{jags}}
#'
#' @references
#' Gelman A, Rubin DB. Inference from iterative simulation using multiple sequences. \emph{Stat Sci} 1992;\bold{7}:457--472. [\doi{10.1214/ss/1177011136}]
#'
#' @examples
#' data("nma.baker2009.RData")
#'
#' # Perform a random-effects network meta-analysis
#' res1 <- run.model(data = nma.baker2009, measure = "OR", model = "RE", assumption = "IDE-ARM", heter.prior = list("halfnormal", 0, 1), mean.misspar = 0, var.misspar = 1, D = 1, n.chains = 3, n.iter = 10000, n.burnin = 1000, n.thin = 1)
#'
#' # Run separate random-effects pairwise meta-analyses
#' run.series.meta(data = nma.baker2009, full = res1, n.chains = 3, n.iter = 10000, n.burnin = 1000, n.thin = 1)
#'
#' @export
run.series.meta <- function(data, full, n.chains, n.iter, n.burnin, n.thin) {


  measure <- full$measure
  model <- full$model
  assumption <- full$assumption
  heter.prior <- full$heter.prior
  mean.misspar <- full$mean.misspar
  var.misspar <- full$var.misspar

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
    list(EM = EM, tau = tau)
  } else {
    list(EM = EM)
  }

  return(return.results)
}


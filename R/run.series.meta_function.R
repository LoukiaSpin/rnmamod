#' Perform a series of Bayesian pairwise meta-analyses
#'
#' @description Performs a Bayesian pairwise meta-analysis for each pairwise
#'   comparison with at least two trials in the network.
#'
#' @param full An object of S3 class \code{\link{run_model}}. See 'Value' in
#'   \code{\link{run_model}}.
#' @param n_chains Integer specifying the number of chains for the MCMC
#'   sampling; an argument of the \code{\link[R2jags:jags]{jags}} function of
#'   the R-package \href{https://CRAN.R-project.org/package=R2jags}{R2jags}.
#'   The default argument is 2.
#' @param n_iter Positive integer specifying the number of Markov chains for the
#'   MCMC sampling; an argument of the \code{\link[R2jags:jags]{jags}} function
#'   of the R-package \href{https://CRAN.R-project.org/package=R2jags}{R2jags}.
#'   The default argument is 10000.
#' @param n_burnin Positive integer specifying the number of iterations to
#'   discard at the beginning of the MCMC sampling; an argument of the
#'   \code{\link[R2jags:jags]{jags}} function of the R-package
#'   \href{https://CRAN.R-project.org/package=R2jags}{R2jags}.
#'   The default argument is 1000.
#' @param n_thin Positive integer specifying the thinning rate for the MCMC
#'   sampling; an argument of the \code{\link[R2jags:jags]{jags}} function of
#'   the R-package \href{https://CRAN.R-project.org/package=R2jags}{R2jags}.
#'   The default argument is 1.
#'
#' @return An R2jags output on the summaries of the posterior distribution, and
#'   the Gelman-Rubin convergence diagnostic (Gelman et al., 1992) of the
#'   following monitored parameters:
#'   \item{EM}{The summary effect estimate (according to the argument
#'   \code{measure} defined in \code{\link{run_model}}) of each observed
#'   pairwise comparison with at least two trials in the network.}
#'   \item{tau}{The between-trial standard deviation for pairwise comparisons
#'   with at least two trials, when the random-effects model has been
#'   specified.}
#'   \item{single}{A binary vector that indicates the comparisons in \code{EM}
#'   with one trial.}
#'
#' @details \code{run_series_meta} inherits the arguments \code{data},
#'   \code{measure}, \code{model}, \code{assumption}, \code{heter_prior},
#'   \code{mean_misspar}, and \code{var_misspar} from \code{\link{run_model}}
#'   (now contained in the argument \code{full}). This prevents specifying
#'   a different Bayesian model from that considered in \code{\link{run_model}}.
#'   Therefore, the user needs first to apply \code{\link{run_model}}, and then
#'   use \code{run_series_meta} (see 'Examples').
#'
#'   For a binary outcome, when \code{measure} is "RR" (relative risk) or "RD"
#'   (risk difference) in \code{\link{run_model}}, \code{run_series_meta}
#'   currently performs a series of pairwise meta-analysis using the odds ratio
#'   as effect measure for being the \strong{base-case} effect measure in
#'   \code{\link{run_model}} for a binary outcome (see also 'Details' in
#'   \code{\link{run_model}}).
#'
#'   \code{run_series_meta} runs a series of Bayesian pairwise meta-analyses
#'   in \code{JAGS}. The progress of the simulation appears on the R console.
#'   The number of times the function is used is also printed on the console
#'   (in red) and is equal to the number of  observed pairwise comparisons
#'   in the network (see 'Examples').
#'
#'   The output of \code{run_series_meta} is not end-user-ready. The
#'   \code{\link{series_meta_plot}} function inherits the output of
#'   \code{run_series_meta} as an S3 object and processes it further to provide
#'   an end-user-ready output.
#'
#'   \code{run_series_meta} can be used only for a network of interventions.
#'   In the case of two interventions, the execution of the function will
#'   be stopped and an error message will be printed on the R console.
#'
#' @author {Loukia M. Spineli}
#'
#' @seealso \code{\link[R2jags:jags]{jags}},
#'   \code{\link{run_model}}, \code{\link{series_meta_plot}}
#'
#' @references
#' Gelman A, Rubin DB. Inference from iterative simulation using multiple
#' sequences. \emph{Stat Sci} 1992;\bold{7}:457--472.
#'
#' @examples
#' data("nma.dogliotti2014")
#'
#' # Show the first six trials of the dataset (one-trial-per-row format)
#' head(nma.dogliotti2014)
#'
#' # Read results from 'run_model' (using the default arguments)
#' res <- readRDS(system.file('extdata/res_dogliotti.rds', package = 'rnmamod'))
#'
#' \donttest{
#' # Run separate random-effects pairwise meta-analyses
#' # Note: Ideally, set 'n_iter' to 10000 and 'n_burnin' to 1000
#' run_series_meta(full = res,
#'                 n_chains = 3,
#'                 n_iter = 1000,
#'                 n_burnin = 100,
#'                 n_thin = 1)
#' }
#'
#' @export
run_series_meta <- function(full, n_chains, n_iter, n_burnin, n_thin) {

  data <- full$data
  measure <- if (is.element(full$measure, c("RR", "RD"))) {
    "OR"
  } else {
    full$measure
  }
  model <- full$model
  assumption <- full$assumption
  heter_prior0 <- if (model == "FE") {
    list(NA, NA, NA)
  } else {
    as.list(full$heter_prior)
  }
  heter_prior0[[3]] <- if (model == "RE" & heter_prior0[[3]] == 1) {
    "halfnormal"
  } else if (model == "RE" & heter_prior0[[3]] == 2) {
    "uniform"
  } else if (model == "RE" & heter_prior0[[3]] == 3) {
    "lognormal"
  } else if (model == "RE" & heter_prior0[[3]] == 4) {
    "logt"
  } else if (model == "FE") {
    NA
  }
  heterog_prior <- if (model == "FE") {
    list(NA, NA, NA)
  } else if (model == "RE") {
    list(heter_prior0[[3]], heter_prior0[[1]], heter_prior0[[2]])
  }
  mean_misspar <- full$mean_misspar
  var_misspar <- full$var_misspar

  # Prepare the dataset for the R2jags
  item <- data_preparation(data, measure)
  if (item$nt < 3) {
    stop("This function is *not* relevant for a pairwise meta-analysis.",
         call. = FALSE)
  }

  # For a continuous outcome
  if (is.element(measure, c("MD", "SMD", "ROM"))) {
    # Turn into contrast-level data ('netmeta')
    pairwise_observed <-
      pairwise(as.list(item$t),
               mean = as.list(item$y0),
               sd = as.list(item$sd0),
               n = as.list(item$N),
               data = cbind(item$t, item$y0, item$sd0, item$N),
               studlab = 1:item$ns)[, c(3:5, 7, 10, 8, 11, 6, 9)]
    colnames(pairwise_observed) <- c("study",
                                      "arm1",
                                      "arm2",
                                      "y1",
                                      "y2",
                                      "sd1",
                                      "sd2",
                                      "n1",
                                      "n2")

    # Maintain MOD and merge with 'pairwise_observed'
    pairwise_mod <- pairwise(as.list(item$t),
                              mean = as.list(item$y0),
                              sd = as.list(item$sd0),
                              n = as.list(item$m),
                              data = cbind(item$t, item$y0, item$sd0, item$m),
                              studlab = 1:item$ns)[, c(6, 9)]
    colnames(pairwise_mod) <- c("m1", "m2")

  } else {
    # Turn into contrast-level data ('netmeta')
    pairwise_observed <-
      pairwise(as.list(item$t),
               event = as.list(item$r),
               n = as.list(item$N),
               data = cbind(item$t, item$r, item$N),
               studlab = 1:item$ns)[, c(3:6, 8, 7, 9)]
    colnames(pairwise_observed) <- c("study",
                                      "arm1",
                                      "arm2",
                                      "r1",
                                      "r2",
                                      "n1",
                                      "n2")

    # Maintain MOD and merge with 'pairwise_observed'
    pairwise_mod <- pairwise(as.list(item$t),
                             event = as.list(item$m),
                             n = as.list(item$N),
                             data = cbind(item$t, item$m, item$N),
                             studlab = 1:item$ns)[, c(6, 8)]
    colnames(pairwise_mod) <- c("m1", "m2")
  }

  # The dataset for the analysis
  pairwise_data <- data.frame(pairwise_observed, pairwise_mod)
  # Control arm
  pairwise_data$t1 <- rep(1, dim(pairwise_data)[1])
  # Experimental arm
  pairwise_data$t2 <- rep(2, dim(pairwise_data)[1])

  # A function to extract numbers from a character.
  # Source: http://stla.github.io/stlapblog/posts/numextract.html
  numextract <- function(string) {
    unlist(regmatches(string, gregexpr("[[:digit:]]+\\.*[[:digit:]]*", string)))
  }

  # Observed comparisons in the network
  comp <- as.data.frame(
    table(paste0(pairwise_data$arm1, "vs", pairwise_data$arm2)))
  colnames(comp) <- c("comparison", "frequency")

  # Indicate all observed comparisons
  obs_comp <- matrix(as.numeric(numextract(comp[, 1])),
                     nrow = dim(comp)[1],
                     ncol = 2,
                     byrow = TRUE)
  n_obs_comp <- dim(obs_comp)[1]

  # Indicate comparisons with one trial
  #keep_comp0 <- subset(comp, frequency > 1)
  #keep_comp <- matrix(as.numeric(numextract(keep_comp0[, 1])),
  #                    nrow = dim(keep_comp0)[1],
  #                    ncol = 2,
  #                    byrow = TRUE)
  #n_comp <- dim(keep_comp)[1]
  single <- ifelse (comp[, 2] < 2, 1, 0) # 1:yes, 0:no

  # Run each random-effects pairwise meta-analysis
  meta <- list()
  for (i in 1:n_obs_comp) {
    message(paste(i, "out of", n_obs_comp, "observed comparisons"))
    # 'D' and 'base_risk' do not matter in pairwise meta-analysis
    meta[[i]] <-
      run_model(data =
                  pairwise_data[pairwise_data$arm1 == obs_comp[i, 1] &
                                  pairwise_data$arm2 == obs_comp[i, 2], ],
                measure,
                model,
                assumption,
                heter_prior = heterog_prior,
                mean_misspar,
                var_misspar,
                D = 1,
                ref = 1,
                base_risk = 0.5,
                n_chains,
                n_iter,
                n_burnin,
                n_thin)
  }

  EM <- data.frame(obs_comp,
                   do.call(rbind,
                           lapply(1:n_obs_comp, function(i) meta[[i]]$EM)))
  colnames(EM) <- c("t1", "t2", "mean", "sd", "2.5%", "25%", "50%", "75%",
                    "97.5%", "Rhat", "n.eff")
  rownames(EM) <- NULL

  if (model == "RE") {
    tau0 <- data.frame(obs_comp,
                       do.call(rbind,
                               lapply(1:n_obs_comp, function(i) meta[[i]]$tau)))
    colnames(tau0) <- c("t1", "t2", "mean", "sd", "2.5%", "25%", "50%", "75%",
                        "97.5%", "Rhat", "n.eff")
    rownames(tau0) <- NULL
    tau <- subset(tau0, single == 0)
  } else {
    tau <- NA
  }

  # Return results based on the model
  return_results <- if (model == "RE") {
    list(EM = EM,
         tau = tau,
         single = single,
         n_chains = n_chains,
         n_iter = n_iter,
         n_burnin = n_burnin,
         n_thin = n_thin,
         measure = measure)
  } else {
    list(EM = EM,
         single = single,
         n_chains = n_chains,
         n_iter = n_iter,
         n_burnin = n_burnin,
         n_thin = n_thin,
         measure = measure)
  }

  return(return_results)
}

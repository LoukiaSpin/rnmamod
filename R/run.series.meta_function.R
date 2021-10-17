#' A series of Bayesian pairwise meta-analyses from a network of interventions
#'
#' @description Performs separate Bayesian pairwise meta-analysis for each
#'   pairwise comparison with at least two trials in the network.
#'
#' @param full An object of S3 class \code{\link{run_model}}. See 'Value' in
#'   \code{\link{run_model}}.
#' @param n_chains Integer specifying the number of chains for the MCMC
#'   sampling; an argument of the \code{\link[R2jags]{jags}} function of the
#'   R-package \href{https://CRAN.R-project.org/package=R2jags}{R2jags}.
#'   The default argument is 2.
#' @param n_iter Positive integer specifying the number of Markov chains for the
#'   MCMC sampling; an argument of the \code{\link[R2jags]{jags}} function of
#'   the R-package \href{https://CRAN.R-project.org/package=R2jags}{R2jags}.
#'   The default argument is 10000.
#' @param n_burnin Positive integer specifying the number of iterations to
#'   discard at the beginning of the MCMC sampling; an argument of the
#'   \code{\link[R2jags]{jags}} function of the R-package
#'   \href{https://CRAN.R-project.org/package=R2jags}{R2jags}.
#'   The default argument is 1000.
#' @param n_thin Positive integer specifying the thinning rate for the MCMC
#'   sampling; an argument of the \code{\link[R2jags]{jags}} function of the
#'   R-package \href{https://CRAN.R-project.org/package=R2jags}{R2jags}.
#'   The default argument is 1.
#'
#' @return An R2jags output on the summaries of the posterior distribution, and
#'   the Gelman-Rubin convergence diagnostic of the following monitored
#'   parameters:
#'   \tabular{ll}{
#'    \code{EM} \tab The summary effect estimate of each pairwise comparison
#'    with at least two trials observed in the network.\cr
#'    \tab \cr
#'    \code{tau} \tab The between-trial standard deviation for pairwise
#'    comparisons with at least two trials, when the random-effects model has
#'    been specified.\cr
#'   }
#'
#' @details \code{run_series_meta} inherits the arguments \code{data},
#'   \code{measure}, \code{model}, \code{assumption}, \code{heter_prior},
#'   \code{mean_misspar}, and \code{var_misspar} from \code{run_model}, now
#'   contained in the argument \code{full}. This prevents misspecifying the
#'   Bayesian model that has been considered in \code{run_model}. Instead,
#'   these arguments are contained in the argument \code{full} of the function.
#'   Therefore, the user needs first to apply \code{run_model}, and then use
#'   \code{run_series_meta} (see, 'Examples').
#'
#'   \code{run_series_meta} performs Bayesian pairwise meta-analysis in
#'   \code{JAGS}. The progress of the simulation appears in the R console. The
#'   number of times a pairwise meta-analysis is preformed is also printed in
#'   the console (in red) and is equal to the number of pairwise comparisons
#'   available in the network (see 'Examples').
#'
#'   The output of \code{run_series_meta} is not end-user-ready. The
#'   \code{series_meta_plot} function uses the output of \code{run_series_meta}
#'   as an S3 object and processes it further to provide an end-user-ready
#'   output.
#'
#'   \code{run_series_meta} can be used only for a network of interventions.
#'
#' @author {Loukia M. Spineli}
#'
#' @seealso \code{\link{run_model}}, \code{\link{series_meta_plot}},
#'   \code{\link[R2jags]{jags}}
#'
#' @references
#' Gelman A, Rubin DB. Inference from iterative simulation using multiple
#' sequences. \emph{Stat Sci} 1992;\bold{7}:457--472.
#' [\doi{10.1214/ss/1177011136}]
#'
#' @examples
#' data("nma.dogliotti2014")
#'
#' # Show the first six trials of the dataset (one-trial-per-row format)
#' head(nma.dogliotti2014)
#' #            study t1 t2 t3   r1   r2 r3  m1  m2 m3   n1   n2 n3
#' #     BAATAF, 1990  1  7 NA  195  188 NA   0  21 NA  208  212 NA
#' #     SPINAF, 1992  1  7 NA  186  172 NA  56  81 NA  265  260 NA
#' #    SPAF-II, 1994  2  7 NA  480  478 NA  23  38 NA  545  555 NA
#' #      PATAF, 1999  2  7 NA  243   86 NA  54  42 NA  319  131 NA
#' # ACTIVE (W), 2006  3  7 NA 2825 3089 NA 410 223 NA 3335 3371 NA
#' #       JAST, 2006  1  2 NA  338  313 NA  89  96 NA  445  426 NA
#'
#' \dontrun{
#' # Perform a random-effects network meta-analysis
#' res <- run_model(data = nma.dogliotti2014,
#'                  measure = "OR",
#'                  model = "RE",
#'                  assumption = "IDE-ARM",
#'                  heter_prior = list("halfnormal", 0, 1),
#'                  mean_misspar = c(0, 0),
#'                  var_misspar = 1,
#'                  D = 1,
#'                  n_chains = 3,
#'                  n_iter = 10000,
#'                  n_burnin = 1000,
#'                  n_thin = 1)
#'
#' # Run separate random-effects pairwise meta-analyses
#' run_series_meta(full = res,
#'                 n_chains = 3,
#'                 n_iter = 10000,
#'                 n_burnin = 1000,
#'                 n_thin = 1)
#' }
#' @export
run_series_meta <- function(full, n_chains, n_iter, n_burnin, n_thin) {

  data <- full$data
  measure <- full$measure
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

  # Turn off warning when variables in the 'data.jag' are not used
  options(warn = -1)

  # Prepare the dataset for the R2jags
  item <- data_preparation(data, measure)
  if (item$nt < 3) {
    stop("This function is *not* relevant for a pairwise meta-analysis",
         call. = F)
  }

  # For a continuous outcome
  if (is.element(measure, c("MD", "SMD", "ROM"))) {
    # Turn into contrast-level data ('netmeta')
    pairwise_observed <-
      pairwise(as.list(item$t),
               mean = as.list(item$y0),
               sd = as.list(item$sd0),
               n = as.list(item$N),
               data = data,
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
                             data = data,
                             studlab = 1:item$ns)[, c(6, 9)]
    colnames(pairwise_mod) <- c("m1", "m2")
  } else {
    # Turn into contrast-level data ('netmeta')
    pairwise_observed <-
      pairwise(as.list(item$t),
               event = as.list(item$r),
               n = as.list(item$N),
               data = data,
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
                             data = data,
                             studlab = 1:item$ns)[, c(6, 8)]
    colnames(pairwise_mod) <- c("m1", "m2")
  }

  # The dataset for the analysis
  pairwise_data <- data.frame(pairwise_observed, pairwise_mod)
  pairwise_data$t1 <- rep(1, dim(pairwise_data)[1])
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

  # Keep comparisons withat least two trials
  keep_comp0 <- subset(comp, frequency > 1)
  keep_comp <- matrix(as.numeric(numextract(keep_comp0[, 1])),
                      nrow = dim(keep_comp0)[1],
                      ncol = 2,
                      byrow = T)
  n_comp <- dim(keep_comp)[1]

  # Run each random-effects pairwise meta-analysis
  meta <- list()
  for (i in 1:n_comp) {
    message(paste(i, "out of", n_comp, "observed comparisons"))
    # 'D' does not matter in pairwise meta-analysis
    meta[[i]] <-
      run_model(data =
                  pairwise_data[pairwise_data$arm1 == keep_comp[i, 1] &
                                  pairwise_data$arm2 == keep_comp[i, 2], ],
                measure,
                model,
                covar_assumption = "NO",
                assumption,
                heter_prior = heterog_prior,
                mean_misspar,
                var_misspar,
                D = 1,
                n_chains,
                n_iter,
                n_burnin,
                n_thin)
  }

  EM <- data.frame(keep_comp,
                   do.call(rbind, lapply(1:n_comp, function(i) meta[[i]]$EM)))
  colnames(EM) <- c("t1", "t2", "mean", "sd", "2.5%", "25%", "50%", "75%",
                    "97.5%", "Rhat", "n.eff")
  rownames(EM) <- NULL

  if (model == "RE") {
    tau <- data.frame(keep_comp,
                      do.call(rbind,
                              lapply(1:n_comp, function(i) meta[[i]]$tau)))
    colnames(tau) <- c("t1", "t2", "median", "sd", "2.5%", "25%", "50%", "75%",
                       "97.5%", "Rhat", "n.eff")
    rownames(tau) <- NULL
  } else {
    tau <- NA
  }

  # Return results based on the model
  return_results <- if (model == "RE") {
    list(EM = EM,
         tau = tau)
  } else {
    list(EM = EM)
  }

  return(return_results)
}

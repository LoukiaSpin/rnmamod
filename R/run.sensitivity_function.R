#' Perform sensitivity analysis for missing participant outcome data
#'
#' @description Performs a sensitivity analysis by applying pairwise
#'   meta-analysis or network meta-analysis for a series of different scenarios
#'   about the informative missingness parameter.
#'
#' @param full An object of S3 class \code{\link{run_model}}. See 'Value' in
#'   \code{\link{run_model}}.
#' @param assumption Character string indicating the structure of the
#'   informative missingness parameter. Set \code{assumption} equal to one of
#'   the following two: \code{"HIE-ARM"}, or \code{"IDE-ARM"} (see 'Details').
#'   The default argument is \code{"IDE-ARM"}. The abbreviations \code{"IDE"},
#'   and \code{"HIE"} stand for identical, and hierarchical, respectively.
#' @param mean_scenarios A vector with numeric values for the mean of the normal
#'   distribution of the informative missingness parameter (see 'Details').
#'   The vector should have a length equal to 5 or larger.
#'   The missing-at-random (MAR) assumption should be the median of the vector,
#'   so that the same number of informative scenarios appear before and after
#'   the MAR. The default scenarios are c(-log(3), -log(2), log(0.9999), log(2),
#'   log(3)) and c(-2, -1, 0, 1, 2) for binary and continuous outcome data,
#'   respectively.
#' @param var_misspar A positive non-zero number for the variance of the normal
#'   distribution of the informative missingness parameter. When the
#'   \code{measure} (defined in \code{\link{run_model}}) is \code{"OR"},
#'   \code{"MD"}, or \code{"SMD"} the default argument is 1. When the
#'   \code{measure} is \code{"ROM"}, the default argument is 0.04.
#' @param n_chains Integer specifying the number of chains for the MCMC
#'   sampling; an argument of the \code{\link[R2jags:jags]{jags}} function of
#'   the R-package \href{https://CRAN.R-project.org/package=R2jags}{R2jags}.
#'   The default argument is 2.
#' @param n_iter Integer specifying the number of Markov chains for the MCMC
#'   sampling; an argument of the \code{\link[R2jags:jags]{jags}} function of
#'   the R-package \href{https://CRAN.R-project.org/package=R2jags}{R2jags}.
#'   The default argument is 10000.
#' @param n_burnin Integer specifying the number of iterations to discard at the
#'   beginning of the MCMC sampling; an argument of the
#'   \code{\link[R2jags:jags]{jags}} function of the R-package
#'   \href{https://CRAN.R-project.org/package=R2jags}{R2jags}.
#'   The default argument is 1000.
#' @param n_thin Integer specifying the thinning rate for the MCMC sampling; an
#'   argument of the \code{\link[R2jags:jags]{jags}} function of the R-package
#'   \href{https://CRAN.R-project.org/package=R2jags}{R2jags}.
#'   The default argument is 1.
#' @param inits A list with the initial values for the parameters; an argument
#'   of the \code{\link[R2jags:jags]{jags}} function of the R-package
#'   \href{https://CRAN.R-project.org/package=R2jags}{R2jags}.
#'   The default argument is \code{NULL}, and JAGS generates the initial values.
#'
#' @return A list of R2jags outputs on the summaries of the posterior
#'   distribution, and the Gelman-Rubin convergence diagnostic
#'   (Gelman et al., 1992) of the following monitored parameters for a
#'   random-effects pairwise meta-analysis:
#'   \item{EM}{The estimated summary effect measure (according to the
#'    argument \code{measure} defined in \code{\link{run_model}}).}
#'   \item{EM_pred}{The predicted summary effect measure (according to the
#'    argument \code{measure} defined in \code{\link{run_model}}). This element
#'    does not appear in the case of a fixed-effect meta-analysis.}
#'   \item{EM_LOR}{The estimated summary odd ratio in the logarithmic scale when
#'   \code{measure = "RR"} or \code{measure = "RD"}.}
#'   \item{EM_LOR_pred}{The predicted summary odd ratio in the logarithmic scale
#'    when \code{measure = "RR"} or \code{measure = "RD"}. This element does
#'    not appear in the case of a fixed-effect meta-analysis.}
#'   \item{tau}{The between-trial standard deviation. This element does
#'    not appear in the case of a fixed-effect meta-analysis.}
#'
#'   In a random-effects network meta-analysis, \code{EM} refer to all possible
#'   pairwise comparisons of interventions in the network. Furthermore,
#'   \code{tau} is typically assumed to be common for all observed comparisons
#'   in the network.
#'
#' @details The model runs in \code{JAGS} and the progress of the simulation
#'   appears on the R console. The number of times \code{run_sensitivity} is
#'   used appears on the R console as a text in red and it equals the
#'   \bold{number of scenarios} defined as \emph{the square of the length of the
#'   vector specified in \code{mean_scenarios}} (see 'Examples').
#'   The output of \code{run_sensitivity} is used as an S3 object by other
#'   functions of the package to be processed further and provide an
#'   end-user-ready output.
#'
#'   In the case of pairwise meta-analysis, \code{EM} and \code{tau} are
#'   estimated as many times as the number of scenarios considered.
#'   In the case of network meta-analysis, each possible pairwise comparison is
#'   estimated as many times as the number of scenarios considered.
#'
#'   The informative missingness parameter is assumed to differ only across the
#'   interventions of the dataset. Therefore, the user can specify the
#'   informative missingness parameter to be arm-specific \emph{and} identical
#'   (\code{assumption = "IDE-ARM"}), or arm-specific \emph{and}
#'   hierarchical (\code{assumption = "HIE-ARM"}) (Spineli et al., 2021).
#'
#'   The length of the vector specified in argument \code{mean_scenarios} should
#'   be equal to or more than 5 (a positive odd integer) to allow for an
#'   adequate number of scenarios. It is important that the number
#'   corresponding to the MAR assumption is the middle of the numbers in the
#'   vector specified in argument \code{mean_scenarios}. The \bold{MAR
#'   assumption} constitutes the \bold{primary analysis}.
#'   Under the informative missingness difference of means parameter (relevant
#'   for the raw and standardised mean diffenre), the MAR assumption equals 0.
#'   Under the informative missingness odds ratio parameter (IMOR; relevant for
#'   the odds ratio) and the informative missingness ratio of means (IMRoM;
#'   relevant for the ratio of means) parameter, the MAR assumption equals 1;
#'   however, both parameters are analysed in the logarithmic scale. We advise
#'   using the value \code{0.999} rather than \code{1} in \code{mean_scenarios}
#'   for the IMOR and IMRoM parameters; otherwise, the execution of the function
#'   will be stopped and the error 'Invalid parent values' will be printed on
#'   the R console.
#'
#'   Currently, there are no empirically-based prior distributions for the
#'   informative missingness parameters. The users may refer to Spineli (2019),
#'   Mavridis et al. (2015), Turner et al. (2015), and White et al. (2008) to
#'   determine \code{mean_scenarios} for an informative missingness mechanism
#'   and select a proper value for \code{var_misspar}.
#'
#'   \code{run_sensitivity} inherits the arguments \code{data},
#'   \code{measure}, \code{model}, \code{heter_prior}, \code{D}, \code{indic},
#'   \code{base_risk}, and \code{ref} from \code{\link{run_model}}
#'   (now contained in the argument \code{full}). This prevents specifying a
#'   different Bayesian model from that considered in the primary analysis
#'   (via \code{\link{run_model}})--an exception in the \code{assumption}
#'   argument as it is restricted to only two character strings. Therefore, the
#'   user needs first to apply \code{\link{run_model}}, and then use
#'   \code{run_sensitivity} (see 'Examples').
#'
#'   The \code{run_sensitivity} function also returns the arguments
#'   \code{measure}, \code{scenarios}, \code{D}, \code{heter}, \code{n_chains},
#'   \code{n_iter}, \code{n_burnin}, and \code{n_thin} as specified by the user
#'   to be inherited by other relevant functions of the package.
#'
#'   The model is updated until convergence using the
#'   \code{\link[R2jags:autojags]{autojags}} function of the R-package
#'   \href{https://CRAN.R-project.org/package=R2jags}{R2jags} with 2 updates and
#'   number of iterations and thinning equal to \code{n_iter} and \code{n_thin},
#'   respectively.
#'
#'   \code{run_sensitivity} can be used only when missing participant
#'   outcome data have been extracted for at least one trial. Otherwise, the
#'   execution of the function will be stopped and an error message will be
#'   printed on the R console.
#'
#' @author {Loukia M. Spineli}
#'
#' @seealso \code{\link[R2jags:autojags]{autojags}},
#'   \code{\link[R2jags:jags]{jags}}, \code{\link{run_model}}
#'
#' @references
#' Gelman, A, Rubin, DB. Inference from iterative simulation using multiple
#' sequences. \emph{Stat Sci} 1992;\bold{7}(4):457--72.
#' doi: 10.1214/ss/1177011136
#'
#' Mavridis D, White IR, Higgins JP, Cipriani A, Salanti G. Allowing for
#' uncertainty due to missing continuous outcome data in pairwise and network
#' meta-analysis. \emph{Stat Med} 2015;\bold{34}(5):721--41.
#' doi: 10.1002/sim.6365
#'
#' Spineli LM, Kalyvas C, Papadimitropoulou K. Quantifying the robustness of
#' primary analysis results: A case study on missing outcome data in pairwise
#' and network meta-analysis.
#' \emph{Res Synth Methods} 2021;\bold{12}(4):475--90. doi: 10.1002/jrsm.1478
#'
#' Spineli LM. An empirical comparison of Bayesian modelling strategies for
#' missing binary outcome data in network meta-analysis.
#' \emph{BMC Med Res Methodol} 2019;\bold{19}(1):86.
#' doi: 10.1186/s12874-019-0731-y
#'
#' Turner NL, Dias S, Ades AE, Welton NJ. A Bayesian framework to account for
#' uncertainty due to missing binary outcome data in pairwise meta-analysis.
#' \emph{Stat Med} 2015;\bold{34}(12):2062--80. doi: 10.1002/sim.6475
#'
#' White IR, Higgins JP, Wood AM. Allowing for uncertainty due to missing data
#' in meta-analysis--part 1: two-stage methods.
#' \emph{Stat Med} 2008;\bold{27}(5):711--27. doi: 10.1002/sim.3008
#'
#' @examples
#' data("pma.taylor2004")
#'
#' # Read results from 'run_model' (using the default arguments)
#' res <- readRDS(system.file('extdata/res_taylor.rds', package = 'rnmamod'))
#'
#' \donttest{
#' # Perform the sensitivity analysis (default arguments)
#' # Note: Ideally, set 'n_iter' to 10000 and 'n_burnin' to 1000
#' run_sensitivity(full = res,
#'                 assumption = "IDE-ARM",
#'                 var_misspar = 1,
#'                 n_chains = 3,
#'                 n_iter = 1000,
#'                 n_burnin = 100,
#'                 n_thin = 5)
#' }
#'
#' @export
run_sensitivity <- function(full,
                            assumption,
                            mean_scenarios,
                            var_misspar,
                            n_chains,
                            n_iter,
                            n_burnin,
                            n_thin,
                            inits = NULL) {


  if (!inherits(full, "run_model") || is.null(full)) {
    stop("'full' must be an object of S3 class 'run_model'.",
         call. = FALSE)
  }

  data <- full$data
  measure <- full$measure
  model <- full$model
  heterog_prior <- full$heter_prior
  D <- full$D
  ref <- full$ref
  indic <- full$indic
  base_risk <- full$base_risk
  assumption <- if (missing(assumption)) {
    "IDE-ARM"
  } else if(!is.element(assumption, c("IDE-ARM", "HIE-ARM"))) {
    stop("Insert 'IDE-ARM', or, 'HIE-ARM'", call. = FALSE)
  } else {
    assumption
  }

  # Prepare the dataset for the R2jags
  item <- data_preparation(data, measure)

  if (min(na.omit(unlist(item$I))) == 0 & max(na.omit(unlist(item$I))) == 0) {
    aa <- "Missing participant outcome data have *not* been collected."
    stop(paste(aa, "This function cannot be used."), call. = FALSE)
    return(NA)
  } else if (min(na.omit(unlist(item$I))) == 0 &
             max(na.omit(unlist(item$I))) == 1) {
    aa <- "Missing participant outcome data have been collected *partially*."
    stop(paste(aa, "This function cannot be used."), call. = FALSE)
    return(NA)
  }

  # Scenarios for missingness mechanism in an intervention (PMID: 30223064)
  mean_scenarios <- if (missing(mean_scenarios) &
                        is.element(measure, c("MD", "SMD"))) {
    aa <- "The default scenarios were considered:"
    message(paste(aa, "c(-2, -1, 0, 1, 2)."))
    c(-2, -1, 0, 1, 2)
  } else if (missing(mean_scenarios) & is.element(measure,
                                                  c("OR", "RR", "RD", "ROM"))) {
    aa <- "The default scenarios were considered:"
    message(paste(aa, "c(-log(3), -log(2), log(0.9999), log(2), log(3))."))
    c(-log(3), -log(2), log(0.9999), log(2), log(3))
  } else if (length(mean_scenarios) < 5) {
    stop("The argument 'mean_scenarios' must have a length of at least 5.",
         call. = FALSE)
  } else {
    mean_scenarios
  }
  var_misspar <- if (missing(var_misspar) &
                     is.element(measure, c("OR", "RR", "RD", "MD", "SMD"))) {
    1
  } else if (missing(var_misspar) & measure == "ROM") {
    0.2^2
  } else if (var_misspar <= 0) {
    stop("The argument 'var_misspar' must be a positive non-zero number.",
         call. = FALSE)
  } else {
    var_misspar
  }
  ref_base <- if (is.element(measure, c("OR", "RR", "RD"))) {
    baseline_model(base_risk,
                   n_chains,
                   n_iter,
                   n_burnin,
                   n_thin)$ref_base
  } else if (!is.element(measure, c("OR", "RR", "RD"))) {
    NA
  }
  n_chains <- if (missing(n_chains)) {
    2
  } else if (n_chains < 1) {
    stop("The argument 'n_chains' must be a positive integer.", call. = FALSE)
  } else {
    n_chains
  }
  n_iter <- if (missing(n_iter)) {
    10000
  } else if (n_iter < 1) {
    stop("The argument 'n_iter' must be a positive integer.", call. = FALSE)
  } else {
    n_iter
  }
  n_burnin <- if (missing(n_burnin)) {
    1000
  } else if (n_burnin < 1) {
    stop("The argument 'n_burnin' must be a positive integer.", call. = FALSE)
  } else {
    n_burnin
  }
  n_thin <- if (missing(n_thin)) {
    1
  } else if (n_thin < 1) {
    stop("The argument 'n_thin' must be a positive integer.", call. = FALSE)
  } else {
    n_thin
  }
  inits <- if (is.null(inits)) {
    message("JAGS generates initial values for the parameters.")
    NULL
  } else {
    inits
  }

  # A 2x2 matrix of 25 reference-specific scenarios (PMID: 30223064)
  mean_misspar <- as.matrix(cbind(rep(mean_scenarios,
                                      each = length(mean_scenarios)),
                                  rep(mean_scenarios, length(mean_scenarios))))

  # Prepare parameters for JAGS
  jagsfit0 <- jagsfit <- data_jag <- list()

  ## Parameters to save
  param_jags <- if (model == "RE") {
    c("EM", "EM.pred", "tau")
  } else {
    c("EM")
  }

  param_jags <- if (is.element(measure, c("RR", "RD"))) {
    append(param_jags, "EM.LOR", "EM.LOR.pred")
  } else {
    param_jags
  }

  # Calculate time needed for all models
  for (i in seq_len(length(mean_misspar[, 1]))) {
    data_jag[[i]] <- list("m" = item$m,
                          "N" = item$N,
                          "t" = item$t,
                          "na" = item$na,
                          "nt" = item$nt,
                          "ns" = item$ns,
                          "ref" = ref,
                          "indic" = indic,
                          "I" = item$I,
                          "meand.phi" = mean_misspar[i, ],
                          "precd.phi" = 1 / var_misspar,
                          "D" = D,
                          "cov_value" = 0,
                          "beta.n" = rep(0, item$nt),
                          "beta" = rep(0, item$nt),
                          "wgt.value" = rep(1, item$ns))

    if (is.element(measure, c("MD", "SMD", "ROM"))) {
      data_jag[[i]] <- append(data_jag[[i]],
                              list("y.o" = item$y0,
                                   "se.o" = item$se0))
    } else if (is.element(measure, c("OR", "RR", "RD"))) {
      data_jag[[i]] <- append(data_jag[[i]], list("r" = item$r,
                                                  "ref_base" = ref_base))
    }

    data_jag[[i]] <- if (model == "RE") {
      append(data_jag[[i]], list("heter.prior" = heterog_prior))
    } else {
      data_jag[[i]]
    }

    message(paste(i, "out of", length(mean_misspar[, 1]), "total scenarios"))
    suppressWarnings({
    jagsfit0[[i]] <- jags(data = data_jag[[i]],
                          inits = inits,
                          parameters.to.save = param_jags,
                          model.file =
                            textConnection(prepare_model(measure,
                                                         model,
                                                         covar_assumption ="NO",
                                                         assumption,
                                                         trans_wgt = "no")),
                          n.chains = n_chains,
                          n.iter = n_iter,
                          n.burnin = n_burnin,
                          n.thin = n_thin)

    # Update until convergence is necessary
    message(paste("Updating model for scenario", i, "until convergence"))
    jagsfit[[i]] <- autojags(jagsfit0[[i]],
                             n.iter = n_iter,
                             n.thin = n_thin,
                             n.update = 2)
    })
  }

  # Obtain the posterior distribution of the necessary model paramters
  get_results <- list()
  for (i in seq_len(length(mean_misspar[, 1]))) {
    get_results[[i]] <- as.data.frame(t(jagsfit[[i]]$BUGSoutput$summary))
  }


  EM <- do.call(rbind,
                lapply(
                  seq_len(length(mean_misspar[, 1])),
                  function(i) t(get_results[[i]])[
                    startsWith(rownames(t(get_results[[i]])), "EM["), ]))
  EM_pred <- do.call(rbind,
                     lapply(seq_len(length(mean_misspar[, 1])),
                            function(i) t(get_results[[i]])[
                              startsWith(rownames(t(get_results[[i]])),
                                         "EM.pred["), ]))

  if (is.element(measure, c("RR", "RD"))) {
    EM_LOR <- do.call(rbind,
                      lapply(
                        seq_len(length(mean_misspar[, 1])),
                        function(i) t(get_results[[i]])[
                          startsWith(rownames(t(get_results[[i]])),
                                     "EM.LOR["), ]))
    EM_LOR_pred <- do.call(rbind,
                           lapply(seq_len(length(mean_misspar[, 1])),
                                  function(i) t(get_results[[i]])[
                                    startsWith(rownames(t(get_results[[i]])),
                                               "EM.LOR.pred["), ]))
  }

  if (model == "RE") {
    tau <- do.call(rbind,
                   lapply(
                     seq_len(length(mean_misspar[, 1])),
                     function(i) t(get_results[[i]])[
                       startsWith(rownames(t(get_results[[i]])), "tau"), ]))
  }

  # Return results
  results <- if (model == "RE" & !is.element(measure, c("RR", "RD"))) {
    list(EM = EM,
         EM_pred = EM_pred,
         tau = tau,
         measure = measure,
         model = model,
         scenarios = mean_scenarios,
         D = D,
         heter = heterog_prior,
         n_chains = n_chains,
         n_iter = n_iter,
         n_burnin = n_burnin,
         n_thin = n_thin)
  } else if (model == "FE" & !is.element(measure, c("RR", "RD"))) {
    list(EM = EM,
         measure = measure,
         model = model,
         scenarios = mean_scenarios,
         D = D,
         n_chains = n_chains,
         n_iter = n_iter,
         n_burnin = n_burnin,
         n_thin = n_thin)
  } else if (model == "RE" & is.element(measure, c("RR", "RD"))) {
    list(EM = EM,
         EM_LOR = EM_LOR,
         EM_pred = EM_pred,
         EM_LOR_pred = EM_LOR_pred,
         tau = tau,
         measure = measure,
         model = model,
         scenarios = mean_scenarios,
         D = D,
         heter = heterog_prior,
         n_chains = n_chains,
         n_iter = n_iter,
         n_burnin = n_burnin,
         n_thin = n_thin)
  } else if (model == "FE" & is.element(measure, c("RR", "RD"))) {
    list(EM = EM,
         EM_LOR = EM_LOR,
         measure = measure,
         model = model,
         scenarios = mean_scenarios,
         D = D,
         n_chains = n_chains,
         n_iter = n_iter,
         n_burnin = n_burnin,
         n_thin = n_thin)
  }

  class(results) <- "run_sensitivity"

  return(results)
}

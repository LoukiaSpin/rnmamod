#' Perform the unrelated mean effects model
#'
#' @description Performs the unrelated mean effects model of Dias et al. (2013)
#'   that has been refined (Spineli, 2021) and extended to address aggregate
#'   binary and continuous missing participant outcome data via the
#'   pattern-mixture model (Spineli et al. 2021; Spineli, 2019). This model
#'   offers a global evaluation of the plausibility of the consistency
#'   assumption in the network.
#'
#' @param full An object of S3 class \code{\link{run_model}}. See 'Value' in
#'   \code{\link{run_model}}.
#' @param n_chains Positive integer specifying the number of chains for the MCMC
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
#'   \code{measure} defined in \code{\link{run_model}}) for each pairwise
#'   comparison observed in the network.}
#'   \item{dev_o}{The deviance contribution of each trial-arm based on the
#'   observed outcome.}
#'   \item{hat_par}{The fitted outcome at each trial-arm.}
#'   \item{tau}{The between-trial standard deviation (assumed common across the
#'   observed pairwise comparisons) for the whole network, when a random-effects
#'   model has been specified.}
#'   \item{m_tau}{The between-trial standard deviation (assumed common
#'   across the observed pairwise comparisons) for the subset of multi-arm
#'   trials, when a random-effects model has been specified.}
#'
#'   The output also includes the following elements:
#'   \item{leverage_o}{The leverage for the observed outcome at each trial-arm.}
#'   \item{sign_dev_o}{The sign of the difference between observed and fitted
#'   outcome at each trial-arm.}
#'   \item{model_assessment}{A data-frame on the measures of model assessment:
#'   deviance information criterion, number of effective parameters, and total
#'   residual deviance.}
#'   \item{jagsfit}{An object of S3 class \code{\link[R2jags:jags]{jags}} with
#'   the posterior results on all monitored parameters to be used in the
#'   \code{\link{mcmc_diagnostics}} function.}
#'
#'   Furthermore, \code{run_ume} returns a character vector with the pairwise
#'   comparisons observed in the network, \code{obs_comp}, and a character
#'   vector with comparisons between the non-baseline interventions observed in
#'   multi-arm trials only, \code{frail_comp}. Both vectors are used in
#'   \code{\link{ume_plot}} function.
#'
#' @details \code{run_ume} inherits the arguments \code{data},
#'   \code{measure}, \code{model}, \code{assumption}, \code{heter_prior},
#'   \code{mean_misspar}, \code{var_misspar}, and \code{ref} from
#'   \code{\link{run_model}}.
#'   This prevents specifying a different Bayesian model from that considered in
#'   \code{\link{run_model}}.Therefore, the user needs first to apply
#'   \code{\link{run_model}}, and then use \code{run_ume} (see 'Examples').
#'
#'   The \code{run_ume} function also returns the arguments \code{data},
#'   \code{model}, \code{measure}, \code{assumption}, \code{n_chains},
#'   \code{n_iter}, \code{n_burnin}, and \code{n_thin} as specified by the user
#'   to be inherited by other relevant functions of the package.
#'
#'   Initially, \code{run_ume} calls the \code{\link{improved_ume}} function to
#'   identify the \emph{frail comparisons}, that is, comparisons between
#'   non-baseline interventions in multi-arm trials not investigated in any
#'   two-arm or multi-arm trial of the network (Spineli, 2021). The 'original'
#'   model of Dias et al. (2013) omits the frail comparisons from the estimation
#'   process. Consequently, the number of estimated summary effects is less
#'   than those obtained by performing separate pairwise meta-analyses
#'   (see \code{\link{run_series_meta}}).
#'
#'   For a binary outcome, when \code{measure} is "RR" (relative risk) or "RD"
#'   (risk difference) in \code{\link{run_model}}, \code{run_ume} currently
#'   considers the odds ratio as effect measure for being the \strong{base-case}
#'   effect measure in \code{\link{run_model}} for a binary outcome
#'   (see also 'Details' in \code{\link{run_model}}).
#'
#'   \code{run_ume} calls the \code{\link{prepare_ume}} function which contains
#'   the WinBUGS code as written by Dias et al. (2013) for binomial and normal
#'   likelihood to analyse binary and continuous outcome data, respectively.
#'   \code{\link{prepare_ume}} has been extended to incorporate the
#'   pattern-mixture model with informative missingness parameters for binary
#'   and continuous outcome data (see 'Details' in \code{\link{run_model}}).
#'   \code{\link{prepare_ume}} has also been refined to account for the
#'   multi-arm trials by assigning conditional univariate normal distributions
#'   on the underlying trial-specific effect size of comparisons with the
#'   baseline arm of the multi-arm trial (Spineli, 2021).
#'
#'   \code{run_ume} runs Bayesian unrelated mean effects model in \code{JAGS}.
#'   The progress of the simulation appears on the R console.
#'
#'   The output of \code{run_ume} is not end-user-ready. The
#'   \code{\link{ume_plot}} function uses the output of \code{run_ume} as an S3
#'   object and processes it further to provide an end-user-ready output.
#'
#'   \code{run_ume} can be used only for a network of interventions. In the case
#'   of two interventions, the execution of the function will be stopped and an
#'   error message will be printed on the R console.
#'
#' @author {Loukia M. Spineli}
#'
#' @seealso \code{\link[R2jags:jags]{jags}},
#'   \code{\link{prepare_ume}}, \code{\link{run_model}},
#'   \code{\link{run_series_meta}}, \code{\link{ume_plot}}
#'
#' @references
#' Dias S, Welton NJ, Sutton AJ, Caldwell DM, Lu G, Ades AE. Evidence synthesis
#' for decision making 4: inconsistency in networks of evidence based on
#' randomized controlled trials.
#' \emph{Med Decis Making} 2013;\bold{33}(5):641--56.
#' doi: 10.1177/0272989X12455847
#'
#' Gelman A, Rubin DB. Inference from iterative simulation using multiple
#' sequences. \emph{Stat Sci} 1992;\bold{7}(4):457--72.
#' doi: 10.1214/ss/1177011136
#'
#' Spineli LM. A Revised Framework to Evaluate the Consistency Assumption
#' Globally in a Network of Interventions. \emph{Med Decis Making} 2021.
#' doi: 10.1177/0272989X211068005
#'
#' Spineli LM, Kalyvas C, Papadimitropoulou K. Continuous(ly) missing outcome
#' data in network meta-analysis: a one-stage pattern-mixture model approach.
#' \emph{Stat Methods Med Res} 2021;\bold{30}(4):958--75.
#' doi: 10.1177/0962280220983544
#'
#' Spineli LM. An empirical comparison of Bayesian modelling strategies for
#' missing binary outcome data in network meta-analysis.
#' \emph{BMC Med Res Methodol} 2019;\bold{19}(1):86.
#' doi: 10.1186/s12874-019-0731-y
#'
#' @examples
#' data("nma.liu2013")
#'
#' # Read results from 'run_model' (using the default arguments)
#' res <- readRDS(system.file('extdata/res_liu.rds', package = 'rnmamod'))
#'
#' \donttest{
#' # Run random-effects unrelated mean effects model
#' # Note: Ideally, set 'n_iter' to 10000 and 'n_burnin' to 1000
#' run_ume(full = res,
#'         n_chains = 3,
#'         n_iter = 1000,
#'         n_burnin = 100,
#'         n_thin = 1)
#' }
#'
#' @export
run_ume <- function(full, n_iter, n_burnin, n_chains, n_thin) {


  if (full$type != "nma" || is.null(full$type)) {
    stop("'full' must be an object of S3 class 'run_model'.",
         call. = FALSE)
  }

  # Default arguments
  data <- full$data
  measure <- if (is.element(full$measure, c("RR", "RD"))) {
    "OR"
  } else {
    full$measure
  }
  model <- full$model
  assumption <- full$assumption
  heterog_prior <- full$heter_prior
  mean_misspar <- full$mean_misspar
  var_misspar <- full$var_misspar
  ref <- full$ref

  # Prepare the dataset for the R2jags
  item <- data_preparation(data, measure)
  if (item$nt < 3) {
    stop("This function is *not* relevant for a pairwise meta-analysis.",
         call. = FALSE)
  }

  # Default arguments
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

  # Move multi-arm trials at the bottom
  t <- item$t[order(item$na, na.last = TRUE), ]
  m <- item$m[order(item$na, na.last = TRUE), ]
  N <- item$N[order(item$na, na.last = TRUE), ]
  na <- sort(item$na)
  ns <- item$ns

  # A function to extract numbers from a character.
  # Source: http://stla.github.io/stlapblog/posts/Numextract.html
  Numextract <- function(string) {
    unlist(regmatches(string, gregexpr("[[:digit:]]+\\.*[[:digit:]]*", string)))
  }

  # Observed comparisons in the network
  impr_ume <- improved_ume(t, N, ns, na)
  observed_comp0 <- impr_ume$obs_comp
  observed_comp <- matrix(Numextract(observed_comp0[, 1]),
                          nrow = length(observed_comp0[, 1]),
                          ncol = 2,
                          byrow = TRUE)
  t1_obs_com <- as.numeric(as.character(observed_comp[, 1]))
  t2_obs_com <- as.numeric(as.character(observed_comp[, 2]))
  obs_comp <- paste0(t2_obs_com, "vs", t1_obs_com)

  # Keep only comparisons with the baseline intervention
  indic0 <- list()
  for (i in 1:ns) {
    indic0[[i]] <- combn(t(na.omit(t(t[i, ]))), 2)[, 1:(na[i] - 1)]
    }
  indic <- unique(t(do.call(cbind, indic0)))
  t1_indic <- indic[, 1]
  t2_indic <- indic[, 2]
  n_obs <- length(t1_indic)

  # Keep only comparisons with the baseline intervention in multi-arm trials
  ns_multi <- length(na[na > 2])
  if (ns_multi < 1) {
    t1_indic_multi <- 0
    t2_indic_multi <- 0
    connected <- 1
  } else {
    indic_multi0 <- list()
    for (i in (ns - ns_multi + 1):ns) {
      indic_multi0[[i]] <- combn(t(na.omit(t(t[i, ]))), 2)[, 1:(na[i] - 1)]
      }
    indic_multi <- unique(t(do.call(cbind, indic_multi0)))
    t1_indic_multi <- indic_multi[, 1]
    t2_indic_multi <- indic_multi[, 2]
    t_indic_multi <-
      unique(c(t1_indic_multi, t2_indic_multi))[-impr_ume$ref_base]
    t_indic_multi2 <- unique(c(t1_indic_multi, t2_indic_multi))

    # Is the subset of multi-arm trials a connected network?
    multi_network <- pairwise(as.list(t[(ns - ns_multi + 1):ns, ]),
                              event = as.list(N[(ns - ns_multi + 1):ns, ]),
                              n = as.list(N[(ns - ns_multi + 1):ns, ]),
                              data = cbind(t[(ns - ns_multi + 1):ns, ],
                                          N[(ns - ns_multi + 1):ns, ],
                                          N[(ns - ns_multi + 1):ns, ],
                                          N[(ns - ns_multi + 1):ns, ]),
                              studlab = 1:ns_multi)

    connected <-
      netconnection(treat1, treat2, studlab, data = multi_network)$n.subnets

    ## For the case of a disconnected network of multi-arm trials
    if (connected > 1) {
      dist_mat <-
        netconnection(treat1, treat2, studlab, data = multi_network)$D.matrix
      group0 <- apply(dist_mat, 2, function(x) length(which(!is.infinite(x))))
      group <- data.frame("treat" = attributes(group0)$names, "freq" = group0)

      if (length(unique(group$freq)) < 2) {
        find_groups <-
          split(group[, 1],
                sort(rep_len(1:(length(group$freq) / unique(group$freq)),
                             length(group[, 1]))))
        trm2 <- data.frame("t_m1" = as.numeric(group[, 1]),
                           "t_m2" = as.numeric(
                             unlist(
                               lapply(
                                 1:(length(group$freq) / unique(group$freq)),
                                 function(i)
                                   rep(min(as.numeric(find_groups[[i]])),
                                       unique(group$freq))))))
      } else {
        find_groups <- split(group[, 1], rep(seq_len(
          length(unique(group$freq))), unique(group$freq)))
        trm2 <- data.frame("t_m1" = as.numeric(group[, 1]),
                           "t_m2" = as.numeric(
                             unlist(
                               lapply(
                                 seq_len(length(unique(group$freq))),
                                 function(i)
                                   rep(min(as.numeric(find_groups[[i]])),
                                       unique(group$freq)[i])))))
      }

      ref_m <- rep(NA, ns)
      for (i in (ns - ns_multi + 1):ns) {
        ref_m[i] <- trm2[is.element(trm2[, 1], t[i, 1]), 2]
        }
      ref_nbase_multi <- rep(NA, impr_ume$nbase_multi)
      for (i in 1:impr_ume$nbase_multi) {
        ref_nbase_multi[i] <- trm2[is.element(trm2[, 1],
                                              cbind(impr_ume$t1_bn[i],
                                                    impr_ume$t2_bn[i])), 2]
        }
    }
  }

  # Data in list format for R2jags
  data_jag <- list("m" = m,
                   "N" = N,
                   "t" = t,
                   "na" = na,
                   "nt" = item$nt,
                   "ns" = ns,
                   "ref" = ifelse(is.element(assumption,
                                             c("HIE-ARM", "IDE-ARM")),
                                  ref, NA),
                   "I" = item$I[order(item$na, na.last = TRUE), ],
                   "t1" = t1_indic,
                   "t2" = t2_indic,
                   "N.obs" = n_obs)

  if (is.element(measure, c("MD", "SMD", "ROM"))) {
    data_jag <- append(data_jag, list("y.o" =
                                        item$y0[
                                          order(item$na, na.last = TRUE), ],
                                      "se.o" =
                                        item$se0[
                                          order(item$na, na.last = TRUE), ],
                                      "y.m" =
                                      item$y0[
                                        order(item$na, na.last = TRUE), ]))
  } else if (measure == "OR") {
    data_jag <- append(data_jag, list("r" =
                                        item$r[
                                          order(item$na, na.last = TRUE), ],
                                      "r.m" =
                                        item$r[
                                          order(item$na, na.last = TRUE), ]))
  }

  data_jag <- if (max(na) > 2 & !is.null(impr_ume$nbase_multi) & connected > 1) {
    append(data_jag, list("ns.multi" = ns_multi,
                          "t1.bn" = impr_ume$t1_bn,
                          "t2.bn" = impr_ume$t2_bn,
                          "nbase.multi" = impr_ume$nbase_multi,
                          "ref.base" = impr_ume$ref_base,
                          "N.t.m" = length(t_indic_multi),
                          "t.m" = t_indic_multi,
                          "ref.m" = ref_m,
                          "ref.nbase.multi" = ref_nbase_multi,
                          "N.t.m2" = length(t_indic_multi2),
                          "t.m2" = trm2)
           )
    } else if (max(na) > 2 & !is.null(impr_ume$nbase_multi) & connected == 1) {
      append(data_jag, list("ns.multi" = ns_multi,
                            "t1.bn" = impr_ume$t1_bn,
                            "t2.bn" = impr_ume$t2_bn,
                            "nbase.multi" = impr_ume$nbase_multi,
                            "ref.base" = impr_ume$ref_base,
                            "N.t.m" = length(t_indic_multi),
                            "t.m" = t_indic_multi)
             )
    } else if (max(na) < 3 || is.null(impr_ume$nbase_multi)) {
      append(data_jag, list("ns.multi" = ns_multi,
                            "t1.bn" = t1_indic,
                            "t2.bn" = t1_indic,
                            "nbase.multi" = 0,
                            "ref.base" = 1,
                            "N.t.m" = length(2:5),
                            "t.m" = 2:5)
             )
    }

  data_jag <- if (is.element(assumption, "IND-CORR")) {
    append(data_jag, list("M" = ifelse(!is.na(m), mean_misspar, NA),
                          "cov.phi" = 0.5 * var_misspar,
                          "var.phi" = var_misspar))
  } else {
    append(data_jag, list("meand.phi" = mean_misspar,
                          "precd.phi" = 1 / var_misspar))
  }

  data_jag <- if (model == "RE") {
    append(data_jag, list("heter.prior" = heterog_prior))
  } else {
    data_jag
  }

  # Define the nodes to be monitored
  param_jags <- if (model == "RE") {
    c("EM", "dev.o", "resdev.o", "totresdev.o", "tau", "m.tau", "hat.par")
  } else {
    c("EM", "dev.o", "resdev.o", "totresdev.o", "hat.par")
  }

  # Run the Bayesian analysis
  jagsfit <- suppressWarnings({jags(data = data_jag,
                  parameters.to.save = param_jags,
                  model.file = textConnection(prepare_ume(measure,
                                                          model,
                                                          assumption,
                                                          connected)),
                  n.chains = n_chains,
                  n.iter = n_iter,
                  n.burnin = n_burnin,
                  n.thin = n_thin,
                  DIC = FALSE)
  })

  # Turn summary of posterior results (R2jags object) into a data-frame
  # to select model parameters (using 'dplyr')
  get_results <- as.data.frame(t(jagsfit$BUGSoutput$summary))

  # Effect size of all unique pairwise comparisons
  EM <- t(get_results %>% dplyr::select(starts_with("EM[")))

  # Between-trial standard deviation
  tau <- t(get_results %>% dplyr::select(starts_with("tau")))
  # For the subnetwork of multi-arm trials
  m_tau <- t(get_results %>% dplyr::select(starts_with("m.tau")))

  # Trial-arm deviance contribution for observed outcome
  dev_o <- t(get_results %>% dplyr::select(starts_with("dev.o[")))

  # Fitted/predicted outcome
  hat_par <- t(get_results %>% dplyr::select(starts_with("hat.par")))

  # Total residual deviance
  dev <- jagsfit$BUGSoutput$summary["totresdev.o", "mean"]

  # Calculate the deviance at posterior mean of fitted values
  # Turn 'number of observed' into a vector
  # (first column, followed by second column, and so on)
  m_new <- suppressMessages({
    as.vector(na.omit(melt(m)[, 2]))
    })
  N_new <- suppressMessages({
    as.vector(na.omit(melt(N)[, 2]))
    })
  obs <- N_new - m_new

  if (is.element(measure, c("MD", "SMD", "ROM"))) {
    # Turn 'y0', 'se0'into a vector as above
    y0_new <- suppressMessages({
      as.vector(na.omit(melt(item$y0[order(item$na, na.last = TRUE), ])[, 2]))
      })
    se0_new <- suppressMessages({
      as.vector(na.omit(melt(item$se0[order(item$na, na.last = TRUE), ])[, 2]))
      })

    # Deviance at the posterior mean of the fitted mean outcome
    dev_post_o <- (y0_new -
                     as.vector(hat_par[, 1])) *
      (y0_new - as.vector(hat_par[, 1])) * (1 / se0_new^2)

    # Sign of the difference between observed and fitted mean outcome
    sign_dev_o <- sign(y0_new - as.vector(hat_par[, 1]))
  } else {
    # Turn 'r' and number of observed into a vector as above
    r_new <- suppressMessages({
      as.vector(na.omit(melt(item$r[order(item$na, na.last = TRUE), ])[, 2]))
      })

    # Correction for zero events in trial-arm
    r0 <- ifelse(r_new == 0, r_new + 0.01,
                 ifelse(r_new == obs, r_new - 0.01, r_new))

    # Deviance at the posterior mean of the fitted response
    dev_post_o <- 2 *
      (r0 * (log(r0) - log(as.vector(hat_par[, 1]))) +
         (obs - r0) * (log(obs - r0) - log(obs - as.vector(hat_par[, 1]))))

    # Sign of the difference between observed and fitted response
    sign_dev_o <- sign(r0 - as.vector(hat_par[, 1]))
  }

  # Obtain the leverage for observed outcomes
  leverage_o <- as.vector(dev_o[, 1]) - dev_post_o

  # Number of effective parameters
  pD <- dev - sum(dev_post_o)

  # Deviance information criterion
  DIC <- pD + dev

  # Measures of model assessment: DIC, pD, and total residual deviance
  model_assessment <- data.frame(DIC, pD, dev)

  ## Collect the minimum results at common
  results <- if (model == "RE") {
    list(EM = EM,
         dev_o = dev_o,
         hat_par = hat_par,
         leverage_o = leverage_o,
         sign_dev_o = sign_dev_o,
         tau = tau,
         m_tau = m_tau,
         model_assessment = model_assessment,
         obs_comp = obs_comp,
         jagsfit = jagsfit,
         data = data,
         model = model,
         measure = measure,
         assumption = assumption,
         phi = NULL,
         n_chains = n_chains,
         n_iter = n_iter,
         n_burnin = n_burnin,
         n_thin = n_thin,
         type = "ume")
  } else {
    list(EM = EM,
         dev_o = dev_o,
         hat_par = hat_par,
         leverage_o = leverage_o,
         sign_dev_o = sign_dev_o,
         model_assessment = model_assessment,
         obs_comp = obs_comp,
         jagsfit = jagsfit,
         data = data,
         model = model,
         measure = measure,
         assumption = assumption,
         phi = NULL,
         n_chains = n_chains,
         n_iter = n_iter,
         n_burnin = n_burnin,
         n_thin = n_thin,
         type = "ume")
  }

  # Return different list of results according to a condition
  if (is.null(impr_ume$nbase_multi)) {
    return(results)
  } else {
    return(append(results,
                  list(m_tau = m_tau,
                       frail_comp =
                         paste0(impr_ume$t2_bn, "vs", impr_ume$t1_bn))))
  }
}

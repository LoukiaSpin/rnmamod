#' Perform Bayesian pairwise or network meta-regression
#'
#' @description
#'   Performs a one-stage pairwise or network meta-regression while addressing
#'   aggregate binary or continuous missing participant outcome data via the
#'   pattern-mixture model.
#'
#' @param full An object of S3 class \code{\link{run_model}}.
#' See 'Value' in \code{\link{run_model}}.
#' @param covariate A numeric vector or matrix for a trial-specific covariate
#'   that is a potential effect modifier. See 'Details'.
#' @param covar_assumption Character string indicating the structure of the
#'   intervention-by-covariate interaction, as described in
#'   Cooper et al. (2009). Set \code{covar_assumption} equal to
#'   \code{"exchangeable"}, \code{"independent"}, or \code{"common"}.
#' @param cov_value Numeric for the covariate value of interest.
#' @param n_chains Positive integer specifying the number of chains for the
#'   MCMC sampling; an argument of the \code{\link[R2jags:jags]{jags}} function
#'   of the R-package \href{https://CRAN.R-project.org/package=R2jags}{R2jags}.
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
#' @param n_thin Positive integer specifying the thinning rate for the
#'   MCMC sampling; an argument of the \code{\link[R2jags]{jags}} function
#'   of the R-package \href{https://CRAN.R-project.org/package=R2jags}{R2jags}.
#'   The default argument is 1.
#' @param inits A list with the initial values for the parameters; an argument
#'   of the \code{\link[R2jags:jags]{jags}} function of the R-package
#'   \href{https://CRAN.R-project.org/package=R2jags}{R2jags}.
#'   The default argument is \code{NULL}, and JAGS generates the initial values.
#'
#' @return A list of R2jags outputs on the summaries of the posterior
#'   distribution, and the Gelman-Rubin convergence diagnostic
#'   (Gelman et al., 1992) for the following monitored parameters for a
#'   fixed-effect pairwise meta-analysis:
#'   \item{EM}{The estimated summary effect measure (according to the argument
#'   \code{measure} defined in \code{\link{run_model}}).}
#'   \item{beta_all}{The estimated regression coefficient for all possible
#'   pairwise comparisons according to the argument \code{covar_assumption}.}
#'   \item{dev_o}{The deviance contribution of each trial-arm based on the
#'   observed outcome.}
#'   \item{hat_par}{The fitted outcome at each trial-arm.}
#'   \item{phi}{The informative missingness parameter.}
#'
#'   For a fixed-effect network meta-analysis, the output additionally
#'   includes:
#'   \item{SUCRA}{The surface under the cumulative ranking (SUCRA) curve for
#'   each intervention.}
#'   \item{effectiveneness}{The ranking probability of each intervention for
#'   every rank.}
#'
#'   For a random-effects pairwise meta-analysis, the output additionally
#'   includes the following elements:
#'   \item{EM_pred}{The predicted summary effect measure (according to the
#'   argument \code{measure} defined in \code{\link{run_model}}).}
#'   \item{delta}{The estimated trial-specific effect measure (according to the
#'   argument \code{measure} defined in \code{\link{run_model}}).
#'   For a multi-arm trial, we estimate \emph{T-1} effects, where \emph{T}
#'   is the number of interventions in the trial.}
#'   \item{tau}{The between-trial standard deviation.}
#'
#'   In network meta-analysis, \code{EM} and \code{EM_pred} refer to all
#'   possible pairwise comparisons of interventions in the network. Furthermore,
#'   \code{tau} is typically assumed to be common for all observed comparisons
#'   in the network.
#'   For a multi-arm trial, we estimate a total \emph{T-1} of \code{delta} for
#'   comparisons with the baseline intervention of the trial (found in the first
#'   column of the element \bold{t}), with \emph{T} being the number of
#'   interventions in the trial.
#'
#'   Furthermore, the output includes the following elements:
#'   \item{abs_risk}{The adjusted absolute risks for each intervention. This
#'   appears only when \code{measure = "OR"}, \code{measure = "RR"}, or
#'   \code{measure = "RD"}.}
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
#'   The \code{run_metareg} function also returns the arguments \code{data},
#'   \code{measure}, \code{model}, \code{assumption}, \code{covariate},
#'   \code{covar_assumption}, \code{n_chains}, \code{n_iter}, \code{n_burnin},
#'   and \code{n_thin} to be inherited by other relevant functions of the
#'   package.
#'
#' @details \code{run_metareg} inherits the arguments \code{data},
#'   \code{measure}, \code{model}, \code{assumption}, \code{heter_prior},
#'   \code{mean_misspar}, \code{var_misspar}, \code{D}, \code{ref},
#'   \code{indic}, and \code{base_risk} from \code{\link{run_model}}
#'   (now contained in the argument \code{full}). This prevents specifying a
#'   different Bayesian model from that considered in \code{\link{run_model}}.
#'   Therefore, the user needs first to apply \code{\link{run_model}}, and then
#'   use \code{run_metareg} (see 'Examples').
#'
#'   The model runs in \code{JAGS} and the progress of the simulation appears on
#'   the R console. The output of \code{run_metareg} is used as an S3 object by
#'   other functions of the package to be processed further and provide an
#'   end-user-ready output. The model is updated until convergence using the
#'   \code{\link[R2jags:autojags]{autojags}} function of the R-package
#'   \href{https://CRAN.R-project.org/package=R2jags}{R2jags} with 2 updates and
#'   number of iterations and thinning equal to \code{n_iter} and \code{n_thin},
#'   respectively.
#'
#'   The models described in Spineli et al. (2021), and Spineli (2019) have
#'   been extended to incorporate one \emph{study-level covariate} variable
#'   following the assumptions of Cooper et al. (2009) for the structure of the
#'   intervention-by-covariate interaction. The covariate can be either a
#'   numeric vector or matrix with columns equal to the maximum number of arms
#'   in the dataset.
#'
#' @author {Loukia M. Spineli}
#'
#' @seealso \code{\link[R2jags:autojags]{autojags}}, \code{\link[R2jags]{jags}},
#'   \code{\link{run_model}}
#'
#' @references
#' Cooper NJ, Sutton AJ, Morris D, Ades AE, Welton NJ. Addressing between-study
#' heterogeneity and inconsistency in mixed treatment comparisons: Application
#' to stroke prevention treatments in individuals with non-rheumatic atrial
#' fibrillation. \emph{Stat Med} 2009;\bold{28}(14):1861--81.
#' doi: 10.1002/sim.3594
#'
#' Gelman A, Rubin DB. Inference from iterative simulation using multiple
#' sequences. \emph{Stat Sci} 1992;\bold{7}(4):457--72.
#' doi: 10.1214/ss/1177011136
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
#' data("nma.baker2009")
#'
#' # Read results from 'run_model' (using the default arguments)
#' res <- readRDS(system.file('extdata/res_baker.rds', package = 'rnmamod'))
#'
#' # Publication year
#' pub_year <- c(1996, 1998, 1999, 2000, 2000, 2001, rep(2002, 5), 2003, 2003,
#'               rep(2005, 4), 2006, 2006, 2007, 2007)
#'
#' \donttest{
#' # Perform a random-effects network meta-regression (exchangeable structure)
#' # Note: Ideally, set 'n_iter' to 10000 and 'n_burnin' to 1000
#' run_metareg(full = res,
#'             covariate = pub_year,
#'             covar_assumption = "exchangeable",
#'             cov_value = 2007,
#'             n_chains = 3,
#'             n_iter = 1000,
#'             n_burnin = 100,
#'             n_thin = 1)
#' }
#'
#' @export
run_metareg <- function(full,
                        covariate,
                        covar_assumption,
                        cov_value,
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
  assumption <- full$assumption
  heterog_prior <- full$heter_prior
  mean_misspar <- full$mean_misspar
  var_misspar <- full$var_misspar
  D <- full$D
  ref <- full$ref
  indic <- full$indic
  base_risk <- full$base_risk

  # Prepare the dataset for the R2jags
  item <- data_preparation(data, measure)

  # Missing and default arguments
  covariate <- if (missing(covariate)) {
    stop("The argument 'covariate' needs to be defined.", call. = FALSE)
  } else {
    covariate
  }
  covar_assumption <- if (missing(covar_assumption)) {
    stop("The argument 'covar_assumption' needs to be defined.",
         call. = FALSE)
  } else if (!is.element(covar_assumption,
                         c("exchangeable", "independent", "common"))) {
    stop("Insert 'exchangeable', 'independent', or 'common'.",
         call. = FALSE)
  } else {
    covar_assumption
  }
  cov_val <- if (missing(cov_value)) {
    stop("The argument 'cov_value' needs to be defined.", call. = FALSE)
  } else if (length(unique(covariate)) < 3) {
    cov_value
  } else if (length(unique(covariate)) > 2) {
    cov_value - mean(covariate)
  }
  if (length(unique(covariate)) < 3 &
      !is.element(cov_value, covariate)) {
    aa <- "The argument 'cov_value' is out of the value"
    stop(paste(aa, "range of the analysed covariate."), call. = FALSE)
  } else if (length(unique(covariate)) > 2 &
             (cov_value < min(covariate) |
              cov_value > max(covariate))) {
    aa <- "The argument 'cov_value' is out of the value"
    stop(paste(aa, "range of the analysed covariate."), call. = FALSE)
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

  # Data in list format for R2jags
  data_jag <- list("m" = item$m,
                   "N" = item$N,
                   "t" = item$t,
                   "na" = item$na,
                   "nt" = item$nt,
                   "ns" = item$ns,
                   "ref" = ref,
                   "I" = item$I,
                   "indic" = indic,
                   "D" = D,
                   "cov_value" = cov_val,
                   "wgt.value" = rep(1, item$ns))

  data_jag <- if (!is.element(measure, c("OR", "RR", "RD"))) {
    append(data_jag, list("y.o" = item$y0,
                          "se.o" = item$se0))
  } else if (is.element(measure, c("OR", "RR", "RD"))) {
    append(data_jag, list("r" = item$r,
                          "ref_base" = ref_base))
  }

  data_jag <- if (length(unique(covariate)) > 2) {
    append(data_jag, list("cov.vector" = covariate - mean(covariate),
                          "cov.matrix" = matrix(0,
                                                nrow = item$ns,
                                                ncol = max(item$na))))
  } else if (length(unique(covariate)) < 3) {
    append(data_jag, list("cov.vector" = covariate,
                          "cov.matrix" = matrix(0,
                                                nrow = item$ns,
                                                ncol = max(item$na))))
  } else if (!is.vector(covariate)) {
    append(data_jag, list("cov.vector" = rep(0, item$ns),
                          "cov.matrix" = covariate))
  }

  data_jag <- if (is.element(assumption, "IND-CORR")) {
    append(data_jag, list("M" = ifelse(!is.na(item$m), mean_misspar, NA),
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

  param_jags <- c("delta",
                  "EM",
                  "EM.pred",
                  "tau",
                  "beta",
                  "beta.all",
                  "SUCRA",
                  "effectiveness",
                  "dev.o",
                  "totresdev.o",
                  "hat.par")

  param_jags <- if (is.element(assumption,
                               c("HIE-COMMON", "HIE-TRIAL", "HIE-ARM"))) {
    append(param_jags, "mean.phi")
  } else {
    append(param_jags, "phi")
  }

  param_jags <- if (model == "RE") {
    param_jags
  } else {
    param_jags[!is.element(param_jags,
                           c("EM.pred", "tau", "delta"))]
  }

  param_jags <- if (is.element(measure, c("OR", "RR", "RD"))) {
    append(param_jags, "abs_risk")
  } else {
    param_jags
  }

  # Run the Bayesian analysis
  message("Running the model ...")
  suppressWarnings({
  jagsfit0 <- jags(data = data_jag,
                   inits = inits,
                   parameters.to.save = param_jags,
                   model.file = textConnection(prepare_model(measure,
                                                             model,
                                                             covar_assumption,
                                                             assumption)),
                   n.chains = n_chains,
                   n.iter = n_iter,
                   n.burnin = n_burnin,
                   n.thin = n_thin)

  # Update until convergence is necessary
  message("... Updating the model until convergence")
  jagsfit <- autojags(jagsfit0, n.iter = n_iter, n.thin = n_thin, n.update = 2)
  })

  # Turn R2jags object into a data-frame
  get_results <- as.data.frame(t(jagsfit$BUGSoutput$summary))

  # Effect size of all unique pairwise comparisons
  EM <- t(get_results)[startsWith(rownames(t(get_results)), "EM["), ]

  # Predictive effects of all unique pairwise comparisons
  EM_pred <- t(get_results)[startsWith(rownames(t(get_results)), "EM.pred["), ]

  # Unique absolute risks for all interventions (only binary data)
  abs_risk <-
    t(get_results)[startsWith(rownames(t(get_results)), "abs_risk["), ]

  # Between-trial standard deviation
  tau <- t(get_results)[startsWith(rownames(t(get_results)), "tau"), ]

  # Regression coefficient for all comparisons with the reference intervention
  beta <-
    t(get_results)[startsWith(rownames(t(get_results)), "beta[") |
                     startsWith(rownames(t(get_results)), "beta"), ]

  # Regression coefficient for all unique pairwise comparisons
  # (not applicable for 'common' covar_assumption)
  beta_all <-
    t(get_results)[startsWith(rownames(t(get_results)), "beta.all["), ]

  # SUrface under the Cumulative RAnking curve values
  SUCRA <- t(get_results)[startsWith(rownames(t(get_results)), "SUCRA"), ]

  # Within-trial effects size
  delta <- t(get_results)[startsWith(rownames(t(get_results)), "delta") &
                            !endsWith(rownames(t(get_results)), ",1]"), ]

  # Ranking probability of each intervention for every rank
  effectiveness <-
    t(get_results)[startsWith(rownames(t(get_results)), "effectiveness"), ]

  # Estimated missingness parameter
  phi <- if (min(na.omit(unlist(item$I))) == 1 &
             max(na.omit(unlist(item$I))) == 1) {
    t(get_results)[startsWith(rownames(t(get_results)), "phi") |
                     startsWith(rownames(t(get_results)), "mean.phi") |
                     startsWith(rownames(t(get_results)), "mean.phi[") |
                     startsWith(rownames(t(get_results)), "phi["), ]
  } else if (min(na.omit(unlist(item$I))) == 0 &
             max(na.omit(unlist(item$I))) == 1) {
    t(get_results)[startsWith(rownames(t(get_results)), "phi["), ]
  } else if (min(na.omit(unlist(item$I))) == 0 &
             max(na.omit(unlist(item$I))) == 0) {
    NULL
  }

  # Trial-arm deviance contribution for observed outcome
  dev_o <- t(get_results)[startsWith(rownames(t(get_results)), "dev.o"), ]

  # Fitted/predicted outcome
  hat_par <- t(get_results)[startsWith(rownames(t(get_results)), "hat.par"), ]

  # Total residual deviance
  dev <- jagsfit$BUGSoutput$summary["totresdev.o", "mean"]

  # Calculate the deviance at posterior mean of fitted values
  # Turn 'N' and 'm' into a vector (first column, followed by second, etc)
  m_new <- suppressMessages({
    as.vector(na.omit(melt(item$m)[, 2]))
    })
  N_new <- suppressMessages({
    as.vector(na.omit(melt(item$N)[, 2]))
    })
  obs <- N_new - m_new

  if (is.element(measure, c("MD", "SMD", "ROM"))) {
    # Turn 'y0', 'se0' into a vector as above
    y0_new <- suppressMessages({
      as.vector(na.omit(melt(item$y0)[, 2]))
      })
    se0_new <- suppressMessages({
      as.vector(na.omit(melt(item$se0)[, 2]))
      })

    # Deviance contribution at the posterior mean of the fitted mean outcome
    dev_post_o <- (y0_new -
                     as.vector(hat_par[, 1])) *
      (y0_new - as.vector(hat_par[, 1])) * (1 / se0_new^2)

    # Sign of the difference between observed and fitted mean outcome
    sign_dev_o <- sign(y0_new - as.vector(hat_par[, 1]))
  } else {
    # Turn 'r' and number of observed into a vector as above
    r_new <- suppressMessages({
      as.vector(na.omit(melt(item$r)[, 2]))
      })

    # Correction for zero events in trial-arm
    r0 <- ifelse(r_new == 0, r_new + 0.01,
                 ifelse(r_new == obs, r_new - 0.01, r_new))

    # Deviance contribution at the posterior mean of the fitted response
    dev_post_o <- 2 * (r0 * (log(r0) -
                               log(as.vector(hat_par[, 1]))) +
                         (obs - r0) * (log(obs - r0) -
                                         log(obs - as.vector(hat_par[, 1]))))

    # Sign of the difference between observed and fitted response
    sign_dev_o <- sign(r0 - as.vector(hat_par[, 1]))
  }

  # Obtain the leverage for observed outcomes
  leverage_o <- as.vector(dev_o[, 1]) - dev_post_o

  # Number of effective parameters
  pD <- dev - sum(dev_post_o)

  # Deviance information criterion
  DIC <- pD + dev

  # A data-frame on the measures of model assessment:
  # DIC, pD, and total residual deviance
  model_assessment <- data.frame(DIC, pD, dev)

  # Return a list of results
  if (model == "RE") {
    ma_results <- list(EM = EM,
                       EM_pred = EM_pred,
                       tau = tau,
                       delta = delta,
                       beta_all = beta_all,
                       dev_o = dev_o,
                       hat_par = hat_par,
                       leverage_o = leverage_o,
                       sign_dev_o = sign_dev_o,
                       phi = phi,
                       model_assessment = model_assessment,
                       measure = measure,
                       model = model,
                       assumption = assumption,
                       covariate = covariate,
                       covar_assumption = covar_assumption,
                       cov_value = cov_value,
                       jagsfit = jagsfit,
                       data = data,
                       n_chains = n_chains,
                       n_iter = n_iter,
                       n_burnin = n_burnin,
                       n_thin = n_thin,
                       abs_risk = abs_risk)
    nma_results <- append(ma_results, list(SUCRA = SUCRA,
                                           effectiveness = effectiveness,
                                           D = D))
  } else if (model == "FE") {
    ma_results <- list(EM = EM,
                       beta_all = beta_all,
                       dev_o = dev_o,
                       hat_par = hat_par,
                       leverage_o = leverage_o,
                       sign_dev_o = sign_dev_o,
                       phi = phi,
                       model_assessment = model_assessment,
                       measure = measure,
                       model = model,
                       assumption = assumption,
                       covariate = covariate,
                       covar_assumption = covar_assumption,
                       cov_value = cov_value,
                       jagsfit = jagsfit,
                       data = data,
                       n_chains = n_chains,
                       n_iter = n_iter,
                       n_burnin = n_burnin,
                       n_thin = n_thin,
                       abs_risk = abs_risk)
    nma_results <- append(ma_results, list(SUCRA = SUCRA,
                                           effectiveness = effectiveness,
                                           D = D))
  }

  class(ma_results) <- class(nma_results) <- "run_metareg"

  ifelse(item$nt > 2, return(nma_results), return(ma_results))
}

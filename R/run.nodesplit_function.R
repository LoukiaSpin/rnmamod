#' Perform the node-splitting approach
#'
#' @description
#'   Performs the Bayesian node-splitting approach of Dias et al. (2010)
#'   extended to address aggregate binary and continuous missing participant
#'   outcome data via the pattern-mixture model (Spineli et al., 2021;
#'   Spineli, 2019). This model offers a local evaluation of the
#'   plausibility of the consistency assumption in the network
#'   (Dias et al., 2010).
#'
#' @param full An object of S3 class \code{\link{run_model}}.
#'   See 'Value' in \code{\link{run_model}}.
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
#'   MCMC sampling; an argument of the \code{\link[R2jags:jags]{jags}} function
#'   of the R-package \href{https://CRAN.R-project.org/package=R2jags}{R2jags}.
#'   The default argument is 1.
#'
#' @return An R2jags output on the summaries of the posterior distribution,
#'   and the Gelman-Rubin convergence diagnostic of the
#'   following monitored parameters:
#'   \item{direct}{The summary effect measure (according to the argument
#'   \code{measure} defined in \code{\link{run_model}}) of each split node based
#'   on the corresponding trials.}
#'   \item{indirect}{The indirect summary effect measure (according to the
#'   argument \code{measure} defined in \code{\link{run_model}}) of each split
#'   node based on the remaining network after removing (splitting) the
#'   corresponding node.}
#'   \item{diff}{The inconsistency parameter for each split node defined as the
#'   difference between the direct and indirect effect of the corresponding
#'   split node.}
#'   \item{tau}{The between-trial standard deviation after each split node, when
#'   the random-effects model has been specified.}
#'
#'   Furthermore, the output includes the following element:
#'   \item{model_assessment}{A data-frame on the measures of model assessment
#'   after each split node: deviance information criterion,
#'   total residual deviance, and number of effective parameters.}
#'
#' @details \code{run_nodesplit} inherits the arguments \code{data},
#'   \code{measure}, \code{model}, \code{assumption}, \code{heter_prior},
#'   \code{mean_misspar}, \code{var_misspar}, \code{ref}, and \code{indic} from
#'   \code{\link{run_model}} (now contained in the argument \code{full}).
#'   This prevents specifying a different Bayesian model from that considered
#'   in \code{\link{run_model}}. Therefore, the user needs first to apply
#'   \code{\link{run_model}}, and then use \code{run_nodesplit}
#'   (see 'Examples').
#'
#'   For a binary outcome, when \code{measure} is "RR" (relative risk) or "RD"
#'   (risk difference) in \code{\link{run_model}}, \code{run_nodesplit}
#'   currently performs node-splitting using the odds ratio as effect measure
#'   for being the \strong{base-case} effect measure in \code{\link{run_model}}
#'   for a binary outcome (see also 'Details' in \code{\link{run_model}}).
#'
#'   To perform the Bayesian node-splitting approach, the
#'   \code{\link{prepare_nodesplit}} function is called which contains the
#'   WinBUGS code as written by Dias et al. (2010) for binomial and normal
#'   likelihood to analyse binary and continuous outcome data, respectively.
#'   \code{\link{prepare_nodesplit}} has been extended to incorporate the
#'   pattern-mixture model with informative missingness parameters for binary
#'   and continuous outcome data (see 'Details' in \code{\link{run_model}}).
#'
#'   \code{run_nodesplit} runs the Bayesian node-splitting approach in
#'   \code{JAGS}. The progress of the simulation appears on the R console. The
#'   number of times \code{run_nodesplit} is used appears on the R console as a
#'   text in red and it equals the number of split nodes (see 'Examples').
#'   If there are no split nodes in the network, the execution of the function
#'   will be stopped and an error message will be printed on the R console.
#'
#'   \code{run_nodesplit} uses the
#'   \code{\link[gemtc:mtc.nodesplit.comparisons]{mtc.nodesplit.comparisons}}
#'   function of the R-package
#'   \href{https://CRAN.R-project.org/package=gemtc}{gemtc}
#'   to obtain automatically the nodes to split based on the decision rule of
#'   van Valkenhoef et al. (2016).
#'   \code{run_nodesplit} uses the option (1) in van Valkenhoef et al. (2016)
#'   to parameterise multi-arm trials that contain the node-to-split.
#'   In contrast,
#'   \code{\link[gemtc:mtc.nodesplit.comparisons]{mtc.nodesplit.comparisons}}
#'   uses the option (3) in van Valkenhoef et al. (2016).
#'   Option (1) keeps the baseline arm of the node-to-split in the corresponding
#'   multi-arms. Option (3) excludes both arms of the node-to-split from the
#'   corresponding multi-arm trials.
#'
#'   The output of \code{run_nodesplit} is not end-user-ready.
#'   The \code{\link{nodesplit_plot}} function inherits the output of
#'   \code{run_nodesplit} as an S3 object and processes it further to provide an
#'   end-user-ready output.
#'
#'   \code{run_nodesplit} can be used only for a network of interventions.
#'   In the case of two interventions, the execution of the function will
#'   be stopped and an error message will be printed on the R console.
#'
#' @author {Loukia M. Spineli}
#'
#' @seealso \code{\link[R2jags:jags]{jags}},
#'   \code{\link[gemtc:mtc.nodesplit.comparisons]{mtc.nodesplit.comparisons}},
#'   \code{\link{nodesplit_plot}}, \code{\link{prepare_nodesplit}},
#'   \code{\link{run_model}}
#'
#' @references
#' Dias S, Welton NJ, Caldwell DM, Ades AE. Checking consistency in mixed
#' treatment comparison meta-analysis.
#' \emph{Stat Med} 2010;\bold{29}(7-8):932--44.
#' doi: 10.1002/sim.3767
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
#' van Valkenhoef G, Dias S, Ades AE, Welton NJ. Automated generation of
#' node-splitting models for assessment of inconsistency in network
#' meta-analysis. \emph{Res Synth Methods} 2016;\bold{7}(1):80--93.
#' doi: 10.1002/jrsm.1167
#'
#' @examples
#' data("nma.baker2009")
#'
#' # Read results from 'run_model' (using the default arguments)
#' res <- readRDS(system.file('extdata/res_baker.rds', package = 'rnmamod'))
#'
#' \donttest{
#' # Run random-effects node-splitting approach
#' # Note: Ideally, set 'n_iter' to 10000 and 'n_burnin' to 1000
#' run_nodesplit(full = res,
#'               n_chains = 3,
#'               n_iter = 1000,
#'               n_burnin = 100,
#'               n_thin = 1)
#' }
#'
#' @export
run_nodesplit <- function(full,
                          n_chains,
                          n_iter,
                          n_burnin,
                          n_thin) {

  if (full$type != "nma" || is.null(full$type)) {
    stop("'full' must be an object of S3 class 'run_meta'.",
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
  indic <- full$indic
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

  # Prepare the dataset for the R2jags
  item <- data_preparation(data, measure)
  if (item$nt < 3) {
    stop("This function is *not* relevant for a pairwise meta-analysis.",
         call. = FALSE)
  }

  if (is.element(measure, c("MD", "SMD", "ROM"))) {
    # Rename columns to agree with gemtc
    names(item$y0) <- paste0("y..", seq_len(max(item$na)), ".")
    names(item$se0) <- paste0("se..", seq_len(max(item$na)), ".")
    names(item$N) <- paste0("n..", seq_len(max(item$na)), ".")
    names(item$t) <- paste0("t..", seq_len(max(item$na)), ".")
    na.. <- item$na

    # Convert to one-arm-per-row as required in GeMTC
    transform <- mtc.data.studyrow(cbind(item$t,
                                         item$y0,
                                         item$se0,
                                         item$N,
                                         na..),
                                   armVars = c("treatment" = "t",
                                               "mean" = "y",
                                               "std.error" = "se",
                                               "sampleSize" = "n"),
                                   nArmsVar = "na")
  } else {
    # Rename columns to agree with gemtc
    names(item$r) <- paste0("r..", seq_len(max(item$na)), ".")
    names(item$N) <- paste0("n..", seq_len(max(item$na)), ".")
    names(item$t) <- paste0("t..", seq_len(max(item$na)), ".")
    na.. <- item$na

    # Convert to one-arm-per-row as required in GeMTC
    transform <- mtc.data.studyrow(cbind(item$t,
                                         item$r,
                                         item$N,
                                         na..),
                                   armVars = c("treatment" = "t",
                                               "response" = "r",
                                               "sampleSize" = "n"),
                                   nArmsVar = "na")
  }

  # Detect the nodes to split (GeMTC functions)
  transform$treatment <- as.numeric(transform$treatment)
  splitting <- mtc.nodesplit.comparisons(mtc.network(transform))
  colnames(splitting) <- NULL
  rownames(splitting) <- NULL

  if (dim(splitting)[1] < 1) {
    stop("There is *no* loop to evaluate", call. = FALSE)
  } else {

    # Define node to split: AB=(1,2)
    pair <- if (dim(splitting)[1] == 1) {
      t(apply(as.matrix(splitting, ncol = 2), 2, as.numeric))
    } else if (dim(splitting)[1] > 1) {
      t(apply(apply(as.matrix(splitting, ncol = 2), 2, as.numeric), 1, sort))
    }

    # Parameters to save
    param_jags <- if (model == "RE") {
      c("EM", "direct", "diff", "tau", "totresdev.o", "hat.par")
    } else {
      c("EM", "direct", "diff",  "totresdev.o", "hat.par")
    }

    # Define necessary model components
    jagsfit <- data_jag <- checkPair <- bi <- si <- m <- list()
    checkPair_node <- t_node <- N_node <- m_node <- m
    y_node <- se_node <- r_node <- I_sign <- m

    for (i in seq_len(length(pair[, 1]))) {
      # Calculate split (1 if present node to split) and b (baseline position)
      checkPair[[i]] <- PairXY(as.matrix(item$t), pair[i, ])
      checkPair_node[[i]] <- checkPair[[i]]
      r_node[[i]] <- matrix(nrow = item$ns, ncol = max(na..))
      se_node[[i]] <- y_node[[i]] <- m_node[[i]] <- r_node[[i]]
      N_node[[i]] <- t_node[[i]] <- I_sign[[i]] <- r_node[[i]]

      for (j in seq_len(item$ns)) {
        t_node[[i]][j, 1] <- item$t[j, checkPair[[i]][j, "b"]]
        t_node[[i]][j, 2:max(na..)] <- unlist(item$t[j,
                                                     -checkPair[[i]][j, "b"]])
        N_node[[i]][j, 1] <- item$N[j, checkPair[[i]][j, "b"]]
        N_node[[i]][j, 2:max(na..)] <- unlist(item$N[j,
                                                     -checkPair[[i]][j, "b"]])
        m_node[[i]][j, 1] <- item$m[j, checkPair[[i]][j, "b"]]
        m_node[[i]][j, 2:max(na..)] <- unlist(item$m[j,
                                                     -checkPair[[i]][j, "b"]])

        if (is.element(measure, c("MD", "SMD", "ROM"))) {
          y_node[[i]][j, 1] <- item$y0[j, checkPair[[i]][j, "b"]]
          y_node[[i]][j, 2:max(na..)] <-
            unlist(item$y0[j, -checkPair[[i]][j, "b"]])
          se_node[[i]][j, 1] <- item$se0[j, checkPair[[i]][j, "b"]]
          se_node[[i]][j, 2:max(na..)] <-
            unlist(item$se0[j, -checkPair[[i]][j, "b"]])
        } else {
          r_node[[i]][j, 1] <- item$r[j, checkPair[[i]][j, "b"]]
          r_node[[i]][j, 2:max(na..)] <-
            unlist(item$r[j, -checkPair[[i]][j, "b"]])
        }

        for (k in 2:max(na..)) {
          I_sign[[i]][j, k] <-
            ifelse(t_node[[i]][j, 1] > t_node[[i]][j, k], -1, 1)
        }
      }

      checkPair_node[[i]][, "b"] <-
        ifelse(checkPair[[i]][, "b"] > 1, 1, checkPair[[i]][, "b"])

      # Build vector bi[i] with baseline treatment: t[i, b[i]]
      bi[[i]] <- Basetreat(as.matrix(t_node[[i]]), checkPair_node[[i]][, "b"])

      # Indexes to sweep non-baseline arms only
      m[[i]] <- NonbaseSweep(checkPair_node[[i]], na..)

      # Build matrix si[i,k] with non-baseline treatments: t[i, m[i,k]]
      si[[i]] <- Sweeptreat(as.matrix(t_node[[i]]), m[[i]])

      # Data in list format for R2jags
      data_jag[[i]] <- list("mod" = m_node[[i]], # item$m
                            "N" = N_node[[i]],   # item$N
                            "t" = t_node[[i]],   # item$t
                            "na" = na..,
                            "nt" = item$nt,
                            "ns" = item$ns,
                            "ref" = ref,
                            "indic" = indic,
                            "I" = item$I,
                            "I.sign" = I_sign[[i]],
                            "split" = checkPair_node[[i]][, "split"],
                            "m" = m[[i]],
                            "bi" = bi[[i]],
                            "si" = si[[i]],
                            "pair" = pair[i, ])

      data_jag[[i]]  <- if (is.element(measure, c("MD", "SMD", "ROM"))) {
        append(data_jag[[i]], list("y.o" = y_node[[i]], "se.o" = se_node[[i]]))
      } else if (measure == "OR") {
        append(data_jag[[i]], list("r" = r_node[[i]]))
      }

      data_jag[[i]] <- if (is.element(assumption, "IND-CORR")) {
        append(data_jag[[i]], list("M" = ifelse(!is.na(item$m), mean_misspar,
                                                NA),
                                   "cov.phi" = 0.5 * var_misspar,
                                   "var.phi" = var_misspar))
      } else {
        append(data_jag[[i]], list("meand.phi" = mean_misspar,
                                   "precd.phi" = 1 / var_misspar))
      }

      data_jag[[i]] <- if (model == "RE") {
        append(data_jag[[i]], list("heter.prior" = heterog_prior))
      } else {
        data_jag[[i]]
      }

      # Run the Bayesian analysis
      message(paste(i, "out of", length(pair[, 1]), "split nodes"))
      jagsfit[[i]] <- jags(data = data_jag[[i]],
                           parameters.to.save = param_jags,
                           model.file =
                             textConnection(prepare_nodesplit(measure,
                                                              model,
                                                              assumption)),
                           n.chains = n_chains,
                           n.iter = n_iter,
                           n.burnin = n_burnin,
                           n.thin = n_thin)
    } # Stop loop for 'pair'
  }

  # Indirect effect size for the split node
  EM <- data.frame(pair[, 2],
                   pair[, 1],
                   do.call(rbind,
                           lapply(seq_len(length(pair[, 1])),
                                  function(i)
                                    jagsfit[[i]]$BUGSoutput$summary[
                                      paste0("EM[", pair[i, 2], ",",
                                             pair[i, 1], "]"),
                                      c("mean", "sd", "2.5%", "97.5%", "Rhat",
                                        "n.eff")])))
  colnames(EM) <- c("treat1",
                    "treat2",
                    "mean",
                    "sd",
                    "2.5%",
                    "97.5%",
                    "Rhat",
                    "n.eff")

  # Direct effect size for the split node
  direct <- data.frame(pair[, 2],
                       pair[, 1],
                       do.call(rbind,
                               lapply(seq_len(length(pair[, 1])),
                                      function(i)
                                        jagsfit[[i]]$BUGSoutput$summary[
                                          "direct",
                                          c("mean", "sd", "2.5%", "97.5%",
                                            "Rhat", "n.eff")])))
  colnames(direct) <- c("treat1",
                        "treat2",
                        "mean",
                        "sd",
                        "2.5%",
                        "97.5%",
                        "Rhat",
                        "n.eff")

  # Inconsistency for for split node
  diff <- data.frame(pair[, 2],
                     pair[, 1],
                     do.call(rbind,
                             lapply(seq_len(length(pair[, 1])),
                                    function(i)
                                      jagsfit[[i]]$BUGSoutput$summary[
                                        "diff",
                                        c("mean", "sd", "2.5%", "97.5%", "Rhat",
                                          "n.eff")])))
  colnames(diff) <- c("treat1",
                      "treat2",
                      "mean",
                      "sd",
                      "2.5%",
                      "97.5%",
                      "Rhat",
                      "n.eff")

  # Between-trial variance after node-splitting
  if (model == "RE") {
    tau <- data.frame(pair[, 2],
                      pair[, 1],
                      do.call(rbind,
                              lapply(seq_len(length(pair[, 1])),
                                     function(i)
                                       jagsfit[[i]]$BUGSoutput$summary[
                                         "tau",
                                         c("50%", "sd", "2.5%", "97.5%", "Rhat",
                                           "n.eff")])))
    colnames(tau) <- c("treat1",
                       "treat2",
                       "50%",
                       "sd",
                       "2.5%",
                       "97.5%",
                       "Rhat",
                       "n.eff")
  } else {
    tau <- NA
  }

  obs <- N_new <- m_new <- get_results <- hat_par <- list()
  r0 <- r_new <- se0_new <- y0_new <- dev_post_o <- list()
  dev <- rep(NA, length(pair[, 1]))
  for (i in seq_len(length(pair[, 1]))) {
    get_results[[i]] <- as.data.frame(t(jagsfit[[i]]$BUGSoutput$summary))

    # Total residual deviance
    dev[i] <- t(get_results[[i]] %>% dplyr::select(
      starts_with("totresdev.o")))[, 1]

    # Fitted/predicted number of observed data (hat.par")
    hat_par[[i]] <- t(get_results[[i]] %>% dplyr::select(
      starts_with("hat.par[")))

    # Calculate the deviance contribution at posterior mean of fitted values
    # Turn 'N' and 'm' into a vector (first column, followed by second, etc)
    m_new[[i]] <- suppressMessages({
      as.vector(na.omit(melt(m_node[[i]])[, 3]))
      })
    N_new[[i]] <- suppressMessages({
      as.vector(na.omit(melt(N_node[[i]])[, 3]))
      })
    obs[[i]] <- N_new[[i]] - m_new[[i]]

    if (is.element(measure, c("MD", "SMD", "ROM"))) {
      # Turn 'y0', 'se0'into a vector as above
      y0_new[[i]] <- suppressMessages({
        as.vector(na.omit(melt(y_node[[i]])[, 3]))
        })
      se0_new[[i]] <- suppressMessages({
        as.vector(na.omit(melt(se_node[[i]])[, 3]))
        })
      # Deviance contribution at the posterior mean of the fitted mean outcome
      dev_post_o[[i]] <- (y0_new[[i]] -
                            as.vector(hat_par[[i]][, 1])) *
        (y0_new[[i]] - as.vector(hat_par[[i]][, 1])) * (1 / (se0_new[[i]]^2))
    } else {
      # Turn 'r' and number of observed into a vector as above
      r_new[[i]] <- suppressMessages({
        as.vector(na.omit(melt(r_node[[i]])[, 3]))
        })
      # Correction for zero events in trial-arm
      r0[[i]] <- ifelse(r_new[[i]] == 0, r_new[[i]] + 0.01,
                        ifelse(r_new[[i]] == obs[[i]], r_new[[i]] - 0.01,
                               r_new[[i]]))
      # Deviance contribution at the posterior mean of the fitted mean outcome
      dev_post_o[[i]] <- 2 *
        (r0[[i]] * (log(r0[[i]]) -
                      log(as.vector(hat_par[[i]][, 1]))) +
           (obs[[i]] - r0[[i]]) * (log(obs[[i]] - r0[[i]]) -
                                     log(obs[[i]] -
                                           as.vector(hat_par[[i]][, 1]))))
    }
  }

  # Number of effective parameters
  pD <- as.vector(do.call(rbind,
                          lapply(seq_len(length(pair[, 1])),
                                 function(i) dev[i] - sum(dev_post_o[[i]]))))

  # Deviance information criterion
  DIC <- as.vector(do.call(rbind,
                           lapply(seq_len(length(pair[, 1])),
                                  function(i) pD[i] + dev[i])))

  # A data-frame on the measures of model assessment:
  # DIC, pD, and total residual deviance
  model_assessment <- data.frame(pair[, 2],
                                 pair[, 1],
                                 unlist(dev),
                                 DIC,
                                 pD)
  colnames(model_assessment) <- c("treat1",
                                  "treat2",
                                  "deviance",
                                  "DIC",
                                  "pD")

  results <- if (model == "RE") {
    list(direct = direct,
         indirect = EM,
         diff = diff,
         tau = tau,
         model_assessment = model_assessment,
         n_chains = n_chains,
         n_iter = n_iter,
         n_burnin = n_burnin,
         n_thin = n_thin,
         type = "node")
  } else {
    list(direct = direct,
         indirect = EM,
         diff = diff,
         model_assessment = model_assessment,
         n_chains = n_chains,
         n_iter = n_iter,
         n_burnin = n_burnin,
         n_thin = n_thin,
         type = "node")
  }

  return(results)
}

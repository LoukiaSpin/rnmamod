#' A function to perform Bayesian node-splitting approach
#'
#' @description This function performs a Bayesian network meta-analysis based on the node-splitting approach of Dias et al. (2010) extended to address aggregate binary and continuous participant outcome data via the pattern-mixture model
#'   (Spineli, 2019; Spineli et al., 2021). This model offers a local evaluation of the plausibility of the consistency assumption in the network (Dias et al. (2010)).
#'
#' @param full An object of S3 class \code{\link{run.model}}. See 'Value' in \code{\link{run.model}}.
#' @param n.chains Positive integer specifying the number of chains for the MCMC sampling; an argument of the \code{\link[R2jags]{jags}} function of the R-package \href{https://CRAN.R-project.org/package=R2jags}{R2jags}.
#'   The default argument is 2.
#' @param n.iter Positive integer specifying the number of Markov chains for the MCMC sampling; an argument of the \code{\link[R2jags]{jags}} function of the R-package \href{https://CRAN.R-project.org/package=R2jags}{R2jags}.
#'   The default argument is 10000.
#' @param n.burnin Positive integer specifying the number of iterations to discard at the beginning of the MCMC sampling; an argument of the \code{\link[R2jags]{jags}} function of the R-package \href{https://CRAN.R-project.org/package=R2jags}{R2jags}.
#'   The default argument is 1000.
#' @param n.thin Positive integer specifying the thinning rate for the MCMC sampling; an argument of the \code{\link[R2jags]{jags}} function of the R-package \href{https://CRAN.R-project.org/package=R2jags}{R2jags}.
#'   The default argument is 1.
#'
#' @return An R2jags output on the summaries of the posterior distribution, and the Gelman-Rubin convergence diagnostic of the following monitored parameters:
#' \tabular{ll}{
#'  \code{direct} \tab The (summary) direct effect estimate of each split node based on the corresponding trials.\cr
#'  \tab \cr
#'  \code{indirect} \tab The indirect effect estimate of each split node based on the remaining network after removing (splitting) the corresponding node.\cr
#'  \tab \cr
#'  \code{diff} \tab The inconsistency parameter for each split node defined as the difference between the direct and indirect effect of the corresponding split node.\cr
#'  \tab \cr
#'  \code{tau} \tab The between-trial standard deviation after each split node, when the random-effects model has been specified.\cr
#' }
#'
#' Furthermore, the output includes the following element:
#' \tabular{ll}{
#'  \code{model.assessment} \tab A data-frame on the measures of model assessment after each split node: deviance information criterion, total residual deviance, and number of effective parameters.\cr
#' }
#'
#' @details \code{run.nodesplit} does not contain the arguments \code{data}, \code{measure}, \code{model}, \code{assumption}, \code{heter.prior}, \code{mean.misspar}, and \code{var.misspar} that are found in \code{run.model}.
#'   This is to prevent misspecifying the Bayesian model as it would make the comparison of the consistency model (via \code{run.model}) with the node-splitting approach meaningless.
#'   Instead, these arguments are contained in the argument \code{full} of the function. Therefore, the user needs first to apply \code{run.model}, and then use \code{run.nodesplit} (see, 'Examples').
#'
#'   To perform the Bayesian node-splitting approach, the \code{prepare.nodesplit} function is called which contains the WinBUGS code as written by Dias et al. (2010) for binomial and normal likelihood to analyse binary and continuous outcome data, respectively.
#'   \code{prepare.nodesplit} has been extended to incorporate the pattern-mixture model with informative missingness parameters for binary and continuous outcome data (see, 'Details' in \code{run.model}).
#'
#'   \code{run.nodesplit} runs Bayesian network meta-analysis based on the node-splitting approach in \code{JAGS}. The progress of the simulation appears in the R console. The number of times \code{run.nodesplit} is used
#'   appears in the R console as a text in red and it equals the number of split nodes (see 'Examples'). If there are no split nodes in the network, the execution of the function will be stopped and an error message will be printed in the R console.
#'
#'   The output of \code{run.nodesplit} is not end-user-ready. The \code{nodesplit.plot} function uses the output of \code{run.nodesplit} as an S3 object and processes it further to provide an end-user-ready output.
#'
#'   \code{run.nodesplit} can be used only for a network of interventions. In the case of two interventions, the execution of the function will be stopped and an error message will be printed in the R console.
#'
#' @author {Loukia M. Spineli}
#'
#' @seealso \code{\link{run.model}}, \code{\link{nodesplit.plot}}, \code{\link{prepare.nodesplit}}, \code{\link[R2jags]{jags}}
#'
#' @references
#' Spineli LM, Kalyvas C, Papadimitropoulou K. Continuous(ly) missing outcome data in network meta-analysis: a one-stage pattern-mixture model approach. \emph{Stat Methods Med Res} 2021. [\doi{10.1177/0962280220983544}]
#'
#' Spineli LM. An empirical comparison of Bayesian modelling strategies for missing binary outcome data in network meta-analysis. \emph{BMC Med Res Methodol} 2019;\bold{19}(1):86. [\doi{10.1186/s12874-019-0731-y}]
#'
#' Dias S, Welton NJ, Caldwell DM, Ades AE. Checking consistency in mixed treatment comparison meta-analysis. \emph{Stat Med} 2010;\bold{29}(7-8):932--44. [\doi{10.1002/sim.3767}]
#'
#' Gelman A, Rubin DB. Inference from iterative simulation using multiple sequences. \emph{Stat Sci} 1992;\bold{7}:457--472. [\doi{10.1214/ss/1177011136}]
#'
#' @examples
#' data("nma.baker2009")
#'
#' \dontrun{
#' # Perform a random-effects network meta-analysis
#' res <- run.model(data = nma.baker2009,
#'                  measure = "OR",
#'                  model = "RE",
#'                  assumption = "IDE-ARM",
#'                  heter.prior = list("halfnormal", 0, 1),
#'                  mean.misspar = c(0, 0),
#'                  var.misspar = 1,
#'                  D = 1,
#'                  n.chains = 3,
#'                  n.iter = 10000,
#'                  n.burnin = 1000,
#'                  n.thin = 1)
#'
#' # Run random-effects network meta-analysis with node-splitting approach
#' run.nodesplit(full = res, n.chains = 3, n.iter = 10000, n.burnin = 1000, n.thin = 1)
#' }
#'
#' @export
run.nodesplit <- function(full, n.chains, n.iter, n.burnin, n.thin){


  data <- full$data
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


  ## Default arguments
  n.chains <- ifelse(missing(n.chains), 2, n.chains)
  n.iter <- ifelse(missing(n.iter), 10000, n.iter)
  n.burnin <- ifelse(missing(n.burnin), 1000, n.burnin)
  n.thin <- ifelse(missing(n.thin), 1, n.thin)


  if (is.element(measure, c("MD", "SMD", "ROM"))) {

    ## Rename columns to agree with gemtc
    names(item$y0) <- paste0("y..", 1:length(item$y0[1, ]), ".")
    names(item$se0) <- paste0("se..", 1:length(item$se0[1, ]), ".")
    names(item$N) <- paste0("n..", 1:length(item$N[1, ]), ".")
    names(item$t) <- paste0("t..", 1:length(item$t[1, ]), ".")
    na.. <- item$na

    ## Convert one-study-per-row data to one-arm-per-row as required in GeMTC
    transform <- mtc.data.studyrow(cbind(item$t, item$y0, item$se0, item$N, na..), armVars = c('treatment'= 't', 'mean'='y', 'std.error'='se', 'sampleSize'='n'), nArmsVar='na')

  } else {

    ## Rename columns to agree with gemtc
    names(item$r) <- paste0("r..", 1:length(item$r[1, ]), ".")
    names(item$N) <- paste0("n..", 1:length(item$N[1, ]), ".")
    names(item$t) <- paste0("t..", 1:length(item$t[1, ]), ".")
    na.. <- item$na

    ## Convert one-study-per-row data to one-arm-per-row as required in GeMTC
    transform <- mtc.data.studyrow(cbind(item$t, item$r, item$N, na..), armVars = c('treatment'= 't', 'response'='r', 'sampleSize'='n'), nArmsVar='na')
  }


  ## Detect the nodes to split
  transform$treatment <- as.numeric(transform$treatment)
  splitting <- mtc.nodesplit.comparisons(mtc.network(transform))
  colnames(splitting) <- NULL
  rownames(splitting) <- NULL

  if(dim(splitting)[1] < 1) {

    stop("There is no loop to evaluate", call. = F)

  } else {

    ## Define node to split: AB=(1,2)
    pair <- if (dim(splitting)[1] == 1) {
      t(apply(as.matrix(splitting, ncol = 2), 2, as.numeric))
    } else if (dim(splitting)[1] > 1) {
      t(apply(apply(as.matrix(splitting, ncol = 2), 2, as.numeric), 1, sort))
    }



    ## Parameters to save
    param.jags <- if (model == "RE") {
      c("EM", "direct", "diff", "tau", "totresdev.o", "hat.par")
    } else {
      c("EM", "direct", "diff",  "totresdev.o", "hat.par")
    }


    ## Define necessary model components
    jagsfit <- data.jag <- checkPair <- bi <- si <- m <- list()
    checkPair.new <- t.new <- N.new <- m.new <- y.new <- se.new <- r.new <- I.sign <- list()


    for (i in 1:length(pair[, 1])) {

      ## Calculate split (1 if node to split is present) and b (baseline position)
      checkPair.new[[i]] <- checkPair[[i]] <- PairXY(as.matrix(item$t), pair[i, ])

      r.new[[i]] <- se.new[[i]] <- y.new[[i]] <- m.new[[i]] <- N.new[[i]] <- t.new[[i]] <- item$t
      I.sign[[i]] <- matrix(nrow = item$ns, ncol = max(na..))
      for(j in 1:item$ns){
        t.new[[i]][j, 1] <- item$t[j, checkPair[[i]][j,"b"]]
        t.new[[i]][j, 2:max(na..)] <- unique(item$t[j, -checkPair[[i]][j,"b"]])
        N.new[[i]][j, 1] <- item$N[j, checkPair[[i]][j,"b"]]
        N.new[[i]][j, 2:max(na..)] <- item$N[j, -checkPair[[i]][j,"b"]]
        m.new[[i]][j, 1] <- item$m[j, checkPair[[i]][j,"b"]]
        m.new[[i]][j, 2:max(na..)] <- item$m[j, -checkPair[[i]][j,"b"]]

        if (is.element(measure, c("MD", "SMD", "ROM"))) {
          y.new[[i]][j, 1] <- item$y0[j, checkPair[[i]][j,"b"]]
          y.new[[i]][j, 2:max(na..)] <- item$y0[j, -checkPair[[i]][j,"b"]]
          se.new[[i]][j, 1] <- item$se0[j, checkPair[[i]][j,"b"]]
          se.new[[i]][j, 2:max(na..)] <- item$se0[j, -checkPair[[i]][j,"b"]]
        } else {
          r.new[[i]][j, 1] <- item$r[j, checkPair[[i]][j,"b"]]
          r.new[[i]][j, 2:max(na..)] <- item$r[j, -checkPair[[i]][j,"b"]]
        }

        for(k in 2:max(na..)){
          I.sign[[i]][j, k] <- ifelse(t.new[[i]][j, 1] > t.new[[i]][j, k], -1, 1)
        }
      }

      checkPair.new[[i]][,"b"] <- ifelse(checkPair[[i]][,"b"] > 1, 1, checkPair[[i]][,"b"])

      ## Build vector bi[i] with baseline treatment: t[i, b[i]]
      #bi[[i]] <- Basetreat(as.matrix(item$t), checkPair[[i]][,"b"])
      bi[[i]] <- Basetreat(as.matrix(t.new[[i]]), checkPair.new[[i]][,"b"])

      ## Indexes to sweep non-baseline arms only
      #m[[i]] <- NonbaseSweep(checkPair[[i]], na..)
      m[[i]] <- NonbaseSweep(checkPair.new[[i]], na..)

      ## Build matrix si[i,k] with non-baseline treatments: t[i, m[i,k]]
      #si[[i]] <- Sweeptreat(as.matrix(item$t), m[[i]])
      si[[i]] <- Sweeptreat(as.matrix(t.new[[i]]), m[[i]])


      ## Data in list format for R2jags
      data.jag[[i]] <- list("mod" = m.new[[i]], # item$m
                            "N" = N.new[[i]],   # item$N
                            "t" = t.new[[i]],   # item$t
                            "na" = na..,
                            "nt" = item$nt,
                            "ns" = item$ns,
                            "ref" = item$ref,
                            "I" = item$I,
                            "I.sign" = I.sign[[i]],
                            "M" = ifelse(!is.na(item$m), mean.misspar, NA),
                            "cov.phi" = 0.5*var.misspar,
                            "var.phi" = var.misspar,
                            "meand.phi" = mean.misspar,
                            "precd.phi" = 1/var.misspar,
                            "split" = checkPair.new[[i]][, "split"], # checkPair[[i]][, "split"]
                            "m" = m[[i]],
                            "bi" = bi[[i]],
                            "si" = si[[i]],
                            "pair" = pair[i, ],
                            "heter.prior" = heter.prior)


      if (is.element(measure, c("MD", "SMD", "ROM"))) {
        data.jag[[i]]  <- append(data.jag[[i]] , list("y.o" = y.new[[i]], "se.o" = se.new[[i]])) # list("y.o" = item$y0, "se.o" = item$se0)
      } else if (measure == "OR") {
        data.jag[[i]]  <- append(data.jag[[i]] , list("r" = r.new[[i]])) #list("r" = item$r)
      }


      ## Run the Bayesian analysis
      message(paste(i, "out of", length(pair[, 1]), "split nodes"))
      jagsfit[[i]] <- jags(data = data.jag[[i]],
                           parameters.to.save = param.jags,
                           model.file = textConnection(prepare.nodesplit(measure, model, assumption)),
                           n.chains = n.chains,
                           n.iter = n.iter,
                           n.burnin = n.burnin,
                           n.thin = n.thin)

    } # Stop loop for 'pair'
  }


  ## Obtain the posterior distribution of the necessary model parameters (node: 'treat1' versus 'treat2')
  # Indirect effect size for the split node
  EM <- data.frame(pair[, 2], pair[, 1], do.call(rbind,lapply(1:length(pair[, 1]), function(i) jagsfit[[i]]$BUGSoutput$summary[paste0("EM[", pair[i, 2], ",", pair[i, 1], "]"), c("mean", "sd", "2.5%", "97.5%", "Rhat", "n.eff")])))
  colnames(EM) <- c("treat1", "treat2", "mean", "sd", "2.5%", "97.5%", "Rhat", "n.eff")

  # Direct effect size for the split node
  direct <- data.frame(pair[, 2], pair[, 1], do.call(rbind,lapply(1:length(pair[, 1]), function(i) jagsfit[[i]]$BUGSoutput$summary["direct", c("mean", "sd", "2.5%", "97.5%", "Rhat", "n.eff")])))
  colnames(direct) <- c("treat1", "treat2", "mean", "sd", "2.5%", "97.5%", "Rhat", "n.eff")

  # Inconsistency for for split node
  diff <- data.frame(pair[, 2], pair[, 1], do.call(rbind,lapply(1:length(pair[, 1]), function(i) jagsfit[[i]]$BUGSoutput$summary["diff", c("mean", "sd", "2.5%", "97.5%", "Rhat", "n.eff")])))
  colnames(diff) <- c("treat1", "treat2", "mean", "sd", "2.5%", "97.5%", "Rhat", "n.eff")

  # Between-trial variance after node-splitting
  if (model == "RE") {
    tau <- data.frame(pair[, 2], pair[, 1], do.call(rbind,lapply(1:length(pair[, 1]), function(i) jagsfit[[i]]$BUGSoutput$summary["tau", c("50%", "sd", "2.5%", "97.5%", "Rhat", "n.eff")])))
    colnames(tau) <- c("treat1", "treat2", "50%", "sd", "2.5%", "97.5%", "Rhat", "n.eff")
  } else {
    tau <- NA
  }


  getResults <- hat.par <- list()
  dev <- rep(NA, length(pair[, 1]))
  for (i in 1:length(pair[, 1])) {
    getResults[[i]] <- as.data.frame(t(jagsfit[[i]]$BUGSoutput$summary))

    # Total residual deviance
    dev[i] <- t(getResults[[i]] %>% dplyr::select(starts_with("totresdev.o")))[, 1]

    # Fitted/predicted number of observed data (hat.par")
    hat.par[[i]] <- t(getResults[[i]] %>% dplyr::select(starts_with("hat.par[")))
  }



  ## Calculate the deviance at posterior mean of fitted values
  # Turn 'number of observed' and 'm' into a vector (first column, followed by second column, and so on)
  m.new <- suppressMessages({as.vector(na.omit(melt(item$m)[, 2]))})
  N.new <- suppressMessages({as.vector(na.omit(melt(item$N)[, 2]))})
  obs <- N.new - m.new

  if (is.element(measure, c("MD", "SMD", "ROM"))) {

    # Turn 'y0', 'se0'into a vector (first column, followed by second column, and so on)
    y0.new <- suppressMessages({as.vector(na.omit(melt(item$y0)[, 2]))})
    se0.new <- suppressMessages({as.vector(na.omit(melt(item$se0)[, 2]))})

    dev.post.o <- list()
    for (i in 1:length(pair[, 1])) {
      # Deviance at the posterior mean of the fitted mean outcome
      dev.post.o[[i]] <- (y0.new - as.vector(hat.par[[i]][, 1]))*(y0.new - as.vector(hat.par[[i]][, 1]))*(1/se0.new^2)
    }

  } else {

    # Turn 'r' and number of observed into a vector (first column, followed by second column, and so on)
    r.new <- suppressMessages({as.vector(na.omit(melt(item$r)[, 2]))})

    # Correction for zero events in trial-arm
    r0 <- ifelse(r.new == 0, r.new + 0.01, ifelse(r.new == obs, r.new - 0.01, r.new))

    dev.post.o <- list()
    for (i in 1:length(pair[, 1])) {
      # Deviance at the posterior mean of the fitted mean outcome
      dev.post.o[[i]] <- 2*(r0*(log(r0) - log(as.vector(hat.par[[i]][, 1]))) + (obs - r0)*(log(obs - r0) - log(obs - as.vector(hat.par[[i]][, 1]))))
    }
  }


  # Number of effective parameters
  pD <- do.call(rbind, lapply(1:length(pair[, 1]), function(i) dev[[i]] - sum(dev.post.o[[i]])))


  # Deviance information criterion
  DIC <- do.call(rbind, lapply(1:length(pair[, 1]), function(i) pD[[i]] + dev[[i]]))

  # A data-frame on the measures of model assessment: DIC, pD, and total residual deviance
  model.assessment <- data.frame(pair[, 2], pair[, 1], DIC, unlist(dev), pD)
  colnames(model.assessment) <- c("treat1", "treat2", "DIC", "deviance", "pD")

  results <- if (model == "RE") {
    list(direct = direct,
         indirect = EM,
         diff = diff,
         tau = tau,
         model.assessment = model.assessment)
  } else {
    list(direct = direct,
         indirect = EM,
         diff = diff,
         model.assessment = model.assessment)
  }

  return(results)
}




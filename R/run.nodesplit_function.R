#' A function to perform Bayesian node-splitting approach for aggregate binary or continuous outcomes
#'
#' @export
run.nodesplit <- function(data, measure, model, assumption, heter.prior, mean.misspar, var.misspar, n.chains, n.iter, n.burnin, n.thin){


  ## Turn off warning when variables in the 'data.jag' are not used
  options(warn = -1)


  ## Prepare the dataset for the R2jags
  item <- data.preparation(data, measure)


  ## Default arguments
  model <- if (missing(model)) {
    "RE"
  } else if (!is.element(model, c("RE", "FE"))) {
    stop("Insert 'RE', or 'FE'")
  } else {
    model
  }
  assumption <- if (missing(assumption)) {
    "IDE-ARM"
  } else if (!is.element(assumption,  c("IDE-ARM", "IDE-TRIAL", "IDE-COMMON", "HIE-ARM", "HIE-TRIAL", "HIE-COMMON", "IND-CORR", "IND-UNCORR"))) {
    stop("Insert 'IDE-ARM', 'IDE-TRIAL', 'IDE-COMMON', 'HIE-ARM', 'HIE-TRIAL', 'HIE-COMMON', 'IND-CORR', or 'IND-UNCORR'")
  } else {
    assumption
  }
  mean.misspar <- missingness.param.prior(assumption, mean.misspar)
  heter.prior <- heterogeneity.param.prior(measure, model, heter.prior)
  var.misspar <- ifelse(missing(var.misspar) & (is.element(measure, c("OR", "MD", "SMD"))), 1, ifelse(missing(var.misspar) & measure == "ROM", 0.2^2, var.misspar))
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
    transform <- mtc.data.studyrow(cbind(item$t, item$y0, item$se0, item$N, item$na), armVars = c('treatment'= 't', 'mean'='y', 'std.error'='se', 'sampleSize'='n'), nArmsVar='na')

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

    suppressMessages({
      message("Called from: rnmamod::run.nodesplit")
    })

  } else {

    ## Define node to split: AB=(1,2)
    pair <- t(apply(apply(as.matrix(splitting, ncol = 2), 2, as.numeric), 1, sort))


    ## Parameters to save
    param.jags <- if (model == "RE") {
      c("EM", "direct", "diff", "tau", "totresdev.o", "hat.par")
    } else {
      c("EM", "direct", "diff",  "totresdev.o", "hat.par")
    }


    ## Define necessary model components
    jagsfit <- data.jag <- checkPair <- bi <- si <- m <- list()


    for (i in 1:length(pair[, 1])) {

      ## Calculate split (1 if node to split is present) and b (baseline position)
      checkPair[[i]] <- PairXY(as.matrix(item$t), pair[i, ])

      ## Build vector bi[i] with baseline treatment: t[i, b[i]]
      bi[[i]] <- Basetreat(as.matrix(item$t), checkPair[[i]][,"b"])

      ## Indexes to sweep non-baseline arms only
      m[[i]] <- NonbaseSweep(checkPair[[i]], na..)

      ## Build matrix si[i,k] with non-baseline treatments: t[i, m[i,k]]
      si[[i]] <- Sweeptreat(as.matrix(item$t), m[[i]])


      ## Data in list format for R2jags
      data.jag[[i]] <- list("mod" = item$m,
                            "N" = item$N,
                            "t" = item$t,
                            "na" = na..,
                            "nt" = item$nt,
                            "ns" = item$ns,
                            "ref" = item$ref,
                            "I" = item$I,
                            "M" = ifelse(!is.na(item$m), mean.misspar, NA),
                            "cov.phi" = 0.5*var.misspar,
                            "var.phi" = var.misspar,
                            "meand.phi" = mean.misspar,
                            "precd.phi" = 1/var.misspar,
                            "split" = checkPair[[i]][, "split"],
                            "m" = m[[i]],
                            "bi" = bi[[i]],
                            "si" = si[[i]],
                            "pair" = pair[i, ],
                            "heter.prior" = heter.prior)


      if (is.element(measure, c("MD", "SMD", "ROM"))) {
        data.jag[[i]]  <- append(data.jag[[i]] , list("y.o" = item$y0, "se.o" = item$se0))
      } else if (measure == "OR") {
        data.jag[[i]]  <- append(data.jag[[i]] , list("r" = item$r))
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
         model.assessment = model.assessment,
         measure = measure,
         model = model)
  } else {
    list(direct = direct,
         indirect = EM,
         diff = diff,
         model.assessment = model.assessment,
         measure = measure,
         model = model)
  }

  return(results)
}




#' A function to perform Bayesian node-splitting approach for aggregate binary or continuous outcomes
#'
#' @export
run.nodesplit <- function(data, measure, model, assumption, heter.prior, mean.misspar, var.misspar, n.chains, n.iter, n.burnin, n.thin){


  options(warn = -1)

  ## Default arguments
  measure <- if (missing(measure)) {
    stop("The 'measure' needs to be defined")
  } else if (measure != "MD" & measure != "SMD" & measure != "ROM" & measure != "OR") {
    stop("Insert 'MD', 'SMD', 'ROM', or 'OR'")
  } else {
    measure
  }
  model <- if (missing(model)) {
    "RE"
  } else if (model != "RE" & model != "FE") {
    stop("Insert 'RE', or 'FE'")
  } else {
    model
  }
  assumption <- ifelse(missing(assumption), "IDE-ARM", assumption)
  heter.prior <- heterogeneity.param.prior(measure, model, heter.prior)
  var.misspar <- ifelse(missing(var.misspar) & (measure == "OR" || measure == "MD"|| measure == "SMD"), 1, ifelse(missing(var.misspar) & measure == "ROM", 0.2^2, var.misspar))
  n.chains <- ifelse(missing(n.chains), 2, n.chains)
  n.iter <- ifelse(missing(n.iter), 10000, n.iter)
  n.burnin <- ifelse(missing(n.burnin), 1000, n.burnin)
  n.thin <- ifelse(missing(n.thin), 1, n.thin)


  if(measure == "MD" || measure == "SMD"|| measure == "ROM"){

    ## Continuous: arm-level, wide-format dataset
    y.obs <- data %>% dplyr::select(starts_with("y"))                               # Observed mean value in each arm of every trial
    sd.obs <- data %>% dplyr::select(starts_with("sd"))                             # Observed standard deviation in each arm of every trial
    mod0 <- data %>% dplyr::select(starts_with("m"))                                # Number of missing participants in each arm of every trial
    c <- data %>% dplyr::select(starts_with("c"))                                   # Number of completers in each arm of every trial
    se.obs <- sd.obs/sqrt(c)                                                        # Observed standard error in each arm of every trial
    rand <- mod0 + c                                                                # Number of randomised participants in each arm of every trial
    treat <- data %>% dplyr::select(starts_with("t"))                               # Intervention studied in each arm of every trial
    na.. <- apply(treat, 1, function(x) length(which(!is.na(x))))                   # Number of interventions investigated in every trial per network
    nt <- length(table(as.matrix(treat)))                                           # Total number of interventions per network
    ns <- length(y.obs[, 1])                                                        # Total number of included trials per network
    ref <- 1                                                                        # The first intervention (t1 = 1) is the reference of the network


    ## Order by 'id of t1' < 'id of t1'
    y0 <- se0 <- mod <- N <- t <- t0 <- treat
    for (i in 1:ns) {
      t0[i, ] <- order(treat[i, ], na.last = T)
      y0[i, ] <- y.obs[i, order(t0[i, ], na.last = T)]
      se0[i, ] <- se.obs[i, order(t0[i, ], na.last = T)]
      mod[i, ] <- mod0[i, order(t0[i, ], na.last = T)]
      N[i, ] <- rand[i, order(t0[i, ], na.last = T)]
      t[i, ] <- sort(treat[i, ], na.last = T)
    }



    ## Condition regarding the specification of the prior mean ('mean.misspar') for the missingness parameter
    if(missing(mean.misspar)) {

      mean.misspar <- rep(0, 2)

    } else if(!missing(mean.misspar) & (assumption == "HIE-ARM" || assumption == "IDE-ARM" ) & !is.null(dim(mean.misspar))) {

      mean.misspar <- as.vector(mean.misspar)

    } else if(!missing(mean.misspar) & (assumption == "HIE-ARM" || assumption == "IDE-ARM" ) & is.null(dim(mean.misspar))) {

      mean.misspar <- rep(mean.misspar, 2)

    } else {

      mean.misspar <- mean.misspar

    }


    ## Rename columns to agree with gemtc
    names(y0) <- paste0("y..",1:length(y0[1, ]),".")
    names(se0) <- paste0("se..",1:length(se0[1, ]),".")
    names(N) <- paste0("n..",1:length(N[1, ]),".")
    names(t) <- paste0("t..",1:length(t[1, ]),".")


    ## Convert one-study-per-row data to one-arm-per-row as required in GeMTC
    transform <- mtc.data.studyrow(cbind(t, y0, se0, N, na..), armVars = c('treatment'= 't', 'mean'='y', 'std.error'='se', 'sampleSize'='n'), nArmsVar='na')
    transform$treatment <- as.numeric(transform$treatment)

  } else {


    ## Binary: arm-level, wide-format dataset
    (event <- data %>% dplyr::select(starts_with("r")))                             # Number of observed events in each arm of every trial
    (mod0 <- data %>% dplyr::select(starts_with("m")))                              # Number of missing participants in each arm of every trial
    (rand <- data %>% dplyr::select(starts_with("n")))                              # Number randomised participants in each arm of every trial
    (treat <- data %>% dplyr::select(starts_with("t")))                             # Intervention studied in each arm of every trial
    na.. <- apply(treat, 1, function(x) length(which(!is.na(x))))                   # Number of interventions investigated in every trial per network
    nt <- length(table(as.matrix(treat)))                                           # Total number of interventions per network
    ns <- length(event[, 1])                                                        # Total number of included trials per network
    ref <- 1                                                                        # The first intervention (t1 = 1) is the reference of the network



    ## Order by 'id of t1' < 'id of t1'
    r <- mod <- N <- t <- t0 <- treat
    for (i in 1:ns) {
      t0[i, ] <- order(treat[i, ], na.last = T)
      r[i, ] <- event[i, order(t0[i, ], na.last = T)]
      mod[i, ] <- mod0[i, order(t0[i, ], na.last = T)]
      N[i, ] <- rand[i, order(t0[i, ], na.last = T)]
      t[i, ] <- sort(treat[i, ], na.last = T)
    }


    ## Condition regarding the specification of the prior mean ('mean.misspar') for the missingness parameter
    if(missing(mean.misspar)) {

      mean.misspar <- rep(0.0001, 2)

    } else if(!missing(mean.misspar) & (assumption == "HIE-ARM" || assumption == "IDE-ARM" ) & !is.null(dim(mean.misspar))) {

      mean.misspar <- as.vector(mean.misspar)
      mean.misspar[1] <- ifelse(mean.misspar[1] == 0, 0.0001, mean.misspar[1])
      mean.misspar[2] <- ifelse(mean.misspar[2] == 0, 0.0001, mean.misspar[2])

    } else if(!missing(mean.misspar) & (assumption == "HIE-ARM" || assumption == "IDE-ARM" ) & is.null(dim(mean.misspar))) {

      mean.misspar <- rep(ifelse(mean.misspar == 0, 0.0001, mean.misspar), 2)

    } else if(!missing(mean.misspar) & (assumption != "HIE-ARM" || assumption != "IDE-ARM")) {

      mean.misspar <- ifelse(mean.misspar == 0, 0.0001, mean.misspar)

    }


    ## Rename columns to agree with gemtc
    names(r) <- paste0("r..",1:length(r[1, ]),".")
    names(N) <- paste0("n..",1:length(N[1, ]),".")
    names(t) <- paste0("t..",1:length(t[1, ]),".")


    ## Convert one-study-per-row data to one-arm-per-row as required in GeMTC
    transform <- mtc.data.studyrow(cbind(t, r, N, na..), armVars = c('treatment'= 't', 'response'='r', 'sampleSize'='n'), nArmsVar='na')
    transform$treatment <- as.numeric(transform$treatment)

  }


  M <- ifelse(!is.na(mod), mean.misspar, NA)   # Vector of the mean value of the normal distribution of the informative missingness parameter as the number of arms in trial i (independent structure)
  prec.misspar <- 1/var.misspar
  psi.misspar <- sqrt(var.misspar)           # the lower bound of the uniform prior distribution for the prior standard deviation of the missingness parameter (hierarchical structure)
  cov.misspar <- 0.5*var.misspar             # covariance of pair of missingness parameters in a trial (independent structure)



  ## Detect the nodes to split
  splitting <- mtc.nodesplit.comparisons(mtc.network(transform))
  colnames(splitting) <- NULL
  rownames(splitting) <- NULL

  if(dim(splitting)[1] < 1) {

    stop("There is no loop to evaluate", call. = F)

    suppressMessages({
      message("Called from: netmodr::run.nodesplit")
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


    for(i in 1:length(pair[, 1])){


      ## Calculate split (1 if node to split is present) and b (baseline position)
      checkPair[[i]] <- PairXY(as.matrix(t), pair[i, ])

      ## Build vector bi[i] with baseline treatment: t[i, b[i]]
      bi[[i]] <- Basetreat(as.matrix(t), checkPair[[i]][,"b"])

      ## Indexes to sweep non-baseline arms only
      m[[i]] <- NonbaseSweep(checkPair[[i]], na..)

      ## Build matrix si[i,k] with non-baseline treatments: t[i, m[i,k]]
      si[[i]] <- Sweeptreat(as.matrix(t), m[[i]])


      data.jag[[i]] <- list("mod" = mod,
                            "N" = N,
                            "t" = t,
                            "na" = na..,
                            "nt" = nt,
                            "ns" = ns,
                            "ref" = ref,
                            "meand.phi" = mean.misspar,
                            "precd.phi" = prec.misspar,
                            "split" = checkPair[[i]][, "split"],
                            "m" = m[[i]],
                            "bi" = bi[[i]],
                            "si" = si[[i]],
                            "pair" = pair[i, ],
                            "heter.prior" = heter.prior)


      if (measure == "MD" || measure == "SMD" || measure == "ROM") {
        data.jag[[i]] <- append(data.jag[[i]], list("y.o" = y0, "se.o" = se0))
      } else if (measure == "OR") {
        data.jag[[i]] <- append(data.jag[[i]], list("r" = r))
      }


      ## Run the Bayesian analysis
      message(paste(i, "out of", length(pair[, 1]), "split nodes"))
      jagsfit[[i]] <- jags(data = data.jag[[i]],
                           parameters.to.save = param.jags,
                           model.file = textConnection(prepare.nodesplit(measure, model, assumption)),
                           n.chains = n.chains,
                           n.iter = n.iter,
                           n.burnin = n.burnin,
                           n.thin = n.thin,
                           DIC = T)

    } # Stop loop for 'pair'


  }



  ## Obtain the posterior distribution of the necessary model parameters (node: 'treat1' versus 'treat2'
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
  m.new <- suppressMessages({as.vector(na.omit(melt(mod)[, 2]))})
  N.new <- suppressMessages({as.vector(na.omit(melt(N)[, 2]))})
  obs <- N.new - m.new

  if (measure == "MD" || measure == "SMD"|| measure == "ROM") {

    # Turn 'y0', 'se0'into a vector (first column, followed by second column, and so on)
    y0.new <- suppressMessages({as.vector(na.omit(melt(y0)[, 2]))})
    se0.new <- suppressMessages({as.vector(na.omit(melt(se0)[, 2]))})

    dev.post.o <- list()
    for (i in 1:length(pair[, 1])) {
      # Deviance at the posterior mean of the fitted mean outcome
      dev.post.o[[i]] <- (y0.new - as.vector(hat.par[[i]][, 1]))*(y0.new - as.vector(hat.par[[i]][, 1]))*(1/se0.new^2)
    }

  } else {

    # Turn 'r' and number of observed into a vector (first column, followed by second column, and so on)
    r.new <- suppressMessages({as.vector(na.omit(melt(r)[, 2]))})

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




#' A function to perform Bayesian node-splitting approach for aggregate binary or continuous outcomes
#'
#' @export
run.nodesplit <- function(data, measure, assumption, heter.prior, net.ref, mean.misspar, var.misspar, n.chains, n.iter, n.burnin, n.thin){


  ## Default arguments
  assumption <- ifelse(missing(assumption), "IDE-ARM", assumption)
  var.misspar <- ifelse(missing(var.misspar) & (measure == "OR" || measure == "MD"|| measure == "SMD"), 1, ifelse(missing(var.misspar) & measure == "ROM", 0.2^2, var.misspar))
  n.chains <- ifelse(missing(n.chains), 2, n.chains)
  n.iter <- ifelse(missing(n.iter), 10000, n.iter)
  n.burnin <- ifelse(missing(n.burnin), 1000, n.burnin)
  n.thin <- ifelse(missing(n.thin), 1, n.thin)


  if(measure == "MD" || measure == "SMD"|| measure == "ROM"){

    ## Continuous: arm-level, wide-format dataset
    (y.obs <- data %>% dplyr::select(starts_with("y")))                             # Observed mean value in each arm of every trial
    (sd.obs <- data %>% dplyr::select(starts_with("sd")))                           # Observed standard deviation in each arm of every trial
    (mod0 <- data %>% dplyr::select(starts_with("m")))                              # Number of missing participants in each arm of every trial
    (c <- data %>% dplyr::select(starts_with("c")))                                 # Number of completers in each arm of every trial
    (se.obs <- sd.obs/sqrt(c))                                                      # Observed standard error in each arm of every trial
    (rand <- mod0 + c)                                                              # Number of randomised participants in each arm of every trial
    (treat <- data %>% dplyr::select(starts_with("t")))                             # Intervention studied in each arm of every trial
    na.. <- apply(treat, 1, function(x) length(which(!is.na(x))))                   # Number of interventions investigated in every trial per network
    nt <- length(table(as.matrix(treat)))                                           # Total number of interventions per network
    ns <- length(y.obs[, 1])                                                        # Total number of included trials per network
    ref <- ifelse(missing(net.ref), which.max(table(as.matrix(treat))), net.ref)    # Reference intervention per network. If not specify, the most frequently appeared intervention in the network is selected
    # Trial-specific observed pooled standard deviation
    (sigma <- sqrt(apply((sd.obs^2)*(c - 1), 1, sum, na.rm = T)/(apply(c, 1, sum, na.rm = T) - na..)))


    ## Order by 'id of t1' < 'id of t1'
    y0 <- se0 <- mod <- N <- t <- t0 <- treat
    for(i in 1:ns){

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


    ## Information for the prior distribution on the missingness parameter (IMDOM or logIMROM)
    M <- ifelse(!is.na(y0), mean.misspar, NA)  # Vector of the mean value of the normal distribution of the informative missingness parameter as the number of arms in trial i (independent structure)
    prec.misspar <- 1/var.misspar
    psi.misspar <- sqrt(var.misspar)           # the lower bound of the uniform prior distribution for the prior standard deviation of the missingness parameter (hierarchical structure)
    cov.misspar <- 0.5*var.misspar             # covariance of pair of missingness parameters in a trial (independent structure)



    ## Specification of the prior distribution for the between-trial parameter
    if (heter.prior[[1]] == "halfnormal") {

      heter.prior <- as.numeric(c(0, heter.prior[[3]], 1))

    } else if (heter.prior[[1]] == "uniform") {

      heter.prior <- as.numeric(c(0, heter.prior[[3]], 2))

    } else if (measure == "SMD" & heter.prior[[1]] == "logt") {

      heter.prior <- as.numeric(c(heter.prior[[2]], heter.prior[[3]], 3))

    } else if (measure != "SMD" & heter.prior[[1]] == "logt") {

      stop("There are currently no empirically-based prior distributions for MD and ROM. Choose a half-normal or a uniform prior distribution, instead")

    }



    ## Rename columns to agree with gemtc
    names(y0) <- paste0("y..",1:length(y0[1, ]),".")
    names(se0) <- paste0("se..",1:length(se0[1, ]),".")
    names(N) <- paste0("n..",1:length(N[1, ]),".")
    names(t) <- paste0("t..",1:length(t[1, ]),".")


    ## Convert one-study-per-row data to one-arm-per-row as required in GeMTC
    transform <- mtc.data.studyrow(cbind(t, y0, se0, N, na..), armVars = c('treatment'= 't', 'mean'='y', 'std.error'='se', 'sampleSize'='n'), nArmsVar='na')
    transform$treatment <- as.numeric(transform$treatment)


    ## Detect the nodes to split
    splitting <- mtc.nodesplit.comparisons(mtc.network(transform))
    colnames(splitting) <- NULL
    rownames(splitting) <- NULL

    if(dim(splitting)[1] < 1) {

      stop("There is no loop to evaluate", call. = FALSE)

      suppressMessages({
        message("Called from: netmodr::run.nodesplit")

      })

    } else {

      pair <- t(apply(apply(as.matrix(splitting, ncol = 2), 2, as.numeric), 1, sort))


      ## Define necessary model components
      jagsfit <- data.jag <- checkPair <- bi <- si <- m <- list()


      ## Parameters to save
      param.jags <- c("EM", "direct", "diff", "tau")


      for(i in 1:length(pair[, 1])){


        ## Calculate split (1 if node to split is present) and b (baseline position)
        checkPair[[i]] <- PairXY(as.matrix(t), pair[i, ])

        ## Build vector bi[i] with baseline treatment: t[i, b[i]]
        bi[[i]] <- Basetreat(as.matrix(t), checkPair[[i]][,"b"])

        ## Indexes to sweep non-baseline arms only
        m[[i]] <- NonbaseSweep(checkPair[[i]], na..)

        ## Build matrix si[i,k] with non-baseline treatments: t[i, m[i,k]]
        si[[i]] <- Sweeptreat(as.matrix(t), m[[i]])


        # Under the Independent structure with or without SMD as effect measure
        if (measure == "SMD" & assumption != "IND-CORR") {

          data.jag[[i]] <- list("y.o" = y0, "se.o" = se0, "mod" = mod, "N" = N, "t" = t, "na" = na.., "nt" = nt, "ns" = ns, "ref" = ref, "sigma" = sigma, "meand.phi" = mean.misspar, "precd.phi" = prec.misspar,
                                "split" = checkPair[[i]][, "split"], "m" = m[[i]], "bi" = bi[[i]], "si" = si[[i]], "pair" = pair[i, ], "heter.prior" = heter.prior)

        } else if (measure == "SMD" & assumption == "IND-CORR"){

          data.jag[[i]] <- list("y.o" = y0, "se.o" = se0, "mod" = mod, "N" = N, "t" = t, "na" = na.., "nt" = nt, "ns" = ns, "ref" = ref, "sigma" = sigma, "M" = M, "cov.phi" = cov.misspar, "var.phi" = var.misspar,
                                "split" = checkPair[[i]][, "split"], "m" = m[[i]], "bi" = bi[[i]], "si" = si[[i]], "pair" = pair[i, ], "heter.prior" = heter.prior)

        } else if (measure != "SMD" & assumption == "IND-CORR") {

          data.jag[[i]] <- list("y.o" = y0, "se.o" = se0, "mod" = mod, "N" = N, "t" = t, "na" = na.., "nt" = nt, "ns" = ns, "ref" = ref, "M" = M, "cov.phi" = cov.misspar, "var.phi" = var.misspar,
                                "split" = checkPair[[i]][, "split"], "m" = m[[i]], "bi" = bi[[i]], "si" = si[[i]], "pair" = pair[i, ], "heter.prior" = heter.prior)

        } else {

          data.jag[[i]] <- list("y.o" = y0, "se.o" = se0, "mod" = mod, "N" = N, "t" = t, "na" = na.., "nt" = nt, "ns" = ns, "ref" = ref, "meand.phi" = mean.misspar, "precd.phi" = prec.misspar,
                                "split" = checkPair[[i]][, "split"], "m" = m[[i]], "bi" = bi[[i]], "si" = si[[i]], "pair" = pair[i, ], "heter.prior" = heter.prior)

        } # Stop if-statement for 'measure' and 'assumption'


        ## Run the Bayesian analysis
        jagsfit[[i]] <- jags(data = data.jag[[i]], parameters.to.save = param.jags, model.file = paste0("./model/node-splitting/RE-Node-Splitting_", measure, "_Pattern-mixture_", assumption, ".txt"),
                             n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin, DIC = T)

      }  # Stop loop for 'pair'


    }

    ## Define node to split: AB=(1,2)


  } else {


    ## Binary: arm-level, wide-format dataset
    (event <- data %>% dplyr::select(starts_with("r")))                             # Number of observed events in each arm of every trial
    (mod0 <- data %>% dplyr::select(starts_with("m")))                              # Number of missing participants in each arm of every trial
    (rand <- data %>% dplyr::select(starts_with("n")))                              # Number randomised participants in each arm of every trial
    (treat <- data %>% dplyr::select(starts_with("t")))                             # Intervention studied in each arm of every trial
    na.. <- apply(treat, 1, function(x) length(which(!is.na(x))))                   # Number of interventions investigated in every trial per network
    nt <- length(table(as.matrix(treat)))                                           # Total number of interventions per network
    ns <- length(event[, 1])                                                        # Total number of included trials per network
    ref <- ifelse(missing(net.ref), which.max(table(as.matrix(treat))), net.ref)    # Reference intervention per network. If not specify, the most frequently appeared intervention in the network is selected



    ## Order by 'id of t1' < 'id of t1'
    r <- mod <- N <- t <- t0 <- treat
    for(i in 1:ns){

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


    M <- ifelse(!is.na(r), mean.misspar, NA)   # Vector of the mean value of the normal distribution of the informative missingness parameter as the number of arms in trial i (independent structure)
    prec.misspar <- 1/var.misspar
    psi.misspar <- sqrt(var.misspar)           # the lower bound of the uniform prior distribution for the prior standard deviation of the missingness parameter (hierarchical structure)
    cov.misspar <- 0.5*var.misspar             # covariance of pair of missingness parameters in a trial (independent structure)



    ## Specification of the prior distribution for the between-trial parameter
    if (heter.prior[[1]] == "halfnormal") {

      heter.prior <- as.numeric(c(0, heter.prior[[3]], 1))

    } else if (heter.prior[[1]] == "uniform") {

      heter.prior <- as.numeric(c(0, heter.prior[[3]], 2))

    } else if (heter.prior[[1]] == "lognormal")  {

      heter.prior <- as.numeric(c(heter.prior[[2]], heter.prior[[3]], 3))

    }



    ## Rename columns to agree with gemtc
    names(r) <- paste0("r..",1:length(r[1, ]),".")
    names(N) <- paste0("n..",1:length(N[1, ]),".")
    names(t) <- paste0("t..",1:length(t[1, ]),".")


    ## Convert one-study-per-row data to one-arm-per-row as required in GeMTC
    transform <- mtc.data.studyrow(cbind(t, r, N, na..), armVars = c('treatment'= 't', 'response'='r', 'sampleSize'='n'), nArmsVar='na')
    transform$treatment <- as.numeric(transform$treatment)


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
      param.jags <- c("EM", "direct", "diff", "tau")


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


        ## Condition for the Independent structure
        if (assumption != "IND-CORR") {

          data.jag[[i]] <- list("r" = r, "mod" = mod, "N" = N, "t" = t, "na" = na.., "nt" = nt, "ns" = ns, "ref" = ref, "meand.phi" = mean.misspar, "precd.phi" = prec.misspar,
                                "split" = checkPair[[i]][, "split"], "m" = m[[i]], "bi" = bi[[i]], "si" = si[[i]], "pair" = pair[i, ], "heter.prior" = heter.prior)

        } else {

          data.jag[[i]] <- list("r" = r, "mod" = mod, "N" = N, "t" = t, "na" = na.., "nt" = nt, "ns" = ns, "ref" = ref, "M" = M, "cov.phi" = cov.misspar, "var.phi" = var.misspar,
                                "split" = checkPair[[i]][, "split"], "m" = m[[i]], "bi" = bi[[i]], "si" = si[[i]], "pair" = pair[i, ], "heter.prior" = heter.prior)

        } # Stop if-statement for 'assumption'


        ## Run the Bayesian analysis
        jagsfit[[i]] <- jags(data = data.jag[[i]], parameters.to.save = param.jags, model.file = paste0("./model/node-splitting/RE-Node-Splitting_", measure, "_Pattern-mixture_", assumption, ".txt"),
                             n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin, DIC = T)

      } # Stop loop for 'pair'


    }

  }


  ## Obtain the posterior distribution of the necessary model parameters (node: 'treat1' versus 'treat2')
  EM <- data.frame(pair[, 2], pair[, 1], do.call(rbind,lapply(1:length(pair[, 1]), function(i) jagsfit[[i]]$BUGSoutput$summary[paste0("EM[", pair[i, 2], ",", pair[i, 1], "]"), c("mean", "sd", "2.5%", "97.5%", "Rhat", "n.eff")])))
  colnames(EM) <- c("treat1", "treat2", "mean", "sd", "2.5%", "97.5%", "Rhat", "n.eff")
  direct <- data.frame(pair[, 2], pair[, 1], do.call(rbind,lapply(1:length(pair[, 1]), function(i) jagsfit[[i]]$BUGSoutput$summary["direct", c("mean", "sd", "2.5%", "97.5%", "Rhat", "n.eff")])))
  colnames(direct) <- c("treat1", "treat2", "mean", "sd", "2.5%", "97.5%", "Rhat", "n.eff")
  diff <- data.frame(pair[, 2], pair[, 1], do.call(rbind,lapply(1:length(pair[, 1]), function(i) jagsfit[[i]]$BUGSoutput$summary["diff", c("mean", "sd", "2.5%", "97.5%", "Rhat", "n.eff")])))
  colnames(diff) <- c("treat1", "treat2", "mean", "sd", "2.5%", "97.5%", "Rhat", "n.eff")
  tau <- data.frame(pair[, 2], pair[, 1], do.call(rbind,lapply(1:length(pair[, 1]), function(i) jagsfit[[i]]$BUGSoutput$summary["tau", c("50%", "sd", "2.5%", "97.5%", "Rhat", "n.eff")])))
  colnames(tau) <- c("treat1", "treat2", "50%", "sd", "2.5%", "97.5%", "Rhat", "n.eff")
  dev <- do.call(rbind,lapply(1:length(pair[, 1]), function(i) jagsfit[[i]]$BUGSoutput$summary["deviance", "mean"]))
  DIC <- do.call(rbind,lapply(1:length(pair[, 1]), function(i) jagsfit[[i]]$BUGSoutput$DIC))
  pD <- do.call(rbind,lapply(1:length(pair[, 1]), function(i) jagsfit[[i]]$BUGSoutput$pD))
  model.assessment <- data.frame(pair[, 2], pair[, 1], DIC, dev, pD)
  colnames(model.assessment) <- c("treat1", "treat2", "DIC", "deviance", "pD")

  return(list(EM = EM, direct = direct, diff = diff, tau = tau, model.assessment = model.assessment))


}



## Necessary arguments (with explanations)
# data0: an arm-level, wide-format dataset where each line is a trial. The file should be a delimited text-file (.txt).
#        tX, intervention identifier in arm X (X = 1, 2, 3, 4); rX, number of events in arm X; mX, number of missing outcome data in arm X,
#        nX, number of randomised participant in arm X; id.NMA, identifier for each network (necessary even when a single network is analysed);
#        Keep.NMA, indicator of whether a network is eligible or not for the node-splitting approach (necessary even when a single network is analysed).
#        A network without closed loops or with closed loops that are informed by multi-arm trials exclusively are not eligible to apply the node-splitting approach;
#        id.Node, identifier for each eligible network to apply the node-splitting approach (necessary even when a single network is analysed).
# tau2.priors: elected predictive distributions for between-trial variance tailored to the outcome and intervention-comparison type.
#              It should be a delimited text file (.txt) where the first column in the identifier for each network (necessary even when a single network is analysed),
#              the second and the third columns are the mean and the standard deviation of the log-normal distribution, respectively.
# meta.model: a character to indicate whether random-effects or fixed-effect network meta-analysis will be applied.
#             Possible values are "RE", and "FE" to indicate a random-effects and a fixed-effect network meta-analysis, respectively.
# model: a character to indicate the prior structure for the informative missingness odd ratio (IMOR).
#        Possible values are "identical", "hierarchical", "independent".
# var.logimor: a numerical value for the prior variance of the IMOR in the logarithim value. Plausible values include [0.25, 1, and 4].
# n.chains: number of Markov chains (default: 3 - https://www.rdocumentation.org/packages/R2jags/versions/0.5-7/topics/jags).
# n.iter: number of total iterations per chain (including burn in; default: 2000 - https://www.rdocumentation.org/packages/R2jags/versions/0.5-7/topics/jags).
# n.burnin: length of burn in, i.e. number of iterations to discard at the beginning (https://www.rdocumentation.org/packages/R2jags/versions/0.5-7/topics/jags).
# n.thin: thinning rate. Must be a positive integer (https://www.rdocumentation.org/packages/R2jags/versions/0.5-7/topics/jags).
# NOTE: in networks with many closed loops, it takes time before the function provides results.

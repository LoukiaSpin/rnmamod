#' A function to perform Bayesian unrelated mean effects model
#'
#' @description This function performs a Bayesian network meta-analysis based on the unrelated mean effects model of Dias et al. (2013a) extended to address aggregate binary and continuous participant outcome data via the pattern-mixture model
#'   (Spineli, 2019; Spineli et al., 2021). This model offers a global evaluation of the plausibility of the consistency assumption in the network (Dias et al. (2013b)).
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
#'  \code{EM} \tab The summary effect estimate for each pairwise comparison observed in the network.\cr
#'  \tab \cr
#'  \code{dev.o} \tab The deviance contribution of each trial-arm based on the observed outcome.\cr
#'  \tab \cr
#'  \code{hat.par} \tab The fitted outcome at each trial-arm.\cr
#'  \tab \cr
#'  \code{tau} \tab The between-trial standard deviation (assumed common across the observed pairwise comparisons) for the whole network, when a random-effects model has been specified.\cr
#'  \tab \cr
#'  \code{m.tau} \tab The between-trial standard deviation (assumed common across the observed pairwise comparisons) for the subset of multi-arm trials, when a random-effects model has been specified.\cr
#' }
#'
#' The output also includes the following elements - the first three resulting from relevant monitored parameters:
#' \tabular{ll}{
#'  \code{leverage.o} \tab The leverage for the observed outcome at each trial-arm.\cr
#'  \tab \cr
#'  \code{sign.dev.o} \tab The sign of the difference between observed and fitted outcome at each trial-arm.\cr
#'  \tab \cr
#'  \code{model.assessment} \tab A data-frame on the measures of model assessment: deviance information criterion, number of effective parameters, and total residual deviance.\cr
#'  \tab \cr
#'  \code{jagsfit} \tab An object of S3 class \code{\link[R2jags]{jags}} with the posterior results on all monitored parameters to be used in the \code{mcmc.diagnostics} function.\cr
#' }
#' Furthermore, \code{run.UME} returns a character vector with the pairwise comparisons observed in the network, \code{obs.comp},
#' and a character vector with comparisons between non-baseline interventions observed in multi-arm trials only, \code{frail.comp}. Both vectors are used in \code{UME.plot} function.
#'
#' @details \code{run.UME} does not contain the arguments \code{data}, \code{measure}, \code{model}, \code{assumption}, \code{heter.prior}, \code{mean.misspar}, and \code{var.misspar} that are found in \code{run.model}.
#'   This is to prevent misspecifying the Bayesian model as it would make the comparison of the consistency model (via \code{run.model}) with the unrelated mean effects model meaningless.
#'   Instead, these arguments are contained in the argument \code{full} of the function. Therefore, the user needs first to apply \code{run.model}, and then use \code{run.UME} (see, 'Examples').
#'
#'   Initially, \code{run.UME} calls the \code{improved.UME} function to identify the \emph{frail comparisons}, that is, comparisons between non-baseline interventions in multi-arm trials -- not investigated in any two-arm trial of the network (Spineli, 2021).
#'   The 'original' model of Dias et al. (2013a) omits the frail comparisons from the estimation process. Consequently, the number of estimated summary effect sizes is less than those obtained by performing separate pairwise meta-analyses (see \code{run.series.meta}).
#'
#'   Then, to perform the Bayesian unrelated effect measures model, \code{run.UME} calls the \code{prepare.UME} function which contains the WinBUGS code as written by Dias et al. (2013a) for binomial and normal likelihood to analyse binary and continuous outcome data, respectively.
#'   \code{prepare.UME} has been extended to incorporate the pattern-mixture model with informative missingness parameters for binary and continuous outcome data (see, 'Details' in \code{run.model}).
#'   \code{prepare.UME} has been also 'upgraded' to account for the multi-arm trials by assigning conditional univariate normal distributions on the basic parameters of these trials, that is, effect parameters between non-baseline and baseline arms (Dias et al., 2013b; Spineli, 2021).
#'
#'   \code{run.UME} runs Bayesian unrelated mean effects model in \code{JAGS}. The progress of the simulation appears in the R console.
#'
#'   The output of \code{run.UME} is not end-user-ready. The \code{UME.plot} function uses the output of \code{run.UME} as an S3 object and processes it further to provide an end-user-ready output.
#'
#'   \code{run.UME} can be used only for a network of interventions. In the case of two interventions, the execution of the function will be stopped and an error message will be printed in the R console.
#'
#' @author {Loukia M. Spineli}
#'
#' @seealso \code{\link{run.model}}, \code{\link{UME.plot}}, \code{\link{prepare.UME}}, \code{\link{run.series.meta}}, \code{\link[R2jags]{jags}}
#'
#' @references
#' Spineli LM. A novel framework to evaluate the consistency assumption globally in a network of interventions. \emph{submitted} 2021.
#'
#' Spineli LM, Kalyvas C, Papadimitropoulou K. Continuous(ly) missing outcome data in network meta-analysis: a one-stage pattern-mixture model approach. \emph{Stat Methods Med Res} 2021. [\doi{10.1177/0962280220983544}]
#'
#' Spineli LM. An empirical comparison of Bayesian modelling strategies for missing binary outcome data in network meta-analysis. \emph{BMC Med Res Methodol} 2019;\bold{19}(1):86. [\doi{10.1186/s12874-019-0731-y}]
#'
#' Dias S, Welton NJ, Sutton AJ, Caldwell DM, Lu G, Ades AE. Evidence synthesis for decision making 4: inconsistency in networks of evidence based on randomized controlled trials. \emph{Med Decis Making} 2013a;\bold{33}(5):641--56. [\doi{10.1177/0272989X12455847}]
#'
#' Dias S, Sutton AJ, Ades AE, Welton NJ. Evidence synthesis for decision making 2: a generalized linear modeling framework for pairwise and network meta-analysis of randomized controlled trials. \emph{Med Decis Making} 2013b;\bold{33}(5):607--617. [\doi{10.1177/0272989X12458724}]
#'
#' Gelman A, Rubin DB. Inference from iterative simulation using multiple sequences. \emph{Stat Sci} 1992;\bold{7}:457--472. [\doi{10.1214/ss/1177011136}]
#'
#' @examples
#' data("nma.liu2013")
#'
#' \dontrun{
#' # Perform a random-effects network meta-analysis (consistency model)
#' res <- run.model(data = nma.liu2013,
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
#' # Run random-effects unrelated mean effects model
#' run.UME(full = res, n.chains = 3, n.iter = 10000, n.burnin = 1000, n.thin = 1)
#' }
#'
#' @export
run.UME <- function(full, n.iter, n.burnin, n.chains, n.thin) {


  ## Turn off warning when variables in the 'data.jag' are not used
  options(warn = -1)


  ## Default arguments
  data <- full$data
  measure <- full$measure
  model <- full$model
  assumption <- full$assumption
  heter.prior <- full$heter.prior
  mean.misspar <- full$mean.misspar
  var.misspar <- full$var.misspar


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


  ## Move multi-arm trials at the bottom
  t <- item$t[order(item$na, na.last = T), ]
  m <- item$m[order(item$na, na.last = T), ]
  N <- item$N[order(item$na, na.last = T), ]
  na <- sort(item$na)
  ns <- item$ns


  ## Unique comparisons with the baseline intervention
  # A function to extract numbers from a character. Source: http://stla.github.io/stlapblog/posts/Numextract.html
  Numextract <- function(string){
    unlist(regmatches(string,gregexpr("[[:digit:]]+\\.*[[:digit:]]*",string)))
  }


  ## Observed comparisons in the network
  impr.UME <- improved.UME(t, N, ns, na)
  observed.comp0 <- impr.UME$obs.comp
  observed.comp <- matrix(Numextract(observed.comp0[, 1]), nrow = length(observed.comp0[, 1]), ncol = 2, byrow = T)
  t1.obs.com <- as.numeric(as.character(observed.comp[, 1]))
  t2.obs.com <- as.numeric(as.character(observed.comp[, 2]))
  obs.comp <- paste0(t2.obs.com, "vs", t1.obs.com)


  ## Keep only comparisons with the baseline intervention
  indic0 <- list()
  for (i in 1:ns) {indic0[[i]] <- combn(t(na.omit(t(t[i, ]))), 2)[, 1:(na[i] - 1)]}
  indic <- unique(t(do.call(cbind, indic0)))
  t1.indic <- indic[, 1]
  t2.indic <- indic[, 2]
  N.obs <- length(t1.indic)


  ## Keep only comparisons with the baseline intervention in multi-arm trials
  ns.multi <- length(na[na > 2])
  if (ns.multi < 1) {
    N.obs.multi <- 0
    t1.indic.multi <- 0
    t2.indic.multi <- 0
  } else {
    indic.multi0 <- list()
    for (i in (ns - ns.multi + 1):ns) {indic.multi0[[i]] <- combn(t(na.omit(t(t[i, ]))), 2)[, 1:(na[i] - 1)]}
    indic.multi <- unique(t(do.call(cbind, indic.multi0)))
    t1.indic.multi <- indic.multi[, 1]
    t2.indic.multi <- indic.multi[, 2]
    t.indic.multi <- unique(c(t1.indic.multi, t2.indic.multi))[-impr.UME$ref.base]
    t.indic.multi2 <- unique(c(t1.indic.multi, t2.indic.multi))
    N.obs.multi <- length(t1.indic.multi)

    ## Is the subset of multi-arm trials a connected network?
    multi.network <- pairwise(as.list(t[(ns - ns.multi + 1):ns, ]), mean = as.list(N[(ns - ns.multi + 1):ns, ]), sd = as.list(N[(ns - ns.multi + 1):ns, ]), n = as.list(N[(ns - ns.multi + 1):ns, ]), data = cbind(t[(ns - ns.multi + 1):ns, ], N[(ns - ns.multi + 1):ns, ], N[(ns - ns.multi + 1):ns, ], N[(ns - ns.multi + 1):ns, ]), studlab = 1:ns.multi)
    connected <- netconnection(treat1, treat2, studlab, data = multi.network)$n.subnets

    ## For the case of a disconnected network of multi-arm trials
    if (connected > 1) {
      dist.mat <- netconnection(treat1, treat2, studlab, data = multi.network)$D.matrix
      group0 <- apply(dist.mat, 2, function(x) length(which(!is.infinite(x))))
      group <- data.frame("treat" = attributes(group0)$names, "freq" = group0)

      if (length(unique(group$freq)) < 2) {
        find.groups <- split(group[, 1], sort(rep_len(1:(length(group$freq)/unique(group$freq)), length(group[, 1]))))
        t.m2 <- data.frame("t.m1" = as.numeric(group[, 1]), "t.m2" = as.numeric(unlist(lapply(1:(length(group$freq)/unique(group$freq)), function(i) rep(min(as.numeric(find.groups[[i]])), unique(group$freq))))))
        #find.groups <- split(group[, 1], sort(rep_len(1:unique(group$freq), length(group[, 1]))))
        #t.m2 <- data.frame("t.m1" = as.numeric(group[, 1]), "t.m2" = as.numeric(unlist(lapply(1:unique(group$freq), function(i) rep(min(as.numeric(find.groups[[i]])),  unique(group$freq)[i])))))
      } else {
        find.groups <- split(group[, 1], rep(1:length(unique(group$freq)), unique(group$freq)))
        t.m2 <- data.frame("t.m1" = as.numeric(group[, 1]), "t.m2" = as.numeric(unlist(lapply(1:length(unique(group$freq)), function(i) rep(min(as.numeric(find.groups[[i]])),  unique(group$freq)[i])))))
      }

      ref.m <- rep(NA, ns)
      for (i in (ns - ns.multi + 1):ns) {ref.m[i] <- t.m2[is.element(t.m2[, 1], t[i, 1]), 2]}

      ref.nbase.multi <- rep(NA, impr.UME$nbase.multi)
      for (i in 1:impr.UME$nbase.multi) {ref.nbase.multi[i] <- t.m2[is.element(t.m2[, 1], cbind(impr.UME$t1.bn[i], impr.UME$t2.bn[i]) ), 2]}
    }
  }


  ## Data in list format for R2jags
  data.jag <- list("m" = m,
                   "N" = N,
                   "t" = t,
                   "na" = na,
                   "nt" = item$nt,
                   "ns" = ns,
                   "ref" = ifelse(is.element(assumption, c("HIE-ARM", "IDE-ARM")), item$ref, NA),
                   "I" = item$I[order(item$na, na.last = T), ],
                   "M" = ifelse(!is.na(m), mean.misspar, NA),
                   "cov.phi" = 0.5*var.misspar,
                   "var.phi" = var.misspar,
                   "meand.phi" = mean.misspar,
                   "precd.phi" = 1/var.misspar,
                   "heter.prior" = heter.prior,
                   "t1" = t1.indic,
                   "t2" = t2.indic,
                   "N.obs" = N.obs)


  if (is.element(measure, c("MD", "SMD", "ROM"))) {
    data.jag <- append(data.jag, list("y" = item$y0[order(item$na, na.last = T), ], "se" = item$se0[order(item$na, na.last = T), ], "y.m" = item$y0[order(item$na, na.last = T), ]))
  } else if (measure == "OR") {
    data.jag <- append(data.jag, list("r" = item$r[order(item$na, na.last = T), ], "r.m"= item$r[order(item$na, na.last = T), ]))
  }


  data.jag <- if (max(na) > 2 & !is.null(impr.UME$nbase.multi)) {
    append(data.jag, list("ns.multi" = ns.multi,
                          "t1.bn" = impr.UME$t1.bn,
                          "t2.bn" = impr.UME$t2.bn,
                          "nbase.multi" = impr.UME$nbase.multi,
                          "ref.base" = impr.UME$ref.base,                                     # For connected network of multi-arm trials
                          "N.t.m" = length(t.indic.multi),                                    # For connected network of multi-arm trials
                          "t.m" = t.indic.multi,                                              # For connected network of multi-arm trials
                          "ref.m" = if (connected > 1) {ref.m} else {1},                      # For *dis*connected network of multi-arm trials
                          "ref.nbase.multi" = if (connected > 1) {ref.nbase.multi} else {1},  # For *dis*connected network of multi-arm trials
                          "N.t.m2" = ifelse(connected > 1, length(t.indic.multi2), 1),        # For *dis*connected network of multi-arm trials
                          "t.m2" = if (connected > 1) {t.m2} else {t.indic.multi}))           # For *dis*connected network of multi-arm trials
    } else if (max(na) < 3 || is.null(impr.UME$nbase.multi)) {
      append(data.jag, list("ns.multi" = ns.multi,
                            "t1.bn" = t1.indic,            # Pseudo-vector (so that the model runs)
                            "t2.bn" = t1.indic,            # Pseudo-vector
                            "nbase.multi" = 0,             # Pseudo-number
                            "ref.base" = 1,                # Pseudo-number - connected network of multi-arm trials
                            "N.t.m" = length(2:5),         # Pseudo-number - connected network of multi-arm trials
                            "t.m" = 2:5,                   # Pseudo-vector - connected network of multi-arm trials
                            "ref.m" = 1,                   # Pseudo-number - *dis*connected network of multi-arm trials
                            "ref.nbase.multi" = 1,         # Pseudo-number - *dis*connected network of multi-arm trials
                            "N.t.m2" = 1,                  # Pseudo-number - *dis*connected network of multi-arm trials
                            "t.m2" = 2:5))                 # Pseudo-vector - *dis*connected network of multi-arm trials
    }


  ## Define the nodes to be monitored
  param.jags <- if (model == "RE") {
    c("EM", "dev.o", "resdev.o", "totresdev.o", "tau", "m.tau", "hat.par")
  } else {
    c("EM", "dev.o", "resdev.o", "totresdev.o", "hat.par")
  }



  ## Run the Bayesian analysis
  jagsfit <- jags(data = data.jag,
                  parameters.to.save = param.jags,
                  model.file = textConnection(prepare.UME(measure, model, assumption, connected)),
                  n.chains = n.chains,
                  n.iter = n.iter,
                  n.burnin = n.burnin,
                  n.thin = n.thin,
                  DIC = F)


  ## Turn summary of posterior results (R2jags object) into a data-frame to select model parameters (using 'dplyr')
  getResults <- as.data.frame(t(jagsfit$BUGSoutput$summary))

  # Effect size of all unique pairwise comparisons
  EM <- t(getResults %>% dplyr::select(starts_with("EM[")))

  # Between-trial standard deviation
  tau <- t(getResults %>% dplyr::select(starts_with("tau")))
  m.tau <- t(getResults %>% dplyr::select(starts_with("m.tau"))) # For the subnetwork of multi-arm trials

  # Trial-arm deviance contribution for observed outcome
  dev.o <- t(getResults %>% dplyr::select(starts_with("dev.o[")))

  # Fitted/predicted outcome
  hat.par <- t(getResults %>% dplyr::select(starts_with("hat.par")))

  # Total residual deviance
  dev <- jagsfit$BUGSoutput$summary["totresdev.o", "mean"]


  ## Calculate the deviance at posterior mean of fitted values
  # Turn 'number of observed' into a vector (first column, followed by second column, and so on)
  m.new <- suppressMessages({as.vector(na.omit(melt(m)[, 2]))})
  N.new <- suppressMessages({as.vector(na.omit(melt(N)[, 2]))})
  obs <- N.new - m.new

  if (is.element(measure, c("MD", "SMD", "ROM"))) {
    # Turn 'y0', 'se0'into a vector (first column, followed by second column, and so on)
    y0.new <- suppressMessages({as.vector(na.omit(melt(item$y0[order(item$na, na.last = T), ])[, 2]))})
    se0.new <- suppressMessages({as.vector(na.omit(melt(item$se0[order(item$na, na.last = T), ])[, 2]))})

    # Deviance at the posterior mean of the fitted mean outcome
    dev.post.o <- (y0.new - as.vector(hat.par[, 1]))*(y0.new - as.vector(hat.par[, 1]))*(1/se0.new^2)

    # Sign of the difference between observed and fitted mean outcome
    sign.dev.o <- sign(y0.new - as.vector(hat.par[, 1]))

  } else {
    # Turn 'r' and number of observed into a vector (first column, followed by second column, and so on)
    r.new <- suppressMessages({as.vector(na.omit(melt(item$r[order(item$na, na.last = T), ])[, 2]))})

    # Correction for zero events in trial-arm
    r0 <- ifelse(r.new == 0, r.new + 0.01, ifelse(r.new == obs, r.new - 0.01, r.new))

    # Deviance at the posterior mean of the fitted response
    dev.post.o <- 2*(r0*(log(r0) - log(as.vector(hat.par[, 1]))) + (obs - r0)*(log(obs - r0) - log(obs - as.vector(hat.par[, 1]))))

    # Sign of the difference between observed and fitted response
    sign.dev.o <- sign(r0 - as.vector(hat.par[, 1]))
  }


  ## Obtain the leverage for observed outcomes
  leverage.o <- as.vector(dev.o[, 1]) - dev.post.o


  ## Number of effective parameters
  pD <- dev - sum(dev.post.o)


  ## Deviance information criterion
  DIC <- pD + dev


  ## A data-frame on the measures of model assessment: DIC, pD, and total residual deviance
  model.assessment <- data.frame(DIC, pD, dev)


  ## Collect the minimum results at common
  results <- if (model == "RE") {
    list(EM = EM,
         dev.o = dev.o,
         hat.par = hat.par,
         leverage.o = leverage.o,
         sign.dev.o = sign.dev.o,
         tau = tau,
         model.assessment = model.assessment,
         obs.comp = obs.comp,
         jagsfit = jagsfit,
         data = data,
         measure = measure)
  } else {
    list(EM = EM,
         dev.o = dev.o,
         hat.par = hat.par,
         leverage.o = leverage.o,
         sign.dev.o = sign.dev.o,
         model.assessment = model.assessment,
         obs.comp = obs.comp,
         jagsfit = jagsfit,
         data = data,
         measure = measure)
  }


  ## Return different list of results according to a condition
  if (is.null(impr.UME$nbase.multi)) {
    return(results)
  } else {
    return(append(results, list(m.tau = m.tau, frail.comp = paste0(impr.UME$t2.bn, "vs", impr.UME$t1.bn))))
  }
}


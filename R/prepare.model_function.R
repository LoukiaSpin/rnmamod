#' WinBUGS code for Bayesian pairwise or network meta-analysis and
#' meta-regression
#'
#' @description
#'   The WinBUGS code, as written by Dias et al. (2013) to run a one-stage
#'   Bayesian network meta-analysis, extended to incorporate the pattern-mixture
#'   model for binary or continuous missing participant outcome data (Spineli
#'   et al., 2021; Spineli, 2019). The model has been also extended to
#'   incorporate a trial-level covariate to apply meta-regression
#'   (Cooper et al., 2009). In the case of two interventions, the code boils
#'   down to a one-stage Bayesian pairwise meta-analysis with pattern-mixture
#'   model (Turner et al., 2015; Spineli et al, 2021).
#'
#' @param measure Character string indicating the effect measure. For a binary
#'   outcome, the following can be considered: \code{"OR"}, \code{"RR"} or
#'   \code{"RD"} for the odds ratio, relative risk, and risk difference,
#'   respectively. For a continuous outcome, the following can be considered:
#'   \code{"MD"}, \code{"SMD"}, or \code{"ROM"} for mean difference,
#'   standardised mean difference and ratio of means, respectively.
#' @param model Character string indicating the analysis model with values
#'   \code{"RE"}, or \code{"FE"} for the random-effects and fixed-effect model,
#'   respectively. The default argument is \code{"RE"}.
#' @param covar_assumption Character string indicating the structure of the
#'   intervention-by-covariate interaction, as described in
#'   Cooper et al., (2009). Set \code{covar_assumption} equal to one of the
#'   following, when meta-regression is performed: \code{"exchangeable"},
#'   \code{"independent"}, and \code{"common"}. Assign \code{"NO"} to perform
#'   pairwise or network meta-analysis.
#' @param assumption Character string indicating the structure of the
#'   informative missingness parameter. Set \code{assumption} equal to one of
#'   the following: \code{"HIE-COMMON"}, \code{"HIE-TRIAL"}, \code{"HIE-ARM"},
#'   \code{"IDE-COMMON"}, \code{"IDE-TRIAL"}, \code{"IDE-ARM"},
#'   \code{"IND-CORR"}, or \code{"IND-UNCORR"}. The default argument is
#'   \code{"IDE-ARM"}. The abbreviations \code{"IDE"}, \code{"HIE"}, and
#'   \code{"IND"} stand for identical, hierarchical and independent,
#'   respectively. \code{"CORR"} and \code{"UNCORR"} stand for correlated and
#'   uncorrelated, respectively.
#'
#' @return An R character vector object to be passed to \code{\link{run_model}}
#'   and \code{\link{run_metareg}} through the
#'   \code{\link[base:textConnection]{textConnection}} function as the argument
#'   \code{object}.
#'
#' @details \code{prepare_model} creates the model in the JAGS dialect
#'   of the BUGS language. The output of this function constitutes the argument
#'   \code{model.file} of the \code{\link[R2jags:jags]{jags}} function (in the
#'   R-package \href{https://CRAN.R-project.org/package=R2jags}{R2jags}) via the
#'   \code{\link[base:textConnection]{textConnection}} function.
#'
#' @author {Loukia M. Spineli}
#'
#' @seealso \code{\link{run_metareg}}, \code{\link{run_model}},
#'    \code{\link[R2jags:jags]{jags}},
#'    \code{\link[base:textConnection]{textConnection}}
#'
#' @references
#' Cooper NJ, Sutton AJ, Morris D, Ades AE, Welton NJ. Addressing between-study
#' heterogeneity and inconsistency in mixed treatment comparisons: Application
#' to stroke prevention treatments in individuals with non-rheumatic atrial
#' fibrillation. \emph{Stat Med} 2009;\bold{28}(14):1861--81.
#' doi: 10.1002/sim.3594
#'
#' Dias S, Sutton AJ, Ades AE, Welton NJ. Evidence synthesis for decision
#' making 2: a generalized linear modeling framework for pairwise and network
#' meta-analysis of randomized controlled trials. \emph{Med Decis Making}
#' 2013;\bold{33}(5):607--17. doi: 10.1177/0272989X12458724
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
#' Turner NL, Dias S, Ades AE, Welton NJ. A Bayesian framework to account
#' for uncertainty due to missing binary outcome data in pairwise meta-analysis.
#' \emph{Stat Med} 2015;\bold{34}(12):2062--80. doi: 10.1002/sim.6475
#'
#' @export
prepare_model <- function(measure,
                          model,
                          covar_assumption,
                          assumption) {

  stringcode <- "model {
                    for (i in 1:ns) {\n"

  stringcode <- if (model == "RE") {
    paste(stringcode, "delta[i, 1] <- 0
                       w[i, 1] <- 0
                       u[i] ~ dnorm(0, .0001)\n")
  } else {
    paste(stringcode, "u[i] ~ dnorm(0, .0001)\n")
  }

  stringcode <- if (measure == "SMD") {
    paste(stringcode, "theta[i, 1] <- u[i]
                       sigma[i] <- sqrt(sum(nom[i, 1:na[i]])/(sum(c[i, 1:na[i]]) - na[i]))
                       a[i] <- sum(N[i, 1:na[i]] - 1)/2
                       b[i] <- sum(N[i, 1:na[i]] - 1)/(2*sigma[i]*sigma[i])
                       var.pooled[i] ~ dgamma(a[i], b[i])
                       sd.pooled[i] <- sqrt(var.pooled[i])\n")
  } else if (is.element(measure, c("MD", "ROM"))) {
    paste(stringcode, "theta[i, 1] <- u[i]\n")
  } else if (is.element(measure, c("OR", "RR", "RD"))) {
    paste(stringcode, "logit(p[i, 1]) <- u[i]\n")
  }

  stringcode <- paste(stringcode, "for (k in 1:na[i]) { \n")

  stringcode <- if (measure == "SMD") {
    paste(stringcode, "prec.o[i, k] <- pow(se.o[i, k], -2)
                       y.o[i, k] ~ dnorm(theta.o[i, k], prec.o[i, k])
                       c[i, k] <- N[i, k] - m[i, k]
                       sd.obs[i, k] <- se.o[i, k]*sqrt(c[i, k])
                       nom[i, k] <- pow(sd.obs[i, k], 2)*(c[i, k] - 1)\n")
  } else if (is.element(measure, c("MD", "ROM"))) {
    paste(stringcode, "prec.o[i, k] <- pow(se.o[i, k], -2)
                       y.o[i, k] ~ dnorm(theta.o[i, k], prec.o[i, k])\n")
  } else if (is.element(measure, c("OR", "RR", "RD"))) {
    paste(stringcode, "r[i, k] ~ dbin(p_o[i, k], obs[i, k])
                       obs[i, k] <- N[i, k] - m[i, k]\n")
  }

  stringcode <- if (is.element(measure, c("MD", "SMD"))) {
    paste(stringcode, "theta.o[i, k] <- theta[i, k] - phi.m[i, k]*q[i, k]\n")
  } else if (measure == "ROM") {
    paste(stringcode, "theta.o[i, k] <- theta[i, k]/(1 - q[i, k]*(1 - exp(phi.m[i, k])))\n")
  } else if (is.element(measure, c("OR", "RR", "RD"))) {
    paste(stringcode, "p_o[i, k] <- max(0, min(1, ((-((q[i, k] - p[i, k])*(1 - exp(phi.m[i, k])) - 1) - sqrt((pow(((q[i, k] - p[i, k])*(1 - exp(phi.m[i, k])) - 1), 2)) -
                                       ((4*p[i, k])*(1 - q[i, k])*(1 - exp(phi.m[i, k])))))/(2*(1 - q[i, k])*(1 - exp(phi.m[i, k]))))))\n")
  }

  stringcode <- paste(stringcode, "q[i, k] <- q0[i, k]*I[i, k]
                                   m[i, k] ~ dbin(q0[i, k], N[i, k])
                                   q0[i, k] ~ dunif(0, 1)\n")

  stringcode <- if (!is.element(measure, c("OR", "RR", "RD"))) {
    paste(stringcode, "hat.par[i, k] <- theta.o[i, k]
                       dev.o[i, k] <- (y.o[i, k] - theta.o[i, k])*(y.o[i, k] - theta.o[i, k])*prec.o[i, k]
                       }\n")
  } else if (is.element(measure, c("OR", "RR", "RD"))) {
    paste(stringcode, "hat.par[i, k] <- rhat[i, k]
                       rhat[i, k] <- p_o[i, k]*obs[i, k]
                       dev.o[i, k] <- 2*(r[i, k]*(log(r[i, k]) - log(rhat[i, k])) + (obs[i, k] - r[i, k])*(log(obs[i, k] - r[i, k]) - log(obs[i, k] - rhat[i, k])))
                       }\n")
  }

  stringcode <- paste(stringcode, "resdev.o[i] <- sum(dev.o[i, 1:na[i]])
                                   for (k in 2:na[i]) {\n")

  stringcode <- if (measure == "MD") {
    paste(stringcode, "theta[i, k] <- u[i] + delta.star[i, k]\n")
  } else if (measure == "SMD") {
    paste(stringcode, "theta[i, k] <- u[i] + sd.pooled[i]*delta.star[i, k]\n")
  } else if (measure == "ROM") {
    paste(stringcode, "theta[i, k] <- u[i]*exp(delta.star[i, k])\n")
  } else if (is.element(measure, c("OR", "RR", "RD"))) {
    paste(stringcode, "logit(p[i, k]) <- u[i] + delta.star[i, k]\n")
  }

  stringcode <- if (model == "RE") {
    paste(stringcode, "delta.star[i, k] <- delta[i, k] + Beta[i, k]
                       delta[i, k] ~ dnorm(md[i, k], precd[i, k])
                       md[i, k] <- (d[t[i, k]]*indic[i, k] - d[t[i, 1]]*indic[i, 1]) + sw[i, k]
                       precd[i, k] <- (2 * prec * (k - 1)) / k
                       w[i, k] <- delta[i, k] - (d[t[i, k]]*indic[i, k] - d[t[i, 1]]*indic[i, 1])
                       sw[i, k] <- sum(w[i, 1:(k - 1)])/(k - 1)
                       }}\n")
  } else {
    paste(stringcode, "delta.star[i, k] <- (d[t[i, k]]*indic[i, k] - d[t[i, 1]]*indic[i, 1]) + Beta[i, k]
                       }}\n")
  }

  stringcode <- if (covar_assumption == "NO") {
    paste(stringcode, "for (i in 1:ns) {
                         for (k in 2:na[i]) {
                           Beta[i, k] <- 0
                       }}\n")
  } else if (is.element(covar_assumption, c("exchangeable", "independent"))) {
    paste(stringcode, "for (i in 1:ns) {
                         for (k in 2:na[i]) {
                           Beta[i, k] <- (beta[t[i, k]]*indic[i, k] - beta[t[i, 1]]*indic[i, 1])*(cov.vector[i]*(1 - equals(cov.vector[i], 0)) + cov.matrix[i, k]*equals(cov.vector[i], 0))
                       }}\n")
  } else if (covar_assumption == "common") {
    paste(stringcode, "for (i in 1:ns) {
                         for (k in 2:na[i]) {
                           Beta[i, k] <- beta*(cov.vector[i]*(1 - equals(cov.vector[i], 0)) + cov.matrix[i, k]*equals(cov.vector[i], 0))
                       }}\n")
  }

  stringcode <- paste(stringcode, "totresdev.o <- sum(resdev.o[])
                                   d[ref] <- 0
                                   d.n[ref] <- 0
                                   for (t in 1:(ref - 1)) {
                                     d[t] ~ dnorm(0, 0.0001)
                                     d.n[t] <- d[t]*equals(min(t, ref), ref) + d[t]*(-1)*equals(min(t, ref), t) + beta.n[t] * cov_value
                                   }
                                   for (t in (ref + 1):nt) {
                                     d[t] ~ dnorm(0, 0.0001)
                                     d.n[t] <- d[t]*equals(min(t, ref), ref) + d[t]*(-1)*equals(min(t, ref), t) + beta.n[t] * cov_value
                                   }\n")


  stringcode <- if (is.element(measure, c("OR", "RR", "RD"))) {
    paste(stringcode, "mean_logit_base_event <- ref_base[1]
                       prec_logit_base_event <- ref_base[2]*(1 - equals(ref_base[1], ref_base[2])) + equals(ref_base[1], ref_base[2])
                       logit_base_risk ~ dnorm(mean_logit_base_event, prec_logit_base_event)
                       base_risk_logit <- logit_base_risk*(1 - equals(ref_base[1], ref_base[2])) + ref_base[1]*equals(ref_base[1], ref_base[2])
                       base_risk <- exp(base_risk_logit)/(1 + exp(base_risk_logit))
                       for (t in 1:nt) {
                         logit(abs_risk[t]) <- base_risk_logit + d[t] + beta.t[t] * cov_value
                       }\n")
  } else {
    paste(stringcode, " ")
  }

  stringcode <- if (covar_assumption == "exchangeable") {
    paste(stringcode, "beta[ref] <- 0
                       beta.t[ref] <- 0
                       beta.n[ref] <- 0
                       for (t in 1:(ref - 1)) {
                         beta[t] ~ dnorm(mean.B, prec.B)
                         beta.t[t] <- beta[t]
                         beta.n[t] <- beta[t]*equals(min(t, ref), ref) + beta[t]*(-1)*equals(min(t, ref), t)
                       }
                       for (t in (ref + 1):nt) {
                         beta[t] ~ dnorm(mean.B, prec.B)
                         beta.t[t] <- beta[t]
                         beta.n[t] <- beta[t]*equals(min(t, ref), ref) + beta[t]*(-1)*equals(min(t, ref), t)
                       }
                       mean.B ~ dnorm(0, .0001)
                       prec.B <- 1/pow(beta.SD,2)
                       beta.SD ~ dnorm(0, 1)I(0, )
                       for (c in 1:(nt - 1)) {
                         for (k in (c + 1):nt) {
                           beta.all[k, c] <- beta.n[k] - beta.n[c]
                       }}\n")
  } else if (covar_assumption == "independent") {
    paste(stringcode, "beta[ref] <- 0
                       beta.t[ref] <- 0
                       beta.n[ref] <- 0
                       for (t in 1:(ref - 1)) {
                         beta[t] ~ dnorm(0, 0.0001)
                         beta.t[t] <- beta[t]
                         beta.n[t] <- beta[t]*equals(min(t, ref), ref) + beta[t]*(-1)*equals(min(t, ref), t)
                       }
                       for (t in (ref + 1):nt) {
                         beta[t] ~ dnorm(0, 0.0001)
                         beta.t[t] <- beta[t]
                         beta.n[t] <- beta[t]*equals(min(t, ref), ref) + beta[t]*(-1)*equals(min(t, ref), t)
                       }
                       for (c in 1:(nt - 1)) {
                         for (k in (c + 1):nt) {
                           beta.all[k, c] <- beta.n[k] - beta.n[c]
                       }}\n")
  } else if (covar_assumption == "common") {
    paste(stringcode, "beta ~ dnorm(0, 0.0001)
                       beta.t[ref] <- 0
                       beta.n[ref] <- 0
                       for (t in 1:(ref - 1)) {
                         beta.t[t] <- beta
                         beta.n[t] <- beta*equals(min(t, ref), ref) + beta*(-1)*equals(min(t, ref), t)
                       }
                       for (t in (ref + 1):nt) {
                         beta.t[t] <- beta
                         beta.n[t] <- beta*equals(min(t, ref), ref) + beta*(-1)*equals(min(t, ref), t)
                       }
                       for (c in 1:(nt - 1)) {
                         for (k in (c + 1):nt) {
                           beta.all[k, c] <- beta.n[k] - beta.n[c]
                       }}\n")
  } else if (covar_assumption == "NO") {
    paste(stringcode, " ")
  }

  stringcode <- if (assumption == "HIE-ARM") {
    paste(stringcode, "for (i in 1:ns) {
                         for (k in 1:na[i]) {
                           phi.m[i, k] <- phi[i, k]
                           phi[i, k] ~ dnorm(mean.phi[t[i, k]], prec.phi[t[i, k]])
                        }}
                        mean.phi[ref] ~ dnorm(meand.phi[2], precd.phi)
                        prec.phi[ref] <- pow(sd.phi[ref], -2)
                        sd.phi[ref] ~ dunif(0, psi.phi)
                        for (t in 1:(ref - 1)) {
                          mean.phi[t] ~ dnorm(meand.phi[1], precd.phi)
                          prec.phi[t] <- pow(sd.phi[t], -2)
                          sd.phi[t] ~ dunif(0, psi.phi)
                        }
                        for (t in (ref + 1):nt) {
                          mean.phi[t] ~ dnorm(meand.phi[1], precd.phi)
                          prec.phi[t] <- pow(sd.phi[t], -2)
                          sd.phi[t] ~ dunif(0, psi.phi)
                        }
                        psi.phi <- pow(precd.phi, -2)\n")
  } else if (assumption == "HIE-TRIAL") {
    paste(stringcode, "for (i in 1:ns) {
                          for (k in 1:na[i]) {
                            phi.m[i, k] <- phi[i, k]
                            phi[i, k] ~ dnorm(mean.phi[i], prec.phi[i])
                       }}
                       for (i in 1:ns) {
                         mean.phi[i] ~ dnorm(meand.phi, precd.phi)
                         prec.phi[i] <- pow(sd.phi[i], -2)
                         sd.phi[i] ~ dunif(0, psi.phi)
                       }
                       psi.phi <- pow(precd.phi, -2)\n")
  } else if (assumption == "HIE-COMMON") {
    paste(stringcode, "for (i in 1:ns) {
                         for (k in 1:na[i]) {
                           phi.m[i, k] <- phi[i, k]
                           phi[i, k] ~ dnorm(mean.phi, prec.phi)
                       }}
                       mean.phi ~ dnorm(meand.phi, precd.phi)
                       prec.phi <- pow(sd.phi, -2)
                       sd.phi ~ dunif(0, psi.phi)
                       psi.phi <- pow(precd.phi, -2)\n")
  } else if (assumption == "IDE-ARM") {
    paste(stringcode, "for (i in 1:ns) {
                         for (k in 1:na[i]) {
                           phi.m[i, k] <- phi[t[i, k]]
                       }}
                       phi[ref] ~ dnorm(meand.phi[2], precd.phi)
                       for (t in 1:(ref - 1)) {
                         phi[t] ~ dnorm(meand.phi[1], precd.phi)
                       }
                       for (t in (ref + 1):nt) {
                         phi[t] ~ dnorm(meand.phi[1], precd.phi)
                       }\n")
  } else if (assumption == "IDE-TRIAL") {
    paste(stringcode, "for (i in 1:ns) {
                         for (k in 1:na[i]) {
                           phi.m[i, k] <- phi[i]
                       }}
                       for (i in 1:ns) {
                         phi[i] ~ dnorm(meand.phi, precd.phi)
                       }\n")
  } else if (assumption == "IDE-COMMON") {
    paste(stringcode, "for (i in 1:ns) {
                         for (k in 1:na[i]) {
                           phi.m[i, k] <- phi
                       }}
                       phi ~ dnorm(meand.phi, precd.phi)\n")
  } else if (assumption == "IND-CORR") {
    paste(stringcode, "for (i in 1:ns) {
                         for (k in 1:na[i]) {
                           phi.m[i, k] <- phi[i, k]
                           for (l in 1:na[i]) {
                             V[i, k, l] <- cov.phi*(1 - equals(k, l)) + var.phi*equals(k, l)
                         }}
                         Omega[i, 1:na[i], 1:na[i]] <- inverse(V[i, 1:na[i], 1:na[i]])
                         phi[i, 1:na[i]] ~ dmnorm(M[i, 1:na[i]], Omega[i, 1:na[i], 1:na[i]])
                       }\n")
  } else if (assumption == "IND-UNCORR") {
    paste(stringcode, "for (i in 1:ns) {
                         for (k in 1:na[i]) {
                           phi.m[i, k] <- phi[i, k]
                           phi[i, k] ~ dnorm(meand.phi, precd.phi)
                       }}\n")
  }

  stringcode <- if (!is.element(measure, c("RR", "RD"))) {
    paste(stringcode, "sorted <- rank(d.n[])\n")
  } else if (is.element(measure, c("RR", "RD"))) {
    paste(stringcode, "EM.ref[ref] <- 0
                       sorted <- rank(EM.ref[])   # RR or RD
                       sorted.LOR <- rank(d.n[])
                       for (t in 1:nt) {
                         order.LOR[t] <- (nt + 1 - sorted.LOR[t])*equals(D, 1) + sorted.LOR[t]*(1 - equals(D, 1))
                         most.effective.LOR[t] <- equals(order.LOR[t], 1)
                         for (l in 1:nt) {
                           effectiveness.LOR[t, l] <- equals(order.LOR[t], l)
                           cumeffectiveness.LOR[t, l] <- sum(effectiveness.LOR[t, 1:l])
                         }
                         SUCRA.LOR[t] <- sum(cumeffectiveness.LOR[t, 1:(nt - 1)])/(nt - 1)
                       }\n")
  }

  stringcode <- paste(stringcode, "for (t in 1:nt) {
                                     order[t] <- (nt + 1 - sorted[t])*equals(D, 1) + sorted[t]*(1 - equals(D, 1))
                                     most.effective[t] <- equals(order[t], 1)
                                     for (l in 1:nt) {
                                       effectiveness[t, l] <- equals(order[t], l)
                                       cumeffectiveness[t, l] <- sum(effectiveness[t, 1:l])
                                     }
                                     SUCRA[t] <- sum(cumeffectiveness[t, 1:(nt - 1)])/(nt - 1)
                                   }\n")

  stringcode <- if (model == "RE" & !is.element(measure, c("RR", "RD"))) {
    paste(stringcode, "for (c in 1:(nt - 1)) {
                         for (k in (c + 1):nt) {
                           EM[k, c] <- d.n[k] - d.n[c]
                           EM.pred[k, c] ~ dnorm(EM[k, c], prec)
                       }}\n")
  } else if (model == "FE" & !is.element(measure, c("RR", "RD"))) {
    paste(stringcode, "for (c in 1:(nt - 1)) {
                         for (k in (c + 1):nt) {
                           EM[k, c] <- d.n[k] - d.n[c]
                       }}\n")
  } else if (model == "RE" & measure == "RR") {
    paste(stringcode, "for (t in 1:(ref - 1)) {
                         EM.ref.n[t] <- (d[t] + beta.t[t] * cov_value) - log(1 - (1 - exp(d[t] + beta.t[t] * cov_value))*base_risk)
                         EM.ref[t] <- EM.ref.n[t]*equals(min(t, ref), ref) + EM.ref.n[t]*(-1)*equals(min(t, ref), t)
                       }
                       for (t in (ref + 1):nt) {
                         EM.ref.n[t] <- (d[t] + beta.t[t] * cov_value) - log(1 - (1 - exp(d[t] + beta.t[t] * cov_value))*base_risk)
                         EM.ref[t] <- EM.ref.n[t]*equals(min(t, ref), ref) + EM.ref.n[t]*(-1)*equals(min(t, ref), t)
                       }
                       for (c in 1:(nt - 1)) {
                         for (k in (c + 1):nt) {
                           EM.LOR[k, c] <- d.n[k] - d.n[c] # LOR
                           EM[k, c] <- EM.LOR[k, c] - log(1 - abs_risk[c]*(1 - exp(EM.LOR[k, c]))) # LRR
                           EM.pred.LOR[k, c] ~ dnorm(EM.LOR[k, c], prec) # LOR
                           EM.pred[k, c] <- EM.pred.LOR[k, c] - log(1 - abs_risk[c]*(1 - exp(EM.pred.LOR[k, c]))) # LRR
                       }}\n")
  } else if (model == "RE" & measure == "RD") {
    paste(stringcode, "for (t in 1:(ref - 1)) {
                         EM.ref.RR.n[t] <- (d[t] + beta.t[t] * cov_value) - log(1 - (1 - exp(d[t] + beta.t[t] * cov_value))*base_risk)
                         EM.ref.RR[t] <- EM.ref.RR.n[t]*equals(min(t, ref), ref) + EM.ref.RR.n[t]*(-1)*equals(min(t, ref), t)
                         EM.ref.n[t] <- (exp(EM.ref.RR.n[t]) - 1)*base_risk
                         EM.ref[t] <- EM.ref.n[t]*equals(min(t, ref), ref) + EM.ref.n[t]*(-1)*equals(min(t, ref), t)
                       }
                       for (t in (ref + 1):nt) {
                         EM.ref.RR.n[t] <- (d[t] + beta.t[t] * cov_value) - log(1 - (1 - exp(d[t] + beta.t[t] * cov_value))*base_risk)
                         EM.ref.RR[t] <- EM.ref.RR.n[t]*equals(min(t, ref), ref) + EM.ref.RR.n[t]*(-1)*equals(min(t, ref), t)
                         EM.ref.n[t] <- (exp(EM.ref.RR.n[t]) - 1)*base_risk
                         EM.ref[t] <-  EM.ref.n[t]*equals(min(t, ref), ref) + EM.ref.n[t]*(-1)*equals(min(t, ref), t)
                       }
                       for (c in 1:(nt - 1)) {
                         for (k in (c + 1):nt) {
                           EM.LOR[k, c] <- d.n[k] - d.n[c] # LOR
                           EM.LRR[k, c] <- EM.LOR[k, c] - log(1 - abs_risk[c]*(1 - exp(EM.LOR[k, c]))) # LRR
                           EM[k, c] <- abs_risk[c]*(exp(EM.LRR[k, c]) - 1) # RD
                           EM.pred.LOR[k, c] ~ dnorm(EM.LOR[k, c], prec) # LOR
                           EM.pred.LRR[k, c] <- EM.pred.LOR[k, c] - log(1 - abs_risk[c]*(1 - exp(EM.pred.LOR[k, c]))) # LRR
                           EM.pred[k, c] <- abs_risk[c]*(exp(EM.pred.LRR[k, c]) - 1) # RD
                       }}\n")
  } else if (model == "FE" & measure == "RR") {
    paste(stringcode, "for (t in 1:(ref - 1)) {
                         EM.ref.n[t] <- (d[t] + beta.t[t] * cov_value) - log(1 - (1 - exp(d[t] + beta.t[t] * cov_value))*base_risk)
                         EM.ref[t] <- EM.ref.n[t]*equals(min(t, ref), ref) + EM.ref.n[t]*(-1)*equals(min(t, ref), t)
                       }
                       for (t in (ref + 1):nt) {
                         EM.ref.n[t] <- (d[t] + beta.t[t] * cov_value) - log(1 - (1 - exp(d[t] + beta.t[t] * cov_value))*base_risk)
                         EM.ref[t] <- EM.ref.n[t]*equals(min(t, ref), ref) + EM.ref.n[t]*(-1)*equals(min(t, ref), t)
                       for (c in 1:(nt - 1)) {
                         for (k in (c + 1):nt) {
                           EM.LOR[k, c] <- d.n[k] - d.n[c]
                           EM[k, c] <- EM.LOR[k, c] - log(1 - abs_risk[c]*(1 - exp(EM.LOR[k, c]))) # LRR
                        }}\n")
  } else if (model == "FE" & measure == "RD") {
    paste(stringcode, "for (t in 1:(ref - 1)) {
                         EM.ref.RR.n <- (d[t] + beta.t[t] * cov_value) - log(1 - (1 - exp(d[t] + beta.t[t] * cov_value))*base_risk)
                         EM.ref.RR[t] <- EM.ref.RR.n[t]*equals(min(t, ref), ref) + EM.ref.RR.n[t]*(-1)*equals(min(t, ref), t)
                         EM.ref.n[t] <- (exp(EM.ref.RR.n[t] - 1))*base_risk
                         EM.ref[t] <-  EM.ref.n[t]*equals(min(t, ref), ref) + EM.ref.n[t]*(-1)*equals(min(t, ref), t)
                       }
                       for (t in (ref + 1):nt) {
                         EM.ref.RR.n[t] <- (d[t] + beta.t[t] * cov_value) - log(1 - (1 - exp(d[t] + beta.t[t] * cov_value))*base_risk)
                         EM.ref.RR[t] <- EM.ref.RR.n[t]*equals(min(t, ref), ref) + EM.ref.RR.n[t]*(-1)*(equals(min(t, ref), t)
                         EM.ref.n[t] <- (exp(EM.ref.RR.n[t]) - 1)*base_risk
                         EM.ref[t] <-  EM.ref.n[t]*equals(min(t, ref), ref) + EM.ref.n[t]*(-1)*equals(min(t, ref), t)
                       }
                       for (c in 1:(nt - 1)) {
                         for (k in (c + 1):nt) {
                           EM.LOR[k, c] <- d.n[k] - d.n[c]
                           EM.LRR[k, c] <- EM.LOR[k, c] - log(1 - abs_risk[c]*(1 - exp(EM.LOR[k, c])))
                           EM[k, c] <- abs_risk[c]*(exp(EM.LRR[k, c]) - 1)
                        }}\n")
  }

  stringcode <- if (model == "RE") {
    paste(stringcode, "prec <- pow(tau, -2)
                       tau.a ~ dnorm(0, heter.prior[2])I(0, )
                       tau.b ~ dunif(0, heter.prior[2])
                       tau2.c ~ dlnorm(heter.prior[1], heter.prior[2])
                       log.tau2.d ~ dt(heter.prior[1], heter.prior[2], 5)
                       tau <- tau.a*equals(heter.prior[3], 1) + tau.b*equals(heter.prior[3], 2) +
                       pow(tau2, 0.5)*equals(heter.prior[3], 3) + pow(tau2, 0.5)*equals(heter.prior[3], 4)
                       tau2 <- tau2.c*equals(heter.prior[3], 3) + exp(log.tau2.d)*equals(heter.prior[3], 4) \n")
  } else {
    paste(stringcode, " ")
  }

  stringcode <- paste(stringcode, "\n}")

  return(stringcode)
}

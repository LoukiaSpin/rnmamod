#' WinBUGS code for the node-splitting approach
#'
#' @description The WinBUGS code, as written by Dias et al. (2010) to run a
#'   one-stage Bayesian node-splitting model, extended to incorporate the
#'   pattern-mixture model for binary or continuous missing participant outcome
#'   data.
#'
#' @param measure Character string indicating the effect measure with values
#'   \code{"OR"}, \code{"MD"}, \code{"SMD"}, or \code{"ROM"} for the odds ratio,
#'   mean difference, standardised mean difference and ratio of means,
#'   respectively.
#' @param model Character string indicating the analysis model with values
#'   \code{"RE"}, or \code{"FE"} for the random-effects and fixed-effect model,
#'   respectively. The default argument is \code{"RE"}.
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
#' @return An R character vector object to be passed to
#'   \code{\link{run_nodesplit}} through the
#'   \code{\link[base:textConnection]{textConnection}} function as the argument
#'   \code{object}.
#'
#' @details This functions creates the model in the JAGS dialect of the BUGS
#'   language. The output of this function constitutes the argument
#'   \code{model.file} of \code{\link[R2jags:jags]{jags}} (in the R-package
#'   \href{https://CRAN.R-project.org/package=R2jags}{R2jags}) via the
#'   \code{\link[base:textConnection]{textConnection}} function.
#'
#'   The split nodes have been automatically selected via the
#'   \code{\link[gemtc:mtc.nodesplit.comparisons]{mtc.nodesplit.comparisons}}
#'   function of the R-package
#'   \href{https://CRAN.R-project.org/package=gemtc}{gemtc}.
#'   See 'Details' in \code{\link{run_nodesplit}}.
#'
#' @author {Loukia M. Spineli}
#'
#' @seealso \code{\link[R2jags:jags]{jags}},
#'   \code{\link[gemtc:mtc.nodesplit.comparisons]{mtc.nodesplit.comparisons}},
#'   \code{\link{run_nodesplit}},
#'   \code{\link[base:textConnection]{textConnection}}
#'
#' @references
#' Dias S, Welton NJ, Caldwell DM, Ades AE. Checking consistency in mixed
#' treatment comparison meta-analysis.
#' \emph{Stat Med} 2010;\bold{29}(7-8):932--44. \doi{10.1002/sim.3767}
#'
#' Spineli LM, Kalyvas C, Papadimitropoulou K. Continuous(ly) missing outcome
#' data in network meta-analysis: a one-stage pattern-mixture model approach.
#' \emph{Stat Methods Med Res} 2021;\bold{30}(4):958--975.
#' \doi{10.1177/0962280220983544}
#'
#' Spineli LM. An empirical comparison of Bayesian modelling strategies for
#' missing binary outcome data in network meta-analysis.
#' \emph{BMC Med Res Methodol} 2019;\bold{19}(1):86.
#' \doi{10.1186/s12874-019-0731-y}
#'
#' @export
prepare_nodesplit <- function(measure, model, assumption) {

  stringcode <- "model {
                    for (i in 1:ns) {\n"

  stringcode <- if (model == "RE") {
    paste(stringcode, "delta[i, 1] <- 0
                       j[i, 1] <- 0
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
  } else if (measure == "MD" || measure == "ROM") {
    paste(stringcode, "theta[i, 1] <- u[i]\n")
  } else if (measure == "OR") {
    paste(stringcode, "logit(p[i, 1]) <- u[i]\n")
  }

  stringcode <- paste(stringcode, "for (k in 1:na[i]) {\n")

  stringcode <- if (measure == "SMD") {
    paste(stringcode, "prec.o[i, k] <- pow(se.o[i, k], -2)
                       y.o[i, k] ~ dnorm(theta.o[i, k], prec.o[i, k])
                       c[i, k] <- N[i, k] - mod[i, k]
                       sd.obs[i, k] <- se.o[i, k]*sqrt(c[i, k])
                       nom[i, k] <- pow(sd.obs[i, k], 2)*(c[i, k] - 1)\n")
  } else if (measure == "MD" || measure == "ROM") {
    paste(stringcode, "prec.o[i, k] <- pow(se.o[i, k], -2)
                       y.o[i, k] ~ dnorm(theta.o[i, k], prec.o[i, k])\n")
  } else if (measure == "OR") {
    paste(stringcode, "r[i, k] ~ dbin(p_o[i, k], obs[i, k])
                       obs[i, k] <- N[i, k] - mod[i, k]\n")
  }

  stringcode <- if (measure == "MD" || measure == "SMD") {
    paste(stringcode, "theta.o[i, k] <- theta[i, k] - phi.m[i, k]*q[i, k]\n")
  } else if (measure == "ROM") {
    paste(stringcode, "theta.o[i, k] <- theta[i, k]/(1 - q[i, k]*(1 - exp(phi.m[i, k])))\n")
  } else if (measure == "OR") {
    paste(stringcode, "p_o[i, k] <- max(0, min(1, ((-((q[i, k] - p[i, k])*(1 - exp(phi.m[i, k])) - 1) - sqrt((pow(((q[i, k] - p[i, k])*(1 - exp(phi.m[i, k])) - 1), 2)) -
                                       ((4*p[i, k])*(1 - q[i, k])*(1 - exp(phi.m[i, k])))))/(2*(1 - q[i, k])*(1 - exp(phi.m[i, k]))))))\n")
  }

  stringcode <- paste(stringcode, "q[i, k] <- q0[i, k]*I[i, k]
                                   m[i, k] ~ dbin(q0[i, k], N[i, k])
                                   q0[i, k] ~ dunif(0, 1)\n")

  stringcode <- if (measure == "MD" || measure == "SMD" || measure == "ROM") {
    paste(stringcode, "hat.par[i, k] <- theta.o[i, k]
                       dev.o[i, k] <- (y.o[i, k] - theta.o[i, k])*(y.o[i, k] - theta.o[i, k])*prec.o[i, k]\n")
  } else if (measure == "OR") {
    paste(stringcode, "hat.par[i, k] <- rhat[i, k]
                       rhat[i, k] <- p_o[i, k]*obs[i, k]
                       dev.o[i, k] <- 2*(r[i, k]*(log(r[i, k]) - log(rhat[i, k])) + (obs[i, k] - r[i, k])*(log(obs[i, k] - r[i, k]) - log(obs[i, k] - rhat[i, k])))\n")
  }

  stringcode <- paste(stringcode, "index[i, k] <- split[i]*(equals(t[i, k], pair[1]) + equals(t[i, k], pair[2]))
                                   }
                                   resdev.o[i] <- sum(dev.o[i, 1:na[i]])\n")

  stringcode <- paste(stringcode, "for (k in 2:na[i]) {\n")

  stringcode <- if (measure == "MD") {
    paste(stringcode, "theta[i, k] <- u[i] + delta[i, k]\n")
  } else if (measure == "SMD") {
    paste(stringcode, "theta[i, k] <- u[i] + sd.pooled[i]*delta[i, k]\n")
  } else if (measure == "ROM") {
    paste(stringcode, "theta[i, k] <- u[i]*exp(delta[i, k])\n")
  } else if (measure == "OR") {
    paste(stringcode, "logit(p[i, k]) <- u[i] + delta[i, k]\n")
  }

  stringcode <- if (model == "RE") {
    paste(stringcode, "delta[i, k] ~ dnorm(md[i, k], precd[i, k])
                       md[i, k] <- ((d[si[i, k]]*indic[i, k] - d[bi[i]]*indic[i, 1])*I.sign[i, k] + sw[i, k])*(1 - index[i, m[i, k]]) + direct*index[i, m[i, k]]
                       j[i, k] <- k - (equals(1, split[i])*step(k - 3))
                       precd[i, k] <- prec*2*(j[i, k] - 1)/j[i, k]
                       w[i, k] <- ((delta[i, k] - (d[si[i, k]]*indic[i, k] - d[bi[i]]*indic[i, 1]))*I.sign[i, k])*(1 - index[i, k])
                       sw[i, k] <- sum(w[i, 1:(k - 1)])/(j[i, k] - 1)
                       }}\n")
  } else {
    paste(stringcode, "delta[i, k] <- ((d[si[i, k]]*indic[i, k] - d[bi[i]]*indic[i, 1])*I.sign[i, k])*(1 - index[i, m[i, k]]) + direct*index[i, m[i, k]]
                       }}\n")
  }

  stringcode <- paste(stringcode, "totresdev.o <- sum(resdev.o[])
                                   d[ref] <- 0
                                   d.n[ref] <- 0
                                   for (t in 1:(ref - 1)) {
                                     d[t] ~ dnorm(0, 0.0001)
                                     d.n[t] <- d[t]*(1 - (1 - step(t - ref))) + d[t]*(-1)*(1 - step(t - ref))
                                   }
                                   for (t in (ref + 1):nt) {
                                     d[t] ~ dnorm(0, 0.0001)
                                     d.n[t] <- d[t]*(1 - (1 - step(t - ref))) + d[t]*(-1)*(1 - step(t - ref))
                                   }\n")

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

  stringcode <- paste(stringcode, "for (c in 1:(nt - 1)) {
                                     for (k in (c + 1):nt) {
                                       EM[k, c] <- d.n[k] - d.n[c]
                                   }}
                                   direct ~ dnorm(0, .0001)
                                   diff <- direct - EM[pair[2], pair[1]]\n")

  stringcode <- if (model == "RE") {
    paste(stringcode, "prec <- pow(tau, -2)
                       tau.a ~ dnorm(0, heter.prior[2])I(0, )
                       tau.b ~ dunif(0, heter.prior[2])
                       tau2.c ~ dlnorm(heter.prior[1], heter.prior[2])
                       log.tau2.d ~ dt(heter.prior[1], heter.prior[2], 5)
                       tau <- tau.a*equals(heter.prior[3], 1) + tau.b*equals(heter.prior[3], 2) +
                       pow(tau2, 0.5)*equals(heter.prior[3], 3) + pow(tau2, 0.5)*equals(heter.prior[3], 4)
                       tau2 <- tau2.c*equals(heter.prior[3], 3) + exp(log.tau2.d)*equals(heter.prior[3], 4)\n")
  } else {
    paste(stringcode, " ")
  }

  stringcode <- paste(stringcode, "\n}")

  return(stringcode)
}

#' WinBUGS code for the unrelated mean effects model
#'
#' @description The WinBUGS code, as written by Dias et al., (2013) to run a
#'   one-stage Bayesian unrelated mean effects model, extended to incorporate
#'   the pattern-mixture model for binary or continuous missing participant
#'   outcome data.
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
#' @param connected An integer equal to one or larger that indicates the number
#'   of subnetworks.
#'
#' @return An R character vector object to be passed to \code{\link{run_ume}}
#'   through the \code{\link[base]{textConnection}} function as the argument
#'   \code{object}.
#'
#' @details This functions creates the model in the JAGS dialect of the BUGS
#'   language. The output of this function constitutes the argument
#'   \code{model.file} of \code{\link[R2jags]{jags}} via the
#'   \code{\link[base]{textConnection}} function.
#'
#' @author {Loukia M. Spineli}
#'
#' @seealso \code{\link[R2jags]{jags}}, \code{\link{run_ume}},
#'   \code{\link[base]{textConnection}}
#'
#' @references
#' Spineli LM, Kalyvas C, Papadimitropoulou K. Continuous(ly) missing outcome
#' data in network meta-analysis: a one-stage pattern-mixture model approach.
#' \emph{Stat Methods Med Res} 2021. \doi{10.1177/0962280220983544}
#'
#' Spineli LM. An empirical comparison of Bayesian modelling strategies for
#' missing binary outcome data in network meta-analysis.
#' \emph{BMC Med Res Methodol} 2019;\bold{19}(1):86.
#' \doi{10.1186/s12874-019-0731-y}
#'
#' Dias S, Welton NJ, Sutton AJ, Caldwell DM, Lu G, Ades AE. Evidence synthesis
#' for decision making 4: inconsistency in networks of evidence based on
#' randomized controlled trials.
#' \emph{Med Decis Making} 2013;\bold{33}(5):641--56.
#' \doi{10.1177/0272989X12455847}
#'
#' @export
prepare_ume <- function(measure, model, assumption, connected) {

  code <- paste0("model\n{",
                 "\n\tfor (i in 1:ns) {")

  code <- if (model == "RE") {
    paste0(code, "\n\t\tdelta[i, 1] <- 0",
                 "\n\t\tw[i, 1] <- 0",
                 "\n\t\tu[i] ~ dnorm(0, .0001)")
  } else {
    paste0(code, "\n\t\tu[i] ~ dnorm(0, .0001)")
  }

  code <- if (measure == "SMD") {
    paste0(code, "\n\t\ttheta[i, 1] <- u[i]",
                 "\n\t\tsigma[i] <- sqrt(sum(nom[i, 1:na[i]])/(sum(c[i, 1:na[i]]) - na[i]))",
                 "\n\t\ta[i] <- sum(N[i, 1:na[i]] - 1)/2",
                 "\n\t\tb[i] <- sum(N[i, 1:na[i]] - 1)/(2*sigma[i]*sigma[i])",
                 "\n\t\tvar.pooled[i] ~ dgamma(a[i], b[i])",
                "\n\t\tsd.pooled[i] <- sqrt(var.pooled[i])")
  } else if (measure == "MD" || measure == "ROM") {
    paste0(code, "\n\t\ttheta[i, 1] <- u[i]")
  } else if (measure == "OR") {
    paste0(code, "\n\t\tlogit(p[i, 1]) <- u[i]")
  }

  code <- paste0(code, "\n\t\tfor (k in 1:na[i]) {")

  code <- if (measure == "SMD") {
    paste0(code, "\n\t\t\tprec.o[i, k] <- pow(se.o[i, k], -2)",
                 "\n\t\t\ty.o[i, k] ~ dnorm(theta.o[i, k], prec.o[i, k])",
                 "\n\t\t\tc[i, k] <- N[i, k] - m[i, k]",
                 "\n\t\t\tsd.obs[i, k] <- se.o[i, k]*sqrt(c[i, k])",
                 "\n\t\t\tnom[i, k] <- pow(sd.obs[i, k], 2)*(c[i, k] - 1)")
  } else if (measure == "MD" || measure == "ROM") {
    paste0(code, "\n\t\t\tprec.o[i, k] <- pow(se.o[i, k], -2)",
                 "\n\t\t\ty.o[i, k] ~ dnorm(theta.o[i, k], prec.o[i, k])")
  } else if (measure == "OR") {
    paste0(code, "\n\t\t\tr[i, k] ~ dbin(p_o[i, k], obs[i, k])",
                 "\n\t\t\tobs[i, k] <- N[i, k] - m[i, k]")
  }

  code <- if (measure == "MD" || measure == "SMD") {
    paste0(code, "\n\t\t\ttheta.o[i, k] <- theta[i, k] - phi.m[i, k]*q[i, k]")
  } else if (measure == "ROM") {
    paste0(code, "\n\t\t\ttheta.o[i, k] <- theta[i, k]/(1 - q[i, k]*(1 - exp(phi.m[i, k])))")
  } else if (measure == "OR") {
    paste0(code, "\n\t\t\tp_o[i, k] <- max(0, min(1, ((-((q[i, k] - p[i, k])*(1 - exp(phi.m[i, k])) - 1) - sqrt((pow(((q[i, k] - p[i, k])*(1 - exp(phi.m[i, k])) - 1), 2)) -
                                       ((4*p[i, k])*(1 - q[i, k])*(1 - exp(phi.m[i, k])))))/(2*(1 - q[i, k])*(1 - exp(phi.m[i, k]))))))")
  }

  code <- paste0(code, "\n\t\t\tq[i, k] <- q0[i, k]*I[i, k]",
                       "\n\t\t\tm[i, k] ~ dbin(q0[i, k], N[i, k])",
                       "\n\t\t\tq0[i, k] ~ dunif(0, 1)")

  code <- if (measure == "MD" || measure == "SMD" || measure == "ROM") {
    paste0(code, "\n\t\t\that.par[i, k] <- theta.o[i, k]",
                 "\n\t\t\tdev.o[i, k] <- (y.o[i, k] - theta.o[i, k])*(y.o[i, k] - theta.o[i, k])*prec.o[i, k]")
  } else if (measure == "OR") {
    paste0(code, "\n\t\t\that.par[i, k] <- rhat[i, k]",
                 "\n\t\t\trhat[i, k] <- p_o[i, k]*obs[i, k]",
                 "\n\t\t\tdev.o[i, k] <- 2*(r[i, k]*(log(r[i, k]) - log(rhat[i, k])) + (obs[i, k] - r[i, k])*(log(obs[i, k] - r[i, k]) - log(obs[i, k] - rhat[i, k])))")
  }

  code <- paste0(code, "\n\t\t\t}",
                       "\n\t\tresdev.o[i] <- sum(dev.o[i, 1:na[i]])",
                       "\n\t\tfor (k in 2:na[i]) {")

  code <- if (measure == "MD") {
    paste0(code, "\n\t\t\ttheta[i, k] <- u[i] + delta[i, k]")
  } else if (measure == "SMD") {
    paste0(code, "\n\t\t\ttheta[i, k] <- u[i] + sd.pooled[i]*delta[i, k]")
  } else if (measure == "ROM") {
    paste0(code, "\n\t\t\ttheta[i, k] <- u[i]*exp(delta[i, k])")
  } else if (measure == "OR") {
    paste0(code, "\n\t\t\tlogit(p[i, k]) <- u[i] + delta[i, k]")
  }

  code <- if (model == "RE") {
    paste0(code, "\n\t\t\tdelta[i, k] ~ dnorm(md[i, k], precd[i, k])",
                 "\n\t\t\tmd[i, k] <- EM[t[i, k], t[i, 1]] + sw[i, k]",
                 "\n\t\t\tprecd[i, k] <- 2*prec*(k - 1)/k",
                 "\n\t\t\tw[i, k] <- delta[i, k] - EM[t[i, k], t[i, 1]]",
                 "\n\t\t\tsw[i, k] <- sum(w[i, 1:(k - 1)])/(k - 1)")
  } else {
    paste0(code, "\n\t\t\tdelta[i, k] <- EM[t[i, k], t[i, 1]]")
  }

  code <- paste0(code, "\n\t\t\t}}",
                       "\n\ttotresdev.o <- sum(resdev.o[])",
                       "\n\tfor (i in (ns - ns.multi + 1):ns) {")

  code <- if (model == "RE") {
    paste0(code, "\n\t\tdelta.m[i, 1] <- 0",
                 "\n\t\tw.m[i, 1] <- 0",
                 "\n\t\tu.m[i] ~ dnorm(0, .0001)")
  } else {
    paste0(code, "\n\t\tu.m[i] ~ dnorm(0, .0001)")
  }

  code <- if (measure == "MD" || measure == "SMD" || measure == "ROM") {
    paste0(code, "\n\t\ttheta.m[i, 1] <- u.m[i]")
  } else if (measure == "OR") {
    paste0(code, "\n\t\tlogit(p.m[i, 1]) <- u.m[i]")
  }

  code <- paste0(code, "\n\t\tfor (k in 1:na[i]) {")

  code <- if (measure == "MD" || measure == "SMD" || measure == "ROM") {
    paste0(code, "\n\t\t\ty.m[i, k] ~ dnorm(theta.o.m[i, k], prec.o[i, k])")
  } else if (measure == "OR") {
    paste0(code, "\n\t\t\tr.m[i, k] ~ dbin(p_o.m[i, k], obs[i, k])")
  }

  code <- if (measure == "MD" || measure == "SMD") {
    paste0(code, "\n\t\t\ttheta.o.m[i, k] <- theta.m[i, k] - phi.m[i, k]*q[i, k]")
  } else if (measure == "ROM") {
    paste0(code, "\n\t\t\ttheta.o.m[i, k] <- theta.m[i, k]/(1 - q[i, k]*(1 - exp(phi.m[i, k])))")
  } else if (measure == "OR") {
    paste0(code, "\n\t\t\tp_o.m[i, k] <- max(0, min(1, ((-((q[i, k] - p.m[i, k])*(1 - exp(phi.m[i, k])) - 1) - sqrt((pow(((q[i, k] - p.m[i, k])*(1 - exp(phi.m[i, k])) - 1), 2)) -
                                         ((4*p.m[i, k])*(1 - q[i, k])*(1 - exp(phi.m[i, k])))))/(2*(1 - q[i, k])*(1 - exp(phi.m[i, k]))))))")
  }

  code <- paste0(code, "\n\t\t\t}",
                       "\n\t\tfor (k in 2:na[i]) {")

  code <- if (measure == "MD") {
    paste0(code, "\n\t\t\ttheta.m[i, k] <- u.m[i] + delta.m[i, k]")
  } else if (measure == "SMD") {
    paste0(code, "\n\t\t\ttheta.m[i, k] <- u.m[i] + sd.pooled[i]*delta.m[i, k]")
  } else if (measure == "ROM") {
    paste0(code, "\n\t\t\ttheta.m[i, k] <- u.m[i]*exp(delta.m[i, k])")
  } else if (measure == "OR") {
    paste0(code, "\n\t\t\tlogit(p.m[i, k]) <- u.m[i] + delta.m[i, k]")
  }

  code <- if (model == "RE" & connected < 2) {
    paste0(code, "\n\t\t\tdelta.m[i, k] ~ dnorm(md.m[i, k], precd.m[i, k])",
                 "\n\t\t\tprecd.m[i, k] <- 2*(k - 1)*precm/k",
                 "\n\t\t\tsw.m[i, k] <- sum(w.m[i, 1:(k - 1)])/(k - 1)",
                 "\n\t\t\tmd.m[i, k] <- (d.m[t[i, k]] - d.m[t[i, 1]]) + sw.m[i, k]",
                 "\n\t\t\tw.m[i, k] <- delta.m[i, k] - (d.m[t[i, k]] - d.m[t[i, 1]])",
                 "\n\t\t\t}}",
                 "\n\td.m[ref.base] <- 0",
                 "\n\tfor (i in 1:N.t.m) {",
                 "\n\t\td.m[t.m[i]] ~ dnorm(0, .0001)",
                 "\n\t\t}",
                 "\n\tfor (i in 1:nbase.multi) {",
                 "\n\t\tEM[t2.bn[i], t1.bn[i]] <- d.m[t2.bn[i]] - d.m[t1.bn[i]]",
                 "\n\t\t}")
  } else if (model == "RE" & connected >= 2) {
    paste0(code, "\n\t\t\tdelta.m[i, k] ~ dnorm(md.m[i, k], precd.m[i, k])",
                 "\n\t\t\tprecd.m[i, k] <- 2*(k - 1)*precm/k",
                 "\n\t\t\tsw.m[i, k] <- sum(w.m[i, 1:(k - 1)])/(k - 1)",
                 "\n\t\t\tmd.m[i, k] <- (d.multi[t[i, k], ref.m[i]] - d.multi[t[i, 1], ref.m[i]]) + sw.m[i, k]",
                 "\n\t\t\tw.m[i, k] <- delta.m[i, k] - (d.multi[t[i, k], ref.m[i]] - d.multi[t[i, 1], ref.m[i]])",
                 "\n\t\t\t}}",
                 "\n\tfor (i in 1:N.t.m2) {",
                 "\n\t\td.multi[t.m2[i, 1], t.m2[i, 2]] ~ dnorm(0, .0001)",
                 "\n\t\td.m[t.m2[i, 1], t.m2[i, 2]] <- d.multi[t.m2[i, 1], t.m2[i, 2]]*(1 - equals(t.m2[i, 1], t.m2[i, 2]))",
                 "\n\t\t}",
                 "\n\tfor (i in 1:nbase.multi) {",
                 "\n\t\tEM[t2.bn[i], t1.bn[i]] <- d.m[t2.bn[i], ref.nbase.multi[i]] - d.m[t1.bn[i], ref.nbase.multi[i]]",
                 "\n\t\t}")
  } else if (model == "FE" & connected < 2) {
    paste0(code, "\n\t\t\tdelta.m[i, k] <- d.m[t[i, k]] - d.m[t[i, 1]]",
                 "\n\t\t\t}}",
                 "\n\td.m[ref.base] <- 0",
                 "\n\tfor (i in 1:N.t.m) {",
                 "\n\t\td.m[t.m[i]] ~ dnorm(0, .0001)",
                 "\n\t\t}",
                 "\n\tfor (i in 1:nbase.multi) {",
                 "\n\t\tEM[t2.bn[i], t1.bn[i]] <- d.m[t2.bn[i]] - d.m[t1.bn[i]]",
                 "\n\t\t}")
  } else if (model == "FE" & connected >= 2) {
    paste0(code, "\n\t\t\tdelta.m[i, k] <- d.multi[t[i, k], ref.m[i]] - d.multi[t[i, 1], ref.m[i]]",
                 "\n\t\t\t}}",
                 "\n\tfor (i in 1:N.t.m2) {",
                 "\n\t\td.multi[t.m2[i, 1], t.m2[i, 2]] ~ dnorm(0, .0001)",
                 "\n\t\td.m[t.m2[i, 1], t.m2[i, 2]] <- d.multi[t.m2[i, 1], t.m2[i, 2]]*(1 - equals(t.m2[i, 1], t.m2[i, 2]))",
                 "\n\t\t}",
                 "\n\tfor (i in 1:nbase.multi) {",
                 "\n\t\tEM[t2.bn[i], t1.bn[i]] <- d.m[t2.bn[i], ref.nbase.multi[i]] - d.m[t1.bn[i], ref.nbase.multi[i]]",
                 "\n\t\t}")
  }

  code <- paste0(code, "\n\tfor (i in 1:N.obs) {",
                       "\n\t\tEM[t2[i], t1[i]] ~ dnorm(0, .0001)",
                       "\n\t\t}")

  if (assumption == "HIE-ARM") {
    code <- paste0(code, "\n\tfor (i in 1:ns) {",
                         "\n\t\tfor (k in 1:na[i]) {",
                         "\n\t\t\tphi.m[i, k] <- phi[i, k]",
                         "\n\t\t\tphi[i, k] ~ dnorm(mean.phi[t[i, k]], prec.phi[t[i, k]])",
                         "\n\t\t\t}}",
                         "\n\tmean.phi[ref] ~ dnorm(meand.phi[2], precd.phi)",
                         "\n\tprec.phi[ref] <- pow(sd.phi[ref], -2)",
                         "\n\tsd.phi[ref] ~ dunif(0, psi.phi)",
                         "\n\tfor (t in 1:(ref - 1)) {",
                         "\n\t\tmean.phi[t] ~ dnorm(meand.phi[1], precd.phi)",
                         "\n\t\tprec.phi[t] <- pow(sd.phi[t], -2)",
                         "\n\t\tsd.phi[t] ~ dunif(0, psi.phi)",
                         "\n\t\t}",
                         "\n\tfor (t in (ref + 1):nt) {",
                         "\n\t\tmean.phi[t] ~ dnorm(meand.phi[1], precd.phi)",
                         "\n\t\tprec.phi[t] <- pow(sd.phi[t], -2)",
                         "\n\t\tsd.phi[t] ~ dunif(0, psi.phi)",
                         "\n\t\t}",
                         "\n\tpsi.phi <- pow(precd.phi, -2)")
  } else if (assumption == "HIE-TRIAL") {
    code <- paste0(code, "\n\tfor (i in 1:ns) {",
                         "\n\t\tfor (k in 1:na[i]) {",
                         "\n\t\t\tphi.m[i, k] <- phi[i, k]",
                         "\n\t\t\tphi[i, k] ~ dnorm(mean.phi[i], prec.phi[i])",
                         "\n\t\t\t}}",
                         "\n\tfor (i in 1:ns) {",
                         "\n\t\tmean.phi[i] ~ dnorm(meand.phi, precd.phi)",
                         "\n\t\tprec.phi[i] <- pow(sd.phi[i], -2)",
                         "\n\t\tsd.phi[i] ~ dunif(0, psi.phi)",
                         "\n\t\t}",
                         "\n\tpsi.phi <- pow(precd.phi, -2)")
  } else if (assumption == "HIE-COMMON") {
    code <- paste0(code, "\n\tfor (i in 1:ns) {",
                         "\n\t\tfor (k in 1:na[i]) {",
                         "\n\t\t\tphi.m[i, k] <- phi[i, k]",
                         "\n\t\t\tphi[i, k] ~ dnorm(mean.phi, prec.phi)",
                         "\n\t\t\t}}",
                         "\n\tmean.phi ~ dnorm(meand.phi, precd.phi)",
                         "\n\tprec.phi <- pow(sd.phi, -2)",
                         "\n\tsd.phi ~ dunif(0, psi.phi)",
                         "\n\tpsi.phi <- pow(precd.phi, -2)")
  } else if (assumption == "IDE-ARM") {
    code <- paste0(code, "\n\tfor (i in 1:ns) {",
                         "\n\t\tfor (k in 1:na[i]) {",
                         "\n\t\t\tphi.m[i, k] <- phi[t[i, k]]",
                         "\n\t\t\t}}",
                         "\n\tphi[ref] ~ dnorm(meand.phi[2], precd.phi)",
                         "\n\tfor (t in 1:(ref - 1)) {",
                         "\n\t\tphi[t] ~ dnorm(meand.phi[1], precd.phi)",
                         "\n\t\t}",
                         "\n\tfor (t in (ref + 1):nt) {",
                         "\n\t\tphi[t] ~ dnorm(meand.phi[1], precd.phi)",
                         "\n\t\t}")
  } else if (assumption == "IDE-TRIAL") {
    code <- paste0(code, "\n\tfor (i in 1:ns) {",
                         "\n\t\tfor (k in 1:na[i]) {",
                         "\n\t\t\tphi.m[i, k] <- phi[i]",
                         "\n\t\t\t}}",
                         "\n\tfor (i in 1:ns) {",
                         "\n\t\tphi[i] ~ dnorm(meand.phi, precd.phi)",
                         "\n\t\t}")
  } else if (assumption == "IDE-COMMON") {
    code <- paste0(code, "\n\tfor (i in 1:ns) {",
                         "\n\t\tfor (k in 1:na[i]) {",
                         "\n\t\t\tphi.m[i, k] <- phi",
                         "\n\t\t\t}}",
                         "\n\t\tphi ~ dnorm(meand.phi, precd.phi)")
  } else if (assumption == "IND-CORR") {
    code <- paste0(code, "\n\tfor (i in 1:ns) {",
                         "\n\t\tfor (k in 1:na[i]) {",
                         "\n\t\t\tphi.m[i, k] <- phi[i, k]",
                         "\n\t\t\tfor (l in 1:na[i]) {",
                         "\n\t\t\t\tV[i, k, l] <- cov.phi*(1 - equals(k, l)) + var.phi*equals(k, l)",
                         "\n\t\t\t\t}}",
                         "\n\t\tOmega[i, 1:na[i], 1:na[i]] <- inverse(V[i, 1:na[i], 1:na[i]])",
                         "\n\t\tphi[i, 1:na[i]] ~ dmnorm(M[i, 1:na[i]], Omega[i, 1:na[i], 1:na[i]])",
                         "\n\t\t}")
  } else if (assumption == "IND-UNCORR") {
    code <- paste0(code, "\n\tfor (i in 1:ns) {",
                         "\n\t\tfor (k in 1:na[i]) {",
                         "\n\t\t\tphi.m[i, k] <- phi[i, k]",
                         "\n\t\t\tphi[i, k] ~ dnorm(meand.phi, precd.phi)",
                         "\n\t\t\t}}")
  }

  code <- if (model == "RE") {
    paste0(code, "\n\tprec <- pow(tau, -2)",
                 "\n\ttau.a ~ dnorm(0, heter.prior[2])I(0, )",
                 "\n\ttau.b ~ dunif(0, heter.prior[2])",
                 "\n\ttau2.c ~ dlnorm(heter.prior[1], heter.prior[2])",
                 "\n\tlog.tau2.d ~ dt(heter.prior[1], heter.prior[2], 5) ",
                 "\n\ttau <- tau.a*equals(heter.prior[3], 1) + tau.b*equals(heter.prior[3], 2)+
                  pow(tau2, 0.5)*equals(heter.prior[3], 3) + pow(tau2, 0.5)*equals(heter.prior[3], 4)",
                 "\n\ttau2 <- tau2.c*equals(heter.prior[3], 3) + exp(log.tau2.d)*equals(heter.prior[3], 4)",
                 "\n\tprecm <- pow(m.tau, -2)",
                 "\n\tm.tau ~ dnorm(0, 1)I(0, )")
  } else {
    paste0(code, " ")
  }

  code <- paste0(code, "\n}")

  return(code)
}

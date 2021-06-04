prepare.model <- function(measure, model, assumption) {

  code <- paste0("model\n{")

  code <- paste0(code, "\n\tfor (i in 1:ns) {")

  code <- if (model == "RE") {
    paste0(code, "\n\t\tdelta[i, 1] <- 0",
                 "\n\t\tw[i, 1] <- 0",
                 "\n\t\tu[i] ~ dnorm(0, .0001)")
  } else {
    paste0(code, "\n\t\tu[i] ~ dnorm(0, .0001)")
  }

  if (measure == "SMD") {
    code <- paste0(code, "\n\t\ttheta[i, 1] <- u[i]",
                         "\n\t\tsigma[i] <- sqrt(sum(nom[i, 1:na[i]])/(sum(c[i, 1:na[i]]) - na[i]))",
                         "\n\t\ta[i] <- sum(N[i, 1:na[i]] - 1)/2",
                         "\n\t\tb[i] <- sum(N[i, 1:na[i]] - 1)/(2*sigma[i]*sigma[i])",
                         "\n\t\tvar.pooled[i] ~ dgamma(a[i], b[i])",
                         "\n\t\tsd.pooled[i] <- sqrt(var.pooled[i])")
  } else if (measure == "MD" || measure == "ROM") {
    code <- paste0(code, "\n\t\ttheta[i, 1] <- u[i]")
  } else if (measure == "OR") {
    code <- paste0(code, "\n\t\tlogit(p[i, 1]) <- u[i]")
  }

  code <- paste0(code, "\n\t\tfor (k in 1:na[i]) {")

  if (measure == "SMD") {
    code <- paste0(code, "\n\t\t\tprec.o[i, k] <- pow(se.o[i, k], -2)",
                         "\n\t\t\ty.o[i, k] ~ dnorm(theta.o[i, k], prec.o[i, k])",
                         "\n\t\t\tc[i, k] <- N[i, k] - m[i, k]",
                         "\n\t\t\tsd.obs[i, k] <- se.o[i, k]*sqrt(c[i, k])",
                         "\n\t\t\tnom[i, k] <- pow(sd.obs[i, k], 2)*(c[i, k] - 1)")
  } else if (measure == "MD" || measure == "ROM") {
    code <- paste0(code, "\n\t\t\tprec.o[i, k] <- pow(se.o[i, k], -2)",
                         "\n\t\t\ty.o[i, k] ~ dnorm(theta.o[i, k], prec.o[i, k])")
  } else if (measure == "OR") {
    code <- paste0(code, "\n\t\t\tr[i, k] ~ dbin(p_o[i, k], obs[i, k])",
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

  code <- paste0(code, "\n\t\t\tm[i, k] ~ dbin(q[i, k], N[i, k])",
                       "\n\t\t\tq[i, k] ~ dunif(0, 1)",
                       "\n\t\t\tm0[i, k] <- m[i, k] + 0.01*equals(m[i, k], 0)",
                       "\n\t\t\that.m[i, k] <- q[i, k]*N[i, k]",
                       "\n\t\t\tdev.m[i, k] <- 2*(m0[i, k]*(log(m0[i, k]) - log(hat.m[i, k])) +
                        (N[i, k] - m0[i, k])*(log(N[i, k] - m0[i, k]) - log(N[i, k] - hat.m[i, k])))")

  if (measure == "MD" || measure == "SMD" || measure == "ROM") {
    code <- paste0(code, "\n\t\t\that.par[i, k] <- theta.o[i, k]",
                         "\n\t\t\tdev.o[i, k] <- (y.o[i, k] - theta.o[i, k])*(y.o[i, k] - theta.o[i, k])*prec.o[i, k]")
  } else if (measure == "OR") {
    code <- paste0(code, "\n\t\t\that.par[i, k] <- rhat[i, k]",
                         "\n\t\t\trhat[i, k] <- p_o[i, k]*obs[i, k]",
                         "\n\t\t\tdev.o[i, k] <- 2*(r[i, k]*(log(r[i, k]) - log(rhat[i, k])) + (obs[i, k] - r[i, k])*(log(obs[i, k] - r[i, k]) - log(obs[i, k] - rhat[i, k])))")
  }

  code <- paste0(code, "\n\t\t\t}",
                       "\n\t\tresdev.o[i] <- sum(dev.o[i, 1:na[i]])",
                       "\n\t\tresdev.m[i] <- sum(dev.m[i, 1:na[i]])",
                       "\n\t\tfor (k in 2:na[i]) {")

  code <- if (measure == "MD") {
    paste0(code, "\n\t\t\ttheta[i, k] <- u[i] + delta.star[i, k]")
  } else if (measure == "SMD") {
    paste0(code, "\n\t\t\ttheta[i, k] <- u[i] + sd.pooled[i]*delta.star[i, k]")
  } else if (measure == "ROM") {
    paste0(code, "\n\t\t\ttheta[i, k] <- u[i]*exp(delta.star[i, k])")
  } else if (measure == "OR") {
    paste0(code, "\n\t\t\tlogit(p[i, k]) <- u[i] + delta.star[i, k]")
  }

  code <- if (model == "RE") {
    paste0(code, "\n\t\t\tdelta.star[i, k] <- delta[i, k] + (beta[t[i, k]] - beta[t[i, 1]])*eff.mod[i]*(1 - equals(eff.mod[i], 0))",
                 "\n\t\t\tdelta[i, k] ~ dnorm(md[i, k], precd[i, k])",
                 "\n\t\t\tmd[i, k] <- d[t[i, k]] - d[t[i, 1]] + sw[i, k]",
                 "\n\t\t\tprecd[i, k] <- 2*prec*(k - 1)/k",
                 "\n\t\t\tw[i, k] <- delta[i, k] - (d[t[i, k]] - d[t[i, 1]])",
                 "\n\t\t\tsw[i, k] <- sum(w[i, 1:(k - 1)])/(k - 1)")
  } else {
    paste0(code, "\n\t\t\tdelta.star[i, k] <- d[t[i, k]] - d[t[i, 1]] + (beta[t[i, k]] - beta[t[i, 1]])*eff.mod[i]*(1 - equals(eff.mod[i], 0))")
  }

  code <- paste0(code, "\n\t\t\t}}",
                       "\n\ttotresdev.o <- sum(resdev.o[])",
                       "\n\ttotresdev.m <- sum(resdev.m[])",
                       "\n\td[ref] <- 0",
                       "\n\tbeta[ref] <- 0",
                       "\n\tfor (t in 1:(ref - 1)) {",
                       "\n\t\td[t] ~ dnorm(0, 0.0001)",
                       "\n\t\tbeta[t] ~ dnorm(0, 0.0001)",
                       "\n\t\t}",
                       "\n\tfor (t in (ref + 1):nt) {",
                       "\n\t\td[t] ~ dnorm(0, 0.0001)",
                       "\n\t\tbeta[t] ~ dnorm(0, 0.0001)",
                       "\n\t\t}",
                       "\n\tmean.B ~ dnorm(0, .0001)",
                       "\n\tprec.B <- 1/pow(beta.SD,2)",
                       "\n\tbeta.SD ~ dnorm(0, 1)I(0, )")

  code <- if (assumption == "HIE-ARM") {
    paste0(code, "\n\tfor (i in 1:ns) {",
                 "\n\t\tfor(k in 1:na[i]){",
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
    paste0(code, "\n\tfor (i in 1:ns) {",
                 "\n\t\tfor(k in 1:na[i]){",
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
    paste0(code, "\n\tfor (i in 1:ns) {",
                 "\n\t\tfor(k in 1:na[i]){",
                 "\n\t\t\tphi.m[i, k] <- phi[i, k]",
                 "\n\t\t\tphi[i, k] ~ dnorm(mean.phi, prec.phi)",
                 "\n\t\t\t}}",
                 "\n\tmean.phi ~ dnorm(meand.phi, precd.phi)",
                 "\n\tprec.phi <- pow(sd.phi, -2)",
                 "\n\tsd.phi ~ dunif(0, psi.phi)",
                 "\n\tpsi.phi <- pow(precd.phi, -2)")
  } else if (assumption == "IDE-ARM") {
    paste0(code, "\n\tfor (i in 1:ns) {",
                 "\n\t\tfor(k in 1:na[i]){",
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
    paste0(code, "\n\tfor (i in 1:ns) {",
                 "\n\t\tfor(k in 1:na[i]){",
                 "\n\t\t\tphi.m[i, k] <- phi[i]",
                 "\n\t\t\t}}",
                 "\n\tfor (i in 1:ns) {",
                 "\n\t\tphi[i] ~ dnorm(meand.phi, precd.phi)",
                 "\n\t\t}")
  } else if (assumption == "IDE-COMMON") {
    paste0(code, "\n\tfor (i in 1:ns) {",
                 "\n\t\tfor(k in 1:na[i]){",
                 "\n\t\t\tphi.m[i, k] <- phi",
                 "\n\t\t\t}}",
                 "\n\t\tphi ~ dnorm(meand.phi, precd.phi)")
  } else if (assumption == "IND-CORR") {
    paste0(code, "\n\tfor (i in 1:ns) {",
                 "\n\t\tfor (k in 1:na[i]) {",
                 "\n\t\t\tphi.m[i, k] <- phi[i, k]",
                 "\n\t\t\tfor (l in 1:na[i]) {",
                 "\n\t\t\t\tV[i, k, l] <- cov.phi*(1 - equals(k, l)) + var.phi*equals(k, l)",
                 "\n\t\t\t\t}}",
                 "\n\t\tOmega[i, 1:na[i], 1:na[i]] <- inverse(V[i, 1:na[i], 1:na[i]])",
                 "\n\t\tphi[i, 1:na[i]] ~ dmnorm(M[i, 1:na[i]], Omega[i, 1:na[i], 1:na[i]])",
                 "\n\t\t}")
  } else if (assumption == "IND-UNCORR") {
    paste0(code, "\n\tfor (i in 1:ns) {",
                 "\n\t\tfor (k in 1:na[i]) {",
                 "\n\t\t\tphi.m[i, k] <- phi[i, k]",
                 "\n\t\t\tphi[i, k] ~ dnorm(meand.phi, precd.phi)",
                 "\n\t\t\t}}")
  }

  code <- paste0(code, "\n\tsorted <- rank(d[])",
                       "\n\tfor (t in 1:nt) {",
                       "\n\t\torder[t] <- (nt + 1 - sorted[t])*equals(D, 1) + sorted[t]*(1 - equals(D, 1))",
                       "\n\t\tmost.effective[t] <- equals(order[t], 1)",
                       "\n\t\tfor (l in 1:nt) {",
                       "\n\t\t\teffectiveness[t, l] <- equals(order[t], l)",
                       "\n\t\t\tcumeffectiveness[t, l] <- sum(effectiveness[t, 1:l])",
                       "\n\t\t\t}",
                       "\n\t\tSUCRA[t] <- sum(cumeffectiveness[t, 1:(nt - 1)])/(nt - 1)",
                       "\n\t\t}")

  code <- if (model == "RE") {
    paste0(code, "\n\tfor (t in 1:(ref - 1)) {",
                 "\n\t\tEM.ref[t] <- d[t] - d[ref]",
                 "\n\t\tpred.ref[t] ~ dnorm(EM.ref[t], prec)",
                 "\n\t\t}",
                 "\n\tfor (t in (ref + 1):nt) {",
                 "\n\t\tEM.ref[t] <- d[t] - d[ref]",
                 "\n\t\tpred.ref[t] ~ dnorm(EM.ref[t], prec)",
                 "\n\t\t}",
                 "\n\tfor (c in 1:(nt - 1)) {",
                 "\n\t\tfor (k in (c + 1):nt) {",
                 "\n\t\t\tEM[k, c] <- d[k] - d[c]",
                 "\n\t\t\tEM.pred[k, c] ~ dnorm(EM[k, c], prec)")
  } else {
    paste0(code, "\n\tfor (t in 1:(ref - 1)) {",
                 "\n\t\tEM.ref[t] <- d[t] - d[ref]",
                 "\n\t\t}",
                 "\n\tfor (t in (ref + 1):nt) {",
                 "\n\t\tEM.ref[t] <- d[t] - d[ref]",
                 "\n\t\t}",
                 "\n\tfor (c in 1:(nt - 1)) {",
                 "\n\t\tfor (k in (c + 1):nt) {",
                 "\n\t\t\tEM[k, c] <- d[k] - d[c]")
  }

  code <- paste0(code, "\n\t\t\t}}")

  code <- if (model == "RE") {
    paste0(code, "\n\tprec <- pow(tau, -2)",
                 "\n\ttau.a ~ dnorm(0, heter.prior[2])I(0, )",
                 "\n\ttau.b ~ dunif(0, heter.prior[2])",
                "\n\ttau <- tau.a*equals(heter.prior[3], 1) + tau.b*equals(heter.prior[3], 2)")
  } else {
    paste0(code, " ")
  }

  code <- paste0(code, "\n}")

  return(code)
}












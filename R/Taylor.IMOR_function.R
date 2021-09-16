# Function to obtain logORs and SEs using the Taylor-series alongside a pattern-mixture model with IMOR (Strategy 2.1.3)
Taylor.IMOR <- function(data, delta1, delta2, var.delta1, var.delta2, rho){

  for(i in 1:length(data[, 1])){

    # Add 0.5 continuity correction when there is at least on zero cell
    if(data[i, 2] == 0 || data[i, 3] == 0 || data[i, 6] - data[i, 4] - data[i, 2] == 0 || data[i, 7] - data[i, 5] - data[i, 3] == 0){
      data[i, 2] <- data[i, 2] + 0.5
      data[i, 3] <- data[i, 3] + 0.5
      data[i, 6] <- data[i, 6] + 1
      data[i, 7] <- data[i, 7] + 1

    } else {
      data[i, 2] <- data[i, 2]
      data[i, 3] <- data[i, 3]
      data[i, 6] <- data[i, 6]
      data[i, 7] <- data[i, 7]
    }
  }

  # Calculate the probability of event among completers in arm 1 and 2
  p.o1 <- data[, 2]/(data[, 6] - data[, 4]); p.o2 <- data[, 3]/(data[, 7] - data[, 5])

  # Calculate the probability of missing outcome data in arm 1 and 2
  a1 <- data[, 4]/data[, 6]; a2 <- data[, 5]/data[, 7]

  # Calculate the probability of event in randomised sample in arm 1 and 2 via pattern-mixture model
  p.all1 <- (1 - a1)*p.o1 + a1*( (exp(delta1)*p.o1)/(exp(delta1)*p.o1 + 1 - p.o1) ); p.all2 <- (1 - a2)*p.o2 + a2*( (exp(delta2)*p.o2)/(exp(delta2)*p.o2 + 1 - p.o2) )

  # Estimates the odds ratio in the logarithmic scale
  logOR <- log(p.all1/(1 - p.all1)) - log(p.all2/(1 - p.all2))

  #################################
  ## USE OF TAYLOR APPROXIMATION ##
  #################################

  # Derivative of p.all by p.o per arm (first term in Equation 14 in PMID: 17703496)
  A1 <- 1 - a1 + (a1*exp(delta1))/(exp(delta1)*p.o1 + 1 - p.o1)^2; A2 <- 1 - a2 + (a2*exp(delta2))/(exp(delta2)*p.o2 + 1 - p.o2)^2

  # Variance of p.o per arm (second term in Equation 14 in PMID: 17703496)
  B1 <- (p.o1*(1 - p.o1))/(data[, 6] - data[, 4]); B2 <- (p.o2*(1 - p.o2))/(data[, 7] - data[, 5])

  # Derivative of p.all by prob of MOD (i.e. a) per arm (third term in Equation 14 in PMID: 17703496)
  C1 <- (p.o1*(1 - p.o1)*(exp(delta1) - 1))/(exp(delta1)*p.o1 + 1 - p.o1); C2 <- (p.o2*(1 - p.o2)*(exp(delta2) - 1))/(exp(delta2)*p.o2 + 1 - p.o2)

  # Variance of prob of MOD (i.e. a) per arm (fourth term in Equation 14 in PMID: 17703496)
  D1 <- (a1*(1 - a1))/data[, 6]; D2 <- (a2*(1 - a2))/data[, 7]

  # Variance of log odds using delta-method per arm
  E1 <- 1/(p.all1*(1 - p.all1)); E2 <- 1/(p.all2*(1 - p.all2))

  # Derivative of p.all by delta per arm (second Equation after Equation (15) in PMID: 17703496)
  H1 <- (a1*p.o1*(1 - p.o1)*exp(delta1))/(exp(delta1)*p.o1 + 1 - p.o1)^2; H2 <- (a2*p.o2*(1 - p.o2)*exp(delta2))/(exp(delta2)*p.o2 + 1 - p.o2)^2

  # Variance using the observed cases (Equation 13 in PMID: 17703496)
  v.obs <- (A1*A1*B1 + C1*C1*D1)*E1*E1 + (A2*A2*B2 + C2*C2*D2)*E2*E2

  # Variance due to informative missingness (Equation 16 with correlation in PMID: 17703496)
  v.delta <-  H1*H1*var.delta1*E1*E1 + H2*H2*var.delta2*E2*E2 - 2*rho*H1*H2*sqrt(var.delta1)*sqrt(var.delta2)*E1*E2

  # Variance using the randomised sample (Equation 10 in PMID: 17703496)
  v.all <- v.obs + v.delta

  # Include trial-specific adjusted logORs and SEs in the initial dataset
  final <- data.frame(cbind(data, round(logOR, 3),  round(sqrt(v.all), 3)))
  colnames(final) <- c("id", "e1", "e2", "m1", "m2", "n1", "n2", "t1", "t2", "logOR", "SElogOR")

  return(final)
}
## END

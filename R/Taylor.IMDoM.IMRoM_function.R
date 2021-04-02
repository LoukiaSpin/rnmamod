##############################################################################################
#                                                                                            #
#                   Function to obtain MDs and SEs using the Taylor-series                   # 
#               (Pattern-mixture model with IMDoM and IMRoM; PMID: 25393541)                 #
#                                                                                            #
##############################################################################################



Taylor.IMDoM.IMRoM <- function(data, measure, mean.value, var.value, rho){  
  
  
  # Calculate the probability of observing the outcomes in arm 1 (Control) and 2 (Experimental)
  a1 <- data[, 8]/(data[, 8] + data[, 6]); a2 <- data[, 9]/(data[, 9] + data[, 7])
  
  
  # Calculate the adjusted-mean for MOD in the randomised sample of arm 1 and 2
  if(measure == "MD" || measure == "SMD"){
    
    y.all1 <- data[, 2] + mean.value*(1 - a1); y.all2 <- data[, 3] + mean.value*(1 - a2)
    
  } else if(measure == "ROM"){
    
    y.all1 <- data[, 2]*(a1 + exp(mean.value)*(1 - a1)); y.all2 <- data[, 3]*(a2 + exp(mean.value)*(1 - a2))
    
  }

  
  
  # Estimate the adjusted within-trial effect measures (Experimental vs Control)
  if(measure == "MD"){

    MD <- y.all2 - y.all1   
    
  } else if(measure == "SMD"){
    
    nominator1 <- (data[, 8] - 1)*data[, 4]*data[, 4]; nominator2 <- (data[, 9] - 1)*data[, 5]*data[, 5]
    denominator <- (data[, 8] - 1) + (data[, 9] - 1)
    sd.pooled <- sqrt((nominator1 + nominator2)/denominator)
    
    SMD <- (y.all2 - y.all1)/sd.pooled

  } else if(measure == "ROM"){
    
    LROM <- log(y.all2/y.all1)
    
  }

  
  #####################################################################
  ## USE OF TAYLOR APPROXIMATION FOR THE TOTAL WITHIN-TRIAL VARIANCE ##
  #####################################################################
  
  ## Estimating the variance of the mean effect size based on the observed data
  # Derivative of y.all by y.obs per arm 
  if(measure == "MD" || measure == "SMD"){
    
    A1 <- A2 <- 1
    
  } else if(measure == "SMD"){
    
    A1 <- A2 <- 1
    
  } else if(measure == "ROM"){
    
    A1 <-  a1 + exp(mean.value)*(1 - a1);  A2 <-  a2 + exp(mean.value)*(1 - a2)
    
  }

  
  # Variance of y.obs per arm 
  B1 <- (data[, 4]*data[, 4])/data[, 8]; B2 <- (data[, 5]*data[, 5])/data[, 9]   
  

  # Derivative of y.all by prob of MOD (i.e. a) per arm 
  if(measure == "MD" || measure == "SMD"){
    
    C1 <- C2 <- -mean.value
    
  } else if(measure == "ROM"){
    
    C1 <- data[, 2]*(1 - exp(mean.value)); C2 <- data[, 3]*(1 - exp(mean.value))
    
  }
    
  
  # Variance of prob of MOD 
  D1 <- (a1*(1 - a1))/(data[, 8] + data[, 6]); D2 <- (a2*(1 - a2))/(data[, 9] + data[, 7])
  
  
  # Derivative of link function for MD, SMD and LROM per arm
  if(measure == "MD"){
    
    E1 <- E2 <- 1
    
  } else if(measure == "SMD"){
    
    E1 <- E2 <- 1/sd.pooled
    
  } else if(measure == "ROM"){
    
    E1 <- 1/y.all1;  E2 <- 1/y.all2
    
  }

  
  # Variance using the observed cases (Equation 13 in PMID: 17703496)
  v.obs <- (A1*A1*B1 + C1*C1*D1)*E1*E1 + (A2*A2*B2 + C2*C2*D2)*E2*E2
  
  
  ## Estimating the variance of the mean effect size arising from the informative missingness parameter 
  # Derivative of y.all by delta per arm 
  if(measure == "MD" || measure == "SMD"){
    
    H1 <- (1 - a1); H2 <- (1 - a2)
    
  } else if(measure == "ROM"){
    
    H1 <- data[, 2]*exp(mean.value)*(1 - a1); H2 <- data[, 3]*exp(mean.value)*(1 - a2)
    
  }


  # Variance due to informative missingness
  v.delta <-  H1*H1*var.value*E1*E1 + H2*H2*var.value*E2*E2 - 2*rho*H1*H2*E1*E2*sqrt(var.value)*sqrt(var.value)
  
  
  # Variance using the randomised sample 
  v.all <- v.obs + v.delta
  
  
  # Include trial-specific adjusted MDs, SDMs and SEs in the initial dataset
  if(measure == "MD"){
    
    final <- data.frame(cbind(data, round(MD, 3),  round(sqrt(v.all), 3)))
    colnames(final) <- c("id", "mean1", "mean2", "sd1", "sd2", "m1", "m2", "c1", "c2", "t1", "t2", "MD", "SE.MD")
    
  } else if(measure == "SMD"){
    
    final <- data.frame(cbind(data, round(SMD, 3),  round(sqrt(v.all), 3)))
    colnames(final) <- c("id", "mean1", "mean2", "sd1", "sd2", "m1", "m2", "c1", "c2", "t1", "t2", "SMD", "SE.SMD")
    
  } else if(measure == "ROM"){
    
    final <- data.frame(cbind(data, round(LROM, 3),  round(sqrt(v.all), 3)))
    colnames(final) <- c("id", "mean1", "mean2", "sd1", "sd2", "m1", "m2", "c1", "c2", "t1", "t2", "LROM", "SE.LROM")
    
  }

  return(final)
}
## END

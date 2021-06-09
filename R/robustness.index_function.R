#' Robustness index
#'
#' @export
robustness.index <- function(sens, primary.scenar, threshold, nt){


  ES.mat <- sens$EM

  if (missing(threshold) & is.element(sens$measure, "OR")) {
    threshold <- 0.28
    #message("The value 0.28 was assigned on 'threshold' by default")
    message(cat(paste0("\033[0;", col = 32, "m", txt = "The value 0.28 was assigned on 'threshold' by default", "\033[0m", "\n")))
  } else if (missing(threshold) & is.element(sens$measure, c("MD", "SMD", "ROM"))) {
    threshold <- 0.17
    #message("The value 0.17 was assigned on 'threshold' by default")
    message(cat(paste0("\033[0;", col = 32, "m", txt = "The value 0.17 was assigned on 'threshold' by default", "\033[0m", "\n")))
  } else {
    threshold <- threshold
    #message(paste("The value", threshold, "was assigned on 'threshold' for", effect.measure.name(full$measure)))
    message(cat(paste0("\033[0;", col = 32, "m", txt = paste("The value", threshold, "was assigned on 'threshold' for", effect.measure.name(sens$measure)), "\033[0m", "\n")))
  }

  ## Function for the Kullback-Leibler Divergence (comparing two univariate normal distributions)
  KLD.measure.univ <- function(mean.y, sd.y, mean.x, sd.x){

    # x is the 'truth' (e.g. the MAR assumption)
    KLD.xy <- 0.5*(((sd.x/sd.y)^2) + ((mean.y - mean.x)^2)/(sd.y^2) - 1 + 2*log(sd.y/sd.x))

    return(list(KLD.xy = KLD.xy))
  }



  ## Number of scenarios for the sensitivity analysis
  n.scenar <- length(ES.mat[, 1])/(nt*(nt - 1)/2)


  ## A matrix of effect estimates of MCMC standard deviations (or standard errors) for all possible comparisons under each scenario
  sd <- mean <- matrix(NA, nrow = n.scenar, ncol = nt*(nt - 1)/2)

  for(i in 1:n.scenar){

    for(j in 1:(nt*(nt - 1)/2)){

      mean[i, j] <- ES.mat[j + (nt*(nt - 1)/2)*(i - 1), 1]

      sd[i, j] <- ES.mat[j + (nt*(nt - 1)/2)*(i - 1), 2]

    }
  }

  kldxy <- list()

  RI <- nt*(nt - 1)/2

  for(i in 1:(nt*(nt - 1)/2)){ ## We are interested in all possible pairwise comparisons of the network

    kldxy[[i]] <- rep(NA, n.scenar)

    for(j in (1:n.scenar)[-primary.scenar]){


      ## Returns the KLD of informative scenario j when compared with MAR (MAR as 'true') for comparison i
      kldxy[[i]][j] <- KLD.measure.univ(mean[j, i], sd[j, i], mean[primary.scenar, i], sd[primary.scenar, i])[[1]]

    }
    kldxy[[i]][13] <- 0  ## This refers to the primary analysis (here, the MAR assumption)

    ## Returns the Robustness Index of comparison i across all informative scenarios
    RI[i] <- sqrt(round(t(kldxy[[i]][-primary.scenar]) %*% kldxy[[i]][-primary.scenar], 3))

  }

  KLD <- matrix(unlist(kldxy), nrow = (nt*(nt - 1)/2), ncol = n.scenar, byrow = T)


  robust <- ifelse(RI < threshold, "robust", "frail")

  return(list(RI = RI, robust = robust, KLD = KLD, measure = sens$measure, threshold = threshold))
}









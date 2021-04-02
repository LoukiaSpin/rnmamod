#' Preparing the direction variable for the investigation of possible small-study effect
#'
#' @param data
#' @param direction A list refering to the direction of each intervention in the network
#' @param D A binary number for the direction of the outcome. Set \code{D = 1} for a positive outcome and \code{D = 0} for a negative outcome.
#'
#' @export
prepare.small.study.effect <- function(data, measure, var.misspar, direction, D) {

  #data <- data1[, -17]
  #direction <- list(1, "reference")
  #set.seed(123)
  #direction <- list(sample(1:nt), "new-old")
  #var.misspar <- c(1, 0)


  if (measure == "MD" || measure == "SMD"|| measure == "ROM") {

    ## Continuous: arm-level, wide-format dataset
    y.obs <- data %>% dplyr::select(starts_with("y"))                             # Observed mean value in each arm of every trial
    sd.obs <- data %>% dplyr::select(starts_with("sd"))                           # Observed standard deviation in each arm of every trial
    mod <- data %>% dplyr::select(starts_with("m"))                               # Number of missing participants in each arm of every trial
    c0 <- data %>% dplyr::select(starts_with("c"))                                 # Number of completers in each arm of every trial
    rand <- mod + c0                                                               # Number of randomised participants in each arm of every trial
    treat <- data %>% dplyr::select(starts_with("t"))                             # Intervention studied in each arm of every trial
    na <- apply(treat, 1, function(x) length(which(!is.na(x))))                     # Number of interventions investigated in every trial per network
    nt <- length(table(as.matrix(treat)))                                           # Total number of interventions per network
    ns <- length(y.obs[, 1])                                                        # Total number of included trials per network


    ## Order by 'id of t1' < 'id of t1'
    y0 <- sd0 <- m <- c <- N <- t <- t0 <- treat
    for(i in 1:ns){

      t0[i, ] <- order(treat[i, ], na.last = T)
      y0[i, ] <- y.obs[i, order(t0[i, ], na.last = T)]
      sd0[i, ] <- sd.obs[i, order(t0[i, ], na.last = T)]
      m[i, ] <- mod[i, order(t0[i, ], na.last = T)]
      c[i, ] <- c0[i, order(t0[i, ], na.last = T)]
      N[i, ] <- rand[i, order(t0[i, ], na.last = T)]
      t[i, ] <- sort(treat[i, ], na.last = T)
    }


    ## Turn into contrast-level data: one row per possible comparison in each trial ('
    without.MOD <- pairwise(as.list(t), mean = as.list(y0), sd = as.list(sd0), n = as.list(c), data = cbind(t, y0, sd0, c), studlab = 1:ns)[, c(3:5, 7, 10, 8, 11, 6, 9)]
    with.MOD <- pairwise(as.list(t), event = as.list(m), n = as.list(N), data = cbind(t, m, N), studlab = 1:ns)[, c(3:6, 8)]


    ## Prepare the dataset to run the 'Taylor.IMOR' function
    two.stage.dataset <- cbind(without.MOD[, c(1, 4:7)], with.MOD[, 4:5], without.MOD[, c(8:9, 2:3)])


    ## Two-stage pattern-mixture model under MAR
    res.two.stage0 <- Taylor.IMDoM.IMRoM(data = two.stage.dataset, measure = "MD", mean.value = 0, var.value = var.misspar[1], rho = var.misspar[1])[, c(1, 13)]


    ## Indicate the number of possible comparisons in each trial
    arms0 <- arms1 <- list()
    for(i in 1:ns) {
      arms0[[i]] <- 1:dim(combn(na[i], 2))[2]
      arms1[[i]] <- rep(na[i], dim(combn(na[i], 2))[2] )
    }


    ## Obtain the adjusted standard errors per trial in the wide format
    res.two.stage <- as.data.frame(subset(cbind(res.two.stage0, unlist(arms0)), unlist(arms0) <= unlist(arms1) - 1))
    colnames(res.two.stage) <- c("id", "StD.error", "indic")
    SE <- cbind(rep(NA, ns), t(spread(res.two.stage, key = id, value = StD.error))[-1, ])
    Var <- SE*SE
    Prec <- 1/Var


  } else {


    ## Binary: arm-level, wide-format dataset
    event <- data %>% dplyr::select(starts_with("r"))                               # Number of observed events in each arm of every trial
    mod <- data %>% dplyr::select(starts_with("m"))                                 # Number of missing participants in each arm of every trial
    rand <- data %>% dplyr::select(starts_with("n"))                                # Number randomised participants in each arm of every trial
    treat <- data %>% dplyr::select(starts_with("t"))                               # Intervention studied in each arm of every trial
    na <- apply(treat, 1, function(x) length(which(!is.na(x))))                     # Number of interventions investigated in every trial per network
    nt <- length(table(as.matrix(treat)))                                           # Total number of interventions per network
    ns <- length(event[, 1])                                                        # Total number of included trials per network


    ## Order by 'id of t1' < 'id of t1'
    r <- m <- N <- t <- t0 <- treat
    for(i in 1:ns){

      t0[i, ] <- order(treat[i, ], na.last = T)
      r[i, ] <- event[i, order(t0[i, ], na.last = T)]
      m[i, ] <- mod[i, order(t0[i, ], na.last = T)]
      N[i, ] <- rand[i, order(t0[i, ], na.last = T)]
      t[i, ] <- sort(treat[i, ], na.last = T)
    }


    ## Turn into contrast-level data: one row per possible comparison in each trial ('
    without.MOD <- pairwise(as.list(t), event = as.list(r), n = as.list(N), data = cbind(t, r, N), studlab = 1:ns)[, c(3:6, 8, 7, 9)]
    with.MOD <- pairwise(as.list(t), event = as.list(m), n = as.list(N), data = cbind(t, m, N), studlab = 1:ns)[, c(3:6, 8)]


    ## Prepare the dataset to run the 'Taylor.IMOR' function
    two.stage.dataset <- cbind(without.MOD[, c(1, 4:5)], with.MOD[, 4:5], without.MOD[, c(6:7, 2:3)])


    ## Two-stage pattern-mixture model under MAR
    res.two.stage0 <- Taylor.IMOR(data = two.stage.dataset, delta1 = 0, delta2 = 0, var.delta1 = var.misspar[1], var.delta2 = var.misspar[1], rho = var.misspar[1])[, c(1, 11)]


    ## Indicate the number of possible comparisons in each trial
    arms0 <- arms1 <- list()
    for(i in 1:ns) {
      arms0[[i]] <- 1:dim(combn(na[i], 2))[2]
      arms1[[i]] <- rep(na[i], dim(combn(na[i], 2))[2] )
    }


    ## Obtain the adjusted standard errors per trial in the wide format
    res.two.stage <- as.data.frame(subset(cbind(res.two.stage0, unlist(arms0)), unlist(arms0) <= unlist(arms1) - 1))
    colnames(res.two.stage) <- c("id", "StD.error", "indic")
    SE <- cbind(rep(NA, ns), t(spread(res.two.stage, key = id, value = StD.error))[-1, ])
    Var <- SE*SE
    Prec <- 1/Var

  }


  ## Conditions for 'direction' and 'D'
  if (direction[[2]] == "reference" & D == 1) {

    ## Create the dummy variable 'Z'
    Z <- t
    Z[t == direction[[1]]] <- 1
    Z[t != direction[[1]]] <- -1


    ## Create the direction variable 'I'
    I <- cbind(rep(NA, ns),  (Z[, 1] - Z[, -1])*0.5)
    colnames(I) <- paste0("arm", 1:max(na))

  } else if (direction[[2]] == "reference" & D == 0) {

    ## Create the dummy variable 'Z'
    Z <- t
    Z[Z == direction[[1]]] <- -1
    Z[Z != direction[[1]]] <- 1


    ## Create the direction variable 'I'
    I <- cbind(rep(NA, ns),  (Z[, 1] - Z[, -1])*0.5)
    colnames(I) <- paste0("arm", 1:max(na))

  } else if (direction[[2]] == "new-old" & D == 1) {

    ## Re-name the interventions by their order defined in 'direction[[1]]'
    t.new <- t
    for(i in 1:nt){

      t.new[t == sort(na.omit(unique(unlist(t))))[i]] <- direction[[1]][i]
    }


    ## Turn into contrast-level data: one row per possible comparison in each trial ('netmeta')
    wide.format <- pairwise(as.list(t.new), event = as.list(m), n = as.list(N), data = cbind(t.new, m, N), studlab = 1:ns)[, 3:5]


    ## For multi-arm trials, keep comparisons with the baseline arm (i.e. the first arm)
    Z0 <- subset(cbind(wide.format, unlist(arms0)), unlist(arms0) <= unlist(arms1) - 1)[, -4]


    ## Create the dummy variable 'Z'
    Z0[, 3] <- ifelse(Z0[, 2] > Z0[, 3], -1, 1)
    Z0[, 2] <- -Z0[, 3]
    Z1 <- as.data.frame(cbind(Z0[, 1], apply(Z0[, 2:3], 1, mean)))
    colnames(Z1) <- c("studlab", "I")
    Z <- data.frame(Z1, subset(unlist(arms0), unlist(arms0) <= unlist(arms1) - 1))
    colnames(Z) <- c("studlab", "I", "indic")


    ## Create the direction variable 'I'
    I <- cbind(rep(NA, ns),  t(spread(Z, key = studlab, value = I))[-1, ])
    colnames(I) <- paste0("arm", 1:max(na))


  } else if (direction[[2]] == "new-old" & D == 0) {

    ## Re-name the interventions by their order defined in 'direction[[1]]'
    t.new <- t
    for(i in 1:nt){

      t.new[t == sort(na.omit(unique(unlist(t))))[i]] <- direction[[1]][i]
    }


    ## Turn into contrast-level data: one row per possible comparison in each trial ('netmeta')
    wide.format <- pairwise(as.list(t.new), event = as.list(m), n = as.list(N), data = cbind(t.new, m, N), studlab = 1:ns)[, 3:5]


    ## For multi-arm trials, keep comparisons with the baseline arm (i.e. the first arm)
    Z0 <- subset(cbind(wide.format, unlist(arms0)), unlist(arms0) <= unlist(arms1) - 1)[, -4]


    ## Create the dummy variable 'Z'
    Z0[, 3] <- ifelse(Z0[, 2] > Z0[, 3], 1, -1)
    Z0[, 2] <- -Z0[, 3]
    Z1 <- as.data.frame(cbind(Z0[, 1], apply(Z0[, 2:3], 1, mean)))
    colnames(Z1) <- c("studlab", "I")
    Z <- data.frame(Z1, subset(unlist(arms0), unlist(arms0) <= unlist(arms1) - 1))
    colnames(Z) <- c("studlab", "I", "indic")


    ## Create the direction variable 'I'
    I <- cbind(rep(NA, ns),  t(spread(Z, key = studlab, value = I))[-1, ])
    colnames(I) <- paste0("arm", 1:max(na))

  }


  return(list(covariate = I, variance = I*Var, stndard.error = I*SE, precision = I*Prec))

}

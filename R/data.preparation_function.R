#' Prepare the dataaet in the proper format for R2jags
#'
#' @param data The dataset
#' @param measure The effect measure
#'
#' @export
data.preparation <- function(data, measure) {


  measure <- if (missing(measure)) {
    stop("The 'measure' needs to be defined")
  } else if (!is.element(measure, c("MD", "SMD", "ROM", "OR"))) {
    stop("Insert 'MD', 'SMD', 'ROM', or 'OR'")
  } else if ((dim(data %>% dplyr::select(starts_with("r")))[2] > 0) & is.element(measure, c("MD", "SMD", "ROM"))) {
    stop("Inser 'OR' because the outcome data are binary")
  } else if ((dim(data %>% dplyr::select(starts_with("r")))[2] == 0) & is.element(measure, "OR")) {
    stop("Inser 'MD', 'SMD' or 'ROM' because the outcome data are continuous")
  } else {
    measure
  }


  treat <- data %>% dplyr::select(starts_with("t"))                               # Intervention studied in each arm of every trial
  ns <- length(treat[, 1])                                                      # Total number of included trials per network
  na <- apply(treat, 1, function(x) length(which(!is.na(x))))                   # Number of interventions investigated in every trial per network
  nt <- length(table(as.matrix(treat)))                                         # Total number of interventions per network
  ref <- 1                                                                      # The first intervention (t1 = 1) is the reference of the network


  ## For a continuous outcome
  if (is.element(measure, c("MD", "SMD", "ROM"))) {

    ## Continuous: arm-level, wide-format dataset
    y.obs <- data %>% dplyr::select(starts_with("y"))                             # Observed mean value in each arm of every trial
    sd.obs <- data %>% dplyr::select(starts_with("sd"))                           # Observed standard deviation in each arm of every trial
    mod <- data %>% dplyr::select(starts_with("m"))                               # Number of missing participants in each arm of every trial
    c <- data %>% dplyr::select(starts_with("c"))                                 # Number of completers in each arm of every trial
    se.obs <- sd.obs/sqrt(c)                                                      # Observed standard error in each arm of every trial
    rand <- mod + c                                                               # Number of randomised participants in each arm of every trial


    ## Order by 'id of t1' < 'id of t1'
    y0 <- se0 <- m <- N <- t <- t0 <- treat
    for (i in 1:ns) {
      t0[i, ] <- order(treat[i, ], na.last = T)
      y0[i, ] <- y.obs[i, order(t0[i, ], na.last = T)]
      se0[i, ] <- se.obs[i, order(t0[i, ], na.last = T)]
      m[i, ] <- mod[i, order(t0[i, ], na.last = T)]
      N[i, ] <- rand[i, order(t0[i, ], na.last = T)]
      t[i, ] <- sort(treat[i, ], na.last = T)
    }

  } else {


    ## Binary: arm-level, wide-format dataset
    event <- data %>% dplyr::select(starts_with("r"))                               # Number of observed events in each arm of every trial
    mod <- data %>% dplyr::select(starts_with("m"))                                 # Number of missing participants in each arm of every trial
    rand <- data %>% dplyr::select(starts_with("n"))                                # Number randomised participants in each arm of every trial


    ## Order by 'id of t1' < 'id of t2' < 'id of t3', and so on
    r <- m <- N <- t <- t0 <- treat
    for (i in 1:ns) {
      t0[i, ] <- order(treat[i, ], na.last = T)
      r[i, ] <- event[i, order(t0[i, ], na.last = T)]
      m[i, ] <- mod[i, order(t0[i, ], na.last = T)]
      N[i, ] <- rand[i, order(t0[i, ], na.last = T)]
      t[i, ] <- sort(treat[i, ], na.last = T)
    }

  }

  I <- m
  for (i in 1:ns) {
    I[i, ] <- ifelse(is.na(m[i, ]), 0, 1)
  }


  ## Return results
  results <- list(m = m,
                  N = N,
                  t = t,
                  I = I,
                  na = na,
                  nt = nt,
                  ns = ns,
                  ref = ref,
                  measure = measure)

  results <- if (is.element(measure, c("MD", "SMD", "ROM"))) {
    append(results, list(y0 = y0, se0 = se0))
  } else {
    append(results, list(r = r))
  }

  return(results)
}

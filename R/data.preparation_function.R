#' Prepare the dataset in the proper format for R2jags
#'
#' @description This function prepares the dataset in the proper format for R2jags. \code{data.preparation} is found in the \code{run.model} function; thus, the arguments are as specified in the latter.
#'
#' @param data A data-frame of the one-trial-per-row format containing arm-level data of each trial. This format is widely used for BUGS models. See 'Format' in the \code{run.model} function.
#' @param measure Character string indicating the effect measure with values \code{"OR"}, \code{"MD"}, \code{"SMD"}, or \code{"ROM"} for the odds ratio, mean difference,
#'   standardised mean difference and ratio of means, respectively.
#'
#' @return A list of data-frames on the following elements to be passed to \code{\link{run.model}}:
#' \tabular{ll}{
#'  \code{m} \tab The number of missing participant outcome data in each arm of every trial (see 'Details' in the \code{run.model} function).\cr
#'  \tab \cr
#'  \code{N} \tab The number of participants randomised on the assigned intervention in each arm.\cr
#'  \tab \cr
#'  \code{t} \tab The intervention identifier in each arm of every trial.\cr
#'  \tab \cr
#'  \code{I} \tab The pseudo-data-frame \code{I} (see 'Details' in the \code{run.model} function).\cr
#'  \tab \cr
#'  \code{r} \tab The observed number of events of the outcome in each arm of every trial, when the outcome is binary.\cr
#'  \tab \cr
#'  \code{y0} \tab The observed mean value of the outcome in each arm, when the outcome is continuous\cr
#'  \tab \cr
#'  \code{se0} \tab The observed standard deviation of the outcome in each arm, when the outcome is continuous\cr
#' }
#'
#' @author {Loukia M. Spineli}
#'
#' @seealso \code{\link{R2jags}}, \code{\link{run.model}}
#'
#' @export
data.preparation <- function(data, measure) {


  measure <- if (missing(measure)) {
    stop("The argument 'measure' needs to be defined", call. = F)
  } else if ((dim(data %>% dplyr::select(starts_with("r")))[2] > 0) & !is.element(measure, c("MD", "SMD", "ROM", "OR"))) {
    stop("Insert 'OR'", call. = F)
  } else if ((dim(data %>% dplyr::select(starts_with("r")))[2] == 0) & !is.element(measure, c("MD", "SMD", "ROM", "OR"))) {
    stop("Insert 'MD', 'SMD' or 'ROM'", call. = F)
  } else if ((dim(data %>% dplyr::select(starts_with("r")))[2] > 0) & is.element(measure, c("MD", "SMD", "ROM"))) {
    stop("Insert 'OR' because the outcome data are binary", call. = F)
  } else if ((dim(data %>% dplyr::select(starts_with("r")))[2] == 0) & is.element(measure, "OR")) {
    stop("Insert 'MD', 'SMD' or 'ROM' because the outcome data are continuous", call. = F)
  } else {
    measure
  }


  treat <- data %>% dplyr::select(starts_with("t"))                               # Intervention studied in each arm of every trial
  ns <- length(treat[, 1])                                                      # Total number of included trials per network
  na <- apply(treat, 1, function(x) length(which(!is.na(x))))                   # Number of interventions investigated in every trial per network
  nt <- length(table(as.matrix(treat)))                                         # Total number of interventions per network
  ref <- 1                                                                      # The first intervention (t1 = 1) is the reference of the network


  ## When no missing outcome data are collected
  mod <- if (dim(data %>% dplyr::select(starts_with("m")))[2] == 0) {
    message("Missing participant outcome data have *not* been collected")
    as.data.frame(matrix(NA, nrow = nrow(treat), ncol = ncol(treat)))
  } else {
    data %>% dplyr::select(starts_with("m"))                                    # Number of missing participants in each arm of every trial
  }


  ## For a continuous outcome
  if (is.element(measure, c("MD", "SMD", "ROM"))) {

    ## Continuous: arm-level, wide-format dataset
    y.obs <- data %>% dplyr::select(starts_with("y"))                             # Observed mean value in each arm of every trial
    sd.obs <- data %>% dplyr::select(starts_with("sd"))                           # Observed standard deviation in each arm of every trial
    rand <- data %>% dplyr::select(starts_with("n"))                              # Number randomised participants in each arm of every trial
    se.obs <- sd.obs/sqrt(rand - mod)                                             # Observed standard error in each arm of every trial


    ## Order by 'id of t1' < 'id of t1'
    y0 <- sd0 <- se0 <- m <- N <- t <- t0 <- treat
    for (i in 1:ns) {
      t0[i, ] <- order(treat[i, ], na.last = T)
      y0[i, ] <- y.obs[i, order(t0[i, ], na.last = T)]
      sd0[i, ] <- sd.obs[i, order(t0[i, ], na.last = T)]
      se0[i, ] <- se.obs[i, order(t0[i, ], na.last = T)]
      m[i, ] <- mod[i, order(t0[i, ], na.last = T)]
      N[i, ] <- rand[i, order(t0[i, ], na.last = T)]
      t[i, ] <- sort(treat[i, ], na.last = T)
    }

  } else {


    ## Binary: arm-level, wide-format dataset
    event <- data %>% dplyr::select(starts_with("r"))                               # Number of observed events in each arm of every trial
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


  ## Use a dummy variable to indicate which trials report missing outcome data fully (i.e. for each arm)
  ## If a trial reports missing outcome data partially or not at all, insert 'NA'.
  I <- m.new <- m
  for (i in 1:ns) {
    for (k in 1:na[i]) {
      I[i, k] <- if (is.na(m[i, k]) & !is.na(N[i, k])) {
        0
      } else if (is.na(m[i, k]) & is.na(N[i, k])) {
        NA
      } else if (!is.na(m[i, k]) & !is.na(N[i, k])) {
        1
      }

      m.new[i, k] <- if (is.na(m[i, k]) & !is.na(N[i, k])) {
        0
      } else if (is.na(m[i, k]) & is.na(N[i, k])) {
        NA
      } else if (!is.na(m[i, k]) & !is.na(N[i, k])) {
        m[i, k]
      }
    }
  }

  ## Return results
  results <- list(m = m.new,
                  N = N,
                  t = t,
                  I = I,
                  na = na,
                  nt = nt,
                  ns = ns,
                  ref = ref,
                  measure = measure)

  results <- if (is.element(measure, c("MD", "SMD", "ROM"))) {
    append(results, list(y0 = y0, se0 = se0, sd0 = sd0))
  } else {
    append(results, list(r = r))
  }

  return(results)
}

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
#' @details \code{data.preparation} prepares the data for the Bayesian analysis. The data preparation includes the following actions. First, \code{data.preparation} checks
#'   whether the element \strong{m} exists in \code{data}. If this element is missing, \code{data.preparation} creates a pseudo-data-frame for \code{m} that has the zero value for the observed trial-arms, and \code{NA} for the unobserved trial-arms,
#'   and the pseudo-data-frame \code{I} that has the same values with the pseudo-data-frame for \code{m}. If the element \strong{m} exists in \code{data} and has values only for some trials, the pseudo-data-frame for \code{m} has the same values with
#'   \code{m} for the corresponding trial-arms, and the pseudo-data-frame \code{I} has the value one for these trial-arms. Both pseudo-data-frames aim to retain the trials without information on MOD in \code{data} and allow running the Bayesian analysis
#'   when MOD have not been extracted for any trial of the dataset. Second, \code{data.preparation} sorts the interventions across the arms of each trial in an ascending order and correspondingly sorts within each trial the remaining elements in \code{data}
#'   (see also 'Format' in \code{run.model}). Since all Bayesian models in the package consider the first column in \strong{t} to include the control arm for every trial, this sorting ensures that interventions with a lower identifier are consistently treated as the control arm in each trial. This case is relevant in
#'   non-star-shaped networks.
#'
#' @author {Loukia M. Spineli}
#'
#' @seealso \href{https://CRAN.R-project.org/package=R2jags}{R2jags}, \code{\link{run.model}}
#'
#' @export
data.preparation <- function(data, measure) {

  options(warn = -1)



  treat <- if (dim(data %>% select(starts_with("t")))[2] == 0) {
    stop("The information on the individual arms is missing", call. = F)
  } else {
    data %>% select(starts_with("t"))                                           # Intervention studied in each arm of every trial
  }
  ns <- length(treat[, 1])                                                      # Total number of included trials per network
  na <- apply(treat, 1, function(x) length(which(!is.na(x))))                   # Number of interventions investigated in every trial per network
  nt <- length(table(as.matrix(treat)))                                         # Total number of interventions per network
  ref <- 1                                                                      # The first intervention (t1 = 1) is the reference of the network


  measure <- if (missing(measure)) {
    stop("The argument 'measure' needs to be defined", call. = F)
  } else if ((dim(data %>% select(starts_with("r")))[2] > 0) & !is.element(measure, c("MD", "SMD", "ROM", "OR"))) {
    stop("Insert 'OR'", call. = F)
  } else if ((dim(data %>% select(starts_with("r")))[2] == 0) & !is.element(measure, c("MD", "SMD", "ROM", "OR"))) {
    stop("Insert 'MD', 'SMD' or 'ROM'", call. = F)
  } else if ((dim(data %>% select(starts_with("r")))[2] > 0) & is.element(measure, c("MD", "SMD", "ROM"))) {
    stop("Insert 'OR' because the outcome data are binary", call. = F)
  } else if ((dim(data %>% select(starts_with("r")))[2] == 0) & is.element(measure, "OR")) {
    stop("Insert 'MD', 'SMD' or 'ROM' because the outcome data are continuous", call. = F)
  } else {
    measure
  }


  ## When no missing outcome data are collected
  mod <- if (dim(data %>% select(starts_with("m")))[2] == 0) {
    message("Missing participant outcome data have *not* been collected")
    as.data.frame(matrix(NA, nrow = nrow(treat), ncol = ncol(treat)))
  } else {
    data %>% select(starts_with("m"))                                    # Number of missing participants in each arm of every trial
  }


  ## For a continuous outcome
  if (is.element(measure, c("MD", "SMD", "ROM"))) {

    ## Continuous: arm-level, wide-format dataset
    y.obs <- data %>% select(starts_with("y"))                             # Observed mean value in each arm of every trial
    sd.obs <- data %>% select(starts_with("sd"))                           # Observed standard deviation in each arm of every trial
    rand <- data %>% select(starts_with("n"))                              # Number randomised participants in each arm of every trial
    se.obs <- sd.obs/sqrt(rand - mod)                                      # Observed standard error in each arm of every trial

    if ((dim(y.obs)[2] != max(na)) | (dim(sd.obs)[2] != max(na)) | (dim(mod)[2] != max(na)) | (dim(rand)[2] != max(na))) {
      stop("All elements in the argument 'data' (e.g. information on the interventions arms) must have the same dimension", call. = F)
    }

    if ((dim(y.obs)[1] != ns) | (dim(sd.obs)[1] != ns) |(dim(mod)[1] != ns) | (dim(rand)[1] != ns)) {
      stop("The elements in the argument 'data' (e.g. information on the interventions arms) must have the same dimension", call. = F)
    }


    ## Order by 'id of t1' < 'id of t1'
    y0 <- sd0 <- se0 <- m <- N <- t <- treat
    for (i in 1:ns) {
      y0[i, ] <- y.obs[i, order(treat[i, ], na.last = T)]
      sd0[i, ] <- sd.obs[i, order(treat[i, ], na.last = T)]
      se0[i, ] <- se.obs[i, order(treat[i, ], na.last = T)]
      m[i, ] <- mod[i, order(treat[i, ], na.last = T)]
      N[i, ] <- rand[i, order(treat[i, ], na.last = T)]
      t[i, ] <- sort(treat[i, ], na.last = T)
    }

    names(y0) <- paste0("y", 1:length(y0[1, ]))
    names(sd0) <- paste0("sd", 1:length(sd0[1, ]))
    names(se0) <- paste0("se", 1:length(se0[1, ]))
    names(m) <- paste0("m", 1:length(m[1, ]))
    names(N) <- paste0("n", 1:length(N[1, ]))
    names(t) <- paste0("t", 1:length(t[1, ]))


  } else {


    ## Binary: arm-level, wide-format dataset
    event <- data %>% select(starts_with("r"))                               # Number of observed events in each arm of every trial
    rand <- data %>% select(starts_with("n"))                                # Number randomised participants in each arm of every trial


    if ((dim(event)[2] != max(na)) | (dim(mod)[2] != max(na)) | (dim(rand)[2] != max(na))) {
      stop("The elements in the argument 'data' (e.g. information on the interventions arms) must have the same dimension", call. = F)
    }

    if ((dim(event)[1] != ns) | (dim(mod)[1] != ns) | (dim(rand)[1] != ns)) {
      stop("The elements in the argument 'data' (e.g. information on the interventions arms) must have the same dimension", call. = F)
    }


    ## Order by 'id of t1' < 'id of t2' < 'id of t3', and so on
    r <- m <- N <- t <- treat
    for (i in 1:ns) {
      r[i, ] <- event[i, order(treat[i, ], na.last = T)]
      m[i, ] <- mod[i, order(treat[i, ], na.last = T)]
      N[i, ] <- rand[i, order(treat[i, ], na.last = T)]
      t[i, ] <- sort(treat[i, ], na.last = T)
    }

    names(r) <- paste0("r", 1:length(r[1, ]))
    names(m) <- paste0("m", 1:length(m[1, ]))
    names(N) <- paste0("n", 1:length(N[1, ]))
    names(t) <- paste0("t", 1:length(t[1, ]))

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

  names(I) <- paste0("I", 1:length(I[1, ]))

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

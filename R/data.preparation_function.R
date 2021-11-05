#' Prepare the dataset in the proper format for R2jags
#'
#' @description
#'   \code{data_preparation} prepares the dataset in the proper format for
#'   R2jags ans returns alist of elements that \code{\link{run_model}} inherites
#'   via the argument \code{data}.
#'
#' @param data A data-frame of the one-trial-per-row format with arm-level data.
#'   See 'Format' for the specification of the columns.
#' @param measure Character string indicating the effect measure with values
#'   \code{"OR"}, \code{"MD"}, \code{"SMD"}, or \code{"ROM"} for the odds ratio,
#'   mean difference, standardised mean difference and ratio of means,
#'   respectively.
#'
#' @return A list of data-frames on the following elements to be passed
#'   to \code{\link{run_model}}:
#'   \tabular{ll}{
#'    \code{m} \tab The number of missing participant outcome data in each
#'    trial-arm
#'    (see 'Details' in the \code{run_model} function).\cr
#'    \tab \cr
#'    \code{N} \tab The number of randomised participants in each trial-arm.\cr
#'    \tab \cr
#'    \code{t} \tab The intervention identifier in each trial-arm.\cr
#'    \tab \cr
#'    \code{I} \tab The pseudo-data-frame \code{I}
#'    (see 'Details' in \code{run_model}).\cr
#'    \tab \cr
#'    \code{r} \tab The observed number of events of the outcome in each
#'    trial-arm,
#'    when the outcome is binary.\cr
#'    \tab \cr
#'    \code{y0} \tab The observed mean value of the outcome in each trial-arm,
#'    when the outcome is continuous\cr
#'    \tab \cr
#'    \code{se0} \tab The observed standard deviation of the outcome in each
#'    trial-arm,
#'    when the outcome is continuous\cr
#'   }
#'
#' @details \code{data_preparation} prepares the data for the Bayesian analysis.
#'   The data preparation includes the following actions.
#'   First, \code{data_preparation} checks whether the element
#'   \strong{m} exists in the \code{data}. If this element is missing,
#'   \code{data_preparation} creates a pseudo-data-frame for \code{m} that has
#'   zero value for the observed trial-arms, and \code{NA} for the unobserved
#'   trial-arms, and the pseudo-data-frame \code{I} that has the same values
#'   with the pseudo-data-frame for \code{m}. If the element \strong{m} exists
#'   in the \code{data} and has values only for some trials,
#'   the pseudo-data-frame for \code{m} is the same with \code{m} for the
#'   corresponding trial-arms, and the pseudo-data-frame \code{I} has the value
#'   one for these trial-arms. Both pseudo-data-frames aim to retain the trials
#'   without information on missing outcome data in \code{data}. Second,
#'   \code{data_preparation} sorts the interventions across the arms of each
#'   trial in an ascending order and correspondingly the remaining elements in
#'   \code{data} (see 'Format' in \code{run_model}). Since all Bayesian models
#'   in the package consider the first column in \strong{t} as being control
#'   arm for every trial, this sorting ensures that interventions with a lower
#'   identifier are consistently treated as the control arm in each trial.
#'   This case is relevant in non-star-shaped networks.
#'
#' @author {Loukia M. Spineli}
#'
#' @seealso \href{https://CRAN.R-project.org/package=R2jags}{R2jags},
#'   \code{\link{run_model}}
#'
#' @export
data_preparation <- function(data, measure) {

  options(warn = -1)

  # Intervention studied in each arm of every trial
  treat <- if (dim(data %>% select(starts_with("t")))[2] == 0) {
    stop("The information on the individual arms is missing", call. = F)
  } else {
    data %>% select(starts_with("t"))
  }
  # Total number of included trials
  ns <- length(treat[, 1])
  # Number of interventions investigated in every trial
  na <- apply(treat, 1, function(x) length(which(!is.na(x))))
  # Total number of interventions
  nt <- length(table(as.matrix(treat)))
  # The intervention with identifier '1' is the reference of the network
  ref <- 1

  measure <- if (missing(measure)) {
    stop("The argument 'measure' needs to be defined", call. = F)
  } else if ((dim(data %>% select(starts_with("r")))[2] > 0) &
             !is.element(measure, c("MD", "SMD", "ROM", "OR"))) {
    stop("Insert 'OR'", call. = F)
  } else if ((dim(data %>% select(starts_with("r")))[2] == 0) &
             !is.element(measure, c("MD", "SMD", "ROM", "OR"))) {
    stop("Insert 'MD', 'SMD' or 'ROM'", call. = F)
  } else if ((dim(data %>% select(starts_with("r")))[2] > 0) &
             is.element(measure, c("MD", "SMD", "ROM"))) {
    stop("Insert 'OR' for a  binary outcome", call. = F)
  } else if ((dim(data %>% select(starts_with("r")))[2] == 0) &
             is.element(measure, "OR")) {
    stop("Insert 'MD', 'SMD' or 'ROM' for a continuous outcome", call. = F)
  } else {
    measure
  }

  # When no missing outcome data are collected
  mod <- if (dim(data %>% select(starts_with("m")))[2] == 0) {
    message("Missing participant outcome data have *not* been collected")
    as.data.frame(matrix(NA, nrow = nrow(treat), ncol = ncol(treat)))
  } else {
    # Number of missing participants in each arm of every trial
    data %>% select(starts_with("m"))
  }

  # For a continuous outcome
  if (is.element(measure, c("MD", "SMD", "ROM"))) {
    # Observed mean value in each arm of every trial
    y_obs <- data %>% select(starts_with("y"))
    # Observed standard deviation in each arm of every trial
    sd_obs <- data %>% select(starts_with("sd"))
    # Number of randomised participants in each arm of every trial
    rand <- data %>% select(starts_with("n"))
    # Observed standard error in each arm of every trial
    se_obs <- sd_obs / sqrt(rand - mod)

    if ((dim(y_obs)[2] != max(na)) |
        (dim(sd_obs)[2] != max(na)) |
        (dim(mod)[2] != max(na)) |
        (dim(rand)[2] != max(na))) {
      stop("All elements must have the same dimension", call. = F)
    }

    if ((dim(y_obs)[1] != ns) |
        (dim(sd_obs)[1] != ns) |
        (dim(mod)[1] != ns) |
        (dim(rand)[1] != ns)) {
      stop("All elements must have the same dimension", call. = F)
    }

    # Order by 'id of t1' < 'id of t1'
    y0 <- sd0 <- se0 <- m <- N <- t <- treat
    for (i in 1:ns) {
      y0[i, ] <- y_obs[i, order(treat[i, ], na.last = T)]
      sd0[i, ] <- sd_obs[i, order(treat[i, ], na.last = T)]
      se0[i, ] <- se_obs[i, order(treat[i, ], na.last = T)]
      m[i, ] <- mod[i, order(treat[i, ], na.last = T)]
      N[i, ] <- rand[i, order(treat[i, ], na.last = T)]
      t[i, ] <- sort(treat[i, ], na.last = T)
    }

    names(y0) <- paste0("y", seq_len(max(na)))
    names(sd0) <- paste0("sd", seq_len(max(na)))
    names(se0) <- paste0("se", seq_len(max(na)))
    names(m) <- paste0("m", seq_len(max(na)))
    names(N) <- paste0("n", seq_len(max(na)))
    names(t) <- paste0("t", seq_len(max(na)))
  } else {
    # For a binary outcome:
    # Number of observed events in each arm of every trial
    event <- data %>% select(starts_with("r"))
    # Number of randomised participants in each arm of every trial
    rand <- data %>% select(starts_with("n"))

    if ((dim(event)[2] != max(na)) |
        (dim(mod)[2] != max(na)) |
        (dim(rand)[2] != max(na))) {
      stop("All elements must have the same dimension", call. = F)
    }

    if ((dim(event)[1] != ns) |
        (dim(mod)[1] != ns) |
        (dim(rand)[1] != ns)) {
      stop("All elements must have the same dimension", call. = F)
    }

    # Order by 'id of t1' < 'id of t2' < 'id of t3', and so on
    r <- m <- N <- t <- treat
    for (i in 1:ns) {
      r[i, ] <- event[i, order(treat[i, ], na.last = T)]
      m[i, ] <- mod[i, order(treat[i, ], na.last = T)]
      N[i, ] <- rand[i, order(treat[i, ], na.last = T)]
      t[i, ] <- sort(treat[i, ], na.last = T)
    }

    names(r) <- paste0("r", seq_len(max(na)))
    names(m) <- paste0("m", seq_len(max(na)))
    names(N) <- paste0("n", seq_len(max(na)))
    names(t) <- paste0("t", seq_len(max(na)))
  }

  # Indicate the trials with fully reported missing outcome data (for each arm)
  # If a trial reports missing participant outcome data partially or not at all,
  # insert 'NA'.
  I <- m_new <- m
  for (i in 1:ns) {
    for (k in 1:na[i]) {
      I[i, k] <- if (is.na(m[i, k]) & !is.na(N[i, k])) {
        0
      } else if (is.na(m[i, k]) & is.na(N[i, k])) {
        NA
      } else if (!is.na(m[i, k]) & !is.na(N[i, k])) {
        1
      }

      m_new[i, k] <- if (is.na(m[i, k]) & !is.na(N[i, k])) {
        0
      } else if (is.na(m[i, k]) & is.na(N[i, k])) {
        NA
      } else if (!is.na(m[i, k]) & !is.na(N[i, k])) {
        m[i, k]
      }
    }
  }

  names(I) <- paste0("I", seq_len(max(na)))

  # Return results
  results <- list(m = m_new,
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

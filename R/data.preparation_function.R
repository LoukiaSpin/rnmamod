#' Prepare the dataset in the proper format for R2jags
#'
#' @description
#'   \code{data_preparation} prepares the dataset in the proper format for
#'   R2jags and returns a list of elements that \code{\link{run_model}} inherits
#'   via the argument \code{data}.
#'
#' @param data A data-frame of the one-trial-per-row format with arm-level data.
#'   See 'Format' in \code{\link{run_model}}.
#' @param measure Character string indicating the effect measure. For a binary
#'   outcome, the following can be considered: \code{"OR"}, \code{"RR"} or
#'   \code{"RD"} for the odds ratio, relative risk, and risk difference,
#'   respectively. For a continuous outcome, the following can be considered:
#'   \code{"MD"}, \code{"SMD"}, or \code{"ROM"} for mean difference,
#'   standardised mean difference and ratio of means, respectively.
#'
#' @return A list of data-frames on the following elements to be passed
#'   to \code{\link{run_model}}:
#'   \item{m}{The number of missing participant outcome data in each
#'   trial-arm (see 'Details').}
#'   \item{N}{The number of randomised participants in each trial-arm.}
#'   \item{t}{The intervention identifier in each trial-arm.}
#'   \item{I}{A pseudo-data-frame (see 'Details').}
#'   \item{r}{The number of observed events of the outcome in each
#'   trial-arm, when the outcome is binary.}
#'   \item{y0}{The observed mean value of the outcome in each trial-arm,
#'   when the outcome is continuous.}
#'   \item{se0}{The observed standard deviation of the outcome in each
#'   trial-arm, when the outcome is continuous.}
#'
#' @details \code{data_preparation} prepares the data for the Bayesian analysis
#'   (See 'Format' in \code{\link{run_model}}). \code{data_preparation}
#'   creates the pseudo-data-frames \code{m_new}, \code{I}, and \code{m_pseudo}
#'   that have the same dimensions with the element \code{N}. \code{m_new} takes
#'   the zero value for the observed trial-arms with unreported missing
#'   participant outcome data (i.e., \strong{m} equals \code{NA} for the
#'   corresponding trial-arms), the same value with \strong{m} for the observed
#'   trial-arms with reported missing participant outcome data, and \code{NA}
#'   for the unobserved trial-arms. \code{I} is a dummy data-frame and takes the
#'   value one for the observed trial-arms with reported missing participant
#'   outcome data, the zero value for the observed trial-arms with unreported
#'   missing participant outcome data (i.e., \code{m_new} equals zero for the
#'   corresponding trial-arms), and \code{NA} for the unobserved trial-arms.
#'   Thus, \code{I} indicates whether missing participant outcome data have been
#'   collected for the observed trial-arms. If the user has not defined the
#'   element \strong{m} in \code{data_preparation}, \code{m_new} and \code{I}
#'   take the zero value for all observed trial-arms to indicate that no missing
#'   participant outcome data have been collected for the analysed outcome.
#'   \code{I} and \code{m_new} are used from the following functions of the
#'   package: \code{\link{run_model}}, \code{\link{run_metareg}},
#'   \code{\link{prepare_model}}, \code{\link{run_nodesplit}},
#'   \code{\link{prepare_nodesplit}}, \code{\link{run_ume}},
#'   \code{\link{prepare_ume}}, and  \code{\link{run_sensitivity}}.
#'   Lastly, \code{m_pseudo} is a variant of \code{m_new}: it takes the value -1
#'   for the observed trial-arms with unreported missing participant outcome
#'   data (i.e., \strong{m} equals \code{NA} for the corresponding trial-arms),
#'   the same value with \strong{m} for the observed trial-arms with reported
#'   missing participant outcome data, and \code{NA} for the unobserved
#'   trial-arms. It is used in function \code{\link{heatmap_missing_network}} to
#'   calculate and illustrate the percentage of missing participant outcome data
#'   across the observed comparisons and interventions of the network and the
#'   function \code{\link{heatmap_missing_dataset}} to illustrate the trial-arms
#'   with unreported missing participant outcome data. All pseudo-data-frames
#'   aim to retain the trials without information on missing participant outcome
#'   data.
#'
#'   Furthermore, \code{data_preparation} sorts the interventions across
#'   the arms of each trial in an ascending order and correspondingly the
#'   remaining elements in \code{data} (See 'Format' in
#'   \code{\link{run_model}}). \code{data_preparation} considers the first
#'   column in \strong{t} as being the control arm for every trial. Thus,
#'   this sorting ensures that interventions with a lower identifier are
#'   consistently treated as the control arm in each trial. This case is
#'   relevant in non-star-shaped networks.
#'
#' @author {Loukia M. Spineli}
#'
#' @seealso \code{\link{heatmap_missing_dataset}},
#'   \code{\link{heatmap_missing_network}},
#'   \href{https://CRAN.R-project.org/package=R2jags}{R2jags},
#'   \code{\link{run_metareg}}, \code{\link{run_model}},
#'   \code{\link{run_nodesplit}}, \code{\link{run_sensitivity}},
#'   \code{\link{run_ume}}, \code{\link{prepare_model}},
#'   \code{\link{prepare_nodesplit}}, \code{\link{prepare_ume}}
#'
#' @export
data_preparation <- function(data, measure) {

  # Intervention studied in each arm of every trial
  treat <- if (dim(data %>% select(starts_with("t")))[2] == 0) {
    stop("The information on the individual arms is missing", call. = FALSE)
  } else {
    data %>% select(starts_with("t"))
  }
  # Total number of included trials
  ns <- length(treat[, 1])
  # Number of interventions investigated in every trial
  na <- apply(treat, 1, function(x) length(which(!is.na(x))))
  # Total number of interventions
  nt <- length(table(as.matrix(treat)))

  measure <- if (missing(measure)) {
    stop("The argument 'measure' needs to be defined.", call. = FALSE)
  } else if ((dim(data %>% select(starts_with("r")))[2] > 0) &
             !is.element(measure, c("MD", "SMD", "ROM", "OR", "RR", "RD"))) {
    stop("Insert 'OR', 'RR', or 'RD'.", call. = FALSE)
  } else if ((dim(data %>% select(starts_with("r")))[2] == 0) &
             !is.element(measure, c("MD", "SMD", "ROM", "OR", "RR", "RD"))) {
    stop("Insert 'MD', 'SMD' or 'ROM'.", call. = FALSE)
  } else if ((dim(data %>% select(starts_with("r")))[2] > 0) &
             is.element(measure, c("MD", "SMD", "ROM"))) {
    stop("Insert 'OR', 'RR', or 'RD' for a  binary outcome.", call. = FALSE)
  } else if ((dim(data %>% select(starts_with("r")))[2] == 0) &
             is.element(measure, "OR")) {
    stop("Insert 'MD', 'SMD' or 'ROM' for a continuous outcome.", call. = FALSE)
  } else {
    measure
  }

  # When no missing outcome data are collected
  mod <- if (dim(data %>% select(starts_with("m")))[2] == 0) {
    message("Missing participant outcome data have *not* been collected.")
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
      stop("All elements must have the same dimension.", call. = FALSE)
    }

    if ((dim(y_obs)[1] != ns) |
        (dim(sd_obs)[1] != ns) |
        (dim(mod)[1] != ns) |
        (dim(rand)[1] != ns)) {
      stop("All elements must have the same dimension.", call. = FALSE)
    }

    # Order by 'id of t1' < 'id of t1'
    treat_list <- rand_list <- list()
    y_obs_list <- sd_obs_list <- se_obs_list <- mod_list <- treat_list
    for (i in 1:ns) {
      treat_list[[i]] <- treat[i, ]
      y_obs_list[[i]] <- y_obs[i, ]
      sd_obs_list[[i]] <- sd_obs[i, ]
      se_obs_list[[i]] <- se_obs[i, ]
      mod_list[[i]] <- mod[i, ]
      rand_list[[i]] <- rand[i, ]
    }

    y0 <- sd0 <- se0 <- m <- N <- t <- treat
    for (i in 1:ns) {
      y0[i, ] <- unlist(y_obs_list[[i]])[order(unlist(treat_list[[i]]),
                                              na.last = TRUE)]
      sd0[i, ] <- unlist(sd_obs_list[[i]])[order(unlist(treat_list[[i]]),
                                               na.last = TRUE)]
      se0[i, ] <- unlist(se_obs_list[[i]])[order(unlist(treat_list[[i]]),
                                                 na.last = TRUE)]
      m[i, ] <- unlist(mod_list[[i]])[order(unlist(treat_list[[i]]),
                                            na.last = TRUE)]
      N[i, ] <- unlist(rand_list[[i]])[order(unlist(treat_list[[i]]),
                                             na.last = TRUE)]
      t[i, ] <- sort(unlist(treat_list[[i]]), na.last = TRUE)
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
      stop("All elements must have the same dimension.", call. = FALSE)
    }

    if ((dim(event)[1] != ns) |
        (dim(mod)[1] != ns) |
        (dim(rand)[1] != ns)) {
      stop("All elements must have the same dimension.", call. = FALSE)
    }

    # Order by 'id of t1' < 'id of t2' < 'id of t3', and so on
    treat_list <- event_list <- mod_list <- rand_list <- list()
    for (i in 1:ns) {
      treat_list[[i]] <- treat[i, ]
      event_list[[i]] <- event[i, ]
      mod_list[[i]] <- mod[i, ]
      rand_list[[i]] <- rand[i, ]
    }

    r <- m <- N <- t <- treat
    for (i in 1:ns) {
     r[i, ] <- unlist(event_list[[i]])[order(unlist(treat_list[[i]]),
                                             na.last = TRUE)]
     m[i, ] <- unlist(mod_list[[i]])[order(unlist(treat_list[[i]]),
                                           na.last = TRUE)]
     N[i, ] <- unlist(rand_list[[i]])[order(unlist(treat_list[[i]]),
                                            na.last = TRUE)]
     t[i, ] <- sort(unlist(treat_list[[i]]), na.last = TRUE)
    }

    names(r) <- paste0("r", seq_len(max(na)))
    names(m) <- paste0("m", seq_len(max(na)))
    names(N) <- paste0("n", seq_len(max(na)))
    names(t) <- paste0("t", seq_len(max(na)))
  }

  # Indicate the trials with fully reported missing outcome data (for each arm)
  # If a trial reports missing participant outcome data partially or not at all,
  # insert 'NA'.
  I <- m_new <- m_pseudo <- m
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

      #' Matrix that indicates with -1 the trial-arms with non-reported missing
      #' participants. It is used by the following functions:
      #' 'heatmap_missing_dataset' and 'heatmap_missing_network'.
      m_pseudo[i, k] <- if (!is.na(t[i, k]) & is.na(m[i, k])) {
        -1
      } else if (is.na(m[i, k]) & is.na(N[i, k])) {
        NA
      } else if (!is.na(m[i, k]) & !is.na(N[i, k])) {
        m[i, k]
      }
    }
  }
  names(I) <- paste0("I", seq_len(max(na)))

  # Return results
  results <- list(m_pseudo = m_pseudo,
                  m = m_new,
                  N = N,
                  t = t,
                  I = I,
                  na = na,
                  nt = nt,
                  ns = ns,
                  measure = measure)

  results <- if (is.element(measure, c("MD", "SMD", "ROM"))) {
    append(results, list(y0 = y0, se0 = se0, sd0 = sd0))
  } else if (is.element(measure, c("OR", "RR", "RD"))) {
    append(results, list(r = r))
  }

  return(results)
}

#' A function to describe the evidence base
#'
#' @description Calculates the necessary elements to describe the evidence base
#'   for an outcome across the network, the interventions, and observed
#'   comparisons.
#'
#' @param data A data-frame of a one-trial-per-row format containing arm-level
#'   data of each trial. See 'Format' in \code{\link{run_model}}.
#' @param drug_names A vector of labels with the name of the interventions in
#'   the order they appear in the argument \code{data}.
#' @param measure Character string indicating the effect measure. For a binary
#'   outcome, the following can be considered: \code{"OR"}, \code{"RR"} or
#'   \code{"RD"} for the odds ratio, relative risk, and risk difference,
#'   respectively. For a continuous outcome, the following can be considered:
#'   \code{"MD"}, \code{"SMD"}, or \code{"ROM"} for mean difference,
#'   standardised mean difference and ratio of means, respectively.
#' @param save_xls Logical to indicate whether to export the tabulated results
#'   to an 'xlsx' file (via the \code{\link[writexl:write_xlsx]{write_xlsx}}
#'   function of the R-package
#'   \href{https://CRAN.R-project.org/package=writexl}{writexl}) at the working
#'   directory of the user. The default is \code{FALSE} (do not export).
#'
#' @return
#'   \code{describe_network} returns the following data-frames that describe the
#'   evidence base:
#'   \item{network_description}{The number of: interventions, possible
#'   comparisons, direct and indirect comparisons, number of trials in total,
#'   number of two-arm and multi-arm trials, number of randomised participants,
#'   and proportion of participants completing the trial (completers).
#'   When the outcome is binary, the number of trials with at least one zero
#'   event, and the number of trials with all zero events are also
#'   presented.}
#'   \item{table_interventions}{For each intervention, the number of trials,
#'   number of randomised participants, and proportion of completers. When the
#'   outcome is binary, the data-frame presents also the corresponding
#'   proportion of total observed events, the minimum, median and maximum
#'   proportion of observed events across the corresponding trials.}
#'   \item{table_comparisons}{Identical structure to \code{table_interventions}
#'   but for each observed comparison in the network.}
#'
#' @details \code{describe_network} calls \code{\link{data_preparation}} to
#'   facilitate the calculations.
#'
#'   Furthermore, \code{describe_network} exports the data-frames to separate
#'   'xlsx' files (via the \code{\link[writexl:write_xlsx]{write_xlsx}} function
#'   of the R-package
#'   \href{https://CRAN.R-project.org/package=writexl}{writexl}) at the working
#'   directory of the user.
#'
#' @seealso \code{\link{data_preparation}}, \code{\link{run_model}}
#'   \code{\link[writexl:write_xlsx]{write_xlsx}}
#'
#' @author {Loukia M. Spineli}
#'
#' @export
describe_network <- function(data, drug_names, measure, save_xls) {

  drug_names <- if (missing(drug_names)) {
    stop("The argument 'drug_names' has not been defined.", call. = FALSE)
  } else {
    drug_names
  }

  measure <- if (missing(measure)) {
    stop("The argument 'measure' needs to be defined.", call. = FALSE)
  } else if ((dim(data[, startsWith(colnames(data), "r")])[2] > 0) &
             !is.element(measure, c("OR", "RR", "RD"))) {
    stop("Insert 'OR', 'RR', or 'RD' for a  binary outcome.", call. = FALSE)
  } else if ((dim(data[, startsWith(colnames(data), "r")])[2] == 0) &
             !is.element(measure, c("MD", "SMD", "ROM"))) {
    stop("Insert 'MD', 'SMD' or 'ROM' for a  continuous outcome.", call. = FALSE)
  } else {
    measure
  }

  save_xls <- if (missing(save_xls)) {
    FALSE
  } else {
    save_xls
  }

  # Use the 'data_preparation' function
  dat <- data_preparation(data, measure)

  # Number of two-arm trian
  two_arm_ns <- length(which(dat$na < 3))

  # Number of multi-arm trials
  multi_arm_ns <- length(which(dat$na > 2))

  # Number of randomised per intervention
  total_rand_partic_interv <- aggregate(unlist(dat$N),
                                        by = list(unlist(dat$t)), sum)[, 2]

  # Number of completers per intervention
  total_obs_partic_interv <- aggregate(unlist(dat$N - dat$m),
                                       by = list(unlist(dat$t)), sum)[, 2]

  # Proportion of completers per intervention (in %)
  prop_obs_partic_interv <- round(total_obs_partic_interv /
                                  total_rand_partic_interv, 2) * 100

  # Proportion of missing outcome data (MOD) per intervention (in %)
  prop_mod_interv <- 100 - prop_obs_partic_interv

  # Proportion of MOD per trial-arm
  arm_mod <- dat$m / dat$N

  # Minimum proportion of MOD per intervention in % (across the corresponding
  # trials)
  min_mod_interv <- round(aggregate(unlist(arm_mod),
                                    by = list(unlist(dat$t)),
                                    min)[, 2], 2) * 100

  # Median proportion of MOD per intervention in % (across the corresponding
  # trials)
  median_mod_interv <- round(aggregate(unlist(arm_mod),
                                       by = list(unlist(dat$t)),
                                       median)[, 2], 2) * 100

  # Max proportion of MOD per intervention in % (across the corresponding
  # trials)
  max_mod_interv <- round(aggregate(unlist(arm_mod),
                                    by = list(unlist(dat$t)),
                                    max)[, 2], 2) * 100

  # Function to turn wide- to long-format for an element
  log_format <- function (input) {
    if (length(input[1, ]) > 2) {
      long_form0 <- apply(input, 1, function(x) {combn(na.omit(x), 2)})
      long_form <- t(do.call(cbind, long_form0))
    } else {
      long_form <- input
    }
    return(long_form)
  }

  # Turn into long format for MOD
  t_long_form <- log_format(dat$t)
  m_long_form <- log_format(dat$m)
  n_long_form <- log_format(dat$N)
  pair_mod <- data.frame(t_long_form, m_long_form, n_long_form)
  colnames(pair_mod) <- c("t1", "t2", "m1", "m2", "n1", "n2")

  # Name the interventions in each arm
  pair_mod[, 1] <- drug_names[pair_mod$t1]
  pair_mod[, 2] <- drug_names[pair_mod$t2]

  # The comparison between the second and first arms of each trial
  comp <- paste(pair_mod[, "t2"], "vs", pair_mod[, "t1"])

  # Number of direct comparison
  direct_comp <- length(unique(comp))

  # Total number of randomised in the network
  total_rand_network <- sum(apply(dat$N, 1, sum, na.rm = TRUE))

  # Proportion of completers in the network (in %)
  prop_obs_network <-
    round(((total_rand_network - sum(
      apply(dat$m, 1, sum, na.rm = TRUE)
    )) / total_rand_network) * 100, 0)

  # Number of randomised per observed comparison
  total_rand_partic_comp <- aggregate(apply(pair_mod[, c("n1", "n2")], 1, sum),
                                      by = list(comp), sum)[, 2]

  # Number of completers per observed comparison
  total_obs_partic_comp <- aggregate(apply(pair_mod[, c("n1", "n2")] -
                                             pair_mod[, c("m1", "m2")], 1, sum),
                                     by = list(comp), sum)[, 2]

  # Proportion of completers per observed comparison (in %)
  prop_obs_partic_comp <-
    round(total_obs_partic_comp / total_rand_partic_comp, 2) * 100

  # Proportion of MOD per observed comparison (in %)
  prop_mod_comp <- 100 - prop_obs_partic_comp

  # Proportion of MOD per trial-comparison
  trial_mod <- apply(pair_mod[, c("m1", "m2")], 1, sum) /
    apply(pair_mod[, c("n1", "n2")], 1, sum)

  # Minimum proportion of MOD per observed comparison in % (across the
  # corresponding trials)
  min_mod_comp <- round(aggregate(trial_mod,
                                  by = list(comp), min)[, 2], 2) * 100

  # Median proportion of MOD per observed comparison in % (across the
  # corresponding trials)
  median_mod_comp <- round(aggregate(trial_mod,
                                     by = list(comp), median)[, 2], 2) * 100

  # Maximum proportion of MOD per observed comparison in % (across the
  # corresponding trials)
  max_mod_comp <- round(aggregate(trial_mod,
                                  by = list(comp), max)[, 2], 2) * 100

  # Tabulate summary statistics per observed comparison: MOD
  table_comp_mod <- data.frame(as.data.frame(table(comp))[, 1],
                               as.data.frame(table(comp))[, 2],
                               prop_mod_comp,
                               min_mod_comp,
                               median_mod_comp,
                               max_mod_comp)
  colnames(table_comp_mod) <- c("Comparisons",
                                "Total trials",
                                "Missing participants (%)",
                                "Min. missing (%)",
                                "Median missing (%)",
                                "Max. missing (%)")

  # Tabulate summary statistics per intervention: MOD
  table_interv_mod <- data.frame(drug_names,
                                 as.data.frame(table(unlist(dat$t)))[, 2],
                                 prop_mod_interv,
                                 min_mod_interv,
                                 median_mod_interv,
                                 max_mod_interv)
  colnames(table_interv_mod) <- c("Interventions",
                                  "Total trials",
                                  "Missing participants (%)",
                                  "Min. missing (%)",
                                  "Median missing (%)",
                                  "Max. missing (%)")

  if (is.element(measure, c("OR", "RR", "RD"))) {
    # Proportion of observed events per trial-arm
    arm_risk <- dat$r / (dat$N - dat$m)

    # For each trial calculate the number of arms with zero events
    rule <- apply(t(apply(arm_risk, 1, function(x) {
      ifelse(x == 0, 1, 0)
      })), 1,
                  sum, na.rm = TRUE)

    # Number of trials with at least one arm with zero events
    trial_zero_event <- ifelse(length(which(rule > 0)) == 0, 0, which(rule > 0))

    # Number of trials with zero events in *all* arms
    trial_all_zero_event <- ifelse(length(which(rule == dat$na)) == 0, 0,
                                   which(rule == dat$na))

    # Proportion of total events in the network
    prop_event_network <- round(sum(unlist(dat$r), na.rm = TRUE) /
                                  sum(unlist(dat$N) - unlist(dat$m),
                                      na.rm = TRUE),
                                2) * 100

    # Number of events per intervention
    total_event_interv <- aggregate(unlist(dat$r),
                                    by = list(unlist(dat$t)), sum)[, 2]

    # Proportion of observed events per intervention (in %)
    total_risk_interv <- round(total_event_interv /
                                 total_obs_partic_interv, 2) * 100

    # Minimum proportion of observed events per intervention in % (across the
    # corresponding trials)
    min_risk_interv <- round(aggregate(unlist(arm_risk),
                                       by = list(unlist(dat$t)), min)[, 2],
                             2) * 100

    # Median proportion of observed events per intervention in % (across the
    # corresponding trials)
    median_risk_interv <- round(aggregate(unlist(arm_risk),
                                          by = list(unlist(dat$t)),
                                          median)[, 2], 2) * 100

    # Max proportion of observed events per intervention in % (across the
    # corresponding trials)
    max_risk_interv <- round(aggregate(unlist(arm_risk),
                                       by = list(unlist(dat$t)), max)[, 2],
                             2) * 100

    # Turn into long format for events
    r_long_form <- log_format(dat$r)
    pair_bin <- data.frame(t_long_form, r_long_form, n_long_form)
    colnames(pair_bin) <- c("t1", "t2", "r1", "r2", "n1", "n2")

    # Proportion of observed events per trial-comparison
    trial_risk <- apply(pair_bin[, c("r1", "r2")], 1, sum) /
      (apply(pair_bin[, c("n1", "n2")], 1, sum) -
         apply(pair_mod[, c("m1", "m2")], 1, sum))

    # Number of events per observed comparison
    total_event_comp <- aggregate(apply(pair_bin[, c("r1", "r2")], 1, sum),
                                  by = list(comp), sum)[, 2]

    # Proportion of observed events per observed comparison (in %)
    total_risk_comp <- round(total_event_comp / total_obs_partic_comp, 2) * 100

    # Minimum proportion of observed events per observed comparison in % (across
    # the corresponding trials)
    min_risk_comp <- round(aggregate(trial_risk,
                                     by = list(comp), min)[, 2], 2) * 100

    # Median proportion of observed events per observed comparison in % (across
    # the corresponding trials)
    median_risk_comp <- round(aggregate(trial_risk,
                                        by = list(comp), median)[, 2], 2) * 100

    # Maximum proportion of observed events per observed comparison in % (across
    # the corresponding trials)
    max_risk_comp <- round(aggregate(trial_risk,
                                     by = list(comp), max)[, 2], 2) * 100
  }

  # Tabulate summary statistics per intervention
  table_interv <- data.frame(drug_names,
                             as.data.frame(table(unlist(dat$t)))[, 2],
                             total_rand_partic_interv,
                             prop_obs_partic_interv)
  colnames(table_interv) <- c("Interventions",
                              "Total trials",
                              "Total randomised",
                              "Completers (%)")

  if (is.element(measure, c("OR", "RR", "RD"))) {
    table_interv <- cbind(table_interv, data.frame(total_risk_interv,
                                                   min_risk_interv,
                                                   median_risk_interv,
                                                   max_risk_interv))
    colnames(table_interv)[5:8] <- c("Total events (%)",
                                     "Min. events (%)",
                                     "Median events (%)",
                                     "Max. events (%)")
  }

  # Tabulate summary statistics per observed comparison
  table_comp <- data.frame(as.data.frame(table(comp))[, 1],
                           as.data.frame(table(comp))[, 2],
                           total_rand_partic_comp,
                           prop_obs_partic_comp)
  colnames(table_comp) <- c("Comparisons",
                            "Total trials",
                            "Total randomised",
                            "Completers (%)")

  if (is.element(measure, c("OR", "RR", "RD"))) {
    table_comp <- cbind(table_comp, data.frame(total_risk_comp,
                                               min_risk_comp,
                                               median_risk_comp,
                                               max_risk_comp))
    colnames(table_comp)[5:8] <- c("Total events (%)",
                                   "Min. events (%)",
                                   "Median events (%)",
                                   "Max. events (%)")
  }

  # Tabulate general information of the network
  characteristics <- c("Interventions",
                       "Possible comparisons",
                       "Direct comparisons",
                       "Indirect comparisons",
                       "Trials",
                       "Two-arm trials",
                       "Multi-arm trials",
                       "Randomised participants",
                       "Proportion of completers")
  value <- c(length(drug_names),
             dim(combn(length(drug_names), 2))[2],
             direct_comp,
             dim(combn(length(drug_names), 2))[2] - direct_comp,
             dim(data)[1],
             two_arm_ns,
             multi_arm_ns,
             total_rand_network,
             prop_obs_network)

  results <- data.frame(characteristics, value)
  colnames(results) <- c("Characteristic", "Total")

  # Add further information based on the effect measure
  if (measure == "OR") {
    results[10, ] <- rbind("Proportion of observed events",
                           prop_event_network)
    results[11, ] <- rbind("Trials with at least one zero event",
                           trial_zero_event)
    results[12, ] <- rbind("Trials with all zero events",
                           trial_all_zero_event)
  } else {
    results
  }

  # Write the tables as .xlsx
  if (save_xls == TRUE) {
    writexl::write_xlsx(table_interv, "table_interventions.xlsx")
    writexl::write_xlsx(table_comp, "table_comparisons.xlsx")
  }

  # Obtain the element 'm_pseudo'
  na_missing <- data_preparation(data = data, measure = measure)$m_pseudo
  na_missing_trials <- length(which(unlist(na_missing) == -1))

  # Return all the results
  return(list(network_description =
                knitr::kable(results,
                             align = "ll",
                             caption = "Description of the network"),
              table_interventions =
                knitr::kable(table_interv,
                             align = "lccccccc",
                             caption = "Interventions"),
              table_comparisons =
                knitr::kable(table_comp,
                             align = "lccccccc",
                             caption = "Observed comparisons")))

  # Whether there are trials without information on missing participants
  if (na_missing_trials > 0 & na_missing_trials < sum(dat$na)) {
    aa <- "trial-arms without information on"
    bb <- "the number of missing participants."
    message(paste("Note: There are", na_missing_trials, aa, paste0(bb, ".")))
  }
}

#' A function to describe the evidence base
#'
#' @description Calculates the necessary elements to describe the evidence base
#'   for an outcome across the network, the interventions, and observed
#'   comparisons. See also 'Value' in \code{\link{netplot}}.
#'
#' @param data A data-frame of a one-trial-per-row format containing arm-level
#'   data of each trial. See 'Format' in \code{\link{run_model}}.
#' @param drug_names A vector of labels with the name of the interventions in
#'   the order they appear in the argument \code{data}. If \code{drug_names} is
#'   not defined, the order of the interventions as they appear in \code{data}
#'   is used, instead.
#' @param measure Character string indicating the effect measure with values
#'   \code{"OR"}, \code{"MD"}, \code{"SMD"}, or \code{"ROM"} for the odds ratio,
#'   mean difference, standardised mean difference and ratio of means,
#'   respectively.
#'
#' @return A list of scalar results and two data-frames to be passed to
#'   \code{\link{netplot}}. The scalar results include:
#'   \tabular{ll}{
#'    \code{direct_comp} \tab The number of observed comparisons in the
#'    network. \cr
#'    \tab \cr
#'    \code{two_arm_ns} \tab The number of two-arm trials in the network. \cr
#'    \tab \cr
#'    \code{multi_arm_ns} \tab The number of multi-arm trials in the
#'    network. \cr
#'    \tab \cr
#'    \code{total_rand_network} \tab The total number of randomised participants
#'     in the network. \cr
#'    \tab \cr
#'    \code{prop_obs_network} \tab The proportion of participants who completed
#'     the trial. \cr
#'    \tab \cr
#'    \code{prop_event_network} \tab The proportion of observed events in the
#'    network. When the outcome is continuous, this element is omitted. \cr
#'    \tab \cr
#'    \code{trial_zero_event} \tab The number of trials with at least one arm
#'    with zero events. When the outcome is continuous, this element is
#'    omitted. \cr
#'    \tab \cr
#'    \code{trial_all_zero_event} \tab The number of trials with zero events in
#'    all arms. When the outcome is continuous, this element is omitted. \cr
#'   }
#'
#'   The two data-frames include\code{table_interventions} and
#'   \code{table_comparisons}. See 'Value' in \code{\link{netplot}} for these
#'   data-frames.
#'
#' @details \code{describe_network} calls \code{\link{data_preparation}} to
#' facilitate the calculations.
#'
#' @seealso \code{\link{data_preparation}}, \code{\link{netplot}},
#'   \code{\link{run_model}},
#'
#' @author {Loukia M. Spineli}
#'
#' @export
describe_network <- function(data, drug_names, measure) {

  options(warn = -1)

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

  # Turn into long format using the 'pairwise' function (netmeta): MOD
  pair_mod0 <- pairwise(as.list(dat$t),
                        event = as.list(dat$m),
                        n = as.list(dat$N),
                        data = cbind(dat$t, dat$m, dat$N),
                        studlab = 1:dat$ns)[, c(3:6, 8, 7, 9)]
  colnames(pair_mod0) <- c("study", "t1", "t2", "m1", "m2", "n1", "n2")

  # Ensure that t1 < t2 and correspondingly for the other elements
  treat <- treat0 <- pair_mod0[, 2:3]
  miss <- miss0 <- pair_mod0[, 4:5]
  rand <- rand0 <- pair_mod0[, 6:7]
  for (i in seq_len(length(pair_mod0[, 1]))) {
    treat[i, ] <- treat0[i, order(treat0[i, ], na.last = TRUE)]
    miss[i, ] <- miss0[i, order(treat0[i, ], na.last = TRUE)]
    rand[i, ] <- rand0[i, order(treat0[i, ], na.last = TRUE)]
  }
  pair_mod <- data.frame(study = pair_mod0$study, treat, miss, rand)

  # Name the interventions in each arm
  pair_mod[, 2] <- drug_names[pair_mod$t1]
  pair_mod[, 3] <- drug_names[pair_mod$t2]

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

  if (measure == "OR") {
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

    # Turn into long format using the 'pairwise' function (netmeta)
    pair_bin0 <- pairwise(as.list(dat$t),
                          event = as.list(dat$r),
                          n = as.list(dat$N),
                          data = cbind(dat$t, dat$r, dat$N),
                          studlab = 1:dat$ns)[, c(3:6, 8, 7, 9)]
    colnames(pair_bin0) <- c("study", "t1", "t2", "r1", "r2", "n1", "n2")

    # Ensure that t1 < t2 and correspondingly for the other elements
    treat <- treat0 <- pair_bin0[, 2:3]
    resp <- resp0 <- pair_bin0[, 4:5]
    rand <- rand0 <- pair_bin0[, 6:7]
    for (i in seq_len(length(pair_bin0[, 1]))) {
      treat[i, ] <- treat0[i, order(treat0[i, ], na.last = TRUE)]
      resp[i, ] <- resp0[i, order(treat0[i, ], na.last = TRUE)]
      rand[i, ] <- rand0[i, order(treat0[i, ], na.last = TRUE)]
    }
    pair_bin <- data.frame(study = pair_bin0$study, treat, resp, rand)

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

  if (measure == "OR") {
    # Tabulate summary statistics per intervention
    table_interv <- data.frame(drug_names,
                               as.data.frame(table(unlist(dat$t)))[, 2],
                               total_rand_partic_interv,
                               prop_obs_partic_interv,
                               total_risk_interv,
                               min_risk_interv,
                               median_risk_interv,
                               max_risk_interv)
    colnames(table_interv) <- c("Interventions",
                                "Total trials",
                                "Total randomised",
                                "Completers (%)",
                                "Total events (%)",
                                "Min. events (%)",
                                "Median events (%)",
                                "Max. events (%)")

    # Tabulate summary statistics per observed comparison
    table_comp <- data.frame(as.data.frame(table(comp))[, 1],
                             as.data.frame(table(comp))[, 2],
                             total_rand_partic_comp,
                             prop_obs_partic_comp,
                             total_risk_comp,
                             min_risk_comp,
                             median_risk_comp,
                             max_risk_comp)
    colnames(table_comp) <- c("Comparisons",
                              "Total trials",
                              "Total randomised",
                              "Completers (%)",
                              "Total events (%)",
                              "Min. events (%)",
                              "Median events (%)",
                              "Max. events (%)")
  } else {
    # Tabulate summary statistics per intervention
    table_interv <- data.frame(drug_names,
                               as.data.frame(table(unlist(dat$t)))[, 2],
                               total_rand_partic_interv,
                               prop_obs_partic_interv)
    colnames(table_interv) <- c("Interventions",
                                "Total trials",
                                "Total randomised",
                                "Completers (%)")

    # Tabulate summary statistics per observed comparison
    table_comp <- data.frame(as.data.frame(table(comp))[, 1],
                             as.data.frame(table(comp))[, 2],
                             total_rand_partic_comp,
                             prop_obs_partic_comp)
    colnames(table_comp) <- c("Comparisons",
                              "Total trials",
                              "Total randomised",
                              "Completers (%)")
  }

  results <- list(direct_comp = direct_comp,
                  two_arm_ns = two_arm_ns,
                  multi_arm_ns = multi_arm_ns,
                  total_rand_network = total_rand_network,
                  prop_obs_network = prop_obs_network,
                  table_interventions = table_interv,
                  table_comparisons = table_comp)


  if (measure == "OR") {
    results <- append(results,
                      list(prop_event_network = prop_event_network,
                           trial_zero_event = trial_zero_event,
                           trial_all_zero_event = trial_all_zero_event))
  } else {
    results
  }

  return(results)
}

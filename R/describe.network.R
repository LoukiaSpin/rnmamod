#' A function to describe a network of interventions
#'
#' @description This function calculates the necessary elements to describe a network, such as
#'   the number of interventions, trials, randomised participants, and so on. Furthermore, this function provides
#'   summary statistics of the missing participants and the analysed outcome per intervention and observed comparison in a tabulated format.
#'   See 'Value' in the \code{\link[rnmamod]{network.plot}} function for more details.
#'
#' @param data A data-frame of a one-trial-per-row format containing arm-level data of each trial. This format is widely used for BUGS models.
#'   See 'Format' in \code{\link[rnmamod]{run.model}} function for the specification of the columns.
#' @param drug.names A vector of labels with the name of the interventions in the order they appear in the argument \code{data}. If the argument \code{drug.names} is not defined, the order of the interventions
#'   as they appear in \code{data} is used, instead.
#' @param measure Character string indicating the effect measure with values \code{"OR"}, \code{"MD"}, \code{"SMD"}, or \code{"ROM"} for the odds ratio, mean difference,
#'   standardised mean difference and ratio of means, respectively.
#'
#' @return A list of scalar results and four data-frames to be passed to \code{\link{network.plot}}. The scalar results include:
#'   \tabular{ll}{
#'    \code{direct.comp} \tab The number of observed comparisons in the network. \cr
#'    \tab \cr
#'    \code{two.arm.ns} \tab The number of two-arm trials in the network. \cr
#'    \tab \cr
#'    \code{multi.arm.ns} \tab The number of multi-arm trials in the network. \cr
#'    \tab \cr
#'    \code{total.rand.network} \tab The total number of randomised participants in the network. \cr
#'    \tab \cr
#'    \code{prop.obs.network} \tab The proportion of participants who completed the trial. \cr
#'    \tab \cr
#'    \code{prop.event.network} \tab The proportion of observed events in the network. When the outcome is continuous, this element is omitted. \cr
#'    \tab \cr
#'    \code{trial.zero.event} \tab The number of trials with at least one arm with zero events. When the outcome is continuous, this element is omitted. \cr
#'    \tab \cr
#'    \code{trial.all.zero.event} \tab The number of trials with zero events in all arms. When the outcome is continuous, this element is omitted. \cr
#'   }
#'
#'   The four data-frames include \code{Table.interventions.Missing}, \code{Table.comparisons.Missing}, \code{Table.interventions} and \code{Table.comparisons}.
#'   See 'Value' in \code{\link[rnmamod]{network.plot}} that describes these data-frames in detail.
#'
#' @details \code{describe.network} calls the \code{data.preparation} function to facilitate the calculations.
#'
#' @seealso \code{\link[rnmamod]{network.plot}}, \code{\link[rnmamod]{run.model}}, \code{\link[rnmamod]{data.preparation}}
#'
#' @author {Loukia M. Spineli}
#'
#' @export
describe.network <- function(data, drug.names, measure) {

  options(warn = -1)

  # Use the 'data.preparation' function
  dat <- data.preparation(data, measure)

  # Number of two-arm trian
  two.arm.ns <- length(which(dat$na < 3))

  # Number of multi-arm trials
  multi.arm.ns <- length(which(dat$na > 2))

  # Number of randomised per intervention
  total.rand.partic.interv <- aggregate(unlist(dat$N), by = list(unlist(dat$t)), sum)[, 2]

  # Number of completers per intervention
  total.obs.partic.interv <- aggregate(unlist(dat$N - dat$m), by = list(unlist(dat$t)), sum)[, 2]

  # Proportion of completers per intervention (in %)
  prop.obs.partic.interv <- round(total.obs.partic.interv/total.rand.partic.interv, 2)*100

  # Proportion of missing outcome data (MOD) per intervention (in %)
  prop.mod.interv <- 100 - prop.obs.partic.interv

  # Proportion of MOD per trial-arm
  arm.mod <- dat$m/dat$N

  # Minimum proportion of MOD per intervention in % (across the corresponding trials)
  min.mod.interv <- round(aggregate(unlist(arm.mod), by = list(unlist(dat$t)), min)[, 2], 2)*100

  # Median proportion of MOD per intervention in % (across the corresponding trials)
  median.mod.interv <- round(aggregate(unlist(arm.mod), by = list(unlist(dat$t)), median)[, 2], 2)*100

  # Max proportion of MOD per intervention in % (across the corresponding trials)
  max.mod.interv <- round(aggregate(unlist(arm.mod), by = list(unlist(dat$t)), max)[, 2], 2)*100

  # Turn into long format using the 'pairwise' function (netmeta): MOD
  pair.mod <- pairwise(as.list(dat$t),
                       event = as.list(dat$m),
                       n = as.list(dat$N),
                       data = cbind(dat$t, dat$m, dat$N),
                       studlab = 1:dat$ns)[, c(3:6, 8, 7, 9)]
  colnames(pair.mod) <- c("study", "t1", "t2", "m1", "m2", "n1", "n2")

  # Name the interventions in each arm
  pair.mod[, 2] <- drug.names[pair.mod$t1]
  pair.mod[, 3] <- drug.names[pair.mod$t2]

  # The comparison between the second and first arms of each trial
  comp <- paste(pair.mod[, "t2"], "vs", pair.mod[, "t1"])

  # Number of direct comparison
  direct.comp <- length(unique(comp))

  # Total number of randomised in the network
  total.rand.network <- sum(pair.mod[, c("n1", "n2")])

  # Proportion of completers in the network (in %)
  prop.obs.network <- round((sum(pair.mod[, c("n1", "n2")] - pair.mod[, c("m1", "m2")])/total.rand.network), 2)*100

  # Proportion of missing outcome data (MOD) in the network (in %)
  prop.mod.network <- 100 - prop.obs.network

  # Number of randomised per observed comparison
  total.rand.partic.comp <- aggregate(apply(pair.mod[, c("n1", "n2")], 1, sum), by = list(comp), sum)[, 2]

  # Number of completers per observed comparison
  total.obs.partic.comp <- aggregate(apply(pair.mod[, c("n1", "n2")] - pair.mod[, c("m1", "m2")], 1, sum), by = list(comp), sum)[, 2]

  # Proportion of completers per observed comparison (in %)
  prop.obs.partic.comp <- round(total.obs.partic.comp/total.rand.partic.comp, 2)*100

  # Proportion of MOD per observed comparison (in %)
  prop.mod.comp <- 100 - prop.obs.partic.comp

  # Proportion of MOD per trial-comparison
  trial.mod <- apply(pair.mod[, c("m1", "m2")], 1, sum)/apply(pair.mod[, c("n1", "n2")], 1, sum)

  # Minimum proportion of MOD per observed comparison in % (across the corresponding trials)
  min.mod.comp <- round(aggregate(trial.mod, by = list(comp), min)[, 2], 2)*100

  # Median proportion of MOD per observed comparison in % (across the corresponding trials)
  median.mod.comp <- round(aggregate(trial.mod, by = list(comp), median)[, 2], 2)*100

  # Maximum proportion of MOD per observed comparison in % (across the corresponding trials)
  max.mod.comp <- round(aggregate(trial.mod, by = list(comp), max)[, 2], 2)*100

  # Tabulate summary statistics per observed comparison: MOD
  table.comp.mod <- data.frame(as.data.frame(table(comp))[, 1],
                               as.data.frame(table(comp))[, 2],
                               total.rand.partic.comp,
                               prop.obs.partic.comp,
                               prop.mod.comp,
                               min.mod.comp,
                               median.mod.comp,
                               max.mod.comp)
  colnames(table.comp.mod) <- c("Comparisons",
                                "Total trials",
                                "Total randomised",
                                "Completers (%)",
                                "Missing participants (%)",
                                "Min. missing (%)",
                                "Median missing (%)",
                                "Max. missing (%)")

  # Tabulate summary statistics per intervention: MOD
  table.interv.mod <- data.frame(drug.names,
                                 as.data.frame(table(unlist(dat$t)))[, 2],
                                 total.rand.partic.interv,
                                 prop.obs.partic.interv,
                                 prop.mod.interv,
                                 min.mod.interv,
                                 median.mod.interv,
                                 max.mod.interv)
  colnames(table.interv.mod) <- c("Interventions",
                                  "Total trials",
                                  "Total randomised",
                                  "Completers (%)",
                                  "Missing participants (%)",
                                  "Min. missing (%)",
                                  "Median missing (%)",
                                  "Max. missing (%)")

  if (measure == "OR") {
    # Proportion of observed events per trial-arm
    arm.risk <- dat$r/(dat$N - dat$m)

    # For each trial calculate the number of arms with zero events
    rule <- apply(ifelse(arm.risk == 0.0, 1, 0), 1, sum, na.rm = T)

    # Number of trials with at least one arm with zero events
    trial.zero.event <- ifelse(length(which(rule > 0)) == 0, 0, which(rule > 0))

    # Number of trials with zero events in *all* arms
    trial.all.zero.event <- ifelse(length(which(is.element(rule, dat$na) == T)) == 0, 0, which(is.element(rule, dat$na) == T))

    # Total number of events in the network
    prop.event.network <- round(sum(unlist(dat$r), na.rm = T)/sum(unlist(dat$N) - unlist(dat$m), na.rm = T), 2)*100

    # Number of events per intervention
    total.event.interv <- aggregate(unlist(dat$r), by = list(unlist(dat$t)), sum)[, 2]

    # Proportion of observed events per intervention (in %)
    total.risk.interv <- round(total.event.interv/total.obs.partic.interv, 2)*100

    # Minimum proportion of observed events per intervention in % (across the corresponding trials)
    min.risk.interv <- round(aggregate(unlist(arm.risk), by = list(unlist(dat$t)), min)[, 2], 2)*100

    # Median proportion of observed events per intervention in % (across the corresponding trials)
    median.risk.interv <- round(aggregate(unlist(arm.risk), by = list(unlist(dat$t)), median)[, 2], 2)*100

    # Max proportion of observed events per intervention in % (across the corresponding trials)
    max.risk.interv <- round(aggregate(unlist(arm.risk), by = list(unlist(dat$t)), max)[, 2], 2)*100

    # Tabulate summary statistics per intervention
    table.interv.bin <- data.frame(drug.names,
                                   as.data.frame(table(unlist(dat$t)))[, 2],
                                   total.rand.partic.interv,
                                   prop.obs.partic.interv,
                                   prop.mod.interv,
                                   total.risk.interv,
                                   min.risk.interv,
                                   median.risk.interv,
                                   max.risk.interv)
    colnames(table.interv.bin) <- c("Interventions",
                                    "Total trials",
                                    "Total randomised",
                                    "Completers (%)",
                                    "Missing participants (%)",
                                    "Total events (%)",
                                    "Min. events (%)",
                                    "Median events (%)",
                                    "Max. events (%)")

    # Turn into long format using the 'pairwise' function (netmeta): binary outcome
    pair.bin <- pairwise(as.list(dat$t),
                         event = as.list(dat$r),
                         n = as.list(dat$N),
                         data = cbind(dat$t, dat$r, dat$N),
                         studlab = 1:dat$ns)[, c(3:6, 8, 7, 9)]
    colnames(pair.bin) <- c("study", "t1", "t2", "r1", "r2", "n1", "n2")

    # Proportion of observed events per trial-comparison
    trial.risk <- apply(pair.bin[, c("r1", "r2")], 1, sum)/(apply(pair.bin[, c("n1", "n2")], 1, sum) - apply(pair.mod[, c("m1", "m2")], 1, sum))

    # Number of events per observed comparison
    total.event.comp <- aggregate(apply(pair.bin[, c("r1", "r2")], 1, sum), by = list(comp), sum)[, 2]

    # Proportion of observed events per observed comparison (in %)
    total.risk.comp <- round(total.event.comp/total.obs.partic.comp, 2)*100

    # Minimum proportion of observed events per observed comparison in % (across the corresponding trials)
    min.risk.comp <- round(aggregate(trial.risk, by = list(comp), min)[, 2], 2)*100

    # Median proportion of observed events per observed comparison in % (across the corresponding trials)
    median.risk.comp <- round(aggregate(trial.risk, by = list(comp), median)[, 2], 2)*100

    # Maximum proportion of observed events per observed comparison in % (across the corresponding trials)
    max.risk.comp <- round(aggregate(trial.risk, by = list(comp), max)[, 2], 2)*100

    # Tabulate summary statistics per observed comparison
    table.comp.bin <- data.frame(as.data.frame(table(comp))[, 1],
                                 as.data.frame(table(comp))[, 2],
                                 total.rand.partic.comp,
                                 prop.obs.partic.comp,
                                 prop.mod.comp,
                                 total.risk.comp,
                                 min.risk.comp,
                                 median.risk.comp,
                                 max.risk.comp)
    colnames(table.comp.bin) <- c("Comparisons",
                                  "Total trials",
                                  "Total randomised",
                                  "Completers (%)",
                                  "Missing participants (%)",
                                  "Total events (%)",
                                  "Min. events (%)",
                                  "Median events (%)",
                                  "Max. events (%)")
  } else {
    # Minimum t-statistic per intervention (across the corresponding trials)
    min.t.interv <- round(aggregate(unlist(dat$y0)/unlist(dat$se0), by = list(unlist(dat$t)), min)[, 2], 2)

    # Median t-statistic per intervention (across the corresponding trials)
    median.t.interv <- round(aggregate(unlist(dat$y0)/unlist(dat$se0), by = list(unlist(dat$t)), median)[, 2], 2)

    # Maximum t-statistic per intervention (across the corresponding trials)
    max.t.interv <- round(aggregate(unlist(dat$y0)/unlist(dat$se0), by = list(unlist(dat$t)), max)[, 2], 2)

    # Tabulate summary statistics per intervention
    table.interv.con <- data.frame(drug.names,
                                   as.data.frame(table(unlist(dat$t)))[, 2],
                                   total.rand.partic.interv,
                                   prop.obs.partic.interv,
                                   prop.mod.interv,
                                   min.t.interv,
                                   median.t.interv,
                                   max.t.interv)
    colnames(table.interv.con) <- c("Interventions",
                                    "Total trials",
                                    "Total randomised",
                                    "Completers (%)",
                                    "Missing participants (%)",
                                    "Min. t-statistic",
                                    "Median t-statistic",
                                    "Max. t-statistic")

    # Turn into long format using the 'pairwise' function (netmeta): continuous outcome
    pair.con <- pairwise(as.list(dat$t),
                         n = as.list(dat$N),
                         mean = as.list(dat$y0),
                         sd = as.list(dat$se0),
                         data = cbind(dat$t, dat$N, dat$y0, dat$se0),
                         studlab = 1:dat$ns)[, c(3:5, 7, 10, 8, 11, 6, 9)]
    colnames(pair.con) <- c("study", "t1", "t2", "y1", "y2", "se1", "se2","n1", "n2")

    # t-statistic (t2 versus t1) per trial-comparison
    t.stat.trial <- (pair.con$y2 - pair.con$y1)/sqrt((pair.con$se2)^2 + (pair.con$se1)^2)

    # Minimum t-statistic per observed comparison (across the corresponding trials)
    min.t.comp <- round(aggregate(t.stat.trial, by = list(comp), min)[, 2], 2)

    # Median t-statistic per observed comparison (across the corresponding trials)
    median.t.comp <- round(aggregate(t.stat.trial, by = list(comp), median)[, 2], 2)

    # Maximum t-statistic per observed comparison (across the corresponding trials)
    max.t.comp <- round(aggregate(t.stat.trial, by = list(comp), max)[, 2], 2)

    # Tabulate summary statistics per observed comparison
    table.comp.con <- data.frame(as.data.frame(table(comp))[, 1],
                                 as.data.frame(table(comp))[, 2],
                                 total.rand.partic.comp,
                                 prop.obs.partic.comp,
                                 prop.mod.comp,
                                 min.t.comp,
                                 median.t.comp,
                                 max.t.comp)
    colnames(table.comp.con) <- c("Comparisons",
                                  "Total trials",
                                  "Total randomised",
                                  "Completers (%)",
                                  "Missing participants (%)",
                                  "Min. t-statistic",
                                  "Median t-statistic",
                                  "Max. t-statistic")
  }


  results <- list(direct.comp = direct.comp,
                  two.arm.ns = two.arm.ns,
                  multi.arm.ns = multi.arm.ns,
                  total.rand.network = total.rand.network,
                  prop.obs.network = prop.obs.network,
                  Table.interventions.Missing = knitr::kable(table.interv.mod),
                  Table.comparisons.Missing = knitr::kable(table.comp.mod))

  results <- if (measure == "OR") {
    append(results, list(prop.event.network = prop.event.network,
                         trial.zero.event = trial.zero.event,
                         trial.all.zero.event = trial.all.zero.event,
                         Table.interventions = knitr::kable(table.interv.bin),
                         Table.comparisons = knitr::kable(table.comp.bin)))
  } else {
    append(results, list(Table.interventions = knitr::kable(table.interv.con),
                         Table.comparisons = knitr::kable(table.comp.con)))
  }

  return(results)

}

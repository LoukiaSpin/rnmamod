#' Network plot and description of the evidence base
#'
#' @description Illustrates the network plot for one outcome and summarises the
#'   characteristics of the evidence base.
#'
#' @param data A data-frame of a one-trial-per-row format containing arm-level
#'   data of each trial. See 'Format' in \code{\link{run_model}}.
#' @param drug_names A vector of labels with the name of the interventions in
#'   the order they appear in the argument \code{data}.
#' @param save_xls Logical to indicate whether to export the tabulated results
#'   to an 'xlsx' file (via the \code{\link[writexl:write_xlsx]{write_xlsx}}
#'   function of the R-package
#'   \href{https://CRAN.R-project.org/package=writexl}{writexl}) at the working
#'   directory of the user. The default is \code{FALSE} (do not export).
#' @param ... Additional arguments of the
#'   \code{\link[pcnetmeta:nma.networkplot]{nma.networkplot}} function of the
#'   R-package \href{https://CRAN.R-project.org/package=pcnetmeta}{pcnetmeta}.
#'
#' @return A network plot with coloured closed-loops informed by multi-arm
#'   trials. Each node indicates an intervention and each edge an observed
#'   pairwise comparison. The edge thickness is proportional to the number of
#'   trials investigating the corresponding comparison, unless
#'   specified otherwise (see
#'   \code{\link[pcnetmeta:nma.networkplot]{nma.networkplot}} function of the
#'   R-package \href{https://CRAN.R-project.org/package=pcnetmeta}{pcnetmeta}).
#'   The size of the node is weighted by the total number of
#'   trials investigating the corresponding intervention, unless specified
#'   otherwise (see \code{\link[pcnetmeta:nma.networkplot]{nma.networkplot}}
#'   function of the R-package
#'   \href{https://CRAN.R-project.org/package=pcnetmeta}{pcnetmeta}).
#'
#'   \code{netplot} also returns the following data-frames that describe the
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
#' @details \code{netplot} draws the network plot using the
#'   \code{\link[pcnetmeta:nma.networkplot]{nma.networkplot}} function
#'   of the R-package
#'   \href{https://CRAN.R-project.org/package=pcnetmeta}{pcnetmeta}.
#'   The \code{\link[gemtc:mtc.data.studyrow]{mtc.data.studyrow}} function of
#'   the R-package href{https://CRAN.R-project.org/package=gemtc}{gemtc} is
#'   additionally used to convert \code{data} from the required
#'   one-trial-per-row format into the one-arm-per-row format.
#'
#'   Furthermore, \code{netplot} exports the data-frames to separate
#'   'xlsx' files (via the \code{\link[writexl:write_xlsx]{write_xlsx}} function
#'   of the R-package
#'   \href{https://CRAN.R-project.org/package=writexl}{writexl}) at the working
#'   directory of the user.
#'
#' @seealso \code{\link{data_preparation}},
#'   \code{\link[gemtc:mtc.data.studyrow]{mtc.data.studyrow}},
#'   \code{\link[pcnetmeta:nma.networkplot]{nma.networkplot}},
#'   \code{\link{run_model}}
#'   \code{\link[writexl:write_xlsx]{write_xlsx}}
#'
#' @references
#' Lin L, Zhang J, Hodges JS, Chu H. Performing Arm-Based Network Meta-Analysis
#' in R with the pcnetmeta Package.
#' \emph{J Stat Softw} 2017;\bold{80}(5): 1--25. doi: 10.18637/jss.v080.i05
#'
#' van Valkenhoef G, Kuiper J. gemtc: Network Meta-Analysis Using
#' Bayesian Methods. \emph{R package version 1.0-1}. 2021.
#' \url{https://CRAN.R-project.org/package=gemtc}
#'
#' @author {Loukia M. Spineli}
#'
#' @examples
#' data("nma.bottomley2011")
#'
#' # Return the first six trials of the dataset
#' head(nma.bottomley2011)
#'
#' # The names of the interventions in the order they appear in the dataset
#' interv_names <- c("betamethasone dipropionate", "betamethasone valerate",
#'                   "calcipotriol", "calcipotriol plus polytar", "capasal",
#'                   "two-compound formulation gel", "placebo")
#'
#' # Create the network plot
#' netplot(data = nma.bottomley2011,
#'         drug_names = interv_names,
#'         save_xls = FALSE)
#'
#' @export
netplot <- function(data, drug_names, save_xls, ...) {

  save_xls <- if (missing(save_xls)) {
    FALSE
  } else {
    save_xls
  }

  # Obtain dataset
  r <- data %>% dplyr::select(starts_with("r") | starts_with("y"))
  n <- data %>% dplyr::select(starts_with("n"))
  t <- data %>% dplyr::select(starts_with("t"))
  nt <- length(table(as.matrix(t)))
  ns <- length(r[, 1])
  na.. <- rep(0, length(r[, 1]))
  for (i in seq_len(length(r[, 1]))) {
    na..[i] <- table(!is.na(t[i, ]))["TRUE"]
  }

  drug_names <- if (missing(drug_names)) {
    stop("The argument 'drug_names' has not been defined.", call. = FALSE)
  } else {
    drug_names
  }

  # Rename columns to agree with gemtc
  names(r) <- paste0("r..", seq_len(max(na..)), ".")
  names(n) <- paste0("n..", seq_len(max(na..)), ".")
  names(t) <- paste0("t..", seq_len(max(na..)), ".")

  # One row per study arm
  transform <- mtc.data.studyrow(cbind(t, r, n, na..),
                                 armVars = c("treatment" = "t",
                                             "response" = "r",
                                             "sampleSize" = "n"),
                                 nArmsVar = "na")
  transform$treatment1 <- as.numeric(as.character(transform$treatment))

  # Obtain network plot
  network_plot <- if (length(drug_names) > 2) {
    nma.networkplot(study,
                    treatment1,
                    data = transform,
                    trtname = drug_names,
                    multi.show = TRUE, ...)
  } else {
    NA
  }

  if (dim(data %>% dplyr::select(starts_with("r")))[2] > 0) {
    measure <- "OR"
  } else {
    measure <- "MD"
  }

  dat <- suppressMessages({
    describe_network(data, drug_names, measure)
    })
  characteristics <- c("Interventions",
                       "Possible comparisons",
                       "Direct comparisons",
                       "Indirect comparisons",
                       "Trials",
                       "Two-arm trials",
                       "Multi-arm trials",
                       "Randomised participants",
                       "Proportion of completers")
  value <- c(nt,
             dim(combn(nt, 2))[2],
             dat$direct_comp,
             dim(combn(nt, 2))[2] - dat$direct_comp,
             ns,
             dat$two_arm_ns,
             dat$multi_arm_ns,
             dat$total_rand_network,
             dat$prop_obs_network)

  results <- data.frame(characteristics, value)
  colnames(results) <- c("Characteristic", "Total")

  if (measure == "OR") {
    results[10, ] <- rbind("Proportion of observed events",
                           dat$prop_event_network)
    results[11, ] <- rbind("Trials with at least one zero event",
                           dat$trial_zero_event)
    results[12, ] <- rbind("Trials with all zero events",
                           dat$trial_all_zero_event)
  } else {
    results
  }

  # Write the tables as .xlsx
  if (save_xls == TRUE) {
    writexl::write_xlsx(dat$table_interventions, "table_interventions.xlsx")
    writexl::write_xlsx(dat$table_comparisons, "table_comparisons.xlsx")
  }

  # Obtain the element 'm_pseudo'
  na_missing <- data_preparation(data = data, measure = measure)$m_pseudo
  na_missing_trials <- length(which(unlist(na_missing) == -1))

  # Return all the results
  return(list(network_plot = network_plot,
              network_description =
                knitr::kable(results,
                             align = "ll",
                             caption = "Description of the network"),
              table_interventions =
                knitr::kable(dat$table_interventions,
                             align = "lccccccc",
                             caption = "Interventions"),
              table_comparisons =
                knitr::kable(dat$table_comparisons,
                           align = "lccccccc",
                           caption = "Observed comparisons")))

  # Whether there are trials without information on missing participants
  if (na_missing_trials > 0 & na_missing_trials < sum(na..)) {
    aa <- "trial-arms without information on"
    bb <- "the number of missing participants."
    message(paste("Note: There are", na_missing_trials, aa, paste0(bb, ".")))
  }
}

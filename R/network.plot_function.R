#' Create the network plot
#'
#' @description Draws the network plot using the
#'   \code{\link[pcnetmeta]{nma.networkplot}} function
#'   of the R-package
#'   \href{https://CRAN.R-project.org/package=pcnetmeta}{pcnetmeta}.
#'   The \code{\link[gemtc]{mtc.data.studyrow}} function of the R-package
#'   \href{https://CRAN.R-project.org/package=gemtc}{gemtc} is additionally used
#'   to convert \code{data} from the required one-trial-per-row format into the
#'   one-arm-per-row format.
#'
#' @param data A data-frame of a one-trial-per-row format containing arm-level
#'   data of each trial. See 'Format' in \code{\link[rnmamod]{run_model}}
#'   function for the specification of the columns.
#' @param drug_names A vector of labels with the name of the interventions in
#'   the order they appear in the argument \code{data}. If the argument
#'   \code{drug_names} is not defined, the interventions are ordered as they
#'   appear in \code{data}.
#' @param save_xls Logical to indicate whether to export the tabulated results
#'   in Excel 'xlsx' format at the working directory of the user.
#'   The default is \code{FALSE} (do not export in Excel format).
#' @param ... Additional arguments of the
#'   \code{\link[pcnetmeta]{nma.networkplot}} function of the R-package
#'   \href{https://CRAN.R-project.org/package=pcnetmeta}{pcnetmeta}.
#'
#' @return A network plot with coloured closed-loops informed by multi-arm
#'   trials. Each node indicates an intervention and each link an observed
#'   pairwise comparison. The edge thickness is proportional to the number of
#'   direct treatment comparisons, unless specified otherwise
#'   (see \code{\link[pcnetmeta]{nma.networkplot}} function). The size of the
#'   node is weighted by the total number of direct treatment comparisons of the
#'   corresponding treatment, unless specified otherwise
#'   (see \code{\link[pcnetmeta]{nma.networkplot}} function).
#'
#'   \code{netplot} also returns five data-frames that describe characteristics
#'   of the network:
#'   \tabular{ll}{
#'    \code{network_description} \tab The number of: interventions, possible
#'    comparisons, direct and indirect comparisons, number of trials in total,
#'    number of two-arm and multi-arm trials, number of randomised participants,
#'    proportion of participants completing the trial (completers), and
#'    proportion of missing participants. When the outcome is binary, the number
#'    of trials with at least one zero event, and the number of trials with all
#'    zero events are also presented. \cr
#'    \tab \cr
#'    \code{table_interventions} \tab For each intervention, the
#'    number of trials, number of randomised participants, proportion of
#'    completers and of missing participants. When the outcome is binary, the
#'    data-frame presents also the proportion of total observed events, the
#'    minimum, median and maximum proportion of observed events across the
#'    corresponding trials. \cr
#'    \tab \cr
#'    \code{table_comparisons} \tab Identical structure to
#'    \code{table_interventions} but for each observed comparison in the
#'    network. \cr
#'    \tab \cr
#'    \code{table_interventions_missing} \tab The summary results on the
#'    proportion of missing participants for each intervention. See also
#'    'Details'.Identical structure to \code{table_interventions} for the binary
#'    outcome. \cr
#'    \tab \cr
#'    \code{table_comparisons_missing} \tab The summary results on the
#'    proportion of missing participants for each observed comparison in the
#'    network. See also 'Details'. Identical structure to
#'    \code{table_comparisons} for the binary outcome. \cr
#'   }
#'
#' @details The proportion of total missing participants per intervention has
#'   been calculated as the total number of missing participants in an
#'   intervention to the total number of randomised participants in that
#'   intervention. Similarly, for the proportion of total missing participants
#'   per observed comparison,
#'
#'   The median proportion of missing participants in an intervention is the
#'   median of the proportion of missing participants across the trials that
#'   investigated that intervention. Similarly, for the minimum, and the maximum
#'   proportion of missing participants in an intervention. The concept is
#'   similar for the minimum, median, and maximum proportion of missing
#'   participants per observed comparison.
#'
#' @seealso \code{\link[rnmamod]{run_model}},
#'   \href{https://CRAN.R-project.org/package=pcnetmeta}{pcnetmeta} and
#'   \href{https://CRAN.R-project.org/package=gemtc}{gemtc},
#'   \code{\link[rnmamod]{data_preparation}}.
#'
#' @references
#' Lifeng Lin, Jing Zhang, James S. Hodges, Haitao Chu. Performing Arm-Based
#' Network Meta-Analysis in R with the pcnetmeta Package.
#' \emph{J Stat Softw} 2017;\bold{80}(5): 1--25. [\doi{10.18637/jss.v080.i05}]
#'
#' Gert van Valkenhoef and Joel Kuiper. gemtc: Network Meta-Analysis Using
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
#' #           study t1 t2 t3 t4  r1  r2 r3 r4 m1 m2 m3 m4  n1  n2 n3 n4
#' #   Buckley, 2008  1  6 NA NA  67  79 NA NA  2  1 NA NA 110 108 NA NA
#' #    Tyring, 2008  6  7 NA NA  74  12 NA NA  2  0 NA NA 135  42 NA NA
#' # Kragballe, 2009  3  6 NA NA  19 114 NA NA  9  2 NA NA 105 207 NA NA
#' #     Luger, 2008  3  6 NA NA 101 196 NA NA 44  9 NA NA 431 419 NA NA
#' #    Klaber, 1994  2  3 NA NA 175 138 NA NA  2 11 NA NA 234 240 NA NA
#' #   Barrett, 2005  3  4 NA NA  79  79 NA NA 19 18 NA NA 225 236 NA NA
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
    aa <- "The argument 'drug_names' has not been defined."
    bb <- "The intervention ID, as specified in argument 'data' is"
    cc <- "used as intervention names"
    message(cat(paste0("\033[0;", col = 32, "m", aa, " ", bb, " ", cc,
                       "\033[0m", "\n")))
    as.character(1:nt)
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
  network_plot <- nma.networkplot(study,
                                  treatment1,
                                  data = transform,
                                  trtname = drug_names,
                                  multi.show = T, ...)

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
    writexl::write_xlsx(dat$interventions, "table_interventions.xlsx")
    writexl::write_xlsx(dat$comparisons, "table_comparisons.xlsx")
  }

  results <- list(network_plot = network_plot,
                  network_description = knitr::kable(results),
                  table_interventions = dat$interventions,
                  table_comparisons = dat$comparisons)

  return(results)
}

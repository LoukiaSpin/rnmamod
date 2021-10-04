#' Create the network plot
#'
#' @description A function to facilitate drawing the network plot via the \code{\link[pcnetmeta]{nma.networkplot}} function
#'   of the R-package \href{https://CRAN.R-project.org/package=pcnetmeta}{pcnetmeta}. The \code{\link[gemtc]{mtc.data.studyrow}} function
#'   of the R-package \href{https://CRAN.R-project.org/package=gemtc}{gemtc} is also used to convert \code{data} from the required one-trial-per-row format into the one-arm-per-row format.
#'
#' @param data A data-frame of a one-trial-per-row format containing arm-level data of each trial. This format is widely used for BUGS models.
#'   See 'Format' in \code{\link[rnmamod]{run.model}} function for the specification of the columns.
#' @param drug.names A vector of labels with the name of the interventions in the order they appear in the argument \code{data}. If the argument \code{drug.names} is not defined, the order of the interventions
#'   as they appear in \code{data} is used, instead.
#' @param save.xls Logical to indicate whether to export the tabulated results to an Excel 'xlsx' format (via the \code{\link[writexl]{write_xlsx}} function) to the working directory of the user.
#'   The default is \code{FALSE} (do not export to an Excel format).
#' @param ... Additional arguments of the \code{\link[pcnetmeta]{nma.networkplot}} function of the R-package \href{https://CRAN.R-project.org/package=pcnetmeta}{pcnetmeta}.
#'
#' @return A network plot with coloured closed-loops informed by multi-arm trials. Each node refers to the intervention and each link refers to the observed pairwise comparison.
#'   The edges are proportionally to the number of direct treatment comparisons, unless specified otherwise (see \code{\link[pcnetmeta]{nma.networkplot}} function).
#'   the node size is weighted by the total number of direct treatment comparisons of the corresponding treatment, unless specified otherwise (see \code{\link[pcnetmeta]{nma.networkplot}} function).
#'
#'   \code{netplot} also returns five data-frames that describe the network under investigation:
#'   \tabular{ll}{
#'    \code{Network.description} \tab Number of interventions, possible comparisons, direct and indirect comparisons, number of trials in total, number of two-arm and multi-arm trials,
#'    number of randomised participants, proportion of participants completing the trial (completers), and proportion of missing participants. When the outcome is binary, the number of trials with at least one zero event,
#'    and the number of trials with all zero events are also presented. \cr
#'    \tab \cr
#'    \code{Table.interventions} \tab For each interventions, the data-frame presents the number of trials, number of randomised participants, proportion of completers, and thr proportion of missing participants.
#'    When the outcome is binary, the data-frame also presents the proportion of observed events, the minimum, median and maximum proportion of observed events across the corresponding trials.
#'    However, when the outcome is continuous, the data-frame also presents the minimum, median and maximum t-statistic across the corresponding trials. The t-statistic is calculated as the ratio of
#'    the extracted mean outcome, \code{y}, to the standard error, \code{se} (see the \code{data.preparation} function). \cr
#'    \tab \cr
#'    \code{Table.comparisons} \tab The data-frame has the same structure with the \code{Table.interventions}; however, the summary results are illustrated for each observed comparison in the network
#'    When the outcome is continuous, the t-statistic refers to the standardised mean difference defined as the ratio of mean difference to the standard error of mean difference. \cr
#'    \tab \cr
#'    \code{Table.interventions.Missing} \tab The data-frame presents the summary results on the proportion of missing participants for each intervention.
#'    It has the same structure with \code{Table.interventions} for the binary outcome. \cr
#'    \tab \cr
#'    \code{Table.comparisons.Missing} \tab The data-frame presents the summary results on the proportion of missing participants for each observed comparison in the network.
#'    It has the same structure with \code{Table.comparisons} for the binary outcome. \cr
#'   }
#'
#' @seealso \code{\link[rnmamod]{run.model}}, \href{https://CRAN.R-project.org/package=pcnetmeta}{pcnetmeta} and \href{https://CRAN.R-project.org/package=gemtc}{gemtc},
#' \code{\link[rnmamod]{data.preparation}}.
#'
#' @references
#' Lifeng Lin, Jing Zhang, James S. Hodges, Haitao Chu. Performing Arm-Based Network Meta-Analysis in R
#' with the pcnetmeta Package. \emph{J Stat Softw} 2017;\bold{80}(5): 1--25. [\doi{10.18637/jss.v080.i05}]
#'
#' Gert van Valkenhoef and Joel Kuiper. gemtc: Network Meta-Analysis Using Bayesian Methods. \emph{R package version 1.0-1}.
#' 2021. \url{https://CRAN.R-project.org/package=gemtc}
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
#' interv.names <- c("betamethasone dipropionate", "betamethasone valerate", "calcipotriol",
#'                   "calcipotriol plus polytar", "capasal", "two-compound formulation gel",
#'                   "placebo")
#'
#' # Create the network plot
#' netplot(data = nma.bottomley2011, drug.names = interv.names, save.xls = FALSE)
#'
#' @export
netplot <- function(data, drug.names, save.xls, ...){

  save.xls <- if (missing(save.xls)) {
    FALSE
  } else {
    save.xls
  }

  ## Obtain dataset
  r <- data %>% dplyr::select(starts_with("r") | starts_with("y")) # It does not matter whether the outcome is binary or continuous
  n <- data %>% dplyr::select(starts_with("n"))
  t <- data %>% dplyr::select(starts_with("t"))
  nt <- length(table(as.matrix(t)))
  ns <- length(r[, 1])
  na..  <- rep(0, length(r[, 1]))
  for(i in 1:length(r[, 1])){
    na..[i] <- table(!is.na(t[i, ]))["TRUE"]
  }

  drug.names <- if (missing(drug.names)) {
    message(cat(paste0("\033[0;", col = 32, "m", txt = "The argument 'drug.names' has not been defined. The intervention ID, as specified in argument 'data' is used as intervention names", "\033[0m", "\n")))
    as.character(1:nt)
  } else {
    drug.names
  }

  ## Rename columns to agree with gemtc
  names(r) <- paste0("r..", 1:length(r[1, ]), ".")
  names(n) <- paste0("n..", 1:length(n[1, ]), ".")
  names(t) <- paste0("t..", 1:length(t[1, ]), ".")

  ## one row per study arm
  transform <- mtc.data.studyrow(cbind(t, r, n, na..), armVars = c('treatment'= 't', 'response'='r', 'sampleSize'='n'), nArmsVar='na')
  transform$treatment1 <- as.numeric(as.character(transform$treatment))

  ## Obtain network plot
  network.plot <- nma.networkplot(study, treatment1, data = transform, trtname = drug.names, multi.show = T, ...)

  if(dim(data %>% dplyr::select(starts_with("r")))[2] > 0) {
    measure <- "OR"
  } else {
    measure <- "MD"
  }

  dat <- suppressMessages({describe.network(data, drug.names, measure)})

  characteristics <- c("Interventions", "Possible comparisons", "Direct comparisons", "Indirect comparisons",
                       "Trials", "Two-arm trials", "Multi-arm trials", "Randomised participants", "Proportion of completers")
  value <- c(nt,
             dim(combn(nt, 2))[2],
             dat$direct.comp,
             dim(combn(nt, 2))[2] - dat$direct.comp,
             ns,
             dat$two.arm.ns,
             dat$multi.arm.ns,
             dat$total.rand.network,
             dat$prop.obs.network,
             dat$prop.mod.network)

  results <- data.frame(characteristics, value)
  colnames(results) <- c("Characteristic", "Total")

  if (measure == "OR") {
    results[10, ] <- rbind("Proportion of observed events", dat$prop.event.network)
    results[11, ] <- rbind("Trials with at least one zero event", dat$trial.zero.event)
    results[12, ] <- rbind("Trials with all zero events", dat$trial.all.zero.event)
  } else {
    results
  }


  ## Write the tables as .xlsx
  if (save.xls == TRUE) {
    writexl::write_xlsx(dat$Table.interventions, "Table.interventions.xlsx")
    writexl::write_xlsx(dat$Table.comparisons, "Table.comparisons.xlsx")
  }

  results <- list(Network.plot = network.plot,
                  Network.description = knitr::kable(results),
                  Table.interventions = dat$Table.interventions,
                  Table.comparisons   = dat$Table.comparisons)

  if (length(unique(na.omit(unlist(data.preparation(data = data1, measure = "OR")$m)))) > 1) {
    results <- append(results, list(Table.interventions.Missing = dat$Table.interventions.Missing,
                                    Table.comparisons.Missing = dat$Table.comparisons.Missing))
  } else {
    results
  }

  return(results)
}

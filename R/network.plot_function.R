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
#' @param ... Additional arguments of the \code{\link[pcnetmeta]{nma.networkplot}} function of the R-package \href{https://CRAN.R-project.org/package=pcnetmeta}{pcnetmeta}.
#'
#' @return A network plot with coloured closed-loops informed by multi-arm trials. Each node refers to the intervention and each link refers to the observed pairwise comparison.
#'   The edges are proportionally to the number of direct treatment comparisons, unless specified otherwise (see \code{\link[pcnetmeta]{nma.networkplot}} function).
#'   the node size is weighted by the total number of direct treatment comparisons of the corresponding treatment, unless specified otherwise (see \code{\link[pcnetmeta]{nma.networkplot}} function).
#'
#' @seealso \code{\link[rnmamod]{run.model}}, \href{https://CRAN.R-project.org/package=pcnetmeta}{pcnetmeta} and \href{https://CRAN.R-project.org/package=gemtc}{gemtc}.
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
#' netplot(data = nma.bottomley2011, drug.names = interv.names)
#'
#' @export
netplot <- function(data, drug.names, ...){

  options(warn = -1)

  if (missing(drug.names)) {
    message(cat(paste0("\033[0;", col = 32, "m", txt = "The argument 'drug.names' has not been defined. The intervention ID, as specified in argument 'data' is used as intervention names", "\033[0m", "\n")))
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


  ## Rename columns to agree with gemtc
  names(r) <- paste0("r..", 1:length(r[1, ]), ".")
  names(n) <- paste0("n..", 1:length(n[1, ]), ".")
  names(t) <- paste0("t..", 1:length(t[1, ]), ".")


  ## one row per study arm
  transform <- mtc.data.studyrow(cbind(t, r, n, na..), armVars = c('treatment'= 't', 'response'='r', 'sampleSize'='n'), nArmsVar='na')
  transform$treatment1 <- as.numeric(as.character(transform$treatment))


  ## Obtain network plot
  nma.networkplot(study, treatment1, data = transform, trtname = drug.names, multi.show = T, ...)


  if(dim(data %>% dplyr::select(starts_with("r")))[2] > 0) {
    measure <- "OR"
  } else {
    measure <- "MD"
  }

  dat <- describe.network(data, drug.names, measure)

  characteristics <- c("Interventions", "Possible comparisons", "Direct comparisons", "Indirect comparisons",
                       "Trials", "Two-arm trials", "Multi-arm trials", "Randomised participants", "Completers Participants")
  value <- c(nt,
             dim(combn(nt, 2))[2],
             dat$direct.comp,
             dim(combn(nt, 2))[2] - dat$direct.comp,
             ns,
             dat$two.arm.ns,
             dat$multi.arm.ns,
             dat$total.rand.network,
             dat$total.obs.network)

  results <- data.frame(characteristics, value)
  colnames(results) <- c("Characteristic", "Total number")

  if (measure == "OR") {
    results[10, ] <- rbind("Trials with at least one zero event", dat$trial.all.zero.event)
    results[11, ] <- rbind("Trials with all zero events", dat$trial.all.zero.event)
  } else {
    results
  }

  return(list(Network.description = knitr::kable(results),
              Table.interventions = dat$Table.interventions,
              Table.comparisons = dat$Table.comparisons))
}




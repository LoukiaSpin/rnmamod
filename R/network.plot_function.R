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
#' data("nma.baker2009")
#'
#' # Return the first six trials of the dataset
#' head(nma.baker2009)
#' #                 study t1 t2 t3 t4  r1  r2 r3 r4 m1 m2 m3 m4  n1  n2 n3 n4
#' # Llewellyn-Jones, 1996  3  8 NA NA   8   4 NA NA  0  1 NA NA   8   8 NA NA
#' #        Paggiaro, 1998  3  8 NA NA  78  61 NA NA 19 27 NA NA 142 139 NA NA
#' #          Mahler, 1999  6  8 NA NA  98  73 NA NA  9 23 NA NA 135 143 NA NA
#' #        Casaburi, 2000  7  8 NA NA 222 132 NA NA 12 18 NA NA 279 191 NA NA
#' #       van Noord, 2000  6  8 NA NA  29  24 NA NA  7  8 NA NA  47  50 NA NA
#' #         Rennard, 2001  6  8 NA NA  72  65 NA NA 22 29 NA NA 132 135 NA NA
#'
#' # The names of the interventions in the order they appear in the dataset
#' interv.names <- c("budesodine", "budesodine plus formoterol", "fluticasone", "fluticasone plus salmeterol",
#'                   "formoterol", "salmeterol", "tiotropium", "placebo")
#'
#' # Create the network plot
#' netplot(data = nma.baker2009, drug.names = interv.names)
#'
#' @export
netplot <- function(data, drug.names, ...){


  if (missing(drug.names)) {
    message(cat(paste0("\033[0;", col = 32, "m", txt = "The argument 'drug.names' has not been defined. The intervention ID, as specified in argument 'data' is used as intervention names", "\033[0m", "\n")))
  }


  ## Obtain dataset
  r <- data %>% dplyr::select(starts_with("r") | starts_with("y"))
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

}




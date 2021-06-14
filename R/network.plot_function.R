#' Create the network plot
#'
#' @description A function to facilitate drawing the network plot via the \code{\link[pcnetmeta]{nma.networkplot}} function
#'   of the R-package \href{https://CRAN.R-project.org/package=pcnetmeta}{pcnetmeta}. The \code{\link[gemtc]{mtc.data.studyrow}} function
#'   of the R-package \href{https://CRAN.R-project.org/package=gemtc}{gemtc} is also used to convert \code{data} into a one-arm-per-row format.
#'
#' @param data A data-frame of a one-trial-per-row format containing arm-level data of each trial. This format is widely used for BUGS models.
#'   See 'Format' in \code{\link{run.model}} for the specification of the columns.
#' @param drug.names A vector of labels with name of the interventions in the order they appear in the \code{data}.
#' @param ... Additional arguments from the \code{\link[pcnetmeta]{nma.networkplot}} function of the R-package \href{https://CRAN.R-project.org/package=pcnetmeta}{pcnetmeta}.
#'
#' @return A network plot with coloured closed-loops informed by multi-arm trials. Each node refers to the intervention and each link refers to the observed pairwise comparison.
#'
#' @seealso \code{\link[rnmamod]{run.model}}, \href{https://CRAN.R-project.org/package=pcnetmeta}{pcnetmeta} and \href{https://CRAN.R-project.org/package=gemtc}{gemtc}.
#'
#' @references
#' Lifeng Lin, Jing Zhang, James S. Hodges, Haitao Chu. Performing Arm-Based Network Meta-Analysis in R
#' with the pcnetmeta Package. \emph{J Stat Softw} 2017;\bold{80}(5): 1--25. <\doi{10.18637/jss.v080.i05}>
#'
#' Gert van Valkenhoef and Joel Kuiper. gemtc: Network Meta-Analysis Using Bayesian Methods. \emph{R package version 1.0-1}.
#' 2021. \url{https://CRAN.R-project.org/package=gemtc}
#'
#' @author {Loukia M. Spineli}
#'
#' @examples
#' data("cipriani2011")
#'
#' head(cipriani.cont.new)
#'    t1 t2 t3     y1     y2 y3   sd1   sd2 sd3 m1  m2 m3  n1  n2 n3
#' #1  1  2 NA -10.70 -13.31 NA  7.64  7.86  NA  1   6 NA 131 253 NA
#' #2  1  2 NA  -7.19 -12.52 NA 10.57 10.61  NA  3   1 NA 135 137 NA
#' #3  2  7 NA -15.70 -15.65 NA 10.95 10.95  NA  1  10 NA 175 172 NA
#' #4  1  2 NA  -3.35  -8.15 NA 14.03 14.08  NA 10   7 NA 132 130 NA
#' #5  1  2 NA -10.12 -10.80 NA 10.72 10.90  NA  4 138 NA 134 267 NA
#' #6  1 12 NA  21.60  17.50 NA 11.70 12.15  NA 11  10 NA 100  91 NA
#'
#' netplot(data = cipriani.cont.new, drug.names = toupper(letters[1:14]))
#'
#' @export
netplot <- function(data, drug.names, ...){

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
  (transform$treatment1 <- as.numeric(as.character(transform$treatment)))


  ## Obtain network plot
  nma.networkplot(study, treatment1, data = transform, trtname = drug.names, multi.show = T, ...)

}




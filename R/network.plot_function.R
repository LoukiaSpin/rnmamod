#' The network plot
#'
#' @description
#' A function to create the network plot. Each node refers to the intervention and each link refers to the observed pairwise comparison. 
#' Closed loops of interventions informed by multi-arm trials are indicated with different colours. The R-package \code{\link[pcnetmeta]}
#' is used to create the network plot.
#'
#' @param data A data-frame of a one-trial-per-row format containing arm-level data of each trial. This format is widely used for BUGS models.
#' See 'Format' in \code{run.model} for the specification of the columns.
#' @param drug.names A vector of characteristics with name of the interventions as appear in the function \code{run.model}.
#'
#' @return A network plot with coloured closed-loops informed by multi-arm trials.

netplot <- function(data, drug.names, ...){
  
  
  library("pcnetmeta")
  library("gemtc")


  ## Obtain dataset
  r <- data %>% dplyr::select(starts_with("r") | starts_with("y"))
  n <- data %>% dplyr::select(starts_with("n") | starts_with("c"))  
  t <- data %>% dplyr::select(starts_with("t"))
  nt <- length(table(as.matrix(t)))
  ns <- length(r[, 1])
  na..  <- rep(0, length(r[, 1]))
  for(i in 1:length(r[, 1])){
    na..[i] <- table(!is.na(t[i, ]))["TRUE"]
  }


  ## Rename columns to agree with gemtc
  names(r) <- paste0("r..",1:length(r[1, ]),".")
  names(n) <- paste0("n..",1:length(n[1, ]),".")
  names(t) <- paste0("t..",1:length(t[1, ]),".")


  ## one row per study arm
  transform <- mtc.data.studyrow(cbind(t, r, n, na..), armVars = c('treatment'= 't', 'response'='r', 'sampleSize'='n'), nArmsVar='na')
  (transform$treatment1 <- as.numeric(as.character(transform$treatment)))


  ## Obtain network plot
  nma.networkplot(study, treatment1, data = transform, trtname = drug.names, alphabetic = F, title = "", multi.show = T, ...)
  
}




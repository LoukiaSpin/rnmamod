#' The network plot
#'
#' @description
#' A function to create the network plot and illustrate the risk of bias due to missing outcome data in each node and link using different colours.
#' The node refers to the intervention and the link refers to the observed pairwise comparison. The risk of bias due to missing outcome data
#' can be low (the proportion of missing outcome data is up to 5\%), moderate (the proportion of missing outcome data is more than 5\% and up to 20\%),
#' or high (the proportion of missing outcome data is more than 20\%). Green, orange and red represent low, moderate, and high risk of bias.
#' For each node, the risk of bias is determined by the total proportion of missing outcome data in the corresponding intervention, namely, the ratio of
#' the sum of missing outcome data to the sum of randomised in that intervention.
#' For each link, the risk of bias is determined by the total proportion of missing outcome data across the corresponding trials, namely, the ratio of
#' the sum of missing outcome data to the sum of randomised in these trials. The function uses the \code{net.plot} function from the BUGSnet package to
#' colour the nodes and links according to the risk of bias due to missing outcome data. The user can specify to plot the network plot without illustrating
#' the risk of bias due to missing outcome data. In this case, the function uses the \code{nma.networkplot} function from the pcnetmeta package to create
#' the network plot. The fucntion allows the user to incorporate the arguments of \code{net.plot} (when \code{show.bias = TRUE}) and \code{nma.networkplot}
#' (when \code{show.bias = FALSE}).
#'
#' @param data A data-frame of a one-trial-per-row format containing arm-level data of each trial. This format is widely used for BUGS models.
#' See 'Format' in \code{nma.continuous.full.model} for the specification of the columns.
#' @param drug.names A vector of characteristics with name of the interventions as appear in the function \code{run.model}.
#' @param show.bias Indicates whether to show the risk of bias in the nodes and links.
#'
#' @return A network plot with coloured nodes and links to indicate the risk of bias due to missing outcome data.
#'
#' \dontshow{load("./data/NMA Dataset Continuous.RData")}
#' @examples
#' ### Show the data (one-trial-per-row format)
#' (data <- as.data.frame(one.stage.dataset.NMA[[3]]))
#' drug.names <- sapply(1:8, function(x) letters[x])
#'
#' netplot(data = res, drug.names = drug.names)
#'
#' @export
netplot <- function(data, drug.names, show.bias, ...){


  ## Obtain dataset
  m <- data %>% dplyr::select(starts_with("m"))
  n <- data %>% dplyr::select(starts_with("n"))  # Number randomised
  t <- data %>% dplyr::select(starts_with("t"))
  nt <- length(table(as.matrix(t)))
  ns <- length(m[, 1])
  na..  <- rep(0, length(m[, 1]))
  for(i in 1:length(m[, 1])){
    na..[i] <- table(!is.na(t[i, ]))["TRUE"]
  }


  ## Rename columns to agree with gemtc
  names(m) <- paste0("m..",1:length(m[1, ]),".")
  names(n) <- paste0("n..",1:length(n[1, ]),".")
  names(t) <- paste0("t..",1:length(t[1, ]),".")


  ## one row per study arm
  transform0 <- mtc.data.studyrow(cbind(t, m, n, na..), armVars = c('treatment'= 't', 'response'='m', 'sampleSize'='n'), nArmsVar='na')
  (transform0$treatment1 <- as.numeric(as.character(transform0$treatment)))
  (transform0$treatment <- as.numeric(as.character(transform0$treatment)))
  #for(i in sort(unique(transform0$treatment))) {
  #  transform0[transform0$treatment == i, 2] <- drug.names[i]
  #}
  oldvals <- sort(unique(transform0$treatment))
  for(i in 1:length(oldvals)) {
    transform[transform0$treatment == oldvals[i], 2] <- drug.names[i]
  }


  ## Prepare data to use BUGSnet
  transform <- data.prep(arm.data = transform0, varname.t = "treatment", varname.s = "study")


  ## Characteristics of the networks (except for missing outcome data)
  network.char <- net.tab(data = transform, outcome = "response", N = "sampleSize", type.outcome = "binomial", time = NULL)


  ## Turn arm-level to contrast-level dataset
  (pairwise <- pairwise(as.list(t), event = as.list(m), n = as.list(n), data = cbind(t, m, n), studlab = 1:ns)[, c(3:6, 8, 7, 9)])
  colnames(pairwise) <- c("study", "t1", "t2", "m1", "m2", "n1", "n2")


  ## Calculate summary of %MOD in each intervention
  risk.drug <- round(aggregate(unlist(m), by = list(unlist(t)), sum)[, 2]/aggregate(unlist(n), by = list(unlist(t)), sum)[, 2], 2)


  ## Calculate %total MOD per observed comparison
  trial.mod <- apply(pairwise[, 4:5], 1, sum)
  trial.size <- apply(pairwise[, 6:7], 1, sum)
  comp <- paste0(pairwise[, 3], "vs", pairwise[, 2])
  risk.comp <- round(aggregate(trial.mod, by = list(comp), sum)[, 2]/aggregate(trial.size, by = list(comp), sum)[, 2], 2)


  ## Obtain the network plot
  if(show.bias == T) {

    net.plot(transform, node.lab.cex = 1.5, node.scale = 1, edge.scale = 1, label.offset1 = 1.5,
             node.colour = ifelse(risk.drug <= 0.05, "green4", ifelse(risk.drug > 0.20, "brown1", "orange")),
             edge.colour = ifelse(risk.comp <= 0.05, "green4", ifelse(risk.comp > 0.20, "brown1", "orange")), ...)

    #  legend("bottom", legend = c("low", "moderate", "high"), horiz = T, seg.len = 0.2, x.intersp = 0.1,
    #         inset = c(0, -0.1), text.width = c(0, 0.08, 0.05),
    #         bty = "n", xpd = TRUE, cex = 1.2, lty = c(1, 1, 1), lwd = c(2, 2, 2), col = c("green4", "orange", "brown1"))

  } else {

    nma.networkplot(study, treatment1, data = transform0, trtname = drug.names, alphabetic = F, title = "", text.cex = 1.7, multi.show = T, ...)

  }

  return(network.char)

}




#' A function to detect the frail comparisons in multi-arm trials
#'
#' @description \code{improved.UME} detects the frail comparisons in multi-arm trials, that is, comparisons between non-baseline interventions not investigated in any two-arm trial of the network (Spineli, 2021).
#'   The 'original' model of Dias et al. (2013) omits the frail comparisons from the estimation process of the unrelated mean effects model.
#'   Consequently, their posterior distribution coincides with the prior distribution yielding implausible posterior standard deviations.
#'
#' @param t A data-frame of the one-trial-per-row format containing the intervention identifier in each arm of every trial (see 'Details' below, and 'Arguments' in \code{run.model}).
#' @param m A data-frame of the one-trial-per-row format containing the number of missing participant outcome data (MOD) in each arm of every trial (see 'Details' below, and 'Arguments' in \code{run.model}).
#' @param N A data-frame of the one-trial-per-row format containing the number of participants randomised on the assigned intervention in each arm of every trial (see 'Details' below, and 'Arguments' in \code{run.model}).
#' @param ns A scale parameter on the number trials.
#' @param na A vector of length equal to \code{ns} with the number of arms in each trial.
#'
#' @return The output of \code{improved.UME} is a list of elements that are used by the \code{run.UME}:
#' \tabular{ll}{
#'  \code{nbase.multi} \tab A scalar parameter on the number of frail comparisons.\cr
#'  \tab \cr
#'  \code{t1.bn} \tab A vector with numeric values on the first arm of each frail comparison.\cr
#'  \tab \cr
#'  \code{t2.bn} \tab A vector with numeric values on the second arm of each frail comparison.\cr
#'  \tab \cr
#'  \code{base} \tab A vector with numeric values on the baseline intervention of the multi-arm trials that contain the frail comparisons.\cr
#'  \tab \cr
#'  \code{obs.comp} \tab A data-frame that indicates how many two-arms and multi-arm trials have included each pairwise comparison observed in the network.\cr
#' }
#'
#' @details \code{improved.UME} is integrated in the \code{run.UME} function and calls the output of the \code{data.preparation} function after sorting the rows so that multi-arm trials appear at the bottom of the dataset.
#'   When there are no multi-arm trials or no frail comparisons in the network, \code{improved.UME} returns only the element \code{obs.comp} (see, 'Value').
#'
#' @author {Loukia M. Spineli}
#'
#' @seealso \code{\link{run.UME}}, \code{\link{data.preparation}}, \code{\link{run.model}}
#'
#' @references
#' Spineli LM. A novel framework to evaluate the consistency assumption globally in a network of interventions. \emph{submitted} 2021.
#'
#' Dias S, Welton NJ, Sutton AJ, Caldwell DM, Lu G, Ades AE. Evidence synthesis for decision making 4: inconsistency in networks of evidence based on randomized controlled trials. \emph{Med Decis Making} 2013;\bold{33}(5):641--56. [\doi{10.1177/0272989X12455847}]
#'
#' @export
improved.UME <- function(t, m, N, ns, na){


  ## Turn into contrast-level data: one row per possible comparison in each trial ('netmeta')
  wide.format <- pairwise(as.list(t), event = as.list(m), n = as.list(N), data = cbind(t, m, N), studlab = 1:ns)[, c(3:6, 8, 7, 9)]
  colnames(wide.format) <- c("study", "t1", "t2", "m1", "m2", "n1", "n2")


  ## Create a vector with the pairwise comparisons of each row
  wide.format$comp <- paste0(wide.format$t1, "vs", wide.format$t2)


  ## Repeat 'na' (number of arms in each trial) as many times as the number of possible comparisons in each trial
  arms0 <- list()
  for(i in 1:ns) {
    arms0[[i]] <- rep(na[i], dim(combn(na[i], 2))[2] )
  }


  ## Indicate whether a trial is two-arm or multi-arm
  wide.format$arms <- unlist(arms0)
  wide.format[wide.format$arms > 2, 9] <- "multi-arm"
  wide.format[wide.format$arms < 3, 9] <- "two-arm"


  ## The frequency of each observed comparisons in two-arm and multi-arm trials
  tab.comp.arms0 <- xtabs(~ comp + arms, data = wide.format)


  ## Turn 'tab.comp.arms0' into a data-frame
  #if(dim(tab.comp.arms0)[2] == 1 || (length(unique(ifelse(as.matrix(tab.comp.arms0)[, 2] == 0 & as.matrix(tab.comp.arms0)[, 1] != 0, T, F))) == 1)) {
  tab.comp.arms <- if (unique(wide.format$arms == "multi-arm")[1] == T) {
    data.frame(names(tab.comp.arms0[, 1]), tab.comp.arms0[, 1], rep(0, dim(tab.comp.arms0)[1]))
  } else if (unique(wide.format$arms == "two-arm")[2] == T) {
    data.frame(names(tab.comp.arms0[, 1]), rep(0, dim(tab.comp.arms0)[1]), tab.comp.arms0[, 1])
  #} else if(dim(tab.comp.arms0)[2] > 1 & length(unique(ifelse(as.matrix(tab.comp.arms0)[, 2] == 0 & as.matrix(tab.comp.arms0)[, 1] != 0, T, F))) == 2) {
  } else {
    data.frame(names(tab.comp.arms0[, 1]), tab.comp.arms0[, 1], tab.comp.arms0[, 2])
  }
  colnames(tab.comp.arms) <- c("comp", "multi", "two")
  rownames(tab.comp.arms) <- NULL


  ## Keep only those comparisons not studied in two-arm trials
  if (unique((wide.format$arms == "multi-arm")) == TRUE || unique((wide.format$arms == "two-arm")) != TRUE) {
    tab.comp.arms$select <- ifelse(tab.comp.arms$two == 0 & tab.comp.arms$multi != 0, T, F)
    subs <- subset(tab.comp.arms, select == T, select = comp)

    ## Match the selected comparisons with the study id
    pairwise.n0 <- list()
    for(i in 1:length(subs[, 1])) {
      pairwise.n0[[i]] <- wide.format[which(wide.format$comp == unlist(subs)[i]), 1:3]
    }
    pairwise.n1 <- do.call(rbind, pairwise.n0)


    ## When more studies correspond to a comparison, remove the duplicated rows
    #pairwise.n <- pairwise.n1[!duplicated(pairwise.n1[, 2:3]), ]


    ## Sort by the study id in increasing order
    #final0 <- pairwise.n[order(pairwise.n$study), ]
    final0 <- pairwise.n1[order(pairwise.n1$study), ]


    ## Find the unique comparisons with the baseline intervention
    indic0 <- list()
    for(i in 1:ns) {
      indic0[[i]] <- combn(t(na.omit(t(t[i, ]))), 2)[, 1:(na[i] - 1)]
    }
    (indic <- unique(t(do.call(cbind, indic0))))
    t1.indic <- indic[, 1]    # Baseline interventions at the corresponding comparisons
    t2.indic <- indic[, 2]    # Non-baseline interventions at the corresponding comparisons
  #}


  ## Finally, reduce to comparisons between non-baseline interventions
  #if (dim(tab.comp.arms0)[2] == 1 || length(unique(ifelse(as.matrix(tab.comp.arms0)[, 2] == 0 & as.matrix(tab.comp.arms0)[, 1] != 0, T, F))) == 1) {

  #  final <- NA

  #} else if(dim(tab.comp.arms0)[2] > 1 & length(unique(ifelse(as.matrix(tab.comp.arms0)[, 2] == 0 & as.matrix(tab.comp.arms0)[, 1] != 0, T, F))) == 2) {

    pre.final <- final0[!is.element(paste0(final0[, 2], "vs", final0[, 3]), paste0(t1.indic, "vs", t2.indic)), ]

    final <- pre.final[!duplicated(pre.final[, 2:3]), ]

    base <- rep(NA, length(final[, 1]))   # Baseline interventions in the selected trials in 'final'

    for (i in 1:length(final[, 1])) {
       final$base[i] <- unique(t[final$study[i], 1])
    }
  }



  ## Add also the baseline treatment for each selected trial
  #if (dim(tab.comp.arms0)[2] == 1 || length(unique(ifelse(as.matrix(tab.comp.arms0)[, 2] == 0 & as.matrix(tab.comp.arms0)[, 1] != 0, T, F))) == 1) {
  #
  #  base <- nbase.multi <- NA

  #} else if(dim(tab.comp.arms0)[2] > 1 & length(unique(ifelse(as.matrix(tab.comp.arms0)[, 2] == 0 & as.matrix(tab.comp.arms0)[, 1] != 0, T, F))) == 2) {

  #  base <- rep(NA, length(final[, 1]))   # Baseline interventions in the selected trials in 'final'

  #  for(i in 1:length(final[, 1])){
   #   final$base[i] <- unique(t[final$study[i], 1])
   # }

    #nbase.multi <- dim(final[!duplicated(final[, 2:4]), 2:4])[1] # *Unique* non-baseline interventions in the selected trials in 'final'
  #}


  #if (dim(tab.comp.arms0)[2] == 1 || length(unique(ifelse(as.matrix(tab.comp.arms0)[, 2] == 0 & as.matrix(tab.comp.arms0)[, 1] != 0, T, F))) == 1) {

  #  return(list(obs.comp = tab.comp.arms))

  #} else if(dim(tab.comp.arms0)[2] > 1 & length(unique(ifelse(as.matrix(tab.comp.arms0)[, 2] == 0 & as.matrix(tab.comp.arms0)[, 1] != 0, T, F))) == 2) {

    #return(list(nbase.multi = nbase.multi, t1.bn = final$t1, t2.bn = final$t2, base = final$base, obs.comp = tab.comp.arms))

  #}
  if (unique(wide.format$arms == "multi-arm")[1] == T || unique(wide.format$arms == "two-arm")[2] == T) {
    return(list(obs.comp = tab.comp.arms))
  } else {
    return(list(nbase.multi = length(final[, 1]), t1.bn = final$t1, t2.bn = final$t2, ref.base = min(final$base), obs.comp = tab.comp.arms))
  }

}


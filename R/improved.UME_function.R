#' Detect the frail comparisons in multi-arm trials
#'
#' @description Detects the frail comparisons in multi-arm trials, that is,
#'   comparisons between non-baseline interventions not investigated in any
#'   two-arm trial in the network (Spineli, 2021). The 'original' model of
#'   Dias et al. (2013) omits the frail comparisons from the estimation process
#'   of the unrelated mean effects model. Consequently, their posterior
#'   distribution coincides with the prior distribution yielding implausible
#'   posterior standard deviations.
#'
#' @param t A data-frame of the one-trial-per-row format containing the
#'   intervention identifier in each arm of every trial (see 'Details' below,
#'   and 'Format' in \code{\link{run_model}}).
#' @param N A data-frame of the one-trial-per-row format containing the number
#'   of participants randomised to the assigned intervention in each arm of
#'   every trial (see 'Details' below, and 'Format' in \code{\link{run_model}}).
#' @param ns A scale parameter on the number trials.
#' @param na A vector of length equal to \code{ns} with the number of arms in
#'   each trial.
#'
#' @return The output of \code{improved_ume} is a list of elements that are
#'   inherited by \code{\link{run_ume}}:
#'   \item{nbase_multi}{A scalar parameter on the number of frail comparisons.}
#'   \item{t1_bn}{A vector with numeric values referring to the first arm of
#'   each frail comparison.}
#'   \item{t2_bn}{A vector with numeric values referring to the second arm of
#'   each frail comparison.}
#'   \item{ref_base}{A scalar referring to the reference intervention for the
#'   subnetwork of interventions in frail comparisons.}
#'   \item{base}{A vector with numeric values referring to the baseline
#'   intervention of the multi-arm trials that contain the frail comparisons.}
#'   \item{obs_comp}{A data-frame that indicates how many two-arm and
#'   multi-arm trials have included each pairwise comparison observed in the
#'   network.}
#'
#' @details \code{improved_ume} is integrated in \code{\link{run_ume}} and
#'   calls the output of \code{\link{data_preparation}} after sorting the
#'   rows so that multi-arm trials appear at the bottom of the dataset.
#'   When there are no multi-arm trials or no frail comparisons in the network,
#'   \code{improved_ume} returns only the element \code{obs_comp}
#'   (see, 'Value').
#'
#' @author {Loukia M. Spineli}
#'
#' @seealso \code{\link{data_preparation}}, \code{\link{run_model}},
#'   \code{\link{run_ume}}
#'
#' @references
#' Dias S, Welton NJ, Sutton AJ, Caldwell DM, Lu G, Ades AE. Evidence synthesis
#' for decision making 4: inconsistency in networks of evidence based on
#' randomized controlled trials.
#' \emph{Med Decis Making} 2013;\bold{33}(5):641--56.
#' doi: 10.1177/0272989X12455847
#'
#' Spineli LM. A revised framework to evaluate the consistency assumption
#' globally in a network of interventions. \emph{Med Decis Making} 2021.
#' doi: 10.1177/0272989X211068005
#'
#' @export
improved_ume <- function(t, N, ns, na) {

  # Function to turn wide- to long-format for an element
  log_format <- function (input) {
    if (length(input[1, ]) > 2) {
      long_form0 <- apply(input, 1, function(x) {combn(na.omit(x), 2)})
      long_form <- t(do.call(cbind, long_form0))
    } else {
      long_form <- input
    }

    return(long_form)
  }

  # Turn into contrast-level data
  poss_comp <- lapply(na, function(x) {combn(x, 2)})
  len_poss_comp <- unlist(lapply(poss_comp, function(x) {dim(x)[2]}))
  study <- rep(1:ns, len_poss_comp)
  t_long_form <- log_format(t)
  n_long_form <- log_format(N)
  wide_format0 <- data.frame(study, t_long_form, n_long_form, n_long_form)
  colnames(wide_format0) <- c("study", "t1", "t2", "n1", "n2", "n1", "n2")
  #wide_format0 <- pairwise(as.list(t),
  #                         event = as.list(N),
  #                         n = as.list(N),
  #                         data = cbind(t, N, N),
  #                         studlab = 1:ns)[, c(3:6, 8, 7, 9)]
  #colnames(wide_format0) <- c("study", "t1", "t2", "n1", "n2", "n1", "n2")

  # Ensure that t1 < t2 and correspondingly for the other elements
  treat <- treat0 <- wide_format0[, 2:3]
  resp <- resp0 <- wide_format0[, 4:5]
  rand <- rand0 <- wide_format0[, 6:7]
  treat_list <- resp_list <- rand_list <- list()
  for (i in seq_len(length(wide_format0[, 1]))) {
    treat_list[[i]] <- treat0[i, ]
    resp_list[[i]] <- resp0[i, ]
    rand_list[[i]] <- rand0[i, ]
  }
  for (i in seq_len(length(wide_format0[, 1]))) {
    treat[i, ] <- unlist(treat_list[[i]])[order(unlist(treat_list[[i]]),
                                                na.last = TRUE)]
    resp[i, ] <- unlist(resp_list[[i]])[order(unlist(treat_list[[i]]),
                                               na.last = TRUE)]
    rand[i, ] <- unlist(rand_list[[i]])[order(unlist(treat_list[[i]]),
                                              na.last = TRUE)]
  }

  wide_format <- data.frame(study = wide_format0$study, treat, resp, rand)

  # Create a vector with the pairwise comparisons of each row
  wide_format$comp <- paste0(wide_format$t1, "vs", wide_format$t2)

  # Repeat 'na' (number of arms in each trial) as many times as the number of
  # possible comparisons in each trial
  arms0 <- list()
  for (i in 1:ns) {
    arms0[[i]] <- rep(na[i], dim(combn(na[i], 2))[2])
  }

  # Indicate whether a trial is two-arm or multi-arm
  wide_format$arms <- unlist(arms0)
  wide_format[wide_format$arms > 2, 9] <- "multi-arm"
  wide_format[wide_format$arms < 3, 9] <- "two-arm"

  # The frequency of each observed comparisons in two-arm and multi-arm trials
  tab_comp_arms0 <- xtabs(~ comp + arms, data = wide_format)

  ## Turn 'tab_comp_arms0' into a data-frame
  tab_comp_arms <- if (dim(unique(
    as.data.frame(tab_comp_arms0)["arms"]))[1] == 1 &
    levels(unlist(as.data.frame(tab_comp_arms0)["arms"]))[1] == "multi-arm") {

    data.frame(names(tab_comp_arms0[, 1]),
               tab_comp_arms0[, 1],
               rep(0, dim(tab_comp_arms0)[1]))
  } else if (dim(unique(
    as.data.frame(tab_comp_arms0)["arms"]))[1] == 1 &
    levels(unlist(as.data.frame(tab_comp_arms0)["arms"]))[1] == "two-arm") {
    data.frame(names(tab_comp_arms0[, 1]),
               rep(0, dim(tab_comp_arms0)[1]),
               tab_comp_arms0[, 1])
  } else if (dim(unique(as.data.frame(tab_comp_arms0)["arms"]))[1] == 2) {
    data.frame(names(tab_comp_arms0[, 1]),
               tab_comp_arms0[, 1],
               tab_comp_arms0[, 2])
  }
  colnames(tab_comp_arms) <- c("comp", "multi", "two")
  rownames(tab_comp_arms) <- NULL

  # Keep only those comparisons not studied in two-arm trials
  if (dim(unique(
    as.data.frame(tab_comp_arms0)["arms"]))[1] == 1 &
    levels(unlist(as.data.frame(tab_comp_arms0)["arms"]))[1] == "multi-arm" ||
    dim(unique(as.data.frame(tab_comp_arms0)["arms"]))[1] == 2) {
    tab_comp_arms$select <- ifelse(tab_comp_arms$two == 0 &
                                     tab_comp_arms$multi != 0, TRUE, FALSE)
    subs <- subset(tab_comp_arms, select == TRUE, select = comp)

    # Match the selected comparisons with the study id
    pairwise_n0 <- list()
    for (i in seq_len(length(subs[, 1]))) {
      pairwise_n0[[i]] <- wide_format[which(
        wide_format$comp == unlist(subs)[i]), 1:3]
    }
    pairwise_n1 <- do.call(rbind, pairwise_n0)
    final0 <- pairwise_n1[order(pairwise_n1$study), ]

    # Find the unique comparisons with the baseline intervention
    indic0 <- list()
    for (i in 1:ns) {
      indic0[[i]] <- combn(t(na.omit(t(t[i, ]))), 2)[, 1:(na[i] - 1)]
    }
    indic <- unique(t(do.call(cbind, indic0)))
    # Baseline interventions at the corresponding comparisons
    t1_indic <- indic[, 1]
    # Non-baseline interventions at the corresponding comparisons
    t2_indic <- indic[, 2]
    pre_final <- final0[!is.element(paste0(final0[, 2], "vs", final0[, 3]),
                                    paste0(t1_indic, "vs", t2_indic)), ]
    final <- pre_final[!duplicated(pre_final[, 2:3]), ]
    # Baseline interventions in the selected trials in 'final'
    base <- rep(NA, length(final[, 1]))

    for (i in seq_len(length(final[, 1]))) {
       base[i] <- unique(t[final$study[i], 1])
    }
  }

  if (dim(unique(
    as.data.frame(tab_comp_arms0)["arms"]))[1] == 1 &
    levels(unlist(as.data.frame(tab_comp_arms0)["arms"]))[1] == "two-arm") {
    return(list(obs_comp = tab_comp_arms))
  } else {
    return(list(nbase_multi = length(final[, 1]),
                t1_bn = final$t1,
                t2_bn = final$t2,
                ref_base = min(base),
                base = base,
                obs_comp = tab_comp_arms))
  }
}

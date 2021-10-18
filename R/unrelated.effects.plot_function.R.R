#' End-user-ready results for unrelated trials effects model
#'
#' @description Performs the unrelated trial effects model and create a panel of
#'   interval plots on the results of each trial and corresponding pairwise
#'   comparison.
#'
#' @param data A data-frame of a one-trial-per-row format containing arm-level
#'   data of each trial. This format is widely used for BUGS models. See
#'   'Format' in \code{\link[rnmamod]{run_model}} function for the specification
#'   of the columns.
#' @param measure Character string indicating the effect measure with values
#'   \code{"OR"}, \code{"MD"}, \code{"SMD"}, or \code{"ROM"} for the odds ratio,
#'   mean difference, standardised mean difference and ratio of means,
#'   respectively.
#' @param char A data-frame of three columns and number of rows equal to the
#'   number of trials in \code{data}. Each column refers to a
#'   trial-characteristic with nominal elements.
#' @param drug_names A vector of labels with the name of the interventions in
#'   the order they appear in the argument \code{data}. If \code{drug_names} is
#'   not defined, the order of the interventions as they appear in \code{data}
#'   is used, instead.
#' @param trial_names A vector of labels with the name of the trials in the
#'   order they appear in the argument \code{data}. If the argument
#'   \code{drug_names} is not defined, the order of the interventions as they
#'   appear in \code{data} is used, instead.
#' @param mean_value A numeric value for the mean of the normal distribution of
#'   the informative missingness parameter. The same value is considered for all
#'   trial-arms of the dataset. The default argument is 0 and corresponds to the
#'   missing-at-random assumption. For the informative missingness odds ratio
#'   and the informative missingness ratio of means, the mean value is defined
#'   in the logarithmic scale.
#' @param var_value A positive non-zero number for the variance of the normal
#'   distribution of the informative missingness parameter. When the
#'   \code{measure} is \code{"OR"}, \code{"MD"}, or \code{"SMD"} the default
#'   argument is 1; When the \code{measure} is \code{"ROM"} the default argument
#'   is 0.04.
#' @param rho A numeric value in the interval [-1, 1] that indicates the
#'   correlation coefficient between two missingness parameters in a trial.
#'   The same value is considered across all trials of the dataset. The default
#'   argument is 0 and corresponds to uncorrelated missingness parameters.
#' @param save_xls Logical to indicate whether to export the tabulated results
#'   to an Excel 'xlsx' format (via the \code{\link[writexl]{write_xlsx}}
#'   function) to the working directory of the user. The default is \code{FALSE}
#'   (do not export to an Excel format).
#'
#' @return A panel of interval plots for each observed comparison in the
#'   network, when there are up to 21 trials in the \code{data}. Otherwise, an
#'   Excel file with the arm-level data of each trial, the corresponding effect
#'   measure and 95\% confidence interval (lower, upper), the intervention
#'   compared, and the three characteristics (as defined in \code{char}). For
#'   larger datasets, the plot becomes cluttered and it is difficult to identify
#'   the trial names. Hence, exporting the results in an Excel file is an viable
#'   alternative.
#'
#' @details The unrelated trial effects model may be an alternative to network
#'   meta-analysis, when the latter is not deemed appropriate (e.g.,
#'   considerable statistical heterogeneity, or substantial intransitivity). In
#'   the presence of missing participant outcome data (MOD), the effect size and
#'   standard error are adjusted by applying the pattern-mixture model with
#'   Taylor series in trial-arms with MOD (White et al., 2008;
#'   Mavridis et al., 2015). The \code{unrelated_effects_plot} function calls
#'   the \code{taylor_imor} and \code{taylor_continuous} functions (for binary
#'   and continuous outcome, respectively) to employ pattern-mixture model with
#'   Taylor series. The \code{unrelated_effects_plot} function considers the
#'   informative missingness odds ratio in the logarithmic scale for binary
#'   outcome data (White et al., 2008), the informative missingness difference
#'   of means when \code{measure} is \code{"MD"} or \code{"SMD"}, and the
#'   informative missingness ratio of means in the logarithmic scale when
#'   \code{"ROM"} is the effect measure (Mavridis et al., 2015).
#'
#'   The number of interval plots equals the number of observed comparisons in
#'   the network. In each interval plot, the y-axis refers to all trials of the
#'   network and x-axis refers to the selected effect measure. The odds ratio
#'   and ratio of means are calculated in the logarithmic scale but they are
#'   reported in their original scale.
#'
#'   \code{unrelated_effects_plot} depicts all three characteristics for each
#'   trial using different colours, line-types and point-shapes for the
#'   corresponding interval and point estimate. Ideally, each characteristic
#'   should have no more than three categories; otherwise, the plot becomes
#'   cluttered. For now, the \code{unrelated_effects_plot} function uses the
#'   default color palette, line-types and point-shapes.
#'
#' @seealso \code{\link[rnmamod]{run_model}},
#'   \code{\link[rnmamod]{taylor_continuous}},
#'   \code{\link[rnmamod]{taylor_imor}}
#'
#' @references
#' Mavridis D, White IR, Higgins JP, Cipriani A, Salanti G. Allowing for
#' uncertainty due to missing continuous outcome data in pairwise and network
#' meta-analysis. \emph{Stat Med} 2015;\bold{34}(5):721--41.
#' [\doi{10.1002/sim.6365}]
#'
#' White IR, Higgins JP, Wood AM. Allowing for uncertainty due to missing data
#' in meta-analysis--part 1: two-stage methods.
#' \emph{Stat Med} 2008;\bold{27}(5):711--27. [\doi{10.1002/sim.3008}]
#'
#' @author {Loukia M. Spineli}
#'
#' @export
unrelated_effects_plot <- function(data,
                                   measure,
                                   char,
                                   drug_names,
                                   trial_names,
                                   mean_value,
                                   var_value,
                                   rho,
                                   save_xls) {

  item <- data_preparation(data, measure)
  na <- as.vector(do.call(
    rbind, lapply(1:item$ns, function(i) dim(combn(item$na[i], 2))[2])))

  # Default arguments
  char <- if (missing(char)) {
    stop("The argument 'char' needs to be defined", call. = F)
  } else {
    char
  }
  drug_names <- if (missing(drug_names)) {
    message(cat(paste0("\033[0;",
                       col = 32,
                       "m",
                       txt = "The argument 'drug_names' has not been defined.
                       The intervention ID, as specified in 'data' is used as
                       intervention names", "\033[0m", "\n")))
    as.character(1:item$nt)
  } else {
    drug_names
  }
  trial_names <- if (missing(trial_names)) {
    message(cat(paste0("\033[0;",
                       col = 32,
                       "m",
                       txt = "The argument 'trial_names' has not been defined.
                       The trial ID, as specified in the argument 'data' is used
                       as trial names", "\033[0m", "\n")))
    as.character(1:item$ns)
  } else {
    trial_names
  }
  mean_value <- ifelse(missing(mean_value), 0, mean_value)
  var_value <- ifelse(missing(var_value) &
                        (is.element(measure, c("OR", "MD", "SMD"))), 1,
                      ifelse(missing(var_value) & measure == "ROM", 0.2^2,
                             var_value))
  rho <- ifelse(missing(rho), 0, rho)
  save_xls <- if (missing(save_xls)) {
    FALSE
  } else {
    save_xls
  }

  # Turn into contrast-level data ('netmeta')
  if (is.element(measure, c("MD", "SMD", "ROM"))) {
    pairwise_observed0 <-
      pairwise(as.list(item$t),
               mean = as.list(item$y0),
               sd = as.list(item$sd0),
               n = as.list(item$N),
               data = cbind(item$t, item$y0, item$sd0, item$N),
               studlab = 1:item$ns)[, c(3:5, 7, 10, 8, 11, 6, 9)]
    colnames(pairwise_observed0) <- c("study",
                                      "arm1",
                                      "arm2",
                                      "y1",
                                      "y2",
                                      "sd1",
                                      "sd2",
                                      "n1",
                                      "n2")

    pairwise_mod0 <- pairwise(as.list(item$t),
                              mean = as.list(item$y0),
                              sd = as.list(item$sd0),
                              n = as.list(item$m),
                              data = cbind(item$t, item$y0, item$sd0, item$m),
                              studlab = 1:item$ns)[, c(6, 9)]
    colnames(pairwise_mod0) <- c("m1", "m2")

    # Ensure that t1 < t2 and correspondingly for the other elements
    treat <- treat0 <- pairwise_observed0[, 2:3]
    y_mean <- y_mean0 <- pairwise_observed0[, 4:5]
    sd_mean <- sd_mean0 <- pairwise_observed0[, 6:7]
    miss <- miss0 <- pairwise_mod0[, 1:2]
    rand <- rand0 <- pairwise_observed0[, 8:9]
    for (i in seq_len(length(pairwise_observed0[, 1]))) {
      treat[i, ] <- treat0[i, order(treat0[i, ], na.last = T)]
      y_mean[i, ] <- y_mean0[i, order(treat0[i, ], na.last = T)]
      sd_mean[i, ] <- sd_mean0[i, order(treat0[i, ], na.last = T)]
      miss[i, ] <- miss0[i, order(treat0[i, ], na.last = T)]
      rand[i, ] <- rand0[i, order(treat0[i, ], na.last = T)]
    }

    pairwise_observed <- data.frame(study = pairwise_observed0$study,
                                    treat,
                                    y_mean,
                                    sd_mean,
                                    rand)
    pairwise_mod <- miss
    pairwise_data <- data.frame(pairwise_observed[, c(1, 4:7)],
                                pairwise_mod,
                                pairwise_observed[, c(8:9, 2:3)])
  } else {
    pairwise_observed0 <-
      pairwise(as.list(item$t),
               event = as.list(item$r),
               n = as.list(item$N),
               data = cbind(item$t, item$r, item$N),
               studlab = 1:item$ns)[, c(3, 6, 8, 7, 9, 4:5)]
    colnames(pairwise_observed0) <- c("study",
                                      "r1",
                                      "r2",
                                      "n1",
                                      "n2",
                                      "arm1",
                                      "arm2")

    pairwise_mod0 <- pairwise(as.list(item$t),
                              event = as.list(item$m),
                              n = as.list(item$N),
                              data = cbind(item$t, item$m, item$N),
                              studlab = 1:item$ns)[, c(6, 8)]
    colnames(pairwise_mod0) <- c("m1", "m2")

    # Ensure that t1 < t2 and correspondingly for the other elements
    treat <- treat0 <- pairwise_observed0[, 6:7]
    resp <- resp0 <- pairwise_observed0[, 2:3]
    miss <- miss0 <- pairwise_mod0[, 1:2]
    rand <- rand0 <- pairwise_observed0[, 4:5]
    for (i in seq_len(length(pairwise_observed0[, 1]))) {
      treat[i, ] <- treat0[i, order(treat0[i, ], na.last = T)]
      resp[i, ] <- resp0[i, order(treat0[i, ], na.last = T)]
      miss[i, ] <- miss0[i, order(treat0[i, ], na.last = T)]
      rand[i, ] <- rand0[i, order(treat0[i, ], na.last = T)]
    }

    pairwise_observed <- data.frame(study = pairwise_observed0$study,
                                    treat,
                                    resp,
                                    rand)
    pairwise_mod <- miss

    # The dataset to perform the unrelated trial effects model
    pairwise_data <- data.frame(pairwise_observed[, c(1, 4:5)],
                                pairwise_mod,
                                pairwise_observed[, c(6:7, 2:3)])
  }

 if (is.element(measure, c("MD", "SMD", "ROM"))) {
   contrast <- taylor_continuous(pairwise_data,
                                 measure,
                                 mean_value,
                                 var_value,
                                 rho)
  } else {
    contrast <- taylor_imor(pairwise_data,
                            mean_value,
                            var_value,
                            rho)
  }

  # Replace intervention id with their original name
  # All possible comparisons - Treat1 (non-baseline arm)
  for (i in sort(unique(unlist(contrast$t2)))) {
    contrast[contrast$t2 == i, 9] <- drug_names[i]
  }
  # Observed comparisons - Treat2 (baseline arm)
  for (i in sort(unique(unlist(contrast$t1)))) {
    contrast[contrast$t1 == i, 8] <- drug_names[i]
  }

  contrast$lower <- if (!is.element(measure, c("OR", "ROM"))) {
    round(contrast$EM - 1.95 * contrast$se.EM, 2)
  } else {
    round(exp(contrast$EM - 1.95 * contrast$se.EM), 2)
  }
  contrast$upper <- if (!is.element(measure, c("OR", "ROM"))) {
    round(contrast$EM + 1.95 * contrast$se.EM, 2)
  } else {
    round(exp(contrast$EM + 1.95 * contrast$se.EM), 2)
  }
  contrast$EM <-if (!is.element(measure, c("OR", "ROM"))) {
    round(contrast$EM, 2)
  } else {
    round(exp(contrast$EM), 2)
  }
  contrast$studlab <- rep(trial_names, na)
  contrast$comp <- paste(contrast$t2, "versus", contrast$t1)
  contrast$char1 <- rep(char[, 1], na)
  contrast$char2 <- rep(char[, 2], na)
  contrast$char3 <- rep(char[, 3], na)
  table_ute <- contrast[, c(14, 8:9, 2:7, 10:13, 16:18)]

  # Write the table with the EMs from both models as .xlsx
  if (save_xls == TRUE) {
    write_xlsx(contrast[, c(14, 8:9, 2:7, 19, 12, 13, 16:18)],
               paste0("Table Unrelated Trial Effects", ".xlsx"))
  }

  # Present plot or excel under a condition
  results <- if (item$ns < 15) {
    ggplot(contrast,
           aes(x = EM,
               y = studlab,
               xmin = lower,
               xmax = upper,
               color = as.factor(char1),
               linetype = as.factor(char2),
               shape = as.factor(char3))) +
      geom_linerange(size = 2,
                     position = position_dodge(width = 0.5)) +
      geom_vline(xintercept = ifelse(!is.element(measure, c("OR", "ROM")),
                                     0, 1),
                 lty = 2,
                 size = 1.3,
                 col = "grey53") +
      geom_point(size = 3,
                 color = "black") +
      geom_text(aes(x = EM,
                    y = studlab,
                    label = EM),
                color = "black",
                hjust = -0.3,
                vjust = -0.1,
                size = 3.5,
                check_overlap = F,
                parse = F,
                position = position_dodge(width = 0.8),
                inherit.aes = T) +
      facet_wrap(vars(comp), scales = "free_x") +
      labs(y = "",
           x = effect.measure.name(measure),
           color = colnames(char)[1],
           linetype = colnames(char)[2],
           shape = colnames(char)[3]) +
      theme_classic() +
      scale_x_continuous(trans = ifelse(!is.element(measure, c("OR", "ROM")),
                                        "identity", "log10")) +
      theme(axis.title.x = element_text(color = "black", size = 12,
                                        face = "bold"),
            axis.text.x = element_text(color = "black", size = 11),
            axis.text.y = element_text(color = "black", size = 8),
            strip.text = element_text(size = 12),
            legend.position = "bottom",
            legend.text = element_text(color = "black", size = 12),
            legend.title = element_text(color = "black", size = 12,
                                        face = "bold"))
  } else {
    write_xlsx(contrast[, c(14, 8:9, 2:7, 19, 12, 13, 16:18)],
               paste0("Table Unrelated Trial Effects", ".xlsx"))
  }

  return(list(table_unrelated.effects = knitr::kable(table_ute),
              panel_ntervalplots = results))
 }

#' End-user-ready results for unrelated trial effects model
#'
#' @description Performs the unrelated trial effects model (also known as fixed
#'   effects model) and illustrates the results of each trial and corresponding
#'   pairwise comparison.
#'
#' @param data A data-frame of a one-trial-per-row format containing arm-level
#'   data of each trial. See 'Format' in \code{\link{run_model}}.
#' @param measure Character string indicating the effect measure with values
#'   \code{"OR"}, \code{"MD"}, \code{"SMD"}, or \code{"ROM"} for the odds ratio,
#'   mean difference, standardised mean difference and ratio of means,
#'   respectively.
#' @param char A data-frame of three columns and number of rows equal to the
#'   number of trials in \code{data}. Each column refers to a
#'   trial-characteristic with \strong{nominal} elements.
#' @param drug_names A vector of labels with the name of the interventions in
#'   the order they appear in the argument \code{data}. If \code{drug_names} is
#'   not defined, the order of the interventions as they appear in \code{data}
#'   is used, instead.
#' @param trial_names A vector of labels with the name of the trials in the
#'   order they appear in the argument \code{data}. If \code{trial_names} is not
#'   defined, the order of the trials as they appear in \code{data} is used,
#'   instead.
#' @param mean_misspar A numeric value for the mean of the normal distribution
#'   of the informative missingness parameter (see 'Details'). The default
#'   argument is 0 and corresponds to the missing-at-random assumption. The same
#'   value is considered across all trials of the dataset.
#' @param var_misspar A positive non-zero number for the variance of the
#'   normal distribution of the informative missingness parameter.
#'   When the \code{measure} is \code{"OR"}, \code{"MD"}, or \code{"SMD"}
#'   the default argument is 1. When the \code{measure} is \code{"ROM"}
#'   the default argument is 0.04. The same value is considered across all
#'   trials of the dataset.
#' @param rho A numeric value in the interval [-1, 1] that indicates the
#'   correlation coefficient between two informative missingness parameters in
#'   a trial. The same value is considered across all trials of the dataset.
#'   The default argument is 0 and corresponds to uncorrelated missingness
#'   parameters.
#' @param save_xls Logical to indicate whether to export the tabulated results
#'   to an 'xlsx' file (via the \code{\link[writexl:write_xlsx]{write_xlsx}}
#'   function of the R-package
#'   \href{https://CRAN.R-project.org/package=writexl}{writexl}) to the working
#'   directory of the user. The default is \code{FALSE} (do not export).
#'
#' @return A panel of interval plots for each observed comparison in the
#'   network, when there are up to 15 trials in the \code{data}. Otherwise,
#'   \code{unrelated_effects_plot} exports a data-frame to an 'xlsx' file at
#'   the working  directory of the user. This data-frame includes the
#'   \code{data} in the long format, the within-trial effect measure and
#'   95\% confidence interval of the corresponding comparisons, the
#'   interventions compared, and the three characteristics (as defined in
#'   \code{char}).
#'   For datasets with more than 15 trials, the plot becomes cluttered and it is
#'   difficult to identify the trial-names. Hence, exporting the results in an
#'   Excel file is a viable alternative.
#'
#' @details The unrelated trial effects model may be an alternative to network
#'   meta-analysis, when the latter is not deemed appropriate (e.g., there is
#'   considerable statistical heterogeneity, or substantial intransitivity). In
#'   the presence of missing participant outcome data, the effect size and
#'   standard error are adjusted by applying the pattern-mixture model with
#'   Taylor series in trial-arms with reported missing participants (Mavridis et
#'   al., 2015; White et al., 2008). The \code{unrelated_effects_plot} function
#'   calls the \code{\link{taylor_imor}} and \code{\link{taylor_continuous}}
#'   functions (for a binary and continuous outcome, respectively) to employ
#'   pattern-mixture model with Taylor series. The \code{unrelated_effects_plot}
#'   function considers the informative missingness odds ratio in the
#'   logarithmic scale for binary outcome data (White et al., 2008), the
#'   informative missingness difference of means when \code{measure} is
#'   \code{"MD"} or \code{"SMD"}, and the informative missingness ratio of means
#'   in the logarithmic scale when \code{measure} is \code{"ROM"}
#'   (Mavridis et al., 2015).
#'
#'   The number of interval plots equals the number of observed comparisons in
#'   the network. In each interval plot, the y-axis refers to all trials of the
#'   network and x-axis refers to the selected effect measure. The odds ratio
#'   and ratio of means are calculated in the logarithmic scale but they are
#'   reported in their original scale after exponentiation.
#'
#'   \code{unrelated_effects_plot} depicts all three characteristics for each
#'   trial using different colours, line-types and point-shapes for the
#'   corresponding 95\% confidence interval and point estimate. Ideally, each
#'   characteristic should have no more than three categories; otherwise, the
#'   plot becomes cluttered. For now, the \code{unrelated_effects_plot} function
#'   uses the default colour palette, line-types and point-shapes.
#'
#' @seealso \code{\link{run_model}}, \code{\link{taylor_continuous}},
#'   \code{\link{taylor_imor}}, \code{\link[writexl:write_xlsx]{write_xlsx}}
#'
#' @references
#' Mavridis D, White IR, Higgins JP, Cipriani A, Salanti G. Allowing for
#' uncertainty due to missing continuous outcome data in pairwise and network
#' meta-analysis. \emph{Stat Med} 2015;\bold{34}(5):721--41.
#' doi: 10.1002/sim.6365
#'
#' White IR, Higgins JP, Wood AM. Allowing for uncertainty due to missing data
#' in meta-analysis--part 1: two-stage methods.
#' \emph{Stat Med} 2008;\bold{27}(5):711--27. doi: 10.1002/sim.3008
#'
#' @author {Loukia M. Spineli}
#'
#' @export
unrelated_effects_plot <- function(data,
                                   measure,
                                   char,
                                   drug_names,
                                   trial_names,
                                   mean_misspar,
                                   var_misspar,
                                   rho,
                                   save_xls) {

  item <- data_preparation(data, measure)
  na <- as.vector(do.call(
    rbind, lapply(1:item$ns, function(i) dim(combn(item$na[i], 2))[2])))

  # Default arguments
  char <- if (missing(char)) {
    stop("The argument 'char' needs to be defined", call. = FALSE)
  } else {
    char
  }
  drug_names <- if (missing(drug_names)) {
    aa <- "The argument 'drug_names' has not been defined."
    bb <- "The intervention ID, as specified in 'data' is used, instead."
    message(paste(aa, bb))
    as.character(1:item$nt)
  } else {
    drug_names
  }
  trial_names <- if (missing(trial_names)) {
    aa <- "The argument 'trial_names' has not been defined."
    bb <- "The trial ID, as specified in the argument 'data' is used, instead."
    message(paste(aa, bb))
    as.character(1:item$ns)
  } else {
    trial_names
  }
  mean_misspar <- ifelse(missing(mean_misspar), 0, mean_misspar)
  var_misspar <- ifelse(missing(var_misspar) &
                        (is.element(measure, c("OR", "MD", "SMD"))), 1,
                      ifelse(missing(var_misspar) & measure == "ROM", 0.2^2,
                             var_misspar))
  rho <- ifelse(missing(rho), 0, rho)
  save_xls <- if (missing(save_xls)) {
    FALSE
  } else {
    save_xls
  }

  # Turn into contrast-level data ('netmeta')
  if (is.element(measure, c("MD", "SMD", "ROM"))) {
    pairwise_observed <-
      pairwise(as.list(item$t),
               mean = as.list(item$y0),
               sd = as.list(item$sd0),
               n = as.list(item$N),
               data = cbind(item$t, item$y0, item$sd0, item$N),
               studlab = 1:item$ns)[, c(3:5, 7, 10, 8, 11, 6, 9)]
    colnames(pairwise_observed) <- c("study",
                                      "arm1",
                                      "arm2",
                                      "y1",
                                      "y2",
                                      "sd1",
                                      "sd2",
                                      "n1",
                                      "n2")

    pairwise_mod <- pairwise(as.list(item$t),
                              mean = as.list(item$y0),
                              sd = as.list(item$sd0),
                              n = as.list(item$m),
                              data = cbind(item$t, item$y0, item$sd0, item$m),
                              studlab = 1:item$ns)[, c(6, 9)]
    colnames(pairwise_mod) <- c("m1", "m2")

    pairwise_data <- data.frame(pairwise_observed[, c(1, 4:7)],
                                pairwise_mod,
                                pairwise_observed[, c(8:9, 2:3)])
  } else {
    pairwise_observed <-
      pairwise(as.list(item$t),
               event = as.list(item$r),
               n = as.list(item$N),
               data = cbind(item$t, item$r, item$N),
               studlab = 1:item$ns)[, c(3, 6, 8, 7, 9, 4:5)]
    colnames(pairwise_observed) <- c("study",
                                      "r1",
                                      "r2",
                                      "n1",
                                      "n2",
                                      "arm1",
                                      "arm2")

    pairwise_mod <- pairwise(as.list(item$t),
                              event = as.list(item$m),
                              n = as.list(item$N),
                              data = cbind(item$t, item$m, item$N),
                              studlab = 1:item$ns)[, c(6, 8)]
    colnames(pairwise_mod) <- c("m1", "m2")

    # The dataset to perform the unrelated trial effects model
    pairwise_data <- data.frame(pairwise_observed[, 1:3],
                                pairwise_mod,
                                pairwise_observed[, 4:7])
  }

 if (is.element(measure, c("MD", "SMD", "ROM"))) {
   contrast <- taylor_continuous(pairwise_data,
                                 measure,
                                 mean_misspar,
                                 var_misspar,
                                 rho)
  } else {
    contrast <- taylor_imor(pairwise_data,
                            mean_misspar,
                            var_misspar,
                            rho)
  }

  # Replace intervention id with their original name
  # All possible comparisons - Treat1 (non-baseline arm)
  if (!is.element(measure, c("OR", "ROM"))) {
    for (i in sort(unique(unlist(contrast$t2)))) {
      contrast[contrast$t2 == i, 11] <- drug_names[i]
    }
    # Observed comparisons - Treat2 (baseline arm)
    for (i in sort(unique(unlist(contrast$t1)))) {
      contrast[contrast$t1 == i, 10] <- drug_names[i]
    }
  } else {
    for (i in sort(unique(unlist(contrast$t2)))) {
      contrast[contrast$t2 == i, 9] <- drug_names[i]
    }
    # Observed comparisons - Treat2 (baseline arm)
    for (i in sort(unique(unlist(contrast$t1)))) {
      contrast[contrast$t1 == i, 8] <- drug_names[i]
    }
  }

  contrast$lower <- if (!is.element(measure, c("OR", "ROM"))) {
    round(contrast$EM - (1.95 * contrast$se.EM), 2)
  } else {
    round(exp(contrast$EM - (1.95 * contrast$se.EM)), 2)
  }
  contrast$upper <- if (!is.element(measure, c("OR", "ROM"))) {
    round(contrast$EM + (1.95 * contrast$se.EM), 2)
  } else {
    round(exp(contrast$EM + (1.95 * contrast$se.EM)), 2)
  }
  contrast$EM <- if (!is.element(measure, c("OR", "ROM"))) {
    round(contrast$EM, 2)
  } else {
    round(exp(contrast$EM), 2)
  }
  contrast$studlab <- rep(trial_names, na)
  contrast$comp <- paste(contrast$t2, "versus", contrast$t1)
  contrast$char1 <- rep(char[, 1], na)
  contrast$char2 <- rep(char[, 2], na)
  contrast$char3 <- rep(char[, 3], na)
  table_ute <- if (is.element(measure, c("OR", "ROM"))) {
    contrast[, c(14, 8:9, 2:7, 15, 10:13, 16:18)]
  } else {
    contrast[, c(16, 10:11, 2:9, 17, 12:15, 18:20)]
  }


  # Write the table with the EMs from both models as .xlsx
  if (save_xls == TRUE) {
    write_xlsx(table_ute,
               paste0("Table Unrelated Trial Effects", ".xlsx"))
  }

  # Present plot or excel under a condition
  results <- if (item$ns <= 15) {
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
                check_overlap = FALSE,
                parse = FALSE,
                position = position_dodge(width = 0.8),
                inherit.aes = TRUE) +
      facet_wrap(vars(comp), scales = "free_x") +
      labs(y = "",
           x = effect_measure_name(measure, lower = FALSE),
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
   "Plot will *not* be printed for network with more than 15 trials."
  }

  return(list(table_unrelated_effects =
                knitr::kable(table_ute,
                             caption =
                               "Results on unrelated trial effects model"),
              results = results))
 }

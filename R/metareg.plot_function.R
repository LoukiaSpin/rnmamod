#' End-user-ready results for network meta-regression
#'
#' @description Illustrates the effect estimates, predictions and regression
#'   coefficients of comparisons with a specific comparator and also exports
#'   these results to an Excel file.
#'
#' @param full An object of S3 class \code{\link{run_model}}. See 'Value' in
#'   \code{\link{run_model}}.
#' @param reg An object of S3 class \code{\link{run_metareg}}. See 'Value' in
#'   \code{\link{run_metareg}}.
#' @param compar A character to indicate the comparator intervention. It must be
#'   any name found in \code{drug_names}.
#' @param cov_value A list of two elements in the following order: a number
#'   for the covariate value of interest (see 'Arguments' in
#'   \code{\link{run_metareg}}), and a character to indicate the name of
#'   the covariate. See also 'Details'.
#' @param drug_names A vector of labels with the name of the interventions in
#'   the order they appear in the argument \code{data} of
#'   \code{\link{run_model}}. If \code{drug_names} is not defined,
#'   the order of the interventions as they appear in \code{data} is used,
#'   instead.
#' @param save_xls Logical to indicate whether to export the tabulated results
#'   to an 'xlsx' file (via the \code{\link[writexl]{write_xlsx}}
#'   function of the R-package
#'   \href{https://CRAN.R-project.org/package=writexl}{writexl}) to the working
#'   directory of the user. The default is \code{FALSE} (do not export).
#'
#' @return \code{metareg_plot} prints on the R console a message on the most
#'   parsimonious model (if any) based on the DIC (in red text). Furthermore,
#'   the function returns the following list of elements:
#'   \tabular{ll}{
#'    \code{table_estimates} \tab The posterior mean, and 95\% credible interval
#'    of the summary effect measure (according to the argument \code{measure}
#'    defined in \code{\link{run_model}}) for each comparison with the selected
#'    intervention under network meta-analysis and network meta-regression
#'    based on the specified \code{cov_value}.\cr
#'    \tab \cr
#'    \code{table_predictions} \tab The posterior mean, and 95\% predictive
#'    interval of the summary effect measure (according to the argument
#'    \code{measure} defined in \code{\link{run_model}}) for each comparison
#'    with the selected intervention under network meta-analysis and network
#'    meta-regression based on the specified \code{cov_value}.\cr
#'    \tab \cr
#'    \code{table_model_assessment} \tab The DIC, total residual deviance,
#'    number of effective parameters, and the posterior median and 95\% credible
#'    interval of between-trial standard deviation (\emph{tau}) under each model
#'    (Spiegelhalter et al., (2002)). When a fixed-effect model has been
#'    performed, \code{metareg_plot} does not return results on \emph{tau}.\cr
#'    \tab \cr
#'    \code{table_regression_coeffients} \tab The posterior mean and 95\%
#'    credible interval of the regression coefficient(s).\cr
#'    \tab \cr
#'    \code{interval_plots} \tab The panel of forest plots on the estimated and
#'    predicted effect sizes of comparisons with the selected intervention under
#'    network meta-analysis and network meta-regression based on the specified
#'    \code{cov_value}. See 'Details' and 'Value' in
#'    \code{\link{forestplot_metareg}}.\cr
#'    \tab \cr
#'    \code{sucra_scatterplot} \tab A scatterplot of the SUCRA values from the
#'    network meta-analysis against the SUCRA values from the network
#'    meta-regression based on the specified \code{cov_value}. See 'Details'
#'    and 'Value' in \code{\link{scatterplot_sucra}}.
#' }
#'
#' @details The deviance information criterion (DIC) of the network
#'   meta-analysis model is compared with the DIC of the network meta-regression
#'   model. If the difference in DIC exceeds 5, the network meta-regression
#'   model is preferred; if the difference in DIC is less than -5, the network
#'   meta-analysis model is preferred; otherwise, there is little to choose
#'   between the compared models.
#'
#'   When the covariate is binary, specify in the second element of
#'   \code{cov_value} the name of the level for which the output will be
#'   created.
#'
#'   Furthermore, \code{metareg_plot} exports all tabulated results to separate
#'   'xlsx' files (via the \code{\link[writexl]{write_xlsx}} function of the
#'   R-package \href{https://CRAN.R-project.org/package=writexl}{writexl})
#'   to the working directory of the user.
#'
#'   \code{metareg_plot} can be used only for a network of interventions. In the
#'   case of two interventions, the execution of the function will be stopped
#'   and an error message will be printed on the R console.
#'
#' @author {Loukia M. Spineli}
#'
#' @seealso \code{\link{forestplot_metareg}}, \code{\link{run_metareg}},
#'   \code{\link{run_model}}, \code{\link{scatterplot_sucra}},
#'   \href{https://CRAN.R-project.org/package=writexl}{writexl}
#'
#' @references
#' Salanti G, Ades AE, Ioannidis JP. Graphical methods and numerical summaries
#' for presenting results from multiple-treatment meta-analysis: an overview and
#' tutorial. \emph{J Clin Epidemiol} 2011;\bold{64}(2):163--71.
#' \doi{10.1016/j.jclinepi.2010.03.016}
#'
#' Spiegelhalter DJ, Best NG, Carlin BP, van der Linde A. Bayesian measures of
#' model complexity and fit. \emph{J R Stat Soc B} 2002;\bold{64}:583--616.
#' \doi{10.1111/1467-9868.00353}
#'
#' @examples
#' data("nma.baker2009")
#'
#' # Read results from 'run_model' (using the default arguments)
#' res <- readRDS(system.file('extdata/res_baker.rds', package = 'rnmamod'))
#'
#' # Read results from 'run_metareg' (exchangeable structure)
#' reg <- readRDS(system.file('extdata/reg_baker.rds', package = 'rnmamod'))
#'
#' # Publication year as the covariate
#' pub_year <- c(1996, 1998, 1999, 2000, 2000, 2001, rep(2002, 5), 2003, 2003,
#'               rep(2005, 4), 2006, 2006, 2007, 2007)
#'
#' # The names of the interventions in the order they appear in the dataset
#' interv_names <- c("placebo", "budesonide", "budesonide plus formoterol",
#'                   "fluticasone", "fluticasone plus salmeterol",
#'                   "formoterol", "salmeterol", "tiotropium")
#'
#' # Plot the results from both models for all comparisons with salmeterol and
#' # publication year 2000
#' metareg_plot(full = res,
#'              reg = reg,
#'              compar = "salmeterol",
#'              cov_value = list(2000, "publication year"),
#'              drug_names = interv_names)
#'
#' @export
metareg_plot <- function(full,
                         reg,
                         compar,
                         cov_value,
                         drug_names,
                         save_xls) {

  if (length(unique(reg$covariate)) < 3 &
      !is.element(cov_value[[1]], reg$covariate)) {
    aa <- "The first element of the argument 'cov_value' is out of the value"
    stop(paste(aa, "range of the analysed covariate"), call. = FALSE)
  } else if (length(unique(reg$covariate)) > 2 &
             (cov_value[[1]] < min(reg$covariate) |
              cov_value[[1]] > max(reg$covariate))) {
    aa <- "The first element of the argument 'cov_value' is out of the value"
    stop(paste(aa, "range of the analysed covariate"), call. = FALSE)
  }


  if (length(drug_names) < 3) {
    stop("This function is *not* relevant for a pairwise meta-analysis",
         call. = FALSE)
  }

  cov_value <- if (missing(cov_value)) {
    stop("The argument 'cov_value' has not been defined", call. = FALSE)
  } else if (length(cov_value) < 2) {
    aa <- "The argument 'cov_value' must be a list with elements a number and"
    stop(paste(aa, "a character"), call. = FALSE)
  } else if (length(cov_value) == 2) {
    cov_value
  }

  save_xls <- if (missing(save_xls)) {
    FALSE
  } else {
    save_xls
  }

  compar <- if (missing(compar)) {
    stop("The argument 'compar' has not been defined", call. = FALSE)
  } else if (!is.element(compar, drug_names)) {
    stop("The value of the argument 'compar' is not found in the 'drug_names'",
         call. = FALSE)
  } else if (is.element(compar, drug_names)) {
    compar
  }


  covariate <- if (length(unique(reg$covariate)) < 3) {
    unique(reg$covariate)
  } else {
    reg$covariate
  }

  cov_val <- ifelse(length(unique(covariate)) < 3, cov_value[[1]],
                    cov_value[[1]] - mean(covariate))

  drug_names <- if (missing(drug_names)) {
    aa <- "The argument 'drug_names' has not been defined."
    bb <- "The intervention ID, as specified in 'data' is used as"
    cc <- "intervention names"
    message(cat(paste0("\033[0;", col = 32, "m", aa, " ", bb, " ", cc,
                       "\033[0m", "\n")))
    as.character(seq_len(length(full$SUCRA[, 1])))
  } else {
    drug_names
  }

  model <- full$model
  measure <- effect_measure_name(full$measure)

  # Posterior results on the SUCRA value under NMA
  sucra_full <- round(full$SUCRA, 2)
  sucra_full_order <- round(full$SUCRA, 2)[order(sucra_full[, 1],
                                                 decreasing = TRUE), ]

  # Sort the drugs by their NMA-SUCRA in decreasing order
  drug_names_sorted <- drug_names[order(sucra_full[, 1], decreasing = TRUE)]

  # A matrix with all possible comparisons in the network
  poss_pair_comp1 <- data.frame(exp = t(combn(drug_names, 2))[, 2],
                                comp = t(combn(drug_names, 2))[, 1])
  poss_pair_comp2 <- data.frame(exp = t(combn(drug_names, 2))[, 1],
                                comp = t(combn(drug_names, 2))[, 2])
  poss_pair_comp <- rbind(poss_pair_comp1, poss_pair_comp2)

  # Posterior results on NMA for comparisons with the selected intervention
  em_ref00_nma <- cbind(rbind(data.frame(mean = full$EM[, 1],
                                         lower = full$EM[, 3],
                                         upper = full$EM[, 7]),
                              data.frame(mean = full$EM[, 1] * (-1),
                                         lower = full$EM[, 7] * (-1),
                                         upper = full$EM[, 3] * (-1))),
                        poss_pair_comp)
  em_subset_nma <- subset(em_ref00_nma, em_ref00_nma[5] == compar)
  em_ref0_nma <- rbind(em_subset_nma[, 1:3], c(rep(NA, 3)))
  sucra_full_new <- data.frame(sucra_full[, 1], drug_names)[
    order(match(data.frame(sucra_full[, 1], drug_names)[, 2],
                em_subset_nma[, 4])), 1]
  em_ref_nma <- em_ref0_nma[order(sucra_full_new, decreasing = TRUE), ]

  # Posterior mean of regression coefficients for all pairwise comparisons
  if (is.element(reg$covar_assumption, c("exchangeable", "independent"))) {
    beta00 <- cbind(rbind(data.frame(mean = reg$beta_all[, 1],
                                     lower = reg$beta_all[, 3],
                                     upper = reg$beta_all[, 7]),
                          data.frame(mean = reg$beta_all[, 1] * (-1),
                                     lower = reg$beta_all[, 7] * (-1),
                                     upper = reg$beta_all[, 3] * (-1))),
                    poss_pair_comp)
    beta_all_subset <- subset(beta00, beta00[5] == compar)
    beta0 <- rbind(beta_all_subset[, 1:3], c(rep(NA, 3)))
    beta <- beta0[order(sucra_full_new, decreasing = TRUE), ]
    rownames(beta) <- NULL
  } else {
    beta <- reg$beta[1, c(1, 3, 7)]
  }

  # Effect size of all possible pairwise comparisons (NMR)
  par_mean <- as.vector(c(reg$EM[, 1] + reg$beta_all[, 1] * cov_val,
                          (reg$EM[, 1] * (-1)) +
                            (reg$beta_all[, 1] * (-1) * cov_val)))
  par_sd <- as.vector(c(sqrt(((reg$EM[, 2])^2) +
                               ((reg$beta_all[, 2] * cov_val)^2)),
                        sqrt(((reg$EM[, 2])^2) +
                               ((reg$beta_all[, 2] * cov_val)^2))))
  em_ref00_nmr <- cbind(mean = par_mean,
                        lower = par_mean - 1.96 * par_sd,
                        upper = par_mean + 1.96 * par_sd,
                        poss_pair_comp)
  em_subset_nmr <- subset(em_ref00_nmr, em_ref00_nmr[5] == compar)
  em_ref0_nmr <- rbind(em_subset_nmr[, 1:3], c(rep(NA, 3)))
  em_ref_nmr <- em_ref0_nmr[order(sucra_full_new, decreasing = TRUE), ]
  rownames(em_ref_nma) <- rownames(em_ref_nmr) <- NULL

  # Posterior results on the predicted estimates of comparisons with the
  # selected comparator as reference
  if (model == "RE") {
    pred_ref00_nma <- cbind(rbind(data.frame(mean = full$EM_pred[, 1],
                                             lower = full$EM_pred[, 3],
                                             upper = full$EM_pred[, 7]),
                                  data.frame(mean = full$EM_pred[, 1] * (-1),
                                             lower = full$EM_pred[, 7] * (-1),
                                             upper = full$EM_pred[, 3] * (-1))),
                            poss_pair_comp)
    pred_subset_nma <- subset(pred_ref00_nma, pred_ref00_nma[5] == compar)
    pred_ref0_nma <- rbind(pred_subset_nma[, 1:3], c(rep(NA, 3)))
    par_mean <- as.vector(c(reg$EM_pred[, 1] + reg$beta_all[, 1] * cov_val,
                            (reg$EM_pred[, 1] * (-1)) +
                              (reg$beta_all[, 1] * (-1) * cov_val)))
    par_sd <- as.vector(c(sqrt(((reg$EM_pred[, 2])^2) +
                                 ((reg$beta_all[, 2] * cov_val)^2)),
                          sqrt(((reg$EM_pred[, 2])^2) +
                                 ((reg$beta_all[, 2] * cov_val)^2))))
    pred_ref00_nmr <-  cbind(data.frame(mean = par_mean,
                                        lower = par_mean - 1.96 * par_sd,
                                        upper = par_mean + 1.96 * par_sd),
                             poss_pair_comp)
    pred_subset_nmr <- subset(pred_ref00_nmr, pred_ref00_nmr[5] == compar)
    pred_ref0_nmr <- rbind(pred_subset_nmr[, 1:3], c(rep(NA, 3)))

    # Sort by SUCRA in decreasing order and remove the reference intervention
    pred_ref_nma <- pred_ref0_nma[order(sucra_full_new, decreasing = TRUE), ]
    pred_ref_nmr <- pred_ref0_nmr[order(sucra_full_new, decreasing = TRUE), ]
    rownames(pred_ref_nma) <- rownames(pred_ref_nmr) <- NULL
  }

  if (!is.element(measure, c("Odds ratio", "Ratio of means")) & model == "RE") {
    em_ref_nma <- em_ref_nma
    em_ref_nmr <- em_ref_nmr
    pred_ref_nma <- pred_ref_nma
    pred_ref_nmr <- pred_ref_nmr
    beta <- beta

    cri_est_nma <- paste0("(", round(em_ref_nma[, 2], 2), ",", " ",
                          round(em_ref_nma[, 3], 2), ")",
                          ifelse(em_ref_nma[, 2] > 0 |
                                   em_ref_nma[, 3] < 0, "*", " "))
    cri_est_nmr <- paste0("(", round(em_ref_nmr[, 2], 2), ",", " ",
                          round(em_ref_nmr[, 3], 2), ")",
                          ifelse(em_ref_nmr[, 2] > 0 |
                                   em_ref_nmr[, 3] < 0, "*", " "))
    cri_pred_nma <- paste0("(", round(pred_ref_nma[, 2], 2), ",", " ",
                           round(pred_ref_nma[, 3], 2), ")",
                           ifelse(pred_ref_nma[, 2] > 0 |
                                    pred_ref_nma[, 3] < 0, "*", " "))
    cri_pred_nmr <- paste0("(", round(pred_ref_nmr[, 2], 2), ",", " ",
                           round(pred_ref_nmr[, 3], 2), ")",
                           ifelse(pred_ref_nmr[, 2] > 0 |
                                    pred_ref_nmr[, 3] < 0, "*", " "))
    cri_beta <- if (is.element(reg$covar_assumption, c("exchangeable",
                                                       "independent"))) {
      paste0("(", round(beta[, 2], 2), ",", " ", round(beta[, 3], 2), ")",
             ifelse(beta[, 2] > 0 | beta[, 3] < 0, "*", " "))
    } else {
      paste0("(", round(reg$beta[2], 2), ",", " ", round(reg$beta[3], 2), ")",
             ifelse(reg$beta[2] > 0 | reg$beta[3] < 0, "*", " "))
    }

  } else if (is.element(measure, c("Odds ratio", "Ratio of means")) &
             model == "RE") {
    em_ref_nma <- exp(em_ref_nma)
    em_ref_nmr <- exp(em_ref_nmr)
    pred_ref_nma <- exp(pred_ref_nma)
    pred_ref_nmr <- exp(pred_ref_nmr)
    beta <- exp(beta)

    cri_est_nma <- paste0("(", round(em_ref_nma[, 2], 2), ",", " ",
                          round(em_ref_nma[, 3], 2), ")",
                          ifelse(em_ref_nma[, 2] > 1 |
                                   em_ref_nma[, 3] < 1, "*", " "))
    cri_est_nmr <- paste0("(", round(em_ref_nmr[, 2], 2), ",", " ",
                          round(em_ref_nmr[, 3], 2), ")",
                          ifelse(em_ref_nmr[, 2] > 1 |
                                   em_ref_nmr[, 3] < 1, "*", " "))
    cri_pred_nma <- paste0("(", round(pred_ref_nma[, 2], 2), ",", " ",
                           round(pred_ref_nma[, 3], 2), ")",
                           ifelse(pred_ref_nma[, 2] > 1 |
                                    pred_ref_nma[, 3] < 1, "*", " "))
    cri_pred_nmr <- paste0("(", round(pred_ref_nmr[, 2], 2), ",", " ",
                           round(pred_ref_nmr[, 3], 2), ")",
                           ifelse(pred_ref_nmr[, 2] > 1 |
                                    pred_ref_nmr[, 3] < 1, "*", " "))
    cri_beta <- if (is.element(reg$covar_assumption, c("exchangeable",
                                                       "independent"))) {
      paste0("(", round(beta[, 2], 2), ",", " ", round(beta[, 3], 2), ")",
             ifelse(beta[, 2] > 1 | beta[, 3] < 1, "*", " "))
    } else {
      paste0("(", round(beta[2], 2), ",", " ", round(beta[3], 2), ")",
             ifelse(beta[2] > 1 | beta[3] < 1, "*", " "))
    }

  } else if (!is.element(measure, c("Odds ratio", "Ratio of means")) &
             model == "FE") {
    em_ref_nma <- em_ref_nma
    em_ref_nmr <- em_ref_nmr
    beta <- beta

    cri_est_nma <- paste0("(", round(em_ref_nma[, 2], 2), ",", " ",
                          round(em_ref_nma[, 3], 2), ")",
                          ifelse(em_ref_nma[, 2] > 0 |
                                   em_ref_nma[, 3] < 0, "*", " "))
    cri_est_nmr <- paste0("(", round(em_ref_nmr[, 2], 2), ",", " ",
                          round(em_ref_nmr[, 3], 2), ")",
                          ifelse(em_ref_nmr[, 2] > 0 |
                                   em_ref_nmr[, 3] < 0, "*", " "))
    cri_beta <- if (is.element(reg$covar_assumption, c("exchangeable",
                                                       "independent"))) {
      paste0("(", round(beta[, 2], 2), ",", " ", round(beta[, 3], 2), ")",
             ifelse(beta[, 2] > 0 | beta[, 3] < 0, "*", " "))
    } else {
      paste0("(", round(reg$beta[2], 2), ",", " ", round(reg$beta[3], 2), ")",
             ifelse(reg$beta[2] > 0 | reg$beta[3] < 0, "*", " "))
    }
  } else if (is.element(measure, c("Odds ratio", "Ratio of means")) &
             model == "FE") {
    em_ref_nma <- exp(em_ref_nma)
    em_ref_nmr <- exp(em_ref_nmr)
    beta <- exp(beta)

    cri_est_nma <- paste0("(", round(em_ref_nma[, 2], 2), ",", " ",
                          round(em_ref_nma[, 3], 2), ")",
                          ifelse(em_ref_nma[, 2] > 1 |
                                   em_ref_nma[, 3] < 1, "*", " "))
    cri_est_nmr <- paste0("(", round(em_ref_nmr[, 2], 2), ",", " ",
                          round(em_ref_nmr[, 3], 2), ")",
                          ifelse(em_ref_nmr[, 2] > 1 |
                                   em_ref_nmr[, 3] < 1, "*", " "))
    cri_beta <- if (is.element(reg$covar_assumption, c("exchangeable",
                                                       "independent"))) {
      paste0("(", round(beta[, 2], 2), ",", " ", round(beta[, 3], 2), ")",
             ifelse(beta[, 2] > 1 | beta[, 3] < 1, "*", " "))
    } else {
      paste0("(", round(beta[2], 2), ",", " ", round(beta[3], 2), ")",
             ifelse(beta[2] > 1 | beta[3] < 1, "*", " "))
    }

  }

  # Posterior results on between-trial standard deviation under NMA
  tau_nma <- if (model == "RE") {
    round(full$tau, 2)
  }

  # Posterior results on between-trial standard deviation under meta-regression
  tau_nmr <- if (model == "RE") {
    round(reg$tau, 2)
  }

  # Posterior mean of model assessment measures under NMA
  model_assess_nma <- round(full$model_assessment, 2)

  # Posterior mean of model assessment measures under meta-regression
  model_assess_nmr <- round(reg$model_assessment, 2)


  # The 95% CrIs of the between-trial standard deviation under both models
  if (model == "RE") {
    cri_tau_nma <- paste0("(", tau_nma[, 3], ",", " ", tau_nma[, 7], ")")
    cri_tau_nmr <- paste0("(", tau_nmr[, 3], ",", " ", tau_nmr[, 7], ")")
  }

  # Model assessment and between-trial standard deviation for both models
  if (model == "RE") {
    table_model_assess <- data.frame(c("Network meta-analysis",
                                       "Meta-regression"),
                                     rbind(model_assess_nma[c(1, 3, 2)],
                                           model_assess_nmr[c(1, 3, 2)]),
                                     rbind(cbind(tau_nma[, 5],
                                                 tau_nma[, 2],
                                                 cri_tau_nma),
                                           cbind(tau_nmr[, 5],
                                                 tau_nmr[, 2],
                                                 cri_tau_nmr)))
    colnames(table_model_assess) <- c("Analysis",
                                      "DIC", "Mean deviance", "pD",
                                      "Median tau", "SD tau", "95% CrI tau")
  } else {
    table_model_assess <- data.frame(c("Network meta-analysis",
                                       "Meta-regression"),
                                     rbind(model_assess_nma[c(1, 3, 2)],
                                           model_assess_nmr[c(1, 3, 2)]))
    colnames(table_model_assess) <- c("Analysis",
                                      "DIC", "Mean deviance", "pD")
  }
  message(ifelse(model_assess_nma[1] - model_assess_nmr[1] > 5,
                 "NMR preferred when accounting for model fit and complexity",
                 ifelse(
                   model_assess_nma[1] - model_assess_nmr[1] < -5,
                   "NMA preferred when accounting for model fit and complexity",
                        "There is little to choose between the two models")))


  # Tabulate results on comparisons with the reference (both models)
  if (model == "RE") {
    est_both_models <- na.omit(data.frame(drug_names_sorted,
                                          round(em_ref_nma[, 1], 2),
                                          cri_est_nma,
                                          round(em_ref_nmr[, 1], 2),
                                          cri_est_nmr))
    pred_both_models <- na.omit(data.frame(drug_names_sorted,
                                           round(pred_ref_nma[, 1], 2),
                                           cri_pred_nma,
                                           round(pred_ref_nmr[, 1], 2),
                                           cri_pred_nmr))
    colnames(est_both_models) <- colnames(pred_both_models) <-
      c(paste("versus", compar),
        "Mean NMA", "95% CrI NMA", "Mean NMR", "95% CrI NMR")
    rownames(est_both_models) <- rownames(pred_both_models) <- NULL
  } else {
    est_both_models <- na.omit(data.frame(drug_names_sorted,
                                          round(em_ref_nma[, 1], 2),
                                          cri_est_nma,
                                          round(em_ref_nmr[, 1], 2),
                                          cri_est_nmr))
    colnames(est_both_models) <- c(paste("versus", compar),
        "Mean NMA", "95% CrI NMA", "Mean NMR", "95% CrI NMR")
    rownames(est_both_models) <- NULL

  }

  # Results on the regression coefficient
  if (is.element(reg$covar_assumption, c("exchangeable", "independent"))) {
    reg_coeff <- na.omit(data.frame(drug_names_sorted, round(beta[, 1], 2),
                                    cri_beta))
    colnames(reg_coeff) <- c(paste("versus", compar), "Mean beta",
                             "95% CrI beta")
  } else {
    reg_coeff <- data.frame(round(beta[1], 2), cri_beta)
    colnames(reg_coeff) <- c("Mean beta", "95% CrI beta")
  }
  rownames(reg_coeff) <- NULL

  # Forest plots of reference-comparisons on effect estimate
  forest_plots <- forestplot_metareg(full, reg, compar, cov_value, drug_names)
  sucra_scatterplot <- scatterplot_sucra(full, reg, cov_value, drug_names)

  # Write all tables as .xlsx
  if (save_xls == TRUE & model == "RE") {
    write_xlsx(est_both_models, paste0("Table NMA vs NMR_Estimation", ".xlsx"))
    write_xlsx(pred_both_models, paste0("Table NMA vs NMR_Prediction", ".xlsx"))
    write_xlsx(table_model_assess, paste0("Table Model Assessment_NMA vs NMR",
                                          ".xlsx"))
    if (is.element(reg$covar_assumption, c("exchangeable", "independent"))) {
      write_xlsx(reg_coeff, paste0("Table NMA vs NMR_Coefficient", ".xlsx"))
    }
  } else if (save_xls == TRUE & model == "FE") {
    write_xlsx(est_both_models, paste0("Table NMA vs NMR_Estimation", ".xlsx"))
    write_xlsx(table_model_assess, paste0("Table Model Assessment_NMA vs NMR",
                                          ".xlsx"))
    if (is.element(reg$covar_assumption, c("exchangeable", "independent"))) {
      write_xlsx(reg_coeff, paste0("Table NMA vs NMR_Coefficient", ".xlsx"))
    }
  }

  results <- if (model == "RE") {
    list(table_estimates = knitr::kable(est_both_models),
         table_predictions = knitr::kable(pred_both_models),
         table_model_assessment = knitr::kable(table_model_assess),
         table_regression_coeffients = knitr::kable(reg_coeff),
         interval_plots = suppressWarnings(ggarrange(forest_plots)),
         sucra_scatterplot = sucra_scatterplot)
  } else {
    list(table_estimates = knitr::kable(est_both_models),
         table_model_assessment = knitr::kable(table_model_assess),
         table_regression_coeffients = knitr::kable(reg_coeff),
         interval_plots = suppressWarnings(ggarrange(forest_plots)),
         sucra_scatterplot = sucra_scatterplot)
  }

  return(results)
}

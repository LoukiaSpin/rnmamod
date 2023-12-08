#' End-user-ready results for the unrelated mean effects model
#'
#' @description \code{ume_plot} hosts a toolkit of functions that facilitates
#'   the comparison of the consistency model (via \code{\link{run_model}}) with
#'   the unrelated mean effects model (via \code{\link{run_ume}}) regarding the
#'   posterior summaries of the summary effect size for the pairwise comparisons
#'   observed in the network, the between-trial standard deviation (\emph{tau})
#'   and model assessment parameters.
#'
#' @param full An object of S3 class \code{\link{run_model}}. See 'Value' in
#'   \code{\link{run_model}}.
#' @param ume An object of S3 class \code{\link{run_ume}}. See 'Value' in
#'   \code{\link{run_ume}}.
#' @param drug_names A vector of labels with the name of the interventions in
#'   the order they appear in the argument \code{data} of
#'   \code{\link{run_model}}. If \code{drug_names} is not defined,
#'   the order of the interventions as they appear in \code{data} is used,
#'   instead.
#' @param save_xls Logical to indicate whether to export the tabulated results
#'   to an 'xlsx' file (via the \code{\link[writexl:write_xlsx]{write_xlsx}}
#'   function of the R-package
#'   \href{https://CRAN.R-project.org/package=writexl}{writexl}) to the working
#'   directory of the user. The default is \code{FALSE} (do not export).
#'
#' @return \code{ume_plot} prints on the R console a message on the most
#'   parsimonious model (if any) based on the DIC (red text). Then, the function
#'   returns the following list of elements:
#'   \item{table_effect_size}{The posterior median, posterior standard
#'   deviation, and 95\% credible interval of the summary effect size for each
#'   pairwise comparison observed in the network under the consistency model
#'   and the unrelated mean effects model.}
#'   \item{table_model_assessment}{The DIC, number of effective
#'   parameters, and total residual deviance under the consistency model and
#'   the unrelated mean effects model (Spiegelhalter et al., 2002).}
#'   \item{table_tau}{The posterior median and 95\% credible interval of
#'   \emph{tau} under the consistency model and the unrelated mean effects
#'   model. When a fixed-effect model has been performed, \code{ume_plot} does
#'   not return this element.}
#'   \item{scatterplots}{The scatterplot and the Bland-Altman plot on the
#'   posterior mean deviance contribution of the individual data points under
#'   the consistency model and the unrelated mean effects model. See 'Details'
#'   and 'Value' in \code{\link{scatterplots_dev}} and
#'   \code{\link{bland_altman_plot}}, respectively.}
#'   \item{levarage_plots}{The leverage plot under the consistency model
#'   and the unrelated mean effects model, separately. See 'Details' and
#'   'Value' in \code{\link{leverage_plot}}.}
#'   \item{intervalplots}{A panel of interval plots on the summary effect
#'   size under the consistency model and the unrelated mean effects model for
#'   each pairwise comparison observed in the network. See 'Details' and
#'   'Value' in \code{\link{intervalplot_panel_ume}}.}
#'
#' @details The deviance information criterion (DIC) of the consistency model is
#'   compared with the DIC of the unrelated mean effects model
#'   (Dias et al., 2013). If the difference in DIC exceeds 5, the unrelated mean
#'   effects model is preferred. If the difference in DIC is less than -5, the
#'   consistency is preferred; otherwise, there is little to choose between the
#'   compared models.
#'
#'   For a binary outcome, when \code{measure} is "RR" (relative risk) or "RD"
#'   (risk difference) in \code{\link{run_model}}, \code{ume_plot} currently
#'   presents the results from network meta-analysis and unrelated mean effects
#'   in the odds ratio for being the \strong{base-case} effect measure in
#'   \code{\link{run_model}} for a binary outcome (see also 'Details' in
#'   \code{\link{run_model}}).
#'
#'   Furthermore, \code{ume_plot} exports \code{table_effect_size} and
#'   \code{table_model_assessment} to separate 'xlsx' files (via the
#'   \code{\link[writexl:write_xlsx]{write_xlsx}} function) to the working
#'   directory of the user.
#'
#'   \code{ume_plot} can be used only for a network of interventions. In the
#'   case of two interventions, the execution of the function will be stopped
#'   and an error message will be printed on the R console.
#'
#' @author {Loukia M. Spineli}
#'
#' @seealso \code{\link{bland_altman_plot}},
#'   \code{\link{intervalplot_panel_ume}}, \code{\link{leverage_plot}},
#'   \code{\link{run_model}}, \code{\link{run_ume}},
#'   \code{\link[writexl:write_xlsx]{write_xlsx}}
#'
#' @references
#' Dias S, Welton NJ, Sutton AJ, Caldwell DM, Lu G, Ades AE. Evidence synthesis
#' for decision making 4: inconsistency in networks of evidence based on
#' randomized controlled trials.
#' \emph{Med Decis Making} 2013;\bold{33}(5):641--56.
#' doi: 10.1177/0272989X12455847
#'
#' Spiegelhalter DJ, Best NG, Carlin BP, van der Linde A. Bayesian measures of
#' model complexity and fit. \emph{J R Stat Soc B} 2002;\bold{64}(4):583--396.
#' doi: 10.1111/1467-9868.00353
#'
#' @examples
#' data("nma.liu2013")
#'
#' \donttest{
#
#' # Read results from 'run_model' (using the default arguments)
#' res <- readRDS(system.file('extdata/res_liu.rds', package = 'rnmamod'))
#'
#' # Read results from 'run_ume' (using the default arguments)
#' ume <- readRDS(system.file('extdata/ume_liu.rds', package = 'rnmamod'))
#'
#' # The names of the interventions in the order they appear in the dataset
#' interv_names <- c("placebo", "pramipexole", "serotonin-norepinephrine
#'                   reuptake inhibitor", "serotonin reuptake inhibitor",
#'                   "tricyclic antidepressant", "pergolide")
#'
#' # Plot the results from both models
#' ume_plot(full = res,
#'          ume = ume,
#'          drug_names = interv_names)
#' }
#'
#' @export
ume_plot <- function(full, ume, drug_names, save_xls) {


  if (!inherits(full, "run_model") || is.null(full)) {
    stop("'full' must be an object of S3 class 'run_model'.",
         call. = FALSE)
  }

  if (!inherits(ume, "run_ume") || is.null(ume)) {
    stop("'ume' must be an object of S3 class 'run_ume'.",
         call. = FALSE)
  }

  drug_names <- if (missing(drug_names)) {
    aa <- "The argument 'drug_names' has not been defined."
    bb <- "The intervention ID, as specified in 'data' is used, instead."
    message(paste(aa, bb))
    nt <- length(full$SUCRA[, 1])
    as.character(1:nt)
  } else {
    drug_names
  }

  if (length(drug_names) < 3) {
    stop("This function is *not* relevant for a pairwise meta-analysis.",
         call. = FALSE)
  }

  save_xls <- if (missing(save_xls)) {
    FALSE
  } else {
    save_xls
  }

  model <- full$model
  measure <- if (is.element(full$measure, c("RR", "RD"))) {
    "OR"
  } else {
    full$measure
  }

  # Posterior results on the effect estimates under
  em_full <- if (is.element(full$measure, c("RR", "RD"))) {
    full$EM_LOR[, c(5, 2:3, 7)]
  } else {
    full$EM[, c(5, 2:3, 7)]
  }
  em_full[, c(1, 3, 4)] <- if (is.element(measure, c("OR", "ROM"))) {
    exp(em_full[, c(1, 3, 4)])
  } else {
    em_full[, c(1, 3, 4)]
  }

  # Posterior results on the effect estimates under UME
  em_ume <- ume$EM[, c(5, 2:3, 7)]
  em_ume[, c(1, 3, 4)] <- if (is.element(measure, c("OR", "ROM"))) {
    exp(em_ume[, c(1, 3, 4)])
  } else {
    em_ume[, c(1, 3, 4)]
  }

  # Posterior results on between-trial standard deviation under NMA
  tau_full <- if (model == "RE") {
    full$tau
  } else {
    NA
  }

  # Posterior results on between-trial standard deviation under UME
  tau_ume <- if (model == "RE") {
    ume$tau
  } else {
    NA
  }

  # Posterior mean on deviance contribution for observed outcomes under NMA
  dev_o_full <- full$dev_o

  # Posterior mean on deviance contribution for observed outcomes under UME
  dev_o_ume <- ume$dev_o

  # Measures of model assessment: DIC, pD, total residual deviance under NAM
  model_assess_full <- full$model_assessment
  data_points <- model_assess_full[4]

  # Measures of model assessment: DIC, pD, total residual deviance under UME
  model_assess_ume <- ume$model_assessment

  # Observed comparisons in the network
  obs_comp <- ume$obs_comp

  # Possible and observed comparisons (with names)
  possible_comp <- possible_observed_comparisons(drug_names, obs_comp)

  # Keep the effect estimates according to the 'poss.pair.comp.clean' (NMA)
  em_full_clean <- format(round(
    em_full[is.element(possible_comp$poss_comp[, 4], obs_comp), ], 2),
    nsmall = 2)

  # Keep the effect estimates according to the 'poss.pair.comp.clean' (UME)
  em_ume_clean <- format(round(em_ume[, ], 2), nsmall = 2)

  # Keep the 95% credible intervals (CrI) under NMA
  cri_full_clean <- if (is.element(measure, c("OR", "ROM"))) {
    paste0("(", em_full_clean[, 3], ",", " ", em_full_clean[, 4], ")",
           ifelse(as.numeric(em_full_clean[, 3]) > 1 |
                    as.numeric(em_full_clean[, 4]) < 1, "*", " "))
  } else {
    paste0("(", em_full_clean[, 3], ",", " ", em_full_clean[, 4], ")",
           ifelse(as.numeric(em_full_clean[, 3]) > 0 |
                  as.numeric(em_full_clean[, 4]) < 0, "*", " "))
  }

  # Keep the 95% credible intervals (CrI) under UME
  cri_ume_clean <- if (is.element(measure, c("OR", "ROM"))) {
    paste0("(", em_ume_clean[, 3], ",", " ", em_ume_clean[, 4], ")",
           ifelse(as.numeric(em_ume_clean[, 3]) > 1 |
                    as.numeric(em_ume_clean[, 4]) < 1, "*", " "))
  } else {
    paste0("(", em_ume_clean[, 3], ",", " ", em_ume_clean[, 4], ")",
           ifelse(as.numeric(em_ume_clean[, 3]) > 0 |
                    as.numeric(em_ume_clean[, 4]) < 0, "*", " "))
  }

  # Create a data-frame with effect estimates on both models
  em_both <- data.frame(possible_comp$obs_comp[, 4],
                        em_full_clean[, 1:2],
                        cri_full_clean,
                        em_ume_clean[, 1:2],
                        cri_ume_clean)
  colnames(em_both) <- c("Comparison",
                         "Median NMA", "SD NMA",
                         "95% CrI NMA",
                         "Median UME", "SD UME",
                         "95% CrI UME")
  rownames(em_both) <- NULL

  # A data-frame with the measures on model assessment
  model_assessment <- data.frame(t(model_assess_full[-4]), t(model_assess_ume))
  colnames(model_assessment) <- c("Full NMA", "UME model")
  message(ifelse(model_assess_full[1] - model_assess_ume[1] > 5,
                 "UME preferred when accounting for model fit and complexity",
                 ifelse(model_assess_full[1] -
                 model_assess_ume[1] < -5,
                 "NMA preferred when accounting for model fit and complexity.",
                        "There is little to choose between the two models.")))

  # Data-frame  on between-trial standard deviation
  if (model == "RE") {
    between_trial_sd <- rbind(tau_full[c(5, 3, 7)], tau_ume[c(5, 3, 7)])
    colnames(between_trial_sd) <- c("Median", "Lower 95% CrI", "Upper 95% CrI")
    rownames(between_trial_sd) <- c("Full NMA", "UME model")
  } else {
    between_trial_sd <- NA
  }

  # Scatterplot on the deviance contribution of NMA versus UME
  scatterplot_o <- scatterplots_dev(dev_o_full[, 1],
                                    dev_o_ume[, 1],
                                    colour = "#D55E00")

  # Bland-Altman plot on the deviance contribution of consistency versus UME
  bland_altman_observed <- bland_altman_plot(dev_o_full[, 1],
                                             dev_o_ume[, 1],
                                             colour = "#D55E00")

  # Bring together all four plots
  scatterplots <- ggpubr::ggarrange(scatterplot_o,
                                    bland_altman_observed,
                                    nrow = 1,
                                    ncol = 2)

  # Leverage plots
  # Consistency model for observed outcomes
  lever_full_o <- leverage_plot(full, drug_names,
                                title = "Network meta-analysis")

  # UME model for observed outcomes
  lever_ume_o <- leverage_plot(ume, drug_names,
                               title = "Unrelated mean effects model")


  ## Bring together the leverage plots for observed outcome
  lev_plots <- ggpubr::ggarrange(lever_full_o,
                                 lever_ume_o,
                                 nrow = 1, ncol = 2, labels = c("A)",  "B)"))


  intervalplots <- intervalplot_panel_ume(full, ume, drug_names)

  # Write the table with the EMs from both models as .xlsx
  if (save_xls == TRUE) {
    write_xlsx(em_both, paste0("Table NMA vs UME", ".xlsx"))
    write_xlsx(model_assessment,
               paste0("Table Model Assessment_NMA vs UME", ".xlsx"))
  }

  # Return results
  results <- if (model == "RE") {
    list(table_effect_size =
           knitr::kable(em_both,
                        align = "lcccccc",
                        caption = "Estimation for each observed comparison"),
         table_model_assessment =
           knitr::kable(model_assessment,
                        align = "lcc",
                        caption = paste0("Model assessment parameters (",
                                         data_points, " ",
                                         "unconstrained data points)")),
         table_tau =
           knitr::kable(between_trial_sd,
                        align = "lccc",
                        caption = "Between-trial standard deviation"),
         scatterplots = scatterplots,
         leverage_plots = lev_plots,
         intervalplots = intervalplots)
  } else {
    list(table_effect_size =
           knitr::kable(em_both,
                        align = "lcccccc",
                        caption = "Estimation for each observed comparison"),
         table_model_assessment =
           knitr::kable(model_assessment,
                        align = "lcc",
                        caption = paste0("Model assessment parameters (",
                                         data_points, " ",
                                         "unconstrained data points)")),
         scatterplots = scatterplots,
         leverage_plots = lev_plots,
         intervalplots = intervalplots)
  }

  return(results)
}

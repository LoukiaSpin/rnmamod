#' End-user-ready results for series of pairwise meta-analyses
#'
#' @description Facilitates the comparison of the estimated effects and
#'   between-trial standard deviation (where applicable) from the network
#'   meta-analysis versus the results from the pairwise meta-analyses.
#'
#' @param full An object of S3 class \code{\link{run_model}}. See 'Value' in
#'   \code{\link{run_model}}.
#' @param meta An object of S3 class \code{\link{run_series_meta}}. See 'Value'
#'   in \code{\link{run_series_meta}}.
#' @param drug_names A vector of labels with the name of the interventions in
#'   the order they appear in the argument \code{data} of
#'   \code{\link{run_model}}. If the argument \code{drug_names} is not defined,
#'   the interventions are ordered as they appear in \code{data}.
#' @param save_xls Logical to indicate whether to export the tabulated results
#' in Excel 'xlsx' format at the working directory of the user.
#'   The default is \code{FALSE} (do not export in Excel format).
#'
#' @return \code{series_meta_plot} returns a panel of two forest plots: (1) a
#'   forest plot on the posterior mean and 95\% credible interval of the effect
#'   of the observed comparisons under the consistency model and the
#'   corresponding pairwise meta-analyses, and (2) a forest plot on the
#'   posterior median and 95\% credible interval of the between-trial standard
#'   deviation (\eqn{\tau}) for the observed comparisons. The estimated
#'   \eqn{\tau} from the consistency model appears as a rectangle in the forest
#'   plot. When a fixed-effect model has been fitted, only the forest plot on
#'   the estimated effects is shown.
#'
#'   The R console prints the data-frame with the estimated effects of the
#'   observed comparisons under both models and the data-frame on the estimated
#'   (\eqn{\tau})s from the pairwise meta-analyses. Furthermore,
#'   \code{series_meta_plot} exports both data-frames in Excel 'xlsx' format at
#'   the working directory of the user.
#'
#' @details \code{series_meta_plot} can be used only for a network of
#'   interventions. The user can detect any inconsistencies in the estimated
#'   effects from the compared models and explore the gains in their precision
#'   by applying the consistency model. Furthermore, the user can investigate
#'   the plausibility of the common between-trial heterogeneity assumption which
#'   is typically considered in the consistency model.
#'
#' @author {Loukia M. Spineli}
#'
#' @seealso \code{\link{run_model}}, \code{\link{run_series_meta}}
#'
#' @examples
#' data("nma.dogliotti2014")
#'
#' \dontrun{
#' # Perform a random-effects network meta-analysis
#' res <- run_model(data = nma.dogliotti2014,
#'                  measure = "OR",
#'                  model = "RE",
#'                  assumption = "IDE-ARM",
#'                  heter_prior = list("halfnormal", 0, 1),
#'                  mean_misspar = c(0, 0),
#'                  var_misspar = 1,
#'                  D = 1,
#'                  n_chains = 3,
#'                  n_iter = 10000,
#'                  n_burnin = 1000,
#'                  n_thin = 1)
#'
#' # Run separate random-effects pairwise meta-analyses
#' meta <- run_series_meta(full = res,
#'                         n_chains = 3,
#'                         n_iter = 10000,
#'                         n_burnin = 1000,
#'                         n_thin = 1)
#'
#' # The names of the interventions in the order they appear in the dataset
#' interv_names <- c("placebo", "aspirin", "aspirin plus clopidogrel",
#'                   "dabigatran 110 mg", "dabigatran 150 mg", "rivaroxaban",
#'                   "vitamin K antagonist", "apixaban")
#'
#' # Plot the results from both models
#' series_meta_plot(full = res,
#'                  meta = meta,
#'                  drug_names = interv_names)
#' }
#' @export
series_meta_plot <- function(full, meta, drug_names, save_xls) {

  save_xls <- if (missing(save_xls)) {
    FALSE
  } else {
    save_xls
  }

  drug_names <- if (missing(drug_names)) {
    message(cat(paste0("\033[0;",
                       col = 32,
                       "m",
                       txt = "The argument 'drug_names' has not been defined.
                       The intervention ID, as specified in 'data' is used as
                       intervention names", "\033[0m", "\n")))
    nt <- length(full$SUCRA[, 1])
    as.character(1:nt)
  } else {
    drug_names
  }

  if (length(drug_names) < 3) {
    stop("This function is *not* relevant for a pairwise meta-analysis",
         call. = F)
  }

  # Posterior results on the effect estimates under NMA
  em_full0 <- full$EM

  # Posterior results on the effect estimates under separate MAs
  em_meta0 <- meta$EM

  # Posterior results on between-trial standard deviation under NMA
  tau_full <- full$tau

  # Posterior results on between-trial standard deviation under separate MAs
  tau_meta <- meta$tau

  # Possible and observed comparisons
  possible_comp <- possible_observed_comparisons(drug_names,
                                                 obs_comp =
                                                   paste0(meta$EM[, "t2"], "vs",
                                                          meta$EM[, "t1"]))

  # Observed comparisons
  obs_comp <- possible_comp$obs_comp[, 3]

  model <- full$model

  # Keep only the comparisons with at least two trials
  em_full <- em_full0[is.element(possible_comp$poss_comp[, 4], obs_comp),
                      c(1:3, 7)]
  em_full[, c(1, 3:4)] <- if (is.element(full$measure, c("OR", "ROM"))) {
    exp(em_full[, c(1, 3:4)])
  } else if (is.element(full$measure, c("MD", "SMD"))) {
    em_full[, c(1, 3:4)]
  }
  em_full_clean <- format(round(em_full, 2), nsmall = 2)

  # Effect estimate of separate MAs
  em_meta <- round(em_meta0[, c(3:5, 9)], 2)
  em_meta[, c(1, 3:4)] <- if (is.element(full$measure, c("OR", "ROM"))) {
    exp(em_meta[, c(1, 3:4)])
  } else if (is.element(full$measure, c("MD", "SMD"))) {
    em_meta[, c(1, 3:4)]
  }
  em_meta_clean <- format(round(em_meta, 2), nsmall = 2)

  # Between-trial standard deviation of separate MAs
  if (model == "RE") {
    tau_meta_clean <- format(round(tau_meta[, c(7, 4:5, 9)], 2), nsmall = 2)
  } else {
    tau_meta_clean <- NA
  }

  # Credible intervals for comparisons with at least two trials
  cri_full_clean <- if (is.element(full$measure, c("OR", "ROM"))) {
    paste0("(", em_full_clean[, 3], ",", " ", em_full_clean[, 4], ")",
           ifelse(as.numeric(em_full_clean[, 3]) > 1 |
                    as.numeric(em_full_clean[, 4]) < 1, "*", " "))
  } else {
    paste0("(", em_full_clean[, 3], ",", " ", em_full_clean[, 4], ")",
           ifelse(as.numeric(em_full_clean[, 3]) > 0 |
                    as.numeric(em_full_clean[, 4]) < 0, "*", " "))
  }

  # The 95% CrIs of the effect estimate of separate MAs
  cri_meta_clean <- if (is.element(full$measure, c("OR", "ROM"))) {
    paste0("(", em_meta_clean[, 3], ",", " ", em_meta_clean[, 4], ")",
           ifelse(as.numeric(em_meta_clean[, 3]) > 1 |
                    as.numeric(em_meta_clean[, 4]) < 1, "*", " "))
  } else {
    paste0("(", em_meta_clean[, 3], ",", " ", em_meta_clean[, 4], ")",
           ifelse(as.numeric(em_meta_clean[, 3]) > 0 |
                    as.numeric(em_meta_clean[, 4]) < 0, "*", " "))
  }

  # The 95% CrIs of the between-trial standard deviation of separate MAs
  cri_tau_meta <- if (model == "RE") {
    paste0("(", tau_meta_clean[, 3], ",", " ", tau_meta_clean[, 4], ")")
  } else {
    NA
  }

  ## Create a data-frame with effect estimates on both models
  if (model == "RE") {
    em_both <- data.frame(possible_comp$obs_comp[, 4],
                          em_full_clean[, 1:2],
                          cri_full_clean,
                          em_meta_clean[, 1:2],
                          cri_meta_clean,
                          tau_meta_clean[, 1:2],
                          cri_tau_meta)
    colnames(em_both) <- c("Comparison",
                           "Mean NMA", "SD NMA", "95% CrI NMA",
                           "Mean MA", "SD MA", "95% CrI MA",
                           "Median tau", "SD tau", "95% CrI tau")
  } else {
    em_both <- data.frame(possible_comp$obs_comp[, 4],
                          em_full_clean[, 1:2],
                          cri_full_clean,
                          em_meta_clean[, 1:2],
                          cri_meta_clean)
    colnames(em_both) <- c("Comparison",
                           "Mean NMA", "SD NMA", "95% CrI NMA",
                           "Mean MA", "SD MA", "95% CrI MA")
  }
  rownames(em_both) <- NULL

  # Prepare dataset for ggplot2
  # Effect estimate
  prepare <- data.frame(rep(seq_len(length(obs_comp)), 2),
                        rep(possible_comp$obs_comp[, 4], 2),
                        rbind(apply(em_full_clean[, -2], 2, as.numeric),
                              apply(em_meta_clean[, -2], 2, as.numeric)),
                        rep(c("Network meta-analysis",
                              "Pairwise meta-analysis"),
                            each = length(obs_comp)))
  colnames(prepare) <- c("order",
                         "comparison", "mean", "lower", "upper",
                         "analysis")
  rownames(prepare) <- NULL

  # Between-trial standard deviation
  if (model == "RE") {
    prepare_tau <- data.frame(possible_comp$obs_comp[, 4],
                              tau_meta_clean[, -2])
    colnames(prepare_tau) <- c("comparison",
                               "median", "lower", "upper")
  } else {
    prepare_tau <- NA
  }

  # Forest plots of comparisons on effect estimate
  p1 <- ggplot(data = prepare,
               aes(x = as.factor(order),
                   y = mean,
                   ymin = lower,
                   ymax = upper,
                   colour = analysis,
                   group = analysis)) +
          geom_linerange(size = 2,
                         position = position_dodge(width = 0.5)) +
          geom_hline(yintercept =
                       ifelse(!is.element(full$measure, c("OR", "ROM")), 0, 1),
                     lty = 1,
                     size = 1,
                     col = "grey53") +
          geom_point(size = 1.5,
                     colour = "black",
                     stroke = 0.3,
                     position = position_dodge(width = 0.5)) +
          geom_text(aes(x = as.factor(order),
                        y = mean,
                        label = paste0(sprintf("%.2f", mean), " ", "(",
                                       sprintf("%.2f", lower), ",", " ",
                                       sprintf("%.2f", upper), ")"),
                        hjust = 0,
                        vjust = -0.5),
                    color = "black",
                    size = 4.0,
                    position = position_dodge(width = 0.5)) +
          geom_text(aes(x = 0.45,
                        y = ifelse(is.element(full$measure, c("OR", "ROM")),
                                   0.4, -0.2),
                        label = ifelse(full$D == 0, "Favours first arm",
                                       "Favours second arm")),
                    size = 3.5,
                    vjust = 0,
                    hjust = 0,
                    color = "black") +
          geom_text(aes(x = 0.45,
                        y = ifelse(is.element(full$measure, c("OR", "ROM")),
                                   1.2, 0.2),
                        label = ifelse(full$D == 0, "Favours second arm",
                                       "Favours first arm")),
                    size = 3.5,
                    vjust = 0,
                    hjust = 0,
                    color = "black") +
          labs(x = "",
               y = effect_measure_name(full$measure),
               colour = "Analysis") +
          scale_x_discrete(breaks = as.factor(seq_len(length(obs_comp))),
                           labels = prepare$comparison[
                             seq_len(length(obs_comp))]) +
          scale_y_continuous(trans =
                               ifelse(!is.element(full$measure, c("OR", "ROM")),
                                      "identity", "log10")) +
          scale_color_manual(breaks = c("Network meta-analysis",
                                        "Pairwise meta-analysis"),
                             values = c("#009E73", "#D55E00")) +
          coord_flip() +
          theme_classic() +
          theme(axis.text.x = element_text(color = "black", size = 12),
                axis.text.y = element_text(color = "black", size = 12),
                axis.title.x = element_text(color = "black", face = "bold",
                                            size = 12),
                legend.position = "none", legend.justification = c(0.13, 0),
                legend.text =  element_text(color = "black", size = 12),
                legend.title =  element_text(color = "black", face = "bold",
                                             size = 12))

  # Forest plots of comparisons-specific between-trial standard deviation
  p2 <- if (model == "RE") {
    ggplot(data = prepare_tau,
           aes(x = as.factor(seq_len(length(obs_comp))),
               y = as.numeric(median),
               ymin = as.numeric(lower),
               ymax = as.numeric(upper))) +
      geom_rect(aes(xmin = 0,
                    xmax = Inf,
                    ymin = tau_full[3],
                    ymax = tau_full[7]),
                fill = "#1f93ff",
                alpha = 0.1) +
      geom_linerange(size = 2,
                     position = position_dodge(width = 0.5)) +
      geom_hline(yintercept = tau_full[5],
                 lty = 1,
                 size = 1,
                 col = "#006CD1") +
      geom_point(size = 1.5,
                 colour = "white",
                 stroke = 0.3,
                 position = position_dodge(width = 0.5)) +
      geom_text(aes(x = as.factor(seq_len(length(obs_comp))),
                    y = round(as.numeric(median), 2),
                    label = paste0(round(as.numeric(median), 2), " ", "(",
                                   round(as.numeric(lower), 2), ",", " ",
                                   round(as.numeric(upper), 2), ")")),
                color = "black",
                hjust = 0,
                vjust = -0.5,
                size = 4.0,
                check_overlap = F,
                parse = F,
                position = position_dodge(width = 0.8),
                inherit.aes = T) +
      scale_x_discrete(breaks = as.factor(seq_len(length(obs_comp))),
                       labels = prepare_tau$comparison[
                         seq_len(length(obs_comp))]) +
      labs(x = "", y = "Between-trial standard deviation") +
      coord_flip() +
      theme_classic() +
      theme(axis.text.x = element_text(color = "black", size = 12),
            axis.text.y = element_text(color = "black", size = 12),
            axis.title.x = element_text(color = "black", face = "bold",
                                        size = 12))
  } else {
    NA
  }

  # Bring together both forest-plots
  forest_plots <- if (model == "RE") {
    ggarrange(p1, p2, nrow = 1, ncol = 2, labels = c("A)", "B)"),
              common.legend = T, legend = "bottom")
  } else {
    p1
  }

  # Write the table as .xlsx
  if (save_xls == TRUE) {
    write_xlsx(em_both, paste0("Table.NMA.vs.PMA", ".xlsx"))
  }

  return(list(tabulated_results = knitr::kable(em_both),
              forest_plots = forest_plots))
}

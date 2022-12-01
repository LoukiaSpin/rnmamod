#' End-user-ready results for a series of pairwise meta-analyses
#'
#' @description Facilitates the comparison of the consistency model
#'   (via \code{\link{run_model}}) with a series of pairwise meta-analyses
#'   (via \code{\link{run_series_meta}}) regarding the estimated summary effect
#'   sizes and between-trial standard deviation for comparisons with at
#'   least two trials.
#'
#' @param full An object of S3 class \code{\link{run_model}}. See 'Value' in
#'   \code{\link{run_model}}.
#' @param meta An object of S3 class \code{\link{run_series_meta}}. See 'Value'
#'   in \code{\link{run_series_meta}}.
#' @param drug_names A vector of labels with the name of the interventions in
#'   the order they appear in the argument \code{data} of
#'   \code{\link{run_model}}. If \code{drug_names} is not defined,
#'   the order of the interventions as they appear in \code{data} is used,
#'   instead.
#' @param save_xls Logical to indicate whether to export the tabulated results
#'   to an 'xlsx' file (via the \code{\link[writexl:write_xlsx]{write_xlsx}}
#'   function of the R-package
#'   \href{https://CRAN.R-project.org/package=writexl}{writexl}) at the working
#'   directory of the user. The default is \code{FALSE} (do not export).
#'
#' @return The R console prints the data-frame with the estimated summary effect
#'   sizes and between-trial standard deviation of comparisons under both
#'   models. The comparisons have at least two trials. In the case of a
#'   fixed-effect model, the data-frame is printed without the results on the
#'   between-trial standard deviation.
#'
#'   Furthermore, \code{series_meta_plot} exports the data-frame to an 'xlsx'
#'   file at the working  directory of the user.
#'
#'   \code{series_meta_plot} returns a panel of two forest plots: (1) a
#'   forest plot on the posterior mean and 95\% credible interval of the summary
#'   effect size for the observed comparisons from network meta-analysis and the
#'   corresponding pairwise meta-analyses, and (2) a forest plot on the
#'   posterior median and 95\% credible interval of the between-trial standard
#'   deviation for these observed comparisons. The estimated median and 95\%
#'   credible intervals of the between-trial standard deviation from network
#'   meta-analysis appear in the forest plot as a solid and two dotted parallel
#'   blue lines, respectively. The different levels of heterogeneity appear as
#'   green, yellow, orange, and red rectangles to indicate a low, reasonable,
#'   fairly high, and fairly extreme heterogeneity, respectively, following the
#'   classification of Spiegelhalter et al. (2004).
#'   When a fixed-effect model has been fitted, only the forest plot on the
#'   estimated summary effect sizes is shown.
#'
#' @details \code{series_meta_plot} can be used only for a network of
#'   interventions. Otherwise, the execution of the function will be stopped and
#'   an error message will be printed on the R console.
#'
#'   For a binary outcome, when \code{measure} is "RR" (relative risk) or "RD"
#'   (risk difference) in \code{\link{run_model}}, \code{series_meta_plot}
#'   currently presents the results in the odds ratio for being the
#'   \strong{base-case} effect measure in \code{\link{run_model}} for a binary
#'   outcome (see also 'Details' in \code{\link{run_model}}).
#'
#'   The user can detect any inconsistencies in the estimated
#'   effects from the compared models and explore the gains in precision
#'   stemming from applying network meta-analysis. Furthermore, the user can
#'   investigate the plausibility of the common between-trial heterogeneity
#'   assumption which is typically considered in network meta-analysis.
#'
#' @author {Loukia M. Spineli}
#'
#' @seealso \code{\link{run_model}}, \code{\link{run_series_meta}},
#'   \code{\link[writexl:write_xlsx]{write_xlsx}}
#'
#' @references
#' Spiegelhalter DJ, Abrams KR, Myles JP. Bayesian approaches to clinical trials
#' and health-care evaluation. John Wiley and Sons, Chichester, 2004.
#'
#' @examples
#' data("nma.dogliotti2014")
#'
#' # Read results from 'run_model' (using the default arguments)
#' res <- readRDS(system.file('extdata/res_dogliotti.rds', package = 'rnmamod'))
#'
#' # Read results from 'run_series_meta' (using the default arguments)
#' meta <- readRDS(system.file('extdata/meta_dogliotti.rds',
#'                 package = 'rnmamod'))
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
#'
#' @export
series_meta_plot <- function(full, meta, drug_names, save_xls) {

  if (class(full) != "run_model" || is.null(full)) {
    stop("'full' must be an object of S3 class 'run_model'.",
         call. = FALSE)
  }

  if (class(meta) != "run_series_meta" || is.null(meta)) {
    stop("'meta' must be an object of S3 class 'run_series_meta'.",
         call. = FALSE)
  }

  save_xls <- if (missing(save_xls)) {
    FALSE
  } else {
    save_xls
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

  measure <- if (is.element(full$measure, c("RR", "RD"))) {
    "OR"
  } else {
    full$measure
  }

  # Posterior results on the effect estimates under NMA
  em_full0 <- if (is.element(full$measure, c("RR", "RD"))) {
    full$EM_LOR
  } else {
    full$EM
  }

  # Posterior results on the effect estimates under separate MAs
  # Keep only comparisons with at least two trials
  em_meta00 <- meta$EM[meta$single == 0, ]
  rownames(em_meta00) <- seq_len(length(em_meta00[, 1]))

  # Posterior results on between-trial standard deviation under NMA
  tau_full <- full$tau

  # Posterior results on between-trial standard deviation under separate MAs
  tau_meta0 <- meta$tau

  # Possible and observed comparisons
  possible_comp <-
    possible_observed_comparisons(drug_names,
                                  obs_comp =
                                    paste0(
                                      meta$EM[meta$single == 0, "t2"], "vs",
                                      meta$EM[meta$single == 0, "t1"])
                                  )

  # Observed comparisons
  obs_comp <- possible_comp$obs_comp[, 3]

  model <- full$model

  # Keep the same comparisons with those in PMA (they have at least two trials)
  em_full <- em_full0[is.element(possible_comp$poss_comp[, 4], obs_comp),
                      c(1:3, 7)]
  em_full[, c(1, 3:4)] <- if (is.element(measure, c("OR", "ROM"))) {
    exp(em_full[, c(1, 3:4)])
  } else if (is.element(measure, c("MD", "SMD"))) {
    em_full[, c(1, 3:4)]
  }
  em_full_clean <- format(round(em_full, 2), nsmall = 2)

  # Sort comparisons in the order they appear in NMA results
  # Effect estimates
  em_meta01 <- em_meta00[order(em_meta00$t2),] # First, sort by t2
  em_meta0 <- em_meta01[order(em_meta01$t1),]  # Then, sort by t1
  # Between-trial standard deviation tau_meta00
  tau_meta1 <- tau_meta0[order(tau_meta0$t2),] # First, sort by t2
  tau_meta <- tau_meta1[order(tau_meta1$t1),]  # Then, sort by t1

  # Effect estimate of separate MAs
  em_meta <- round(em_meta0[, c(3:5, 9)], 2)
  em_meta[, c(1, 3:4)] <- if (is.element(measure, c("OR", "ROM"))) {
    exp(em_meta[, c(1, 3:4)])
  } else if (is.element(measure, c("MD", "SMD"))) {
    em_meta[, c(1, 3:4)]
  }
  em_meta_clean <- format(round(em_meta, 2), nsmall = 2)

  # Between-trial standard deviation of separate MAs  5, 2:3, 7
  if (model == "RE") {
    tau_meta_clean <- format(round(tau_meta[, c(7, 4:5, 9)], 2), nsmall = 2)
  } else {
    tau_meta_clean <- NA
  }

  # Credible intervals for comparisons with at least two trials
  cri_full_clean <- if (is.element(measure, c("OR", "ROM"))) {
    paste0("(", em_full_clean[, 3], ",", " ", em_full_clean[, 4], ")",
           ifelse(as.numeric(em_full_clean[, 3]) > 1 |
                    as.numeric(em_full_clean[, 4]) < 1, "*", " "))
  } else {
    paste0("(", em_full_clean[, 3], ",", " ", em_full_clean[, 4], ")",
           ifelse(as.numeric(em_full_clean[, 3]) > 0 |
                    as.numeric(em_full_clean[, 4]) < 0, "*", " "))
  }

  # The 95% CrIs of the effect estimate of separate MAs
  cri_meta_clean <- if (is.element(measure, c("OR", "ROM"))) {
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
  add <- ifelse(is.element(measure, c("OR", "ROM")), 1, 4)
  measure2 <- effect_measure_name(measure, lower = FALSE)

  caption <- if (full$D == 0 & is.element(measure, c("OR", "ROM"))) {
    paste(measure2, "< 1, favours the first arm.",
          measure2, "> 1, favours the second arm.")
  } else if (full$D == 1 & is.element(measure, c("OR", "ROM"))) {
    paste(measure2, "< 1, favours the second arm.",
          measure2, "> 1, favours the first arm.")
  } else if (full$D == 0 & !is.element(measure, c("OR", "ROM"))) {
    paste(measure2, "< 0, favours the first arm.",
          measure2, "> 0, favours the second arm.")
  } else if (full$D == 1 & !is.element(measure, c("OR", "ROM"))) {
    paste(measure2, "< 0, favours the second arm.",
          measure2, "> 0, favours the first arm.")
  }

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
                       ifelse(!is.element(measure, c("OR", "ROM")), 0, 1),
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
          scale_x_discrete(breaks = as.factor(seq_len(length(obs_comp))),
                           labels = prepare$comparison[
                             seq_len(length(obs_comp))]) +
          scale_y_continuous(trans =
                               ifelse(!is.element(measure, c("OR", "ROM")),
                                      "identity", "log10")) +
          scale_color_manual(breaks = c("Network meta-analysis",
                                        "Pairwise meta-analysis"),
                             values = c("#009E73", "#D55E00")) +
          labs(x = "",
               y = measure2,
               colour = "Analysis",
               caption = caption) +
          coord_flip() +
          theme_classic() +
          theme(axis.text.x = element_text(color = "black", size = 12),
                axis.text.y = element_text(color = "black", size = 12),
                axis.title.x = element_text(color = "black", face = "bold",
                                            size = 12),
                legend.position = "bottom", legend.justification = c(0.13, 0),
                legend.text =  element_text(color = "black", size = 12),
                legend.title =  element_text(color = "black", face = "bold",
                                             size = 12),
                plot.caption = element_text(hjust = 0.01))

  # Forest plots of comparisons-specific between-trial standard deviation
  p2 <- if (model == "RE") {
   ggplot(data = prepare_tau,
           aes(x = as.factor(seq_len(length(obs_comp))),
               y = as.numeric(median),
               ymin = as.numeric(lower),
               ymax = as.numeric(upper))) +
      geom_rect(aes(xmin = 0, xmax = Inf, ymin = 0, ymax = 0.099,
                    fill = "low"),
                alpha = 0.02) +
      geom_rect(aes(xmin = 0, xmax = Inf, ymin = 0.1, ymax = 0.5,
                    fill = "reasonable"),
                alpha = 0.02) +
      geom_rect(aes(xmin = 0, xmax = Inf, ymin = 0.5, ymax = 1.0,
                    fill = "fairly high"),
                alpha = 0.02) +
      geom_rect(aes(xmin = 0, xmax = Inf, ymin = 1.0, ymax = Inf,
                    fill = "fairly extreme"),
                alpha = 0.02) +
      geom_hline(yintercept = tau_full[5],
                 lty = 1,
                 size = 1,
                 col = "#006CD1") +
      geom_hline(yintercept = tau_full[3],
                 lty = 3,
                 size = 1,
                 col = "#006CD1") +
      geom_hline(yintercept = tau_full[7],
                 lty = 3,
                 size = 1,
                 col = "#006CD1") +
      geom_linerange(size = 2,
                     position = position_dodge(width = 0.5)) +
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
                check_overlap = FALSE,
                parse = FALSE,
                position = position_dodge(width = 0.8),
                inherit.aes = TRUE) +
      scale_x_discrete(breaks = as.factor(seq_len(length(obs_comp))),
                       labels = prepare_tau$comparison[
                         seq_len(length(obs_comp))]) +
      scale_fill_manual(name = "Heterogeneity",
                        values = c("low" = "#009E73",
                                   "reasonable" = "orange",
                                   "fairly high" = "#D55E00",
                                   "fairly extreme" = "red"),
                        limits = c("low", "reasonable", "fairly high",
                                   "fairly extreme")) +
      labs(x = "", y = "Between-trial standard deviation", caption = " ") +
      coord_flip() +
      theme_classic() +
      theme(axis.text.x = element_text(color = "black", size = 12),
            axis.text.y = element_text(color = "black", size = 12),
            axis.title.x = element_text(color = "black", face = "bold",
                                        size = 12),
            legend.position = "bottom",
            legend.text =  element_text(color = "black", size = 12),
            legend.title =  element_text(color = "black", face = "bold",
                                         size = 12))
  }

  # Bring together both forest-plots
  forest_plots <- if (model == "RE") {
    ggarrange(p1, p2 + guides(fill = guide_legend(override.aes =
                                                    list(alpha = 0.4))),
              nrow = 1, ncol = 2, labels = c("A)", "B)"),
              common.legend = FALSE, legend = "bottom")
  } else {
    p1
  }

  # Write the table as .xlsx
  if (save_xls == TRUE) {
    write_xlsx(em_both, paste0("Table.NMA.vs.PMA", ".xlsx"))
  }

  return(list(tabulated_results =
                knitr::kable(em_both,
                             caption =
                               "Estimation for each observed comparison"),
              forest_plots = forest_plots))
}

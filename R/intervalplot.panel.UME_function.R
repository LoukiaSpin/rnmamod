#' A panel of interval plots for the unrelated mean effects model
#'
#' @description Creates a panel of interval plots on the summary effect sizes
#'   under the consistency model and the unrelated mean effects model.
#'   The number of interval plots equals the number of pairwise comparisons
#'   observed in the network.
#'
#' @param full An object of S3 class \code{\link{run_model}}. See 'Value' in
#'   \code{\link{run_model}}.
#' @param ume An object of S3 class \code{\link{run_ume}}. See 'Value' in
#'   \code{\link{run_ume}}.
#' @param drug_names A vector of labels with the name of the interventions in
#'   the order they appear in the argument \code{data} of
#'   \code{\link{run_model}}. If the argument \code{drug_names} is not defined,
#'   the order of the interventions as they appear in \code{data} is used,
#'   instead.
#'
#' @return A panel of interval plots on the posterior mean and 95\% credible
#'   interval of the summary effect size under the consistency model and the
#'   improved unrelated mean effects model (Spineli, 2021) of all pairwise
#'   comparisons observed in the network.
#'
#' @details \code{intervalplot_panel_ume} is integrated in the
#'   \code{\link{ume_plot}} function. The consistency model and the unrelated
#'   mean effects model are abbreviated in the y-axis as 'NMA model' and
#'   'UME model', respectively. The intervals are highlighted with green, when
#'   the corresponding summary effect sizes do not cross the vertical line of no
#'   difference, and red otherwise. Grey panels refer to the frail comparisons
#'   as detected by the \code{\link{improved_ume}} function (see 'Details' in
#'   \code{\link{improved_ume}}).
#'
#'   For a binary outcome, when \code{measure} is "RR" (relative risk) or "RD"
#'   (risk difference) in \code{\link{run_model}}, \code{intervalplot_panel_ume}
#'   currently presents the results in the odds ratio scale.
#'
#' @author {Loukia M. Spineli}
#'
#' @seealso \code{\link{improved_ume}} \code{\link{run_model}},
#'   \code{\link{run_ume}}, \code{\link{ume_plot}}
#'
#' @references
#' Spineli LM. A revised framework to evaluate the consistency assumption
#' globally in a network of interventions. \emph{Med Decis Making} 2021.
#' doi: 10.1177/0272989X211068005
#'
#' @export
intervalplot_panel_ume <- function(full, ume, drug_names) {

  if (full$type != "nma" || is.null(full$type)) {
    stop("'full' must be an object of S3 class 'run_model'.",
         call. = FALSE)
  }

  if (ume$type != "ume" || is.null(ume$type)) {
    stop("'ume' must be an object of S3 class 'run_ume'.",
         call. = FALSE)
  }

  measure <- if (is.element(full$measure, c("RR", "RD"))) {
    "OR"
  } else {
    full$measure
  }
  em_full <- if (is.element(full$measure, c("RR", "RD"))) {
    full$EM_LOR
  } else {
    full$EM
  }
  em_ume <- ume$EM
  obs_comp <- ume$obs_comp
  frail_comp <- ume$frail_comp

  # Possible and observed comparisons
  possible_comp <- possible_observed_comparisons(drug_names, obs_comp)

  # Indicate frail comparisons
  frail_comp_ind <- ifelse(is.element(possible_comp$obs_comp[, 3], frail_comp),
                           "yes", "no")

  # Calculate bounds of 95% CrI
  ume_lower <- em_ume[, 3]
  ume_upper <- em_ume[, 7]

  # Keep only the effect estimates according to the 'poss.pair.comp.clean'
  ume_mean <- round(em_ume[, 1], 2)
  nma_mean <- round(
    em_full[is.element(possible_comp$poss_comp[, 4], obs_comp), 1], 2)
  nma_lower <- round(
    em_full[is.element(possible_comp$poss_comp[, 4], obs_comp), 3], 2)
  nma_upper <- round(
    em_full[is.element(possible_comp$poss_comp[, 4], obs_comp), 7], 2)

  # Indicate statistical significance
  nma_stat_signif <- ifelse(nma_lower > 0 | nma_upper < 0, "strong", "weak")
  ume_stat_signif <- ifelse(ume_lower > 0 | ume_upper < 0, "strong", "weak")

  # Create the data-frame
  data_set <- if (is.element(measure, c("OR", "RR", "ROM"))) {
    data.frame(round(exp(c(nma_mean, ume_mean)), 2),
               round(exp(c(nma_lower, ume_lower)), 2),
               round(exp(c(nma_upper, ume_upper)), 2),
               c(nma_stat_signif, ume_stat_signif),
               rep(possible_comp$obs_comp[, 4], 2),
               rep(c("NMA", "UME"), each = length(obs_comp)),
               frail_comp_ind)
  } else {
    data.frame(c(nma_mean, ume_mean),
               c(nma_lower, ume_lower),
               c(nma_upper, ume_upper),
               c(nma_stat_signif, ume_stat_signif),
               rep(possible_comp$obs_comp[, 4], 2),
               rep(c("NMA", "UME"), each = length(obs_comp)),
               frail_comp_ind)
  }
  colnames(data_set) <- c("mean",
                          "lower",
                          "upper",
                          "stat_sign",
                          "comp",
                          "analysis",
                          "frail")

  # Obtain forestplot
  measure2 <- effect_measure_name(measure, lower = FALSE)
  add <- ifelse(is.element(measure, c("OR", "ROM")), 1, 4)
  caption <- if (full$D == 0 & is.element(measure,
                                          c("OR", "ROM"))) {
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

  ggplot(data = data_set,
         aes(x = as.factor(analysis),
             y = mean,
             ymin = lower,
             ymax = upper,
             colour = stat_sign)) +
    geom_rect(aes(fill = frail),
              xmin = -Inf,
              xmax = Inf,
              ymin = -Inf,
              ymax = Inf,
              alpha = 0.2) +
    geom_linerange(size = 2,
                   position = position_dodge(width = 0.5)) +
    geom_hline(yintercept = ifelse(!is.element(
      measure, c("OR", "RR", "ROM")), 0, 1),
               lty = 1,
               size = 1,
               col = "grey53") +
    geom_point(size = 1.5,
               colour = "white",
               stroke = 0.3,
               position = position_dodge(width = 0.5)) +
    geom_text(aes(x = as.factor(analysis),
                  y = mean,
                  label = paste0(sprintf("%.2f", mean),
                                 " ",
                                 "(",
                                 sprintf("%.2f", lower),
                                 ",",
                                 " ",
                                 sprintf("%.2f", upper), ")")),
              color = "black",
              hjust = 0,
              vjust = -0.5,
              size = 3.3,
              check_overlap = FALSE,
              parse = FALSE,
              position = position_dodge(width = 0.8),
              inherit.aes = TRUE) +
    facet_wrap(vars(factor(comp, levels = unique(data_set$comp))),
               scales = "fixed") +
    scale_fill_manual(breaks = c("yes", "no"),
                      values = c("grey53", "white")) +
    scale_color_manual(breaks = c("strong", "weak"),
                       values = c("#009E73", "#D55E00")) +
    scale_y_continuous(trans = ifelse(
      !is.element(measure, c("OR", "RR", "ROM")), "identity", "log10")) +
    labs(x = "",
         y = measure2,
         colour = "Evidence",
         fill = "",
         caption = caption) +
    coord_flip() +
    theme_classic() +
    guides(fill = "none") +
    theme(axis.text.x = element_text(color = "black", size = 11),
          axis.text.y = element_text(color = "black", size = 11),
          axis.title.x = element_text(color = "black", face = "bold",
                                      size = 11),
          strip.text = element_text(size = 11),
          legend.position = "bottom",
          legend.text = element_text(color = "black", size = 11),
          legend.title = element_text(color = "black", face = "bold",
                                      size = 11),
          plot.caption = element_text(hjust = 0.01))
}

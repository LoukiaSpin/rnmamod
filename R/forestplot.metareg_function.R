#' Comparator-specific forest plot for network meta-regression
#'
#' @description This function illustrates a forest plot of the posterior mean
#'   and 95\% credible and prediction interval of comparisons with the selected
#'   intervention of the network under the network meta-analysis \emph{and}
#'   network meta-regression.
#'
#' @param full An object of S3 class \code{\link{run_model}}. See 'Value' in
#'   \code{\link{run_model}}.
#' @param reg An object of S3 class \code{\link{run_metareg}}. See 'Value' in
#'   \code{\link{run_metareg}}.
#' @param compar A character to indicate the comparator intervention. It must
#'   be any name found in \code{drug_names}.
#' @param cov_value A list of two elements in the following order: a number
#'   for the covariate value of interest (see 'Arguments' in
#'   \code{\link{run_metareg}}), and a character to indicate the name of
#'   the covariate. See also 'Details'.
#' @param drug_names A vector of labels with the name of the interventions in
#'   the order they appear in the argument \code{data} of
#'   \code{\link{run_model}}. If \code{drug_names} is not defined,
#'   the order of the interventions as they appear in \code{data} is used,
#'   instead.
#'
#' @return A panel of two forest plots: (1) a forest plot on the estimated
#'   effect size of comparisons with the selected intervention of the network,
#'   and (2) a forest plot on the predicted effect size of comparisons with the
#'   selected intervention of the network. Both forest plots illustrate the
#'   results from network meta-analysis and network meta-regression using
#'   different colours for the corresponding lines.
#'
#' @details In both plots, the y-axis displays all interventions in the
#'   network; the selected intervention that comprises the \code{compar} is
#'   indicated in the plot with a homonymous label. The numerical results are
#'   displayed above each line.
#'   Odds ratio, relative risks, and ratio of means are reported in the original
#'   scale after exponentiation of the logarithmic scale.
#'
#'   In both plots, the interventions are sorted in the descending order of
#'   their SUCRA values based on the network meta-analysis.
#'
#'   \code{forestplot_metareg} is integrated in \code{\link{metareg_plot}}.
#'
#'   \code{forestplot} can be used only for a network of interventions. In the
#'   case of two interventions, the execution of the function will be stopped
#'   and an error message will be printed on the R console.
#'
#' @author {Loukia M. Spineli}
#'
#' @seealso \code{\link{metareg_plot}}, \code{\link{run_metareg}},
#'   \code{\link{run_model}}
#'
#' @references
#' Salanti G, Ades AE, Ioannidis JP. Graphical methods and numerical summaries
#' for presenting results from multiple-treatment meta-analysis: an overview and
#' tutorial. \emph{J Clin Epidemiol} 2011;\bold{64}(2):163--71.
#' \doi{10.1016/j.jclinepi.2010.03.016}
#'
#' @export
forestplot_metareg <- function(full, reg, compar, cov_value, drug_names) {

  if (full$type != "nma" || is.null(full$type)) {
    stop("The argument 'full' must be an object of S3 class 'run_model'.",
         call. = FALSE)
  }

  if (reg$type != "nmr" || is.null(reg$type)) {
    stop("The argument 'reg' must be an object of S3 class 'run_metareg'.",
         call. = FALSE)
  }

  if (length(unique(reg$covariate)) < 3 &
      !is.element(cov_value[[1]], reg$covariate)) {
    aa <- "The first element of the argument 'cov_value' is out of the value"
    bb <- "range of the analysed covariate."
    stop(paste(aa, bb), call. = FALSE)
  } else if (length(unique(reg$covariate)) > 2 &
             (cov_value[[1]] < min(reg$covariate) |
              cov_value[[1]] > max(reg$covariate))) {
    aa <- "The first element of the argument 'cov_value' is out of the value"
    bb <- "range of the analysed covariate."
    stop(paste(aa, bb), call. = FALSE)
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
  len_drug <- length(drug_names)

  if (length(drug_names) < 3) {
    stop("This function is *not* relevant for a pairwise meta-analysis.",
         call. = FALSE)
  }

  # Sort the drugs by their SUCRA in decreasing order and remove the reference
  # intervention (number 1)
  drug_names_sorted <- drug_names[order(full$SUCRA[, 1], decreasing = TRUE)]

  compar <- if (missing(compar)) {
    stop("The argument 'compar' has not been defined.", call. = FALSE)
  } else if (!is.element(compar, drug_names)) {
    stop("The value of the argument 'compar' is not found in the 'drug_names'.",
         call. = FALSE)
  } else if (is.element(compar, drug_names)) {
    compar
  }

  cov_value <- if (!is.null(reg$beta_all) & missing(cov_value)) {
    stop("The argument 'cov_value' has not been defined.", call. = FALSE)
  } else if (!is.null(reg$beta_all) & length(cov_value) < 2) {
    aa <- "The argument 'cov_value' must be a list with elements a number and"
    stop(paste(aa, "a character."), call. = FALSE)
  } else if (!is.null(reg$beta_all) & length(cov_value) == 2) {
    cov_value
  }

  model <- full$model
  measure <- full$measure
  em_full <- full$EM
  pred_full <- full$EM_pred
  sucra_full <- full$SUCRA

  cov_val <- ifelse(length(unique(reg$covariate)) < 3,
                    cov_value[[1]],
                    cov_value[[1]] - mean(reg$covariate))

  # A matrix with all possible comparisons in the network
  poss_pair_comp1 <- data.frame(exp = t(combn(drug_names, 2))[, 2],
                                comp = t(combn(drug_names, 2))[, 1])
  poss_pair_comp2 <- data.frame(exp = t(combn(drug_names, 2))[, 1],
                                comp = t(combn(drug_names, 2))[, 2])
  poss_pair_comp <- rbind(poss_pair_comp1, poss_pair_comp2)

  # Effect size of all possible pairwise comparisons (NMA)
  em_ref00_nma <- cbind(rbind(data.frame(mean = em_full[, 1],
                                         lower = em_full[, 3],
                                         upper = em_full[, 7]),
                              data.frame(mean = em_full[, 1] * (-1),
                                         lower = em_full[, 7] * (-1),
                                         upper = em_full[, 3] * (-1))),
                        poss_pair_comp)
  em_subset_nma <- subset(em_ref00_nma, em_ref00_nma[5] == compar)
  em_ref0_nma <- rbind(em_subset_nma[, 1:3], c(rep(NA, 3)))
  sucra_new <- data.frame(sucra_full[, 1],
                          drug_names)[order(match(data.frame(sucra_full[, 1],
                                                             drug_names)[, 2],
                                                  em_subset_nma[, 4])), 1]
  em_ref_nma <- em_ref0_nma[order(sucra_new, decreasing = TRUE), ]

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
  em_ref_nmr <- em_ref0_nmr[order(sucra_new, decreasing = TRUE), ]
  rownames(em_ref_nma) <- rownames(em_ref_nmr) <- NULL

  # Posterior results on the predicted estimates of comparisons with the
  # selected comparator as reference
  if (model == "RE") {
    pred_ref00_nma <- cbind(rbind(data.frame(mean = pred_full[, 1],
                                             lower = pred_full[, 3],
                                             upper = pred_full[, 7]),
                                  data.frame(mean = pred_full[, 1] * (-1),
                                             lower = pred_full[, 7] * (-1),
                                             upper = pred_full[, 3] * (-1))),
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
    pred_ref_nma <- pred_ref0_nma[order(sucra_new, decreasing = TRUE), ]
    pred_ref_nmr <- pred_ref0_nmr[order(sucra_new, decreasing = TRUE), ]
    rownames(pred_ref_nma) <- rownames(pred_ref_nmr) <- NULL
  }

  # Create a data-frame with credible and prediction intervals of comparisons
  # with the reference intervention
  if (!is.element(measure, c("OR", "RR", "ROM")) & model == "RE") {
    prepare_em_nma <- data.frame(as.factor(rep(rev(seq_len(len_drug)), 2)),
                                 rep(drug_names_sorted, 2),
                                 round(rbind(em_ref_nma, pred_ref_nma), 2),
                                 rep(c("Estimation", "Prediction"),
                                     each = length(drug_names)))
    prepare_em_nmr <- data.frame(as.factor(rep(rev(seq_len(len_drug)), 2)),
                                 rep(drug_names_sorted, 2),
                                 round(rbind(em_ref_nmr, pred_ref_nmr), 2),
                                 rep(c("Estimation", "Prediction"),
                                     each = length(drug_names)))
    colnames(prepare_em_nma) <- c("order",
                                  "comparison",
                                  "mean", "lower", "upper",
                                  "interval")
    colnames(prepare_em_nmr) <- colnames(prepare_em_nma)
  } else if (is.element(measure, c("OR", "RR", "ROM")) & model == "RE") {
    prepare_em_nma <- data.frame(as.factor(rep(rev(seq_len(len_drug)), 2)),
                                 rep(drug_names_sorted, 2),
                                 round(rbind(exp(em_ref_nma),
                                             exp(pred_ref_nma)), 2),
                                 rep(c("Estimation", "Prediction"),
                                     each = length(drug_names)))
    prepare_em_nmr <- data.frame(as.factor(rep(rev(seq_len(len_drug)), 2)),
                                 rep(drug_names_sorted, 2),
                                 round(rbind(exp(em_ref_nmr),
                                             exp(pred_ref_nmr)), 2),
                                 rep(c("Estimation", "Prediction"),
                                     each = length(drug_names)))
    colnames(prepare_em_nma) <- c("order",
                                  "comparison",
                                  "mean", "lower", "upper",
                                  "interval")
    colnames(prepare_em_nmr) <- colnames(prepare_em_nma)
  } else if (!is.element(measure, c("OR", "RR", "ROM")) & model == "FE") {
    prepare_em_nma <- data.frame(as.factor(rev(seq_len(len_drug))),
                                 drug_names_sorted,
                                 round(em_ref_nma, 2))
    prepare_em_nmr <- data.frame(as.factor(rev(seq_len(len_drug))),
                                 drug_names_sorted,
                                 round(em_ref_nmr, 2))
    colnames(prepare_em_nma) <- c("order",
                                  "comparison",
                                  "mean", "lower", "upper")
    colnames(prepare_em_nmr) <- colnames(prepare_em_nma)
  } else if (is.element(measure, c("OR", "RR", "ROM")) & model == "FE") {
    prepare_em_nma <- data.frame(as.factor(rev(seq_len(len_drug))),
                                 drug_names_sorted,
                                 round(exp(em_ref_nma), 2))
    prepare_em_nmr <- data.frame(as.factor(rev(seq_len(len_drug))),
                                 drug_names_sorted,
                                 round(exp(em_ref_nmr), 2))
    colnames(prepare_em_nma) <- c("order",
                                  "comparison",
                                  "mean", "lower", "upper")
    colnames(prepare_em_nmr) <- colnames(prepare_em_nma)
  }
  prepare_em <- cbind(rbind(prepare_em_nma, prepare_em_nmr),
                      analysis = rep(c("Network meta-analysis",
                                       "Network meta-regression"),
                                     each = length(prepare_em_nma[, 1])))

  # Forest plots on credible/prediction intervals of comparisons with the
  # selected comparator
  subtitle <- paste("Results for", cov_value[[2]],
                   ifelse(length(unique(reg$covariate)) < 3, " ",
                          cov_value[[1]]))

  measure2 <- effect_measure_name(measure, lower = FALSE)
  caption <- if (full$D == 0 & is.element(measure, c("OR", "RR", "ROM"))) {
    paste(measure2, "< 1, favours the first arm.",
          measure2, "> 1, favours", compar)
  } else if (full$D == 1 & is.element(measure, c("OR", "RR", "ROM"))) {
    paste(measure2, "< 1, favours", compar,
          ".", measure2, "> 1, favours the first arm")
  } else if (full$D == 0 & !is.element(measure, c("OR", "RR", "ROM"))) {
    paste(measure2, "< 0, favours the first arm.",
          measure2, "> 0, favours", compar)
  } else if (full$D == 1 & !is.element(measure, c("OR", "RR", "ROM"))) {
    paste(measure2, "< 0, favours", compar,
          ".", measure2, "> 0, favours the first arm")
  }

  forest_plots <- if (model == "RE") {
    ggplot(data = prepare_em,
           aes(x = order,
               y = mean,
               ymin = lower,
               ymax = upper,
               group = analysis,
               colour = analysis)) +
      geom_hline(yintercept = ifelse(!is.element(
        measure, c("OR", "RR", "ROM")), 0, 1),
        lty = 1,
        size = 1,
        col = "grey60") +
      geom_linerange(size = 2,
                     position = position_dodge(width = 0.5)) +
      geom_errorbar(data = prepare_em,
                    aes(x = order,
                        y = mean,
                        ymin = lower,
                        ymax = upper),
                    size = 2,
                    position = position_dodge(width = 0.5),
                    width = 0.0) +
      geom_point(size = 1.5,
                 colour = "white",
                 stroke = 0.3,
                 position = position_dodge(width = 0.5)) +
      geom_text(aes(x = order,
                    y = mean,
                    label = paste0(mean, " ", "(", prepare_em[, 4], ",", " ",
                                   prepare_em[, 5], ")"),
                    hjust = 0,
                    vjust = -0.5),
                color = "black",
                size = 4.0,
                check_overlap = FALSE,
                parse = FALSE,
                position = position_dodge(width = 0.5),
                inherit.aes = TRUE,
                na.rm = TRUE) +
      labs(x = "",
           y = measure2,
           colour = "Analysis",
           subtitle = subtitle,
           caption = caption) +
      facet_wrap(~ interval, ncol = 2, scales = "fixed") +
      scale_x_discrete(breaks = as.factor(seq_len(len_drug)),
                       labels = drug_names_sorted[rev(seq_len(len_drug))]) +
      scale_color_manual(breaks = c("Network meta-analysis",
                                    "Network meta-regression"),
                         values = c("black", "#D55E00")) +
      geom_label(aes(x = unique(order[is.na(mean)]),
                     y = ifelse(!is.element(
                       measure, c("OR", "RR", "ROM")), -0.2, 0.65),
                     hjust = 0,
                     vjust = 1,
                     label = "Comparator intervention"),
                 fill = "beige",
                 colour = "black",
                 fontface = "plain",
                 size = 4) +
      scale_y_continuous(trans = ifelse(!is.element(
        measure, c("OR", "RR", "ROM")), "identity", "log10")) +
      coord_flip() +
      theme_classic() +
      theme(axis.text.x = element_text(color = "black", size = 12),
            axis.text.y = element_text(color = "black", size = 12),
            axis.title.x = element_text(color = "black", face = "bold",
                                        size = 12),
            legend.position = "bottom",
            legend.justification = c(0.43, 0),
            strip.text = element_text(color = "black", face = "bold",
                                      size = 12),
            legend.text =  element_text(color = "black", size = 12),
            legend.title = element_text(color = "black", face = "bold",
                                        size = 12),
            plot.caption = element_text(hjust = 0.01))
  } else {
    ggplot(data = prepare_em,
           aes(x = order,
               y = mean,
               ymin = lower,
               ymax = upper,
               group = analysis,
               colour = analysis)) +
      geom_hline(yintercept = ifelse(!is.element(
        measure, c("OR", "RR", "ROM")), 0, 1),
        lty = 1,
        size = 1,
        col = "grey60") +
      geom_linerange(size = 2,
                     position = position_dodge(width = 0.5)) +
      geom_errorbar(data = prepare_em,
                    aes(x = order,
                        y = mean,
                        ymin = lower,
                        ymax = upper),
                    size = 2,
                    position = position_dodge(width = 0.5),
                    width = 0.0) +
      geom_point(size = 1.5,
                 colour = "white",
                 stroke = 0.3,
                 position = position_dodge(width = 0.5)) +
      geom_text(aes(x = order,
                    y = mean,
                    label = paste0(mean, " ", "(", prepare_em[, 4], ",", " ",
                                   prepare_em[, 5], ")"),
                    hjust = 0,
                    vjust = -0.5),
                color = "black",
                size = 4.0,
                check_overlap = FALSE,
                parse = FALSE,
                position = position_dodge(width = 0.5),
                inherit.aes = TRUE,
                na.rm = TRUE) +
      labs(x = "",
           y = measure2,
           colour = "Analysis",
           subtitle = subtitle,
           caption = caption) +
      scale_x_discrete(breaks = as.factor(seq_len(len_drug)),
                       labels = drug_names_sorted[rev(seq_len(len_drug))]) +
      scale_color_manual(breaks = c("Network meta-analysis",
                                    "Network meta-regression"),
                         values = c("black", "#D55E00")) +
      geom_label(aes(x = unique(order[is.na(mean)]),
                     y = ifelse(!is.element(
                       measure, c("OR", "RR", "ROM")), -0.2, 0.65),
                     hjust = 0,
                     vjust = 1,
                     label = "Comparator intervention"),
                 fill = "beige",
                 colour = "black",
                 fontface = "plain",
                 size = 4) +
      scale_y_continuous(trans = ifelse(!is.element(
        measure, c("OR", "RR", "ROM")), "identity", "log10")) +
      coord_flip() +
      theme_classic() +
      theme(axis.text.x = element_text(color = "black", size = 12),
            axis.text.y = element_text(color = "black", size = 12),
            strip.text = element_text(color = "black", face = "bold",
                                      size = 12),
            axis.title.x = element_text(color = "black", face = "bold",
                                        size = 12),
            legend.position = "bottom",
            legend.justification = c(0.23, 0),
            legend.text =  element_text(color = "black", size = 12),
            legend.title = element_text(color = "black", face = "bold",
                                        size = 12),
            plot.caption = element_text(hjust = 0.01))
  }

 return(forest_plots)
}

#' Comparator-specific forest plot for network meta-analysis
#'
#' @description Provides a forest plot with the posterior mean and 95\% credible
#'   and prediction intervals for comparisons with the selected intervention
#'   (comparator) in the network.
#'
#' @param full An object of S3 class \code{\link{run_model}}. See 'Value' in
#'   \code{\link{run_model}}.
#' @param compar A character to indicate the comparator intervention. It must be
#'   any name found in \code{drug_names}.
#' @param drug_names A vector of labels with the name of the interventions in
#'   the order they appear in the argument \code{data} of
#'   \code{\link{run_model}}.
#'
#' @return A panel of two forest plots: (1) a forest plot on the effect
#'   estimates and predictions of comparisons with the selected intervention in
#'   the network, and (2) a forest plot on the posterior mean and 95\% credible
#'   interval of SUCRA values of the interventions (Salanti et al., 2011).
#'
#' @details The y-axis of the forest plot on the effect sizes displays the
#'   labels of the interventions in the network; the selected intervention that
#'   comprises the \code{compar} argument is annotated in the plot with the
#'   label 'Comparator intervention'.
#'   For each comparison with the selected intervention, the 95\% credible and
#'   prediction intervals are displayed as overlapping lines in different
#'   colours. The corresponding numerical results are displayed above each line:
#'   95\% credible intervals are found in parentheses, and 95\% predictive
#'   intervals are found in brackets. Odds ratios, relative risks, and ratio of
#'   means are reported in the original scale after exponentiation of the
#'   logarithmic scale.
#'
#'   The y-axis for the forest plot on the SUCRA values displays the
#'   labels of the interventions in the network.
#'   The corresponding numerical results are displayed above each line.
#'   Three coloured rectangles appear in the forest plot: a red rectangle for
#'   SUCRA values up to 50\%, a yellow rectangular for SUCRA values between
#'   50\% and 80\%, and a green rectangle for SUCRA values over 80\%.
#'   Interventions falling at the green area are considered as the highest
#'   ranked interventions, whilst interventions falling at the red area are
#'   considered as the lowest ranked interventions.
#'
#'   In both plots, the interventions are sorted in descending order of their
#'   SUCRA values.
#'
#'   \code{forestplot} can be used only for a network of interventions.
#'   Otherwise, the execution of the function will be stopped and an error
#'   message will be printed on the R console.
#'
#' @author {Loukia M. Spineli}
#'
#' @seealso \code{\link{run_model}}
#'
#' @references
#' Salanti G, Ades AE, Ioannidis JP. Graphical methods and numerical summaries
#' for presenting results from multiple-treatment meta-analysis: an overview and
#' tutorial. \emph{J Clin Epidemiol} 2011;\bold{64}(2):163--71.
#' doi: 10.1016/j.jclinepi.2010.03.016
#'
#' @examples
#' data("nma.liu2013")
#'
#' # Show the first six trials of the dataset (one-trial-per-row format)
#' head(nma.liu2013)
#'
#' # Read results from 'run_model' (using the default arguments)
#' res <- readRDS(system.file('extdata/res_liu.rds', package = 'rnmamod'))
#'
#' # The names of the interventions in the order they appear in the dataset
#' interv_names <- c("placebo", "pramipexole", "serotonin-norepinephrine
#' reuptake inhibitor", "serotonin reuptake inhibitor", "tricyclic
#' antidepressant", "pergolide")
#'
#' # Create the forest plot
#' forestplot(full = res,
#'            compar = "placebo",
#'            drug_names = interv_names)
#'
#' @export
forestplot <- function(full, compar, drug_names) {

  if (full$type != "nma" || is.null(full$type)) {
    stop("'full' must be an object of S3 class 'run_model'.",
         call. = FALSE)
  }

  drug_names <- if (missing(drug_names)) {
    stop("The argument 'drug_names' has not been defined.",
         call. = FALSE)
  } else {
    drug_names
  }

  if (length(drug_names) == 2) {
    stop("This function is *not* relevant for a pairwise meta-analysis.",
         call. = FALSE)
  }

  compar <- if (missing(compar)) {
    stop("The argument 'compar' has not been defined.", call. = F)
  } else if (!is.element(compar, drug_names)) {
    stop("The value of 'compar' is not found in the argument 'drug_names'.",
         call. = FALSE)
  } else if (is.element(compar, drug_names)) {
    compar
  }

  # A matrix with all possible comparisons in the network
  poss_pair_comp1 <- data.frame(exp = t(combn(drug_names, 2))[, 2],
                                comp = t(combn(drug_names, 2))[, 1])
  poss_pair_comp2 <- data.frame(exp = t(combn(drug_names, 2))[, 1],
                                comp = t(combn(drug_names, 2))[, 2])
  poss_pair_comp  <- rbind(poss_pair_comp1, poss_pair_comp2)

  measure <- full$measure
  model <- full$model
  sucra <- full$SUCRA
  em_ref00 <- cbind(rbind(data.frame(mean = full$EM[, 5], #1
                                     lower = full$EM[, 3],
                                     upper = full$EM[, 7]),
                          data.frame(mean = full$EM[, 5] * (-1), #1
                                     lower = full$EM[, 7] * (-1),
                                     upper = full$EM[, 3] * (-1))),
                    poss_pair_comp)
  em_subset <- subset(em_ref00, em_ref00[5] == compar)
  em_ref0 <- rbind(em_subset[, 1:3], c(rep(NA, 3)))
  sucra_new <- data.frame(sucra[, 1],
                          drug_names)[order(match(data.frame(sucra[, 1],
                                                             drug_names)[, 2],
                                                  em_subset[, 4])), 1]
  em_ref <- em_ref0[order(sucra_new, decreasing = TRUE), ]
  rownames(em_ref) <- NULL

  # Posterior results on the predicted estimates of comparisons with the
  # selected comparator as reference
  if (model == "RE") {
    pred_ref00 <- cbind(rbind(data.frame(mean = full$EM_pred[, 5], #1
                                         lower = full$EM_pred[, 3],
                                         upper = full$EM_pred[, 7]),
                              data.frame(mean = full$EM_pred[, 5] * (-1), #1
                                         lower = full$EM_pred[, 7] * (-1),
                                         upper = full$EM_pred[, 3] * (-1))),
                        poss_pair_comp)
    pred_subset <- subset(pred_ref00, pred_ref00[5] == compar)
    pred_ref0 <- rbind(pred_subset[, 1:3], c(rep(NA, 3)))

    # Sort by SUCRA in decreasing order and remove the reference intervention
    pred_ref <- pred_ref0[order(sucra_new, decreasing = TRUE), ]
    rownames(pred_ref) <- NULL
  }

  # Sort the drugs by their SUCRA in decreasing order and remove the reference
  # intervention (number 1)
  drug_names_sorted <- drug_names[order(sucra[, 1], decreasing = TRUE)]
  len_drug_names <- length(drug_names_sorted)

  # Create a data-frame with credible and prediction intervals of comparisons
  # with the reference intervention
  if (!is.element(measure, c("OR", "RR", "ROM")) & model == "RE") {
    prepare_em <- data.frame(as.factor(rep(rev(seq_len(len_drug_names)), 2)),
                             rep(drug_names_sorted, 2),
                             round(rbind(em_ref, pred_ref), 2),
                             rep(c("Credible interval", "Prediction interval"),
                                 each = length(drug_names)))
    colnames(prepare_em) <- c("order",
                              "comparison",
                              "mean", "lower", "upper",
                              "interval")
  } else if (is.element(measure, c("OR", "RR", "ROM")) &
             model == "RE") {
    prepare_em <- data.frame(as.factor(rep(rev(seq_len(len_drug_names)), 2)),
                             rep(drug_names_sorted, 2),
                             round(rbind(exp(em_ref), exp(pred_ref)), 2),
                             rep(c("Credible interval", "Prediction interval"),
                                 each = length(drug_names)))
    colnames(prepare_em) <- c("order",
                              "comparison",
                              "mean", "lower", "upper",
                              "interval")
  } else if (!is.element(measure, c("OR", "RR", "ROM")) &
             model == "FE") {
    prepare_em <- data.frame(as.factor(rev(seq_len(len_drug_names))),
                             drug_names_sorted,
                             round(em_ref, 2))
    colnames(prepare_em) <- c("order",
                              "comparison",
                              "mean", "lower", "upper")
  } else if (is.element(measure, c("OR", "RR", "ROM")) &
             model == "FE") {
    prepare_em <- data.frame(as.factor(rev(seq_len(len_drug_names))),
                             drug_names_sorted,
                             round(exp(em_ref), 2))
    colnames(prepare_em) <- c("order",
                              "comparison",
                              "mean", "lower", "upper")
  }

  # Create a data-frame with the SUCRA values
  prepare_sucra <- data.frame(rev(seq_len(len_drug_names)),
                              drug_names[order(sucra[, 1], decreasing = TRUE)],
                              sucra[order(sucra[, 1], decreasing = TRUE),
                                    c(1, 3, 7)])
  colnames(prepare_sucra) <- c("order",
                               "intervention",
                               "mean", "lower", "upper")
  rownames(prepare_sucra) <- NULL

  # Forest plots on credible/prediction intervals of comparisons with the
  # reference
  measure2 <- effect_measure_name(full$measure, lower = FALSE)
  caption <- if (full$D == 0 & is.element(measure, c("OR", "RR", "ROM"))) {
    paste(measure2, "< 1, favours the first arm.",
          measure2, "> 1, favours", paste0(compar, "."))
  } else if (full$D == 1 & is.element(measure, c("OR", "RR", "ROM"))) {
    paste(measure2, "< 1, favours", compar,
          ".", measure2, "> 1, favours the first arm.")
  } else if (full$D == 0 & !is.element(measure, c("OR", "RR", "ROM"))) {
    paste(measure2, "< 0, favours the first arm.",
          measure2, "> 0, favours", paste0(compar, "."))
  } else if (full$D == 1 & !is.element(measure, c("OR", "RR", "ROM"))) {
    paste(measure2, "< 0, favours", compar,
          ".", measure2, "> 0, favours the first arm.")
  }

  p1 <- if (model == "RE") {
    ggplot(data = prepare_em[(length(drug_names_sorted) + 1):
                               (length(drug_names_sorted) * 2), ],
           aes(x = order,
               y = mean,
               ymin = lower,
               ymax = upper,
               colour = interval)) +
      #annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = 0.2,
      #         fill = "grey70") +
      #annotate("rect",xmin = -Inf, xmax = Inf, ymin = 0.2, ymax = 0.7,
      #         fill = "grey85") +
      #annotate("rect",xmin = -Inf, xmax = Inf, ymin = 0.7, ymax = 1.0,
      #         fill = "white") +
      #annotate("rect",xmin = -Inf, xmax = Inf, ymin = 1.0, ymax = 1.5,
      #         fill = "white") +
      #annotate("rect",xmin = -Inf, xmax = Inf, ymin = 1.5, ymax = 5.0,
      #         fill = "grey85") +
      #annotate("rect", xmin = -Inf, xmax = Inf, ymin = 5.0, ymax = Inf,
      #         fill = "grey70") +
      geom_hline(yintercept = ifelse(!is.element(measure,
                                                 c("OR", "RR", "ROM")), 0, 1),
                 lty = 1,
                 size = 1,
                 col = "grey60") +
      geom_linerange(size = 2,
                     position = position_dodge(width = 0.5)) +
      geom_errorbar(data = prepare_em[seq_len(len_drug_names), ],
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
                    label = paste0(mean,
                                   " ",
                                   "(",
                                   prepare_em[seq_len(len_drug_names), 4],
                                   ",",
                                   " ",
                                   prepare_em[seq_len(len_drug_names), 5],
                                   ")",
                                   " ",
                                   "[",
                                   prepare_em[(length(drug_names_sorted) + 1):
                                                (length(drug_names_sorted) * 2),
                                              4],
                                   ",",
                                   " ",
                                   prepare_em[(length(drug_names_sorted) + 1):
                                                (length(drug_names_sorted) * 2),
                                              5],
                                   "]"),
                    hjust = 0,
                    vjust = -0.5),
                color = "black",
                size = 4.0,
                check_overlap = FALSE,
                parse = FALSE,
                position = position_dodge(width = 0.5),
                inherit.aes = TRUE, na.rm = TRUE) +
      scale_x_discrete(breaks = as.factor(seq_len(len_drug_names)),
                       labels = drug_names_sorted[rev(
                         seq_len(len_drug_names))]) +
      scale_y_continuous(trans = ifelse(!is.element(
        measure, c("OR", "RR", "ROM")), "identity", "log10")) +
      scale_color_manual(breaks = c("Credible interval", "Prediction interval"),
                         values = c("black", "#D55E00")) +
      geom_label(aes(x = order[is.na(mean)],
                     y = ifelse(!is.element(measure, c("OR", "RR", "ROM")),
                                0, 1), # -0.2, 0.65
                     hjust = 0,
                     vjust = 1,
                     label = "Comparator intervention"),
                 fill = "beige",
                 colour = "black",
                 fontface = "plain",
                 size = 4) +
      labs(x = "", y = measure2, colour = "Analysis", caption = caption) +
      coord_flip() +
      theme_classic() +
       theme(axis.text.x = element_text(color = "black", size = 12),
             axis.text.y = element_text(color = "black", size = 12),
             axis.title.x = element_text(color = "black", face = "bold",
                                         size = 12),
             legend.position = "bottom",
             legend.text =  element_text(color = "black", size = 12),
             legend.title = element_blank(),
             plot.caption = element_text(hjust = 0.01))
  } else {
    ggplot(data = prepare_em,
           aes(x = order,
               y = mean,
               ymin = lower,
               ymax = upper)) +
      geom_hline(yintercept = ifelse(!is.element(
        measure, c("OR", "RR", "ROM")), 0, 1),
        lty = 2,
        size = 1.3,
        col = "grey53") +
      geom_linerange(size = 2, position = position_dodge(width = 0.5)) +
      geom_errorbar(data = prepare_em[seq_len(len_drug_names), ],
                    aes(x = order,
                        y = mean,
                        ymin = lower,
                        ymax = upper),
                    size = 2, position = position_dodge(width = 0.5),
                    colour = "black", width = 0.0) +
      geom_point(size = 1.5, colour = "white", stroke = 0.3,
                 position = position_dodge(width = 0.5)) +
      geom_text(aes(x = order,
                    y = mean,
                    label = paste0(mean,
                                   " ",
                                   "(",
                                   prepare_em[seq_len(len_drug_names), 4],
                                   ",",
                                   " ",
                                   prepare_em[seq_len(len_drug_names), 5],
                                   ")"),
                    hjust = 0, vjust = -0.5),
                color = "black", size = 4.0, check_overlap = FALSE,
                parse = FALSE,
                position = position_dodge(width = 0.5), inherit.aes = TRUE,
                na.rm = TRUE) +
      scale_x_discrete(breaks = as.factor(seq_len(len_drug_names)),
                       labels = drug_names_sorted[rev(
                         seq_len(len_drug_names))]) +
      geom_label(aes(x = order[is.na(mean)],
                     y = ifelse(!is.element(
                       measure, c("OR", "RR", "ROM")), 0, 1), # -0.2, 0.65
                     hjust = 0, vjust = 1, label = "Comparator intervention"),
                 fill = "beige", colour = "black", fontface = "plain",
                 size = 4) +
      scale_y_continuous(trans = ifelse(!is.element(
        measure, c("OR", "RR", "ROM")), "identity", "log10")) +
      #geom_rect(aes(xmin = 0, xmax = Inf, ymin = 0, ymax = 0.5,
      #              fill = "lowest"),
      #          alpha = 0.02) +
      #geom_rect(aes(xmin = 0, xmax = Inf, ymin = 0.5, ymax = 0.8,
      #              fill = "intermediate"),
      #          alpha = 0.02) +
      scale_fill_manual(name = " ",
                        values = c("lowest" = "white",
                                   "intermediate" = "white",
                                   "highest" = "white")) +
      labs(x = "", y = measure2, caption = caption) +
      coord_flip(clip = "off") +
      theme_classic() +
      theme(axis.text.x = element_text(color = "black", size = 12),
            axis.text.y = element_text(color = "black", size = 12),
            axis.title.x = element_text(color = "black", face = "bold",
                                        size = 12),
            legend.position = "bottom",
            legend.justification = c(0.13, 0), # 0.23
            legend.text =  element_blank(),
            legend.title = element_blank(),
            plot.caption = element_text(hjust = 0.01))
  }


  # Forest plots of SUCRA per intervention
  p2 <- ggplot(data = prepare_sucra[seq_len(len_drug_names), ],
               aes(x = as.factor(order),
                   y = mean,
                   ymin = lower,
                   ymax = upper)) +
    geom_rect(aes(xmin = 0, xmax = Inf, ymin = 0, ymax = 0.5,
                  fill = "lowest"),
              alpha = 0.02) +
    geom_rect(aes(xmin = 0, xmax = Inf, ymin = 0.5, ymax = 0.8,
                  fill = "intermediate"),
              alpha = 0.02) +
    geom_rect(aes(xmin = 0, xmax = Inf, ymin = 0.8, ymax = 1.0,
                  fill = "highest"),
              alpha = 0.02) +
    geom_linerange(size = 2, position = position_dodge(width = 0.5)) +
    geom_point(size = 1.5, colour = "white", stroke = 0.3,
               position = position_dodge(width = 0.5)) +
    geom_text(aes(x = as.factor(order), y = round(mean, 2),
                  label = paste0(round(mean * 100, 0),
                                 " ",
                                 "(",
                                 round(lower * 100, 0),
                                 ",",
                                 " ",
                                 round(upper * 100, 0), ")"),
                        hjust = 0, vjust = -0.5),
              color = "black", size = 4.0, check_overlap = FALSE, parse = FALSE,
              position = position_dodge(width = 0.5), inherit.aes = TRUE) +
    scale_fill_manual(name = "Ranked",
                      values = c("lowest" = "#D55E00",
                                 "intermediate" = "orange",
                                 "highest" = "#009E73")) +
    labs(x = "", y = "Surface under the cumulative ranking curve value",
         caption = " ") +
    scale_x_discrete(breaks = as.factor(seq_len(len_drug_names)),
                     labels = prepare_sucra$intervention[rev(
                       seq_len(len_drug_names))]) +
    scale_y_continuous(labels = percent) +
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

  # Bring together both forest-plots
  forest_plots <- suppressWarnings(
    ggarrange(p1, p2 + guides(fill = guide_legend(override.aes =
                                                    list(alpha = 0.4))),
              nrow = 1, ncol = 2, labels = c("A)", "B)"),
              common.legend = FALSE, legend = "bottom")
  )

  return(forest_plots)
}

#' Forest-plot of comparisons with the reference intervention
#'
#' @export
forestplot.ref <- function(net, drug.names) {


  options(warn = -1)

  ## The results on the following parameters will be used:
  # Effect measure
  measure <- effect.measure.name(net$measure)

  # Analysis model
  model <- net$model

  # Posterior results on the SUCRA of each intervention
  sucra <- net$SUCRA

  # Posterior results on the effect estimates of comparisons with the reference intervention of the network
  EM.ref0 <- rbind(rep(NA, 3), net$EM.ref[, c(1, 3, 7)])

  # Sort by SUCRA in decreasing order and remove the reference intervention (number 1)
  EM.ref <- EM.ref0[order(sucra[, 1], decreasing = T), ]

  # Posterior results on the predicted estimates of comparisons with the reference intervention of the network
  pred.ref0 <- if (model == "RE") {
    rbind(rep(NA, 3), net$pred.ref[, c(1, 3, 7)])
  } else {
    NA
  }

  # Sort by SUCRA in decreasing order and remove the reference intervention (number 1)
  pred.ref <- if (model == "RE") {
    pred.ref0[order(sucra[, 1], decreasing = T), ]
  } else {
    NA
  }


  # Sort the drugs by their SUCRA in decreasing order and remove the reference intervention (number 1)
  drug.names.sorted <- drug.names[order(sucra[, 1], decreasing = T)]



  ## Create a data-frame with credible and predictive intervals of comparisons with the reference intervention
  if (!is.element(measure, c("Odds ratio", "Ratio of means")) & model == "RE") {
    prepare.EM <- data.frame(rep(length(drug.names):1, 2), rep(drug.names.sorted, 2), round(rbind(EM.ref, pred.ref), 2), rep(c("Credible interval", "Predictive interval"), each = length(drug.names)))
    colnames(prepare.EM) <- c("order", "comparison", "mean", "lower", "upper", "interval")
  } else if (is.element(measure, c("Odds ratio", "Ratio of means")) & model == "RE") {
    prepare.EM <- data.frame(rep(length(drug.names):1, 2), rep(drug.names.sorted, 2), round(rbind(exp(EM.ref), exp(pred.ref)), 2), rep(c("Credible interval", "Predictive interval"), each = length(drug.names)))
    colnames(prepare.EM) <- c("order", "comparison", "mean", "lower", "upper", "interval")
  } else if (!is.element(measure, c("Odds ratio", "Ratio of means")) & model == "FE") {
    prepare.EM <- data.frame(length(drug.names):1, drug.names.sorted, round(EM.ref, 2))
    colnames(prepare.EM) <- c("order", "comparison", "mean", "lower", "upper")
  } else if (is.element(measure, c("Odds ratio", "Ratio of means")) & model == "RE") {
    prepare.EM <- data.frame(length(drug.names):1, drug.names.sorted, round(exp(EM.ref), 2))
    colnames(prepare.EM) <- c("order", "comparison", "mean", "lower", "upper")
  }
  rownames(prepare.EM) <- NULL



  ## Create a data-frame with the SUCRA values
  prepare.sucra <- data.frame(length(drug.names):1, drug.names[order(sucra[, 1], decreasing = T)], sucra[order(sucra[, 1], decreasing = T), c(1, 3, 7)])
  colnames(prepare.sucra) <- c("order", "intervention", "mean", "lower", "upper")
  rownames(prepare.sucra) <- NULL



  ## Forest plots on credible/predictive intervals of comparisons with the reference
  p1 <- if (model == "RE") {
    ggplot(data = prepare.EM[(length(drug.names.sorted) + 1):(length(drug.names.sorted)*2), ], aes(x = order, y = mean, ymin = lower, ymax = upper, colour = interval)) +
      geom_linerange(size = 2, position = position_dodge(width = 0.5)) +
      geom_errorbar(data = prepare.EM[1:length(drug.names.sorted), ], aes(x = as.factor(order), y = mean, ymin = lower, ymax = upper), size = 2, position = position_dodge(width = 0.5), width = 0.1) +
      geom_hline(yintercept = ifelse(!is.element(measure, c("Odds ratio", "Ratio of means")), 0, 1), lty = 2, size = 1.3, col = "grey53") +
      geom_point(size = 1.5,  colour = "white", stroke = 0.3, position = position_dodge(width = 0.5)) +
      geom_text(aes(x = order, y = mean, label = paste0(mean, " ", " ", "(", prepare.EM[1:length(drug.names.sorted), 4], ",", " ", prepare.EM[1:length(drug.names.sorted), 5], ")",
                                                        " ", "[", prepare.EM[(length(drug.names.sorted) + 1):(length(drug.names.sorted)*2), 4], ",", " ", prepare.EM[(length(drug.names.sorted) + 1):(length(drug.names.sorted)*2), 5], "]"),
                    hjust = 0, vjust = -0.4), color = "blue", size = 4.0, check_overlap = F, parse = F, position = position_dodge(width = 0.5), inherit.aes = T, na.rm = T) +
      labs(x = "", y = measure, colour = "Analysis") +
      scale_x_discrete(breaks = as.factor(1:length(drug.names.sorted)), labels = drug.names.sorted[length(drug.names.sorted):1]) +
      scale_color_manual(breaks = c("Credible interval", "Predictive interval"), values = c("black", "#D55E00")) +
      geom_label(aes(x = as.factor(order)[is.na(mean)], y = ifelse(!is.element(measure, c("Odds ratio", "Ratio of means")), -0.2, 0.65), hjust = 0, vjust = 1, label = "Reference intervention"), fill = "beige", colour = "black", fontface = "plain", size = 4) +
      scale_y_continuous(trans = ifelse(!is.element(measure, c("Odds ratio", "Ratio of means")), "identity", "log10")) +
      coord_flip(clip = "off") +
      theme_classic() +
      theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12),
            axis.title.x = element_text(color = "black", face = "bold", size = 12), legend.position = "bottom", legend.justification = c(0.23, 0),
            legend.text =  element_text(color = "black", size = 12), legend.title =  element_text(color = "black", face = "bold", size = 12))
  } else {
    ggplot(data = prepare.EM, aes(x = order, y = mean, ymin = lower, ymax = upper)) +
      geom_linerange(size = 2, position = position_dodge(width = 0.5)) +
      geom_errorbar(data = prepare.EM[1:length(drug.names.sorted), ], aes(x = as.factor(order), y = mean, ymin = lower, ymax = upper), size = 2, position = position_dodge(width = 0.5), width = 0.1, colour = "black") +
      geom_hline(yintercept = ifelse(!is.element(measure, c("Odds ratio", "Ratio of means")), 0, 1), lty = 2, size = 1.3, col = "grey53") +
      geom_point(size = 1.5,  colour = "white", stroke = 0.3, position = position_dodge(width = 0.5)) +
      geom_text(aes(x = order, y = mean, label = paste0(mean, " ", " ", "(", prepare.EM[1:length(drug.names.sorted), 4], ",", " ", prepare.EM[1:length(drug.names.sorted), 5], ")"),
                    hjust = 0, vjust = -0.4), color = "blue", size = 4.0, check_overlap = F, parse = F, position = position_dodge(width = 0.5), inherit.aes = T, na.rm = T) +
      labs(x = "", y = measure) +
      scale_x_discrete(breaks = as.factor(1:length(drug.names.sorted)), labels = drug.names.sorted[length(drug.names.sorted):1]) +
      geom_label(aes(x = as.factor(order)[is.na(mean)], y = ifelse(!is.element(measure, c("Odds ratio", "Ratio of means")), -0.2, 0.65), hjust = 0, vjust = 1, label = "Reference intervention"), fill = "beige", colour = "black", fontface = "plain", size = 4) +
      scale_y_continuous(trans = ifelse(!is.element(measure, c("Odds ratio", "Ratio of means")), "identity", "log10")) +
      coord_flip(clip = "off") +
      theme_classic() +
      theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12),
            axis.title.x = element_text(color = "black", face = "bold", size = 12), legend.position = "bottom", legend.justification = c(0.23, 0),
            legend.text =  element_text(color = "black", size = 12), legend.title =  element_text(color = "black", face = "bold", size = 12))
  }




  ## Forest plots of SUCRA per intervention
  p2 <- ggplot(data = prepare.sucra[1:length(drug.names), ], aes(x = as.factor(order), y = mean, ymin = lower, ymax = upper)) +
          geom_linerange(size = 2, position = position_dodge(width = 0.5)) +
          geom_point(size = 1.5,  colour = "white", stroke = 0.3, position = position_dodge(width = 0.5)) +
          geom_text(aes(x = as.factor(order), y = round(mean, 2), label = paste0(round(mean*100, 0), " ", "(", round(lower*100, 0), ",", " ", round(upper*100, 0), ")" ),
                        hjust = 0, vjust = -0.4), color = "black", size = 4.0, check_overlap = F, parse = F, position = position_dodge(width = 0.5), inherit.aes = T) +
          labs(x = "", y = "Surface under the cumulative ranking curve value") +
          scale_x_discrete(breaks = as.factor(1:length(drug.names)), labels = prepare.sucra$intervention[length(drug.names):1]) +
          scale_y_continuous(labels = scales::percent) +
          coord_flip() +
          theme_classic() +
          theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12),
                axis.title.x = element_text(color = "black", face = "bold", size = 12))



  ## Bring together both forest-plots
  forest.plots <- ggarrange(p1, p2, nrow = 1, ncol = 2, labels = c("A)", "B)"), common.legend = T, legend = "bottom")


  return(forest.plots = forest.plots)

}

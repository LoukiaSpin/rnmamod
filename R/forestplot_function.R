#' Forest-plot of comparisons with the selected intervention
#'
#' @description This function illustrates a forest plot of the posterior mean and 95\% credible and predictive interval of comparisons with the selected intervention of the network.
#'
#' @param full An object of S3 class \code{\link{run.model}}. See 'Value' in \code{\link{run.model}}.
#' @param compar A character to indicate the comparator intervention. it must be any name found in \code{drug.names}.
#' @param drug.names A vector of labels with the name of the interventions in the order they appear in the argument \code{data} of \code{\link{run.model}}. If the argument \code{drug.names} is not defined, the order of the interventions
#'   as they appear in \code{data} is used, instead.
#'
#' @return A panel of two forest plots: (1) a forest plot on the effect estimates and predictions of comparisons with the selected intervention of the network, and
#' (2) a forest plot on the posterior mean and 95\% credible interval of SUCRA values of the interventions (Salanti et al., 2011).
#'
#' @details The x-axis in the forest plot of effect sizes displays all interventions in the network; the selected intervention that comprises the \code{compar} is indicated in the plot with a homonymous label.
#'   For each comparison with the selected intervention, the 95\% credible and predictive intervals are displayed as overlapping lines with different colours. When the between-trial variance is very low, these two intervals become indiscernible.
#'   Furthermore, the corresponding numerical results are displayed above each line: 95\% credible intervals are found in parentheses, and 95\% predictive intervals are found in brackets.
#'   Odds ratio and ratio of means are reported in the original scale after exponentiation of the logarithmic scale.
#'
#'   The interventions are sorted in the descending order of their SUCRA values.
#'
#'   \code{forestplot} can be used only for a network of interventions. In the case of two interventions, the execution of the function will be stopped and an error message will be printed in the R console.
#'
#' @author {Loukia M. Spineli}
#'
#' @seealso \code{\link{run.model}}
#'
#' @references
#' Salanti G, Ades AE, Ioannidis JP. Graphical methods and numerical summaries for presenting results from multiple-treatment meta-analysis: an overview and tutorial. \emph{J Clin Epidemiol} 2011;\bold{64}(2):163--71. [\doi{10.1016/j.jclinepi.2010.03.016}]
#'
#' @examples
#' data("nma.liu2013")
#'
#' # Show the first six trials of the dataset (one-trial-per-row format)
#' head(nma.liu2013)
#' #            study t1 t2 t3 r1 r2 r3 m1 m2 m3  n1  n2 n3
#' #    Richard, 2012  1  3  4 15 16 23  6  8  4  39  42 34
#' #     Barone, 2010  1  2 NA 27 38 NA 19 20 NA 152 144 NA
#' # Weinbtraub, 2010  1  3 NA  2  5 NA  6  6 NA  27  28 NA
#' #      Menza, 2009  1  4  5  4  2  9  6  7  5  17  18 17
#' #      Devos, 2008  1  4  5  4  8 11  0  2  1  16  15 17
#' #   Antonini, 2006  4  5 NA 10  8 NA  4  4 NA  16  15 NA
#'
#' \dontrun{
#' # Perform a random-effects network meta-analysis
#' res <- run.model(data = nma.liu2013,
#'                  measure = "OR",
#'                  model = "RE",
#'                  assumption = "IDE-ARM",
#'                  heter.prior = list("halfnormal", 0, 1),
#'                  mean.misspar = c(0, 0),
#'                  var.misspar = 1,
#'                  D = 1,
#'                  n.chains = 3,
#'                  n.iter = 10000,
#'                  n.burnin = 1000,
#'                  n.thin = 1)
#'
#' # The names of the interventions in the order they appear in the dataset
#' interv.names <- c("placebo", "pramipexole", "serotonin–norepinephrine reuptake inhibitor",
#'                   "serotonin reuptake inhibitor", "tricyclic antidepressant", "pergolide")
#'
#' # Create the league heatmap
#' forestplot(full = res, compar = "placebo", drug.names = interv.names)
#' }
#'
#' @export
forestplot <- function(full, compar,  drug.names) {

  options(warn = -1)

  drug.names <- if (missing(drug.names)) {
    message(cat(paste0("\033[0;", col = 32, "m", txt = "The argument 'drug.names' has not been defined. The intervention ID, as specified in 'data' is used as intervention names", "\033[0m", "\n")))
    nt <- length(full$SUCRA[, 1])
    as.character(1:nt)
  } else {
    drug.names
  }


  if(length(drug.names) == 2) {
    stop("This function is *not* relevant for a pairwise meta-analysis", call. = F)
  }


  compar <- if(missing(compar)) {
    stop("The argument 'compar' has not been defined", call. = F)
  } else if(!is.element(compar, drug.names)) {
    stop("The value of 'compar' is not found in the argument 'drug.names'", call. = F)
  } else if(is.element(compar, drug.names)) {
    compar
  }



  ## A matrix with all possible comparisons in the network
  poss.pair.comp1 <- data.frame(exp = t(combn(drug.names, 2))[, 2], comp = t(combn(drug.names, 2))[, 1])
  poss.pair.comp2 <- data.frame(exp = t(combn(drug.names, 2))[, 1], comp = t(combn(drug.names, 2))[, 2])
  poss.pair.comp <- rbind(poss.pair.comp1, poss.pair.comp2)


  measure <- effect.measure.name(full$measure)
  model <- full$model
  sucra <- full$SUCRA
  EM.ref00 <- cbind(rbind(data.frame(mean = full$EM[, 1], lower = full$EM[, 3], upper = full$EM[, 7]),
                          data.frame(mean = full$EM[, 1]*(-1), lower = full$EM[, 7]*(-1), upper = full$EM[, 3]*(-1))),
                    poss.pair.comp)
  EM.subset <- subset(EM.ref00, EM.ref00[5] == compar)

  EM.ref0 <- rbind(EM.subset[, 1:3], c(rep(NA, 3)))

  sucra.new <- data.frame(sucra[, 1], drug.names)[order(match(data.frame(sucra[, 1], drug.names)[, 2], EM.subset[, 4])), 1]

  EM.ref <- EM.ref0[order(sucra.new, decreasing = T), ]
  rownames(EM.ref) <- NULL


  # Posterior results on the predicted estimates of comparisons with the selected comparator as reference
  pred.ref00 <- cbind(rbind(data.frame(mean = full$EM.pred[, 1], lower = full$EM.pred[, 3], upper = full$EM.pred[, 7]),
                           data.frame(mean = full$EM.pred[, 1]*(-1), lower = full$EM.pred[, 7]*(-1), upper = full$EM.pred[, 3]*(-1))),
                      poss.pair.comp)

  pred.subset <- subset(pred.ref00, pred.ref00[5] == compar)

  pred.ref0 <- rbind(pred.subset[, 1:3], c(rep(NA, 3)))


  # Sort by SUCRA in decreasing order and remove the reference intervention (number 1)
  pred.ref <- if (model == "RE") {
    pred.ref0[order(sucra.new, decreasing = T), ]
  } else {
    NA
  }
  rownames(pred.ref) <- NULL


  # Sort the drugs by their SUCRA in decreasing order and remove the reference intervention (number 1)
  drug.names.sorted <- drug.names[order(sucra[, 1], decreasing = T)]


  ## Create a data-frame with credible and predictive intervals of comparisons with the reference intervention
  if (!is.element(measure, c("Odds ratio", "Ratio of means")) & model == "RE") {
    prepare.EM <- data.frame(as.factor(rep(length(drug.names):1, 2)), rep(drug.names.sorted, 2), round(rbind(EM.ref, pred.ref), 2), rep(c("Credible interval", "Predictive interval"), each = length(drug.names)))
    colnames(prepare.EM) <- c("order", "comparison", "mean", "lower", "upper", "interval")
  } else if (is.element(measure, c("Odds ratio", "Ratio of means")) & model == "RE") {
    prepare.EM <- data.frame(as.factor(rep(length(drug.names):1, 2)), rep(drug.names.sorted, 2), round(rbind(exp(EM.ref), exp(pred.ref)), 2), rep(c("Credible interval", "Predictive interval"), each = length(drug.names)))
    colnames(prepare.EM) <- c("order", "comparison", "mean", "lower", "upper", "interval")
  } else if (!is.element(measure, c("Odds ratio", "Ratio of means")) & model == "FE") {
    prepare.EM <- data.frame(as.factor(length(drug.names):1), drug.names.sorted, round(EM.ref, 2))
    colnames(prepare.EM) <- c("order", "comparison", "mean", "lower", "upper")
  } else if (is.element(measure, c("Odds ratio", "Ratio of means")) & model == "RE") {
    prepare.EM <- data.frame(as.factor(length(drug.names):1), drug.names.sorted, round(exp(EM.ref), 2))
    colnames(prepare.EM) <- c("order", "comparison", "mean", "lower", "upper")
  }


  ## Create a data-frame with the SUCRA values
  prepare.sucra <- data.frame(length(drug.names):1, drug.names[order(sucra[, 1], decreasing = T)], sucra[order(sucra[, 1], decreasing = T), c(1, 3, 7)])
  colnames(prepare.sucra) <- c("order", "intervention", "mean", "lower", "upper")
  rownames(prepare.sucra) <- NULL


  ## Forest plots on credible/predictive intervals of comparisons with the reference
  p1 <- if (model == "RE") {
    ggplot(data = prepare.EM[(length(drug.names.sorted) + 1):(length(drug.names.sorted)*2), ], aes(x = order, y = mean, ymin = lower, ymax = upper, colour = interval)) +
      geom_linerange(size = 2, position = position_dodge(width = 0.5)) +
      geom_errorbar(data = prepare.EM[1:length(drug.names.sorted), ], aes(x = order, y = mean, ymin = lower, ymax = upper), size = 2, position = position_dodge(width = 0.5), width = 0.0) +
      geom_hline(yintercept = ifelse(!is.element(measure, c("Odds ratio", "Ratio of means")), 0, 1), lty = 1, size = 1, col = "grey60") +
      geom_point(size = 1.5,  colour = "white", stroke = 0.3, position = position_dodge(width = 0.5)) +
      geom_text(aes(x = order, y = mean, label = paste0(mean, " ", "(", prepare.EM[1:length(drug.names.sorted), 4], ",", " ", prepare.EM[1:length(drug.names.sorted), 5], ")",
                                                        " ", "[", prepare.EM[(length(drug.names.sorted) + 1):(length(drug.names.sorted)*2), 4], ",", " ", prepare.EM[(length(drug.names.sorted) + 1):(length(drug.names.sorted)*2), 5], "]"),
                    hjust = 0, vjust = -0.5), color = "black", size = 4.0, check_overlap = F, parse = F, position = position_dodge(width = 0.5), inherit.aes = T, na.rm = T) +
      geom_text(aes(x = 0.45, y = ifelse(is.element(measure, c("Odds ratio", "Ratio of means")), 0.2, -0.2), label = ifelse(full$D == 0, "Favours first arm", paste("Favours", compar))),
                size = 3.5, vjust = 0, hjust = 0, color = "black") +
      geom_text(aes(x = 0.45, y = ifelse(is.element(measure, c("Odds ratio", "Ratio of means")), 1.2, 0.2), label = ifelse(full$D == 0, paste("Favours", compar), "Favours first arm")),
                size = 3.5, vjust = 0, hjust = 0, color = "black") +
      labs(x = "", y = measure, colour = "Analysis") +
      scale_x_discrete(breaks = as.factor(1:length(drug.names.sorted)), labels = drug.names.sorted[length(drug.names.sorted):1]) +
      scale_color_manual(breaks = c("Credible interval", "Predictive interval"), values = c("black", "#D55E00")) +
      geom_label(aes(x = order[is.na(mean)], y = ifelse(!is.element(measure, c("Odds ratio", "Ratio of means")), -0.2, 0.65), hjust = 0, vjust = 1, label = "Comparator intervention"), fill = "beige", colour = "black", fontface = "plain", size = 4) +
      scale_y_continuous(trans = ifelse(!is.element(measure, c("Odds ratio", "Ratio of means")), "identity", "log10")) +
      coord_flip() +
      theme_classic() +
       theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12),
            axis.title.x = element_text(color = "black", face = "bold", size = 12), legend.position = "bottom", legend.justification = c(0.23, 0),
            legend.text =  element_text(color = "black", size = 12), legend.title = element_text(color = "black", face = "bold", size = 12))
  } else {
    ggplot(data = prepare.EM, aes(x = order, y = mean, ymin = lower, ymax = upper)) +
      geom_linerange(size = 2, position = position_dodge(width = 0.5)) +
      geom_errorbar(data = prepare.EM[1:length(drug.names.sorted), ], aes(x = order, y = mean, ymin = lower, ymax = upper), size = 2, position = position_dodge(width = 0.5), colour = "black", width = 0.0) +
      geom_hline(yintercept = ifelse(!is.element(measure, c("Odds ratio", "Ratio of means")), 0, 1), lty = 2, size = 1.3, col = "grey53") +
      geom_point(size = 1.5,  colour = "white", stroke = 0.3, position = position_dodge(width = 0.5)) +
      geom_text(aes(x = order, y = mean, label = paste0(mean, " ", "(", prepare.EM[1:length(drug.names.sorted), 4], ",", " ", prepare.EM[1:length(drug.names.sorted), 5], ")"),
                    hjust = 0, vjust = -0.5), color = "black", size = 4.0, check_overlap = F, parse = F, position = position_dodge(width = 0.5), inherit.aes = T, na.rm = T) +
      geom_text(aes(x = 0.45, y = ifelse(is.element(measure, c("Odds ratio", "Ratio of means")), 0.2, -0.2), label = ifelse(full$D == 0, "Favours first arm", paste("Favours", compar))),
                size = 3.5, vjust = 0, hjust = 0, color = "black") +
      geom_text(aes(x = 0.45, y = ifelse(is.element(measure, c("Odds ratio", "Ratio of means")), 1.2, 0.2), label = ifelse(full$D == 0, paste("Favours", compar), "Favours first arm")),
                size = 3.5, vjust = 0, hjust = 0, color = "black") +
      labs(x = "", y = measure) +
      scale_x_discrete(breaks = as.factor(1:length(drug.names.sorted)), labels = drug.names.sorted[length(drug.names.sorted):1]) +
      geom_label(aes(x = order[is.na(mean)], y = ifelse(!is.element(measure, c("Odds ratio", "Ratio of means")), -0.2, 0.65), hjust = 0, vjust = 1, label = "Comparator intervention"), fill = "beige", colour = "black", fontface = "plain", size = 4) +
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
                        hjust = 0, vjust = -0.5), color = "black", size = 4.0, check_overlap = F, parse = F, position = position_dodge(width = 0.5), inherit.aes = T) +
    labs(x = "", y = "Surface under the cumulative ranking curve value") +
    scale_x_discrete(breaks = as.factor(1:length(drug.names)), labels = prepare.sucra$intervention[length(drug.names):1]) +
    scale_y_continuous(labels = percent) +
    coord_flip() +
    theme_classic() +
    theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12),
                axis.title.x = element_text(color = "black", face = "bold", size = 12))



  ## Bring together both forest-plots
  forest.plots <- ggarrange(p1, p2, nrow = 1, ncol = 2, labels = c("A)", "B)"), common.legend = T, legend = "bottom")


  return(forest.plots = forest.plots)
}

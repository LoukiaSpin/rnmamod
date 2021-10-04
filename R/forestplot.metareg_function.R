#' Comparator-specific forest-plot on estimated and predicted effect sizes
#'
#' @description This function illustrates a forest plot of the posterior mean and 95\% credible and predictive interval of comparisons with the selected intervention of the network.
#'
#' @param full An object of S3 class \code{\link{run.model}}. See 'Value' in \code{\link{run.model}}.
#' @param reg An object of S3 class \code{\link{run.metareg}}. See 'Value' in \code{\link{run.metareg}}.
#' @param compar A character to indicate the comparator intervention. It must be any name found in \code{drug.names}.
#' @param cov.value A vector of two elements in the following order: a number that corresponds to a value of the covariate considered in \code{\link{run.metareg}},
#'   and a character object to indicate the name of the covariate. See also 'Details'.
#' @param drug.names A vector of labels with the name of the interventions in the order they appear in the argument \code{data} of \code{\link{run.model}}. If the argument \code{drug.names} is not defined, the order of the interventions
#'   as they appear in \code{data} is used, instead.
#'
#' @return A panel of two forest plots: (1) a forest plot on the estimated effect size of comparisons with the selected intervention of the network, and
#'   (2) a forest plot predicted effect size of comparisons with the selected intervention of the network.
#'   Both panels illustrate the results from network meta-analysis and meta-regression using different colours.
#'
#' @details The y-axis displays all interventions in the network; the selected intervention that comprises the \code{compar} is indicated in the plot with a homonymous label.
#'   The numerical results are displayed above each line.
#'   Odds ratio and ratio of means are reported in the original scale after exponentiation of the logarithmic scale.
#'
#'   When the covariate is binary, specify in the second element of \code{cov.value} the name of the level for which the forest plot will be created.
#'
#'   The interventions are sorted in the descending order of their SUCRA values obtain via network meta-analysis.
#'
#'   \code{forestplot} can be used only for a network of interventions. In the case of two interventions, the execution of the function will be stopped and an error message will be printed in the R console.
#'
#' @author {Loukia M. Spineli}
#'
#' @seealso \code{\link{run.model}}, \code{\link{run.metareg}}
#'
#' @references
#' Salanti G, Ades AE, Ioannidis JP. Graphical methods and numerical summaries for presenting results from multiple-treatment meta-analysis: an overview and tutorial. \emph{J Clin Epidemiol} 2011;\bold{64}(2):163--71. [\doi{10.1016/j.jclinepi.2010.03.016}]
#'
#' @export
forestplot.metareg <- function(full, reg, compar, cov.value, drug.names) {

  options(warn = -1)

  if (length(unique(reg$covariate)) < 3 & !is.element(cov.value[1], reg$covariate)) {
    stop("The first element of the argument 'cov.value' is out of the value range of the analysed covariate", call. = F)
  } else if (length(unique(reg$covariate)) > 2 & (cov.value[1] < min(reg$covariate) | cov.value[1] > max(reg$covariate))) {
    stop("The first element of the argument 'cov.value' is out of the value range of the analysed covariate", call. = F)
  }

  drug.names <- if (missing(drug.names)) {
    message(cat(paste0("\033[0;", col = 32, "m", txt = "The argument 'drug.names' has not been defined. The intervention ID, as specified in 'data' is used as intervention names", "\033[0m", "\n")))
    nt <- length(full$SUCRA[, 1])
    as.character(1:nt)
  } else {
    drug.names
  }

  if(length(drug.names) < 3) {
    stop("This function is *not* relevant for a pairwise meta-analysis", call. = F)
  }

  # Sort the drugs by their SUCRA in decreasing order and remove the reference intervention (number 1)
  drug.names.sorted <- drug.names[order(full$SUCRA[, 1], decreasing = T)]

  compar <- if(missing(compar)) {
    stop("The argument 'compar' has not been defined", call. = F)
  } else if(!is.element(compar, drug.names)) {
    stop("The value of the argument 'compar' is not found in the 'drug.names'", call. = F)
  } else if(is.element(compar, drug.names)) {
    compar
  }

  cov.value <- if (!is.null(reg$beta.all) & missing(cov.value)) {
    stop("The argument 'cov.value' has not been defined", call. = F)
  } else if (!is.null(reg$beta.all) & length(cov.value) < 2) {
    stop("The argument 'cov.value' must be a vector with elements a number and a character", call. = F)
  } else if (!is.null(reg$beta.all) & length(cov.value) == 2) {
    cov.value
  }

  measure <- effect.measure.name(full$measure)
  model <- full$model
  cov.val <- ifelse(length(unique(reg$covariate)) < 3, as.numeric(cov.value[1]), as.numeric(cov.value[1]) - mean(reg$covariate))


  ## A matrix with all possible comparisons in the network
  poss.pair.comp1 <- data.frame(exp = t(combn(drug.names, 2))[, 2], comp = t(combn(drug.names, 2))[, 1])
  poss.pair.comp2 <- data.frame(exp = t(combn(drug.names, 2))[, 1], comp = t(combn(drug.names, 2))[, 2])
  poss.pair.comp <- rbind(poss.pair.comp1, poss.pair.comp2)


  ## Effect size of all possible pairwise comparisons (NMA)
  EM.ref00.nma <- cbind(rbind(data.frame(mean = full$EM[, 1], lower = full$EM[, 3], upper = full$EM[, 7]),
                              data.frame(mean = full$EM[, 1]*(-1), lower = full$EM[, 7]*(-1), upper = full$EM[, 3]*(-1))),
                        poss.pair.comp)
  EM.subset.nma <- subset(EM.ref00.nma, EM.ref00.nma[5] == compar)
  EM.ref0.nma <- rbind(EM.subset.nma[, 1:3], c(rep(NA, 3)))
  sucra.new <- data.frame(full$SUCRA[, 1], drug.names)[order(match(data.frame(full$SUCRA[, 1], drug.names)[, 2], EM.subset.nma[, 4])), 1]
  EM.ref.nma <- EM.ref0.nma[order(sucra.new, decreasing = T), ]

  ## Effect size of all possible pairwise comparisons (NMR)
 #if (is.element(reg$covar.assumption, c("exchangeable", "independent"))) {
 # par.mean <- as.vector(c(reg$EM[, 1] + reg$beta.all[, 1]*cov.val,
 #                         (reg$EM[, 1]*(-1)) + (reg$beta.all[, 1]*(-1)*cov.val)))
  # par.sd <- as.vector(c(sqrt(((reg$EM[, 2])^2) + ((reg$beta.all[, 2]*cov.val)^2)),
  #                       sqrt(((reg$EM[, 2])^2) + ((reg$beta.all[, 2]*cov.val)^2))))
 #} else {
#   par.mean <- as.vector(c(reg$EM[, 1], reg$EM[, 1]*(-1)))
#   par.sd <- as.vector(c(reg$EM[, 2], reg$EM[, 2]))
#
#   # Correcting for comparisons with the reference intervention of the network
#   par.mean[1:(length(drug.names) - 1)] <- as.vector(c(reg$EM[1:(length(drug.names) - 1), 1] + reg$beta[1]*cov.val,
#                                                      (reg$EM[1:(length(drug.names) - 1), 1]*(-1)) + (reg$beta[1]*(-1)*cov.val)))
#   par.sd[1:(length(drug.names) - 1)] <- as.vector(c(sqrt(((reg$EM[1:(length(drug.names) - 1), 2])^2) + ((reg$beta[2]*cov.val)^2)),
#                                                     sqrt(((reg$EM[1:(length(drug.names) - 1), 2])^2) + ((reg$beta[2]*cov.val)^2))))
# }
  par.mean <- as.vector(c(reg$EM[, 1] + reg$beta.all[, 1]*cov.val,
                          (reg$EM[, 1]*(-1)) + (reg$beta.all[, 1]*(-1)*cov.val)))
  par.sd <- as.vector(c(sqrt(((reg$EM[, 2])^2) + ((reg$beta.all[, 2]*cov.val)^2)),
                        sqrt(((reg$EM[, 2])^2) + ((reg$beta.all[, 2]*cov.val)^2))))

  EM.ref00.nmr <- cbind(mean = par.mean, lower = par.mean - 1.96*par.sd, upper = par.mean + 1.96*par.sd,
                        poss.pair.comp)
  EM.subset.nmr <- subset(EM.ref00.nmr, EM.ref00.nmr[5] == compar)
  EM.ref0.nmr <- rbind(EM.subset.nmr[, 1:3], c(rep(NA, 3)))
  EM.ref.nmr <- EM.ref0.nmr[order(sucra.new, decreasing = T), ]
  rownames(EM.ref.nma) <- rownames(EM.ref.nmr) <- NULL


  # Posterior results on the predicted estimates of comparisons with the selected comparator as reference
  if (model == "RE") {
    pred.ref00.nma <- cbind(rbind(data.frame(mean = full$EM.pred[, 1], lower = full$EM.pred[, 3], upper = full$EM.pred[, 7]),
                                  data.frame(mean = full$EM.pred[, 1]*(-1), lower = full$EM.pred[, 7]*(-1), upper = full$EM.pred[, 3]*(-1))),
                            poss.pair.comp)
    pred.subset.nma <- subset(pred.ref00.nma, pred.ref00.nma[5] == compar)
    pred.ref0.nma <- rbind(pred.subset.nma[, 1:3], c(rep(NA, 3)))

    #if (is.element(reg$covar.assumption, c("exchangeable", "independent"))) {
    #  par.mean <- as.vector(c(reg$EM.pred[, 1] + reg$beta.all[, 1]*cov.val,
    #                          (reg$EM.pred[, 1]*(-1)) + (reg$beta.all[, 1]*(-1)*cov.val)))
    #  par.sd <- as.vector(c(sqrt(((reg$EM.pred[, 2])^2) + ((reg$beta.all[, 2]*cov.val)^2)),
    #                        sqrt(((reg$EM.pred[, 2])^2) + ((reg$beta.all[, 2]*cov.val)^2))))
    #} else {
    #  par.mean <- as.vector(c(reg$EM.pred[, 1], reg$EM.pred[, 1]*(-1)))
    #  par.sd <- as.vector(c(reg$EM.pred[, 2], reg$EM.pred[, 2]))
    #  # Correcting for comparisons with the reference intervention of the network
    #  par.mean[1:(length(drug.names) - 1)] <- as.vector(c(reg$EM[1:(length(drug.names) - 1), 1] + reg$beta[1]*cov.val,
    #                                                      (reg$EM[1:(length(drug.names) - 1), 1]*(-1)) + (reg$beta[1]*(-1)*cov.val)))
    #  par.sd[1:(length(drug.names) - 1)] <- as.vector(c(sqrt(((reg$EM[1:(length(drug.names) - 1), 2])^2) + ((reg$beta[2]*cov.val)^2)),
    #                                                    sqrt(((reg$EM[1:(length(drug.names) - 1), 2])^2) + ((reg$beta[2]*cov.val)^2))))
    #}

    par.mean <- as.vector(c(reg$EM.pred[, 1] + reg$beta.all[, 1]*cov.val,
                            (reg$EM.pred[, 1]*(-1)) + (reg$beta.all[, 1]*(-1)*cov.val)))
    par.sd <- as.vector(c(sqrt(((reg$EM.pred[, 2])^2) + ((reg$beta.all[, 2]*cov.val)^2)),
                          sqrt(((reg$EM.pred[, 2])^2) + ((reg$beta.all[, 2]*cov.val)^2))))

    pred.ref00.nmr <-  cbind(data.frame(mean = par.mean, lower = par.mean - 1.96*par.sd, upper = par.mean + 1.96*par.sd),
                           poss.pair.comp)
    pred.subset.nmr <- subset(pred.ref00.nmr, pred.ref00.nmr[5] == compar)
    pred.ref0.nmr <- rbind(pred.subset.nmr[, 1:3], c(rep(NA, 3)))
  } else if (model != "RE") {
    pred.ref00.nma <- NA
    pred.ref00.nmr <- NA
  }


  # Sort by SUCRA in decreasing order and remove the reference intervention (number 1)
  if (model == "RE") {
    pred.ref.nma <- pred.ref0.nma[order(sucra.new, decreasing = T), ]
    pred.ref.nmr <- pred.ref0.nmr[order(sucra.new, decreasing = T), ]
  } else {
    NA
  }
  rownames(pred.ref.nma) <- rownames(pred.ref.nmr) <- NULL


  ## Create a data-frame with credible and predictive intervals of comparisons with the reference intervention
  if (!is.element(measure, c("Odds ratio", "Ratio of means")) & model == "RE") {
    prepare.EM.nma <- data.frame(as.factor(rep(length(drug.names):1, 2)),
                                 rep(drug.names.sorted, 2),
                                 round(rbind(EM.ref.nma, pred.ref.nma), 2),
                                 rep(c("Estimation", "Prediction"), each = length(drug.names)))
    prepare.EM.nmr <- data.frame(as.factor(rep(length(drug.names):1, 2)),
                                 rep(drug.names.sorted, 2),
                                 round(rbind(EM.ref.nmr, pred.ref.nmr), 2),
                                 rep(c("Estimation", "Prediction"), each = length(drug.names)))
    colnames(prepare.EM.nma) <- colnames(prepare.EM.nmr) <- c("order", "comparison", "mean", "lower", "upper", "interval")
  } else if (is.element(measure, c("Odds ratio", "Ratio of means")) & model == "RE") {
    prepare.EM.nma <- data.frame(as.factor(rep(length(drug.names):1, 2)),
                                 rep(drug.names.sorted, 2),
                                 round(rbind(exp(EM.ref.nma), exp(pred.ref.nma)), 2),
                                 rep(c("Estimation", "Prediction"), each = length(drug.names)))
    prepare.EM.nmr <- data.frame(as.factor(rep(length(drug.names):1, 2)),
                                 rep(drug.names.sorted, 2),
                                 round(rbind(exp(EM.ref.nmr), exp(pred.ref.nmr)), 2),
                                 rep(c("Estimation", "Prediction"), each = length(drug.names)))
    colnames(prepare.EM.nma) <- colnames(prepare.EM.nmr) <- c("order", "comparison", "mean", "lower", "upper", "interval")
  } else if (!is.element(measure, c("Odds ratio", "Ratio of means")) & model == "FE") {
    prepare.EM.nma <- data.frame(as.factor(length(drug.names):1),
                                 drug.names.sorted,
                                 round(EM.ref.nma, 2))
    prepare.EM.nmr <- data.frame(as.factor(length(drug.names):1),
                                 drug.names.sorted,
                                 round(EM.ref.nmr, 2))
    colnames(prepare.EM.nma) <- colnames(prepare.EM.nmr) <- c("order", "comparison", "mean", "lower", "upper")
  } else if (is.element(measure, c("Odds ratio", "Ratio of means")) & model == "RE") {
    prepare.EM.nma <- data.frame(as.factor(length(drug.names):1),
                                 drug.names.sorted,
                                 round(exp(EM.ref.nma), 2))
    prepare.EM.nmr <- data.frame(as.factor(length(drug.names):1),
                                 drug.names.sorted,
                                 round(exp(EM.ref.nmr), 2))
    colnames(prepare.EM.nma) <- colnames(prepare.EM.nmr) <- c("order", "comparison", "mean", "lower", "upper")
  }

  prepare.EM <- cbind(rbind(prepare.EM.nma, prepare.EM.nmr), analysis = rep(c("Network meta-analysis", "Network meta-regression"), each = length(prepare.EM.nma[, 1])))


  ## Forest plots on credible/predictive intervals of comparisons with the selected comparator
  forest.plots <- if (model == "RE") {
    ggplot(data = prepare.EM, aes(x = order, y = mean, ymin = lower, ymax = upper, group = analysis, colour = analysis)) +
      geom_hline(yintercept = ifelse(!is.element(measure, c("Odds ratio", "Ratio of means")), 0, 1), lty = 1, size = 1, col = "grey60") +
      geom_linerange(size = 2, position = position_dodge(width = 0.5)) +
      geom_errorbar(data = prepare.EM, aes(x = order, y = mean, ymin = lower, ymax = upper), size = 2, position = position_dodge(width = 0.5), width = 0.0) +
      geom_point(size = 1.5, colour = "white", stroke = 0.3, position = position_dodge(width = 0.5)) +
      geom_text(aes(x = order, y = mean, label = paste0(mean, " ", "(", prepare.EM[, 4], ",", " ", prepare.EM[, 5], ")"),
                    hjust = 0, vjust = -0.5), color = "black", size = 4.0, check_overlap = F, parse = F, position = position_dodge(width = 0.5), inherit.aes = T, na.rm = T) +
      geom_text(aes(x = 0.45, y = ifelse(is.element(measure, c("Odds ratio", "Ratio of means")), 0.2, -0.2), label = ifelse(full$D == 0, "Favours first arm", paste("Favours", compar))),
                size = 3.5, vjust = 0, hjust = 0, color = "black") +
      geom_text(aes(x = 0.45, y = ifelse(is.element(measure, c("Odds ratio", "Ratio of means")), 1.2, 0.2), label = ifelse(full$D == 0, paste("Favours", compar), "Favours first arm")),
                size = 3.5, vjust = 0, hjust = 0, color = "black") +
      labs(x = "", y = measure, colour = "Analysis", subtitle = paste("Results for", cov.value[2], ifelse(length(unique(reg$covariate)) < 3, " ", cov.value[1]))) +
      facet_wrap(~ interval, ncol = 2, scales = "fixed") +
      scale_x_discrete(breaks = as.factor(1:length(drug.names.sorted)), labels = drug.names.sorted[length(drug.names.sorted):1]) +
      scale_color_manual(breaks = c("Network meta-analysis", "Network meta-regression"), values = c("black", "#D55E00")) +
      geom_label(aes(x = unique(order[is.na(mean)]), y = ifelse(!is.element(measure, c("Odds ratio", "Ratio of means")), -0.2, 0.65), hjust = 0, vjust = 1, label = "Comparator intervention"), fill = "beige", colour = "black", fontface = "plain", size = 4) +
      scale_y_continuous(trans = ifelse(!is.element(measure, c("Odds ratio", "Ratio of means")), "identity", "log10")) +
      coord_flip() +
      theme_classic() +
      theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12),
            axis.title.x = element_text(color = "black", face = "bold", size = 12), legend.position = "bottom", legend.justification = c(0.23, 0),
            strip.text = element_text(color = "black", face = "bold", size = 12),
            legend.text =  element_text(color = "black", size = 12), legend.title = element_text(color = "black", face = "bold", size = 12))
  } else {
    ggplot(data = prepare.EM, aes(x = order, y = mean, ymin = lower, ymax = upper, group = analysis, colour = analysis)) +
      geom_hline(yintercept = ifelse(!is.element(measure, c("Odds ratio", "Ratio of means")), 0, 1), lty = 1, size = 1, col = "grey60") +
      geom_linerange(size = 2, position = position_dodge(width = 0.5)) +
      geom_errorbar(data = prepare.EM, aes(x = order, y = mean, ymin = lower, ymax = upper), size = 2, position = position_dodge(width = 0.5), width = 0.0) +
      geom_point(size = 1.5, colour = "white", stroke = 0.3, position = position_dodge(width = 0.5)) +
      geom_text(aes(x = order, y = mean, label = paste0(mean, " ", "(", prepare.EM[, 4], ",", " ", prepare.EM[, 5], ")"),
                    hjust = 0, vjust = -0.5), color = "black", size = 4.0, check_overlap = F, parse = F, position = position_dodge(width = 0.5), inherit.aes = T, na.rm = T) +
      geom_text(aes(x = 0.45, y = ifelse(is.element(measure, c("Odds ratio", "Ratio of means")), 0.2, -0.2), label = ifelse(full$D == 0, "Favours first arm", paste("Favours", compar))),
                size = 3.5, vjust = 0, hjust = 0, color = "black") +
      geom_text(aes(x = 0.45, y = ifelse(is.element(measure, c("Odds ratio", "Ratio of means")), 1.2, 0.2), label = ifelse(full$D == 0, paste("Favours", compar), "Favours first arm")),
                size = 3.5, vjust = 0, hjust = 0, color = "black") +
      labs(x = "", y = measure, colour = "Analysis", subtitle = paste("Results for", cov.value[2], ifelse(length(unique(reg$covariate)) < 3, " ", cov.value[1]))) +
      scale_x_discrete(breaks = as.factor(1:length(drug.names.sorted)), labels = drug.names.sorted[length(drug.names.sorted):1]) +
      scale_color_manual(breaks = c("Network meta-analysis", "Network meta-regression"), values = c("black", "#D55E00")) +
      geom_label(aes(x = unique(order[is.na(mean)]), y = ifelse(!is.element(measure, c("Odds ratio", "Ratio of means")), -0.2, 0.65), hjust = 0, vjust = 1, label = "Comparator intervention"), fill = "beige", colour = "black", fontface = "plain", size = 4) +
      scale_y_continuous(trans = ifelse(!is.element(measure, c("Odds ratio", "Ratio of means")), "identity", "log10")) +
      coord_flip() +
      theme_classic() +
      theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12),
            strip.text = element_text(color = "black", face = "bold", size = 12),
            axis.title.x = element_text(color = "black", face = "bold", size = 12), legend.position = "bottom", legend.justification = c(0.23, 0),
            legend.text =  element_text(color = "black", size = 12), legend.title = element_text(color = "black", face = "bold", size = 12))
  }

  return(Forest.plots = forest.plots)
}

#' Plot the results from the meta-regression analysis
#'
#' @param full An object of S3 class \code{\link{run.model}}. See 'Value' in \code{\link{run.model}}.
#' @param metareg An object of S3 class  \code{\link{run.metareg}}. See 'Value' in \code{\link{run.metareg}}.
#' @param compar A character to indicate the comparator intervention. it must be any name found in \code{drug.names}.
#' @param cov.value A vector of two elements in the following order: a number that corresponds to a value of the covariate considered in \code{\link{run.metareg}},
#'   and a character object to indicate the name of the covariate.
#' @param drug.names A vector of labels with the name of the interventions in the order they appear in the argument \code{data} of \code{\link{run.model}}. If the argument \code{drug.names} is not defined, the order of the interventions
#'   as they appear in \code{data} is used, instead.
#' @param save.xls Logical to indicate whether to export the tabulated results to an Excel 'xlsx' format (via the \code{\link[writexl]{write_xlsx}} function) to the working directory of the user.
#'   The default is \code{FALSE} (do not export to an Excel format).
#'
#'
#' @export
metareg.plot <- function(full, reg, compar, cov.value, drug.names, save.xls) {


  options(warn = -1)

  save.xls <- if (missing(save.xls)) {
    FALSE
  } else {
    save.xls
  }

  compar <- if(missing(compar)) {
    stop("The argument 'compar' has not been defined", call. = F)
  } else if(!is.element(compar, drug.names)) {
    stop("The value of 'compar' is not found in the 'drug.names'", call. = F)
  } else if(is.element(compar, drug.names)) {
    compar
  }

  cov.value <- if (missing(cov.value)) {
    stop("The argument 'cov.value' has not been defined", call. = F)
  } else if (length(cov.value) < 2) {
    stop("The argument 'cov.value' must be a vector with elements a number and a character", call. = F)
  } else if (length(cov.value) == 2) {
    cov.value
  }

  covariate <- if (length(unique(reg$covariate)) < 3) {
    unique(reg$covariate)
  } else {
    reg$covariate
  }

  cov.val <- ifelse(length(unique(covariate)) < 3, as.numeric(cov.value[1]), as.numeric(cov.value[1]) - mean(covariate))

  drug.names <- if (missing(drug.names)) {
    message(cat(paste0("\033[0;", col = 32, "m", txt = "The argument 'drug.names' has not been defined. The intervention ID, as specified in 'data' is used as intervention names", "\033[0m", "\n")))
    nt <- as.factor(1:length(full$SUCRA[, 1]))
  } else {
    drug.names
  }

  model <- full$model
  measure <- full$measure

  # Posterior results on the SUCRA value under NMA
  sucra.full <- round(full$SUCRA, 2)
  sucra.full.order <- round(full$SUCRA, 2)[order(sucra.full[, 1], decreasing = T), ]

  # Sort the drugs by their NMA-SUCRA in decreasing order
  drug.names.sorted <- drug.names[order(sucra.full[, 1], decreasing = T)]

  # Posterior results on the SUCRA value under meta-regression
  sucra.meta <- round(reg$SUCRA, 2)
  sucra.meta.order <- round(reg$SUCRA, 2)[order(sucra.full[, 1], decreasing = T), ]

  ## A matrix with all possible comparisons in the network
  poss.pair.comp1 <- data.frame(exp = t(combn(drug.names, 2))[, 2], comp = t(combn(drug.names, 2))[, 1])
  poss.pair.comp2 <- data.frame(exp = t(combn(drug.names, 2))[, 1], comp = t(combn(drug.names, 2))[, 2])
  poss.pair.comp <- rbind(poss.pair.comp1, poss.pair.comp2)

  ## Posterior results on effect size for comparisons with the selected intervention (NMA)
  EM.full00 <- cbind(rbind(data.frame(mean = full$EM[, 1], lower = full$EM[, 3], upper = full$EM[, 7]),
                           data.frame(mean = full$EM[, 1]*(-1), lower = full$EM[, 7]*(-1), upper = full$EM[, 3]*(-1))),
                     poss.pair.comp)

  EM.full.subset <- subset(EM.full00, EM.full00[5] == compar)

  EM.full0 <- rbind(EM.full.subset[, 1:3], c(rep(NA, 3)))

  sucra.full.new <- data.frame(sucra.full[, 1], drug.names)[order(match(data.frame(sucra.full[, 1], drug.names)[, 2], EM.full.subset[, 4])), 1]

  EM.full <- EM.full0[order(sucra.full.new, decreasing = T), ]
  rownames(EM.full) <- NULL

  # Posterior mean of regression coefficients for all unique pairwise comparisons
  beta00 <- cbind(rbind(data.frame(mean = reg$beta.all[, 5], lower = reg$beta.all[, 3], upper = reg$beta.all[, 7]),
                        data.frame(mean = reg$beta.all[, 5]*(-1), lower = reg$beta.all[, 7]*(-1), upper = reg$beta.all[, 3]*(-1))),
                  poss.pair.comp)

  beta.all.subset <- subset(beta00, beta00[5] == compar)

  beta0 <- rbind(beta.all.subset[, 1:3], c(rep(NA, 3)))

  beta <- beta0[order(sucra.full.new, decreasing = T), ]
  rownames(beta) <- NULL

  ## Posterior results on effect size for comparisons with the selected intervention (network meta-regression)
  EM.meta <- rbind(data.frame(mean = EM.full[, 1] + beta[, 1]*cov.val,
                              lower = EM.full[, 2] + beta[, 2]*cov.val,
                              upper = EM.full[, 3] + beta[, 3]*cov.val))

  ## Posterior results on the predicted estimates of comparisons with the selected comparator as reference
  if (model == "RE") {
    pred.ref00.nma <- cbind(rbind(data.frame(mean = full$EM.pred[, 1], lower = full$EM.pred[, 3], upper = full$EM.pred[, 7]),
                                  data.frame(mean = full$EM.pred[, 1]*(-1), lower = full$EM.pred[, 7]*(-1), upper = full$EM.pred[, 3]*(-1))),
                            poss.pair.comp)

    pred.subset.nma <- subset(pred.ref00.nma, pred.ref00.nma[5] == compar)

    pred.ref0.nma <- rbind(pred.subset.nma[, 1:3], c(rep(NA, 3)))
  } else if (model != "RE") {
    pred.ref00.nma <- NA
  }


  # Sort by SUCRA in decreasing order and remove the reference intervention (number 1)
  if (model == "RE") {
    pred.ref.nma <- pred.ref0.nma[order(sucra.full.new, decreasing = T), ]
    pred.ref.nmr <- pred.ref.nma + beta*cov.val
  } else {
    NA
  }
  rownames(pred.ref.nma) <- rownames(pred.ref.nmr) <- NULL


  # Posterior results on between-trial standard deviation under NMA
  tau.full <- if (model == "RE") {
    round(full$tau, 2)
  } else {
    NA
  }

  # Posterior results on between-trial standard deviation under meta-regression
  tau.meta <- if (model == "RE") {
    round(metareg$tau, 2)
  } else {
    NA
  }

  # Posterior mean of model assessment measures under NMA
  model.assess.NMA <- round(full$model.assessment, 2)

  # Posterior mean of model assessment measures under meta-regression
  model.assess.meta <- round(metareg$model.assessment, 2)

  # Effect measure name
  measure <- effect.measure.name(full$measure)


  ## The 95% CrIs of the between-trial standard deviation under NMA and meta-regression
  if (model == "RE") {
    CrI.tau.full <- paste0("(", tau.full[, 3], ",", " ", tau.full[, 7], ")")
    CrI.tau.meta <- paste0("(", tau.meta[, 3], ",", " ", tau.meta[, 7], ")")
  } else {
    CrI.tau.full <- NA
    CrI.tau.meta <- NA
  }


  ## A data-frame on the model assessment results and between-trial standard deviation under NMA & meta-regression
  if (model == "RE") {
    table.model.assess <- data.frame(c("Network meta-analysis", "Meta-regression"),
                                     rbind(model.assess.NMA[c(1, 3, 2)],
                                           model.assess.meta[c(1, 3, 2)]),
                                     rbind(cbind(tau.full[, 1], tau.full[, 2], CrI.tau.full),
                                           cbind(tau.meta[, 1], tau.meta[, 2], CrI.tau.meta)))
    colnames(table.model.assess) <- c("Analysis", "DIC", "Mean deviance", "pD", "Median tau", "SD tau", "95% CrI tau")
  } else {
    table.model.assess <- data.frame(c("Network meta-analysis", "Meta-regression"),
                                     rbind(model.assess.NMA[c(1, 3, 2)],
                                           model.assess.meta[c(1, 3, 2)]))
    colnames(table.model.assess) <- c("Analysis", "DIC", "Mean deviance", "pD")
  }


  ## A data-frame with the effect estimates and regression coefficients of reference-comparisons from both analyses (Sort by NMA-SUCRA in decreasing order)
  if (!is.element(measure, c("Odds ratio", "Ratio of means"))) {
    CrI.full <- paste0("(", round(EM.full[, 2], 2), ",", " ", round(EM.full[, 3], 2), ")", ifelse(EM.full[, 2] > 0 | EM.full[, 3] < 0, "*", " "))
    CrI.meta <- paste0("(", round(EM.meta[, 2], 2), ",", " ", round(EM.meta[, 3], 2), ")", ifelse(EM.meta[, 2] > 0 | EM.meta[, 3] < 0, "*", " "))
    CrI.beta <- if (is.element(metareg$covar.assumption, c("exchangeable", "independent"))) {
      paste0("(", round(beta[, 2], 2), ",", " ", round(beta[, 3], 2), ")", ifelse(beta[, 2] > 0 | beta[, 3] < 0, "*", " "))
    } else {
      paste0("(", round(metareg$beta[3], 2), ",", " ", round(metareg$beta[7], 2), ")", ifelse(metareg$beta[3] > 0 | metareg$beta[7] < 0, "*", " "))
    }

    # Tabulate results on comparisons with the reference (both models)
    EM.both.models <- na.omit(data.frame(drug.names.sorted,
                                         round(EM.full[, 1], 2),
                                         CrI.full,
                                         round(EM.meta[, 1], 2),
                                         CrI.meta))
    colnames(EM.both.models) <- c(paste("versus", compar), "Mean NMA", "95% CrI NMA", "Mean NMR", "95% CrI NMR")
    rownames(EM.both.models) <- NULL

    # Results on the regression coefficient - HERE!
    if (is.element(metareg$covar.assumption, c("exchangeable", "independent")) ) {
      reg.coeff <- na.omit(data.frame(drug.names.sorted, beta[, 1], CrI.beta))
      colnames(reg.coeff) <- c(paste("versus", compar), "Mean beta", "95% CrI beta")
    } else {
      reg.coeff <- data.frame(metareg$beta[1], CrI.beta)
      colnames(reg.coeff) <- c("Mean beta", "95% CrI beta")
    }
    rownames(reg.coeff) <- NULL
  } else if (is.element(measure, c("Odds ratio", "Ratio of means"))) {
    CrI.full <- paste0("(", round(exp(EM.full[, 2]), 2), ",", " ", round(exp(EM.full[, 3]), 2), ")", ifelse(EM.full[, 2] > 0 | EM.full[, 3] < 0, "*", " "))
    CrI.meta <- paste0("(", round(exp(EM.meta[, 2]), 2), ",", " ", round(exp(EM.meta[, 3]), 2), ")", ifelse(EM.meta[, 2] > 0 | EM.meta[, 3] < 0, "*", " "))
    CrI.beta <- if (is.element(metareg$covar.assumption, c("exchangeable", "independent")) ) {
      paste0("(", round(exp(beta[, 2]), 2), ",", " ", round(exp(beta[, 3]), 2), ")", ifelse(beta[, 2] > 0 | beta[, 3] < 0, "*", " "))
    } else {
      paste0("(", round(exp(metareg$beta[3]), 2), ",", " ", round(exp(metareg$beta[7]), 2), ")", ifelse(metareg$beta[3] > 0 | metareg$beta[7] < 0, "*", " "))
    }

    # Tabulate results on comparisons with the reference (both models)
    EM.both.models <- na.omit(data.frame(drug.names.sorted,
                                         round(exp(EM.full[, 1]), 2),
                                         CrI.full,
                                         round(exp(EM.meta[, 1]), 2),
                                         CrI.meta))
    colnames(EM.both.models) <- c(paste("versus", compar), "Mean NMA", "95% CrI NMA", "Mean NMR", "95% CrI NMR")
    rownames(EM.both.models) <- NULL

    # Results on the regression coefficient
    if (is.element(metareg$covar.assumption, c("exchangeable", "independent")) ) {
      reg.coeff <- na.omit(data.frame(drug.names.sorted, round(exp(beta[, 1]), 2), CrI.beta))
      colnames(reg.coeff) <- c(paste("versus", compar), "Mean beta", "95% CrI beta")
    } else {
      reg.coeff <- data.frame(round(exp(metareg$beta[1]), 2), CrI.beta)
      colnames(reg.coeff) <- c("Mean beta", "95% CrI beta")
    }
    rownames(reg.coeff) <- NULL
  }


  ## Create a data-frame with the SUCRA values
  prepare.sucra <- data.frame(rep(length(drug.names):1, 2),
                              rep(drug.names.sorted, 2),
                              round(rbind(sucra.full.order[, c(1, 3, 7)],
                                          sucra.meta.order[, c(1, 3, 7)]), 2),
                              rep(c("Network meta-analysis", "Network meta-regression"), each = length(drug.names)))
  colnames(prepare.sucra) <- c("order", "intervention", "mean", "lower", "upper", "analysis")
  rownames(prepare.sucra) <- NULL


  ## Prepare data for ggplot2 (forest-plot)
  if(!is.element(measure, c("Odds ratio", "Ratio of means"))) {
    prepare.EM <- data.frame(rep(length(drug.names):1, 2),
                             rep(drug.names.sorted, 2),
                             round(rbind(EM.full,
                                         EM.meta), 2),
                             rep(c("Network meta-analysis", "Network meta-regression"), each = length(drug.names)))
  } else if(is.element(measure, c("Odds ratio", "Ratio of means"))) {
    prepare.EM <- data.frame(rep(length(drug.names):1, 2),
                             rep(drug.names.sorted, 2),
                             round(rbind(exp(EM.full),
                                         exp(EM.meta)), 2),
                             rep(c("Network meta-analysis", "Network meta-regression"), each = length(drug.names)))
  }
  colnames(prepare.EM) <- c("order", "comparison", "mean", "lower", "upper", "analysis")
  rownames(prepare.EM) <- NULL


  ## Forest plots of reference-comparisons on effect estimate
  p1 <- ggplot(data = prepare.EM, aes(x = as.factor(order), y = as.numeric(mean), ymin = as.numeric(lower), ymax = as.numeric(upper), colour = analysis, group = analysis)) +
          geom_linerange(size = 2, position = position_dodge(width = 0.5), width = 0.1) +
          geom_hline(yintercept = ifelse(!is.element(measure, c("Odds ratio", "Ratio of means")), 0, 1), lty = 2, size = 1.3, col = "grey53") +
          geom_point(size = 1.5,  colour = "white", stroke = 0.3, position = position_dodge(width = 0.5)) +
          geom_text(aes(x = as.factor(order), y = round(as.numeric(mean), 2), label = paste0(round(as.numeric(mean), 2), " ", "(", round(as.numeric(lower), 2), ",", " ", round(as.numeric(upper), 2), ")"),
                    hjust = 0, vjust = -0.5), color = "black", size = 4.0, check_overlap = F, parse = F, position = position_dodge(width = 0.5), inherit.aes = T) +
          labs(x = "", y = measure, colour = "Analysis") +
          scale_x_discrete(breaks = as.factor(1:length(drug.names.sorted)), labels = drug.names.sorted[length(drug.names.sorted):1]) +
          scale_color_manual(breaks = c("Network meta-analysis", "Network meta-regression"), values = c("#009E73", "black")) +
          geom_label(aes(x = unique(as.factor(order)[is.na(mean)]), y = ifelse(!is.element(measure, c("Odds ratio", "Ratio of means")), 0, 1), hjust = 0, vjust = 1, label = "Comparator intervention"),
                     fill = "beige", colour = "black", fontface = "plain", size = 4) +
          scale_y_continuous(trans = ifelse(!is.element(measure, c("Odds ratio", "Ratio of means")), "identity", "log10")) +
          coord_flip() +
          theme_classic() +
          theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12),
                axis.title.x = element_text(color = "black", face = "bold", size = 12), legend.position = "none",
                legend.text =  element_text(color = "black", size = 12), legend.title =  element_text(color = "black", face = "bold", size = 12))


  ## Forest plots of SUCRA under NMA and meta-regression
  p2 <- ggplot(data = prepare.sucra, aes(x = as.factor(order), y = as.numeric(mean), ymin = as.numeric(lower), ymax = as.numeric(upper), colour = analysis, group = analysis)) +
           geom_linerange(size = 2, position = position_dodge(width = 0.5), width = 0.1) +
           geom_point(size = 1.5,  colour = "white", stroke = 0.3, position = position_dodge(width = 0.5)) +
           geom_text(aes(x = as.factor(order), y = round(as.numeric(mean), 2), label = paste0(round(mean*100, 0), " ", "(", round(lower*100, 0), ",", " ", round(upper*100, 0), ")" ),
                     hjust = 0, vjust = -0.5), color = "black", size = 4.0, check_overlap = F, parse = F, position = position_dodge(width = 0.5), inherit.aes = T) +
           labs(x = "", y = "Surface under the cumulative ranking curve value", colour = "Analysis") +
           scale_x_discrete(breaks = as.factor(1:length(drug.names.sorted)), labels = drug.names.sorted[length(drug.names.sorted):1]) +
           scale_color_manual(breaks = c("Network meta-analysis", "Network meta-regression"), values = c("#009E73", "black")) +
           scale_y_continuous(labels = scales::percent) +
           coord_flip(clip = "off") +
           theme_classic() +
           theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12),
                 axis.title.x = element_text(color = "black", face = "bold", size = 12), legend.position = "none",
                 legend.text =  element_text(color = "black", size = 12), legend.title =  element_text(color = "black", face = "bold", size = 12))


  ## Bring together both forest-plots
  forest.plots <- ggarrange(p1, p2, nrow = 1, ncol = 2, labels = c("A)", "B)"), common.legend = T, legend = "bottom")


  ## Prepare data-frame for scatter-plots/ error-bars
  EM.metareg00 <- cbind(rbind(data.frame(mean = metareg$EM[, 1], lower = metareg$EM[, 3], upper = metareg$EM[, 7]),
                              data.frame(mean = metareg$EM[, 1]*(-1), lower = metareg$EM[, 7]*(-1), upper = metareg$EM[, 3]*(-1))),
                        poss.pair.comp)
  EM.metareg.subset <- subset(EM.metareg00, EM.metareg00[5] == compar)
  EM.metareg0 <- rbind(EM.metareg.subset[, 1:3], c(rep(NA, 3)))
  EM.metareg <- EM.metareg0[order(sucra.full.new, decreasing = T), ]
  rownames(EM.metareg) <- NULL

  covar <- if (length(unique(covariate)) < 3) {
    unique(covariate)
  } else {
    unique(covariate) - mean(unique(covariate))
  }


  if (!is.element(measure, c("Odds ratio", "Ratio of means"))) {
    EM.meta.mean <- na.omit(rep(EM.metareg[, 1], each = length(covar))) + na.omit(rep(beta[, 1], each = length(covar)))*rep(covar, length(drug.names) - 1)
    EM.meta.lower <- na.omit(rep(EM.metareg[, 2], each = length(covar))) + na.omit(rep(beta[, 2], each = length(covar)))*rep(covar, length(drug.names) - 1)
    EM.meta.upper <- na.omit(rep(EM.metareg[, 3], each = length(covar))) + na.omit(rep(beta[, 3], each = length(covar)))*rep(covar, length(drug.names) - 1)
    EM.interc <- rep(na.omit(EM.full[, 1]), each = length(covar))
  } else {
    EM.meta.mean <- round(exp(na.omit(rep(EM.metareg[, 1], each = length(covar)) + (rep(beta[, 1], each = length(covar))*rep(covar, length(drug.names) - 1) ))), 2)
    EM.meta.lower <- round(exp(na.omit(rep(EM.metareg[, 2], each = length(covar)) + (rep(beta[, 2], each = length(covar))*rep(covar, length(drug.names) - 1) ))), 2)
    EM.meta.upper <- round(exp(na.omit(rep(EM.metareg[, 3], each = length(covar)) + (rep(beta[, 3], each = length(covar))*rep(covar, length(drug.names) - 1) ))), 2)
    EM.interc <- rep(round(exp(na.omit(EM.full[, 1])), 2), each = length(covar))
  }

  prepare <- data.frame(paste(rep(drug.names.sorted[drug.names.sorted != compar], each = length(covar)), "versus", compar),
                        rep(unique(covariate), length(drug.names) - 1),
                        EM.meta.mean,
                        EM.meta.lower,
                        EM.meta.upper,
                        EM.interc)
  colnames(prepare) <- c("comparison", "covariate", "mean", "lower", "upper", "intercept")


  ## Keep nodes with mean SUCRA > 0.20
  sucra.thres <- rep(sucra.meta[-1, 1][order(sucra.full[-1, 1], decreasing = T)], length(covariate))
  selection <- subset(prepare, sucra.thres > 0.40)



  ## Create the panel of scatterplots on the regression coefficient of the basic parameters
  p3 <- if(length(drug.names) <= 16) {
    ggplot(data = prepare, aes(x = covariate, y = mean, ymin = lower, ymax = upper) ) +
            geom_errorbar(colour = "#D55E00", size = 1, position = position_dodge(width = 0.5), width = 0.1) +
            geom_point(size = 1.5, colour = "black") +
            geom_line(aes(y = mean), size = 1, color = "#009E73") +
            geom_hline(yintercept = prepare$intercept, lty = 1, size = 1, col = "black") +
            geom_hline(yintercept = ifelse(!is.element(measure, c("Odds ratio", "Ratio of means")), 0, 1), lty = 2, size = 1, col = "grey") +
            facet_wrap(vars(factor(comparison, levels = unique(prepare$comparison))), scales = "free_y") +
            labs(x = cov.value[2], y = measure) +
            scale_y_continuous(trans = ifelse(!is.element(measure, c("Odds ratio", "Ratio of means")), "identity", "log10")) +
            theme_classic() +
            theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12), legend.position = "none",
                  axis.title.x = element_text(color = "black", face = "bold", size = 12), strip.text = element_text(size = 11))
  } else if(length(drug.names) > 16) {
    ggplot(data = selection, aes(x = covariate, y = round(mean, 2), ymin = lower, ymax = upper) ) +
            geom_errorbar(colour = "#D55E00", size = 1, position = position_dodge(width = 0.5), width = 0.1) +
            geom_point(size = 1.5, colour = "black") +
            geom_line(aes(y = mean), size = 1, color = "#009E73") +
            geom_hline(yintercept = intercept, lty = 1, size = 1, col = "black") +
            geom_hline(yintercept = ifelse(!is.element(measure, c("Odds ratio", "Ratio of means")), 0, 1), lty = 2, size = 1.5, col = "grey") +
           # geom_label(aes(x = min(covariate), y = 0, hjust = 0, vjust = 1, label = paste(intercept, ifelse(slope > 0, "+", "-"), abs(slope), "*covariate")), fill = "beige", colour = "black", fontface = "plain", size = 3.1) +
            facet_wrap(vars(factor(comparison, levels = unique(prepare$comparison))), scales = "free_y") +
            labs(x = cov.value[2], y = measure) +
            scale_y_continuous(trans = ifelse(!is.element(measure, c("Odds ratio", "Ratio of means")), "identity", "log10")) +
            theme_classic() +
            theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12), legend.position = "none",
                  axis.title.x = element_text(color = "black", face = "bold", size = 12), strip.text = element_text(size = 11))
  }

  ## Write all tables as .xlsx
  if (save.xls == TRUE) {
    write_xlsx(EM.both.models, paste0("Table NMA vs NMR", ".xlsx"))
    write_xlsx(table.model.assess, paste0("Table MTable Model Assessment_NMA vs NMR", ".xlsx"))

    if (is.element(metareg$covar.assumption, c("exchangeable", "independent")) ) {
      write_xlsx(reg.coeff, paste0("Table NMA vs NMR", ".xlsx"))
    }
  }


  ## Return results
  return(list(Table.effect.size = knitr::kable(EM.both.models),
              Table.model.assessment = knitr::kable(table.model.assess),
              Table.regression.coeffients = knitr::kable(reg.coeff),
              Interval.plots = forest.plots,
              Scatterplots = p3))
}

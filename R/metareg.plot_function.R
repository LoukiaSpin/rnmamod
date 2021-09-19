#' Plot the results from the meta-regression analysis
#'
#' @param covariate.name A factor is the covariate is metric. A factor vector of three elements, if the covariate is binary:
#'   the first element is the name of the covariate, and the last two elements are levels of the covariate (sorted as the reference, and non-reference level)
#'
#' @export
metareg.plot <- function(full, metareg, covariate.name, drug.names, save.xls) {


  options(warn = -1)

  save.xls <- if (missing(save.xls)) {
    FALSE
  } else {
    save.xls
  }


  model <- full$model
  measure <- full$measure
  drug.names <- if (missing(drug.names)) {
    message(cat(paste0("\033[0;", col = 32, "m", txt = "The argument 'drug.names' has not been defined. The intervention ID, as specified in 'data' is used as intervention names", "\033[0m", "\n")))
    nt <- as.factor(1:length(full$SUCRA[, 1]))
  } else {
    drug.names
  }
  # Sort the drugs by their NMA-SUCRA in decreasing order
  drug.names.sorted <- drug.names[order(sucra.full[, 1], decreasing = T)]

  covariate.name <- if (missing(covariate.name)) {
    stop("The argument 'covariate.name' needs to be defined", call. = F)
  } else if (!missing(covariate.name) & (length(unique(metareg$covariate)) > 2) & length(covariate.name) > 1) {
    stop("The argument 'covariate.name' must have one element", call. = F)
  } else if (!missing(covariate.name) & (length(unique(metareg$covariate)) == 2) & (length(covariate.name) < 3 || length(covariate.name) > 3)) {
    stop("The argument 'covariate.name' must be a vector of three elements", call. = F)
  } else {
    covariate.name
  }
  covariate <- if (metareg$covariate) {
    as.factor(metareg$covariate)
  } else {
    metareg$covariate
  }

  # Posterior results on the SUCRA value under NMA
  sucra.full <- round(full$SUCRA, 2)
  sucra.full.order <- round(full$SUCRA, 2)[order(sucra.full[, 1], decreasing = T), ]

  # Posterior results on the SUCRA value under meta-regression
  sucra.meta <- round(metareg$SUCRA, 2)
  sucra.meta.order <- round(metareg$SUCRA, 2)[order(sucra.full[, 1], decreasing = T), ]

  # Posterior results on the effect estimates of comparisons with the reference under NMA (sorted by posterior mean NMA-SUCRA in decreasing order)
  EM.full <- rbind(rep(NA, 9), round(full$EM.ref, 2))[order(sucra.full[, 1], decreasing = T), ]

  # Posterior results on the effect estimates of comparisons with the reference under meta-regression (sorted by posterior mean NMA-SUCRA in decreasing order)
  EM.meta <- rbind(rep(NA, 9), round(metareg$EM.re, 2))[order(sucra.full[, 1], decreasing = T), ]

  # Posterior mean of regression coefficients for the basic parameters (sorted by posterior mean NMA-SUCRA in decreasing order)
  beta <- if (is.element(metareg$covar.assumption, c("exchangeable", "independent")) ) {
    round(metareg$beta, 3)[order(sucra.full[, 1], decreasing = T), ]
  } else {
    round(metareg$beta, 3)
  }

  # Posterior results on the effect estimates of all unique pairwise comparisons (NMA)
  EM.full.all <- round(full$EM[, c(1, 3, 7)], 3)

  # Posterior results on the effect estimates of all unique pairwise comparisons (NMR)
  EM.meta.all <- round(full$EM[, c(1, 3, 7)], 3)

  # Posterior mean of regression coefficients for all unique pairwise comparisons
  beta.all <- if (is.element(metareg$covar.assumption, c("exchangeable", "independent")) ) {
    round(metareg$beta.all[, c(1, 3, 7)], 3)
  }

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
  if (!is.factor(covariate) & (!is.element(measure, c("Odds ratio", "Ratio of means")))) {
    CrI.full <- paste0("(", EM.full[, 3], ",", " ", EM.full[, 7], ")", ifelse(EM.full[, 3] > 0 | EM.full[, 7] < 0, "*", " "))
    CrI.meta <- paste0("(", EM.meta[, 3], ",", " ", EM.meta[, 7], ")", ifelse(EM.meta[, 3] > 0 | EM.meta[, 7] < 0, "*", " "))
    CrI.beta <- if (is.element(metareg$covar.assumption, c("exchangeable", "independent")) ) {
      paste0("(", beta[, 3], ",", " ", beta[, 7], ")", ifelse(beta[, 3] > 0 | beta[, 7] < 0, "*", " "))
    } else {
      paste0("(", beta[3], ",", " ", beta[7], ")", ifelse(beta[3] > 0 | beta[7] < 0, "*", " "))
    }

    # Tabulate results on comparisons with the reference (both models)
    EM.both.models <- na.omit(data.frame(drug.names.sorted,
                                         EM.full[, 1],
                                         CrI.full,
                                         EM.meta[, 1],
                                         CrI.meta))
    colnames(EM.both.models) <- c("Comparison", "Mean NMA", "95% CrI NMA", "Mean NMR", "95% CrI NMR")
    rownames(EM.both.models) <- NULL

    # Results on the regression coefficient
    if (is.element(metareg$covar.assumption, c("exchangeable", "independent")) ) {
      reg.coeff <- na.omit(data.frame(drug.names.sorted, beta[, 1], CrI.beta))
      colnames(reg.coeff) <- c("Comparison", "Mean beta", "95% CrI beta")
    } else {
      reg.coeff <- data.frame(beta[1], CrI.beta)
      colnames(reg.coeff) <- c("Mean beta", "95% CrI beta")
    }
    rownames(reg.coeff) <- NULL
  } else if (!is.factor(covariate) & (is.element(measure, c("Odds ratio", "Ratio of means")))) {
    CrI.full <- paste0("(", round(exp(EM.full[, 3]), 2), ",", " ", round(exp(EM.full[, 7]), 2), ")", ifelse(EM.full[, 3] > 0 | EM.full[, 7] < 0, "*", " "))
    CrI.meta <- paste0("(", round(exp(EM.meta[, 3]), 2), ",", " ", round(exp(EM.meta[, 7]), 2), ")", ifelse(EM.meta[, 3] > 0 | EM.meta[, 7] < 0, "*", " "))
    CrI.beta <- if (is.element(metareg$covar.assumption, c("exchangeable", "independent")) ) {
      paste0("(", round(exp(beta[, 3]), 2), ",", " ", round(exp(beta[, 7]), 2), ")", ifelse(beta[, 3] > 0 | beta[, 7] < 0, "*", " "))
    } else {
      paste0("(", round(exp(beta[3]), 2), ",", " ", round(exp(beta[7]), 2), ")", ifelse(beta[3] > 0 | beta[7] < 0, "*", " "))
    }

    # Tabulate results on comparisons with the reference (both models)
    EM.both.models <- na.omit(data.frame(drug.names.sorted,
                                         round(exp(EM.full[, 1]), 2),
                                         CrI.full,
                                         round(exp(EM.meta[, 1]), 2),
                                         CrI.meta))
    colnames(EM.both.models) <- c("Comparison", "Mean NMA", "95% CrI NMA", "Mean NMR", "95% CrI NMR")
    rownames(EM.both.models) <- NULL

    # Results on the regression coefficient
    if (is.element(metareg$covar.assumption, c("exchangeable", "independent")) ) {
      reg.coeff <- na.omit(data.frame(drug.names.sorted, round(exp(beta[, 1]), 2), CrI.beta))
      colnames(reg.coeff) <- c("Comparison", "Mean beta", "95% CrI beta")
    } else {
      reg.coeff <- data.frame(round(exp(beta[1]), 2), CrI.beta)
      colnames(reg.coeff) <- c("Mean beta", "95% CrI beta")
    }
    rownames(reg.coeff) <- NULL
  } else if (is.factor(covariate) & (!is.element(measure, c("Odds ratio", "Ratio of means")))) {
    EM.meta.nf <- if (is.element(metareg$covar.assumption, c("exchangeable", "independent")) ) {
      EM.meta[, c(1, 3, 7)] + beta[, c(1, 3, 7)]
    } else {
      EM.meta[, c(1, 3, 7)] + matrix(rep(beta[c(1, 3, 7)], length(drug.names)), ncol = 3, byrow = T)
    }
    CrI.full <- paste0("(", EM.full[, 3], ",", " ", EM.full[, 7], ")", ifelse(EM.full[, 3] > 0 | EM.full[, 7] < 0, "*", " "))
    CrI.meta <- paste0("(", EM.meta[, 3], ",", " ", EM.meta[, 7], ")", ifelse(EM.meta[, 3] > 0 | EM.meta[, 7] < 0, "*", " "))
    CrI.meta.nf <- paste0("(", EM.meta.nf[, 2], ",", " ", EM.meta.nf[, 3], ")", ifelse(EM.meta.nf[, 2] > 0 | EM.meta.nf[, 3] < 0, "*", " "))
    CrI.beta <- paste0("(", beta[, 3], ",", " ", beta[, 7], ")", ifelse(beta[, 3] > 0 | beta[, 7] < 0, "*", " "))

    # Tabulate results on comparisons with the reference (both models)
    EM.both.models <- na.omit(data.frame(drug.names.sorted,
                                         round(EM.full[, 1], 2),
                                         CrI.full,
                                         round(EM.meta[, 1], 2),
                                         CrI.meta,
                                         round(EM.meta.nf[, 1], 2),
                                         CrI.meta.nf))
    colnames(EM.both.models) <- c("Comparison", "Mean NMA", "95% CrI NMA",
                                  paste("Mean NMR", covariate.name[2]),
                                  paste("95% CrI NMR", covariate.name[2]),
                                  paste("Mean NMR", covariate.name[3]),
                                  paste("95% CrI NMR", covariate.name[3]))
    rownames(EM.both.models) <- NULL

    # Results on the regression coefficient
    if (is.element(metareg$covar.assumption, c("exchangeable", "independent")) ) {
      reg.coeff <- na.omit(data.frame(drug.names.sorted, beta[, 1], CrI.beta))
      colnames(reg.coeff) <- c("Comparison", "Mean beta", "95% CrI beta")
    } else {
      reg.coeff <- data.frame(beta[1], CrI.beta)
      colnames(reg.coeff) <- c("Mean beta", "95% CrI beta")
    }
    rownames(reg.coeff) <- NULL
  } else if (is.factor(covariate) & (is.element(measure, c("Odds ratio", "Ratio of means")))) {
    EM.meta.nf <- if (is.element(metareg$covar.assumption, c("exchangeable", "independent")) ) {
      round(exp(EM.meta[, c(1, 3, 7)] + beta[, c(1, 3, 7)]), 2)
    } else {
      round(exp(EM.meta[, c(1, 3, 7)] + matrix(rep(beta[c(1, 3, 7)], length(drug.names)), ncol = 3, byrow = T)), 2)
    }
    CrI.full <- paste0("(", round(exp(EM.full[, 3]), 2), ",", " ", round(exp(EM.full[, 7]), 2), ")", ifelse(EM.full[, 3] > 0 | EM.full[, 7] < 0, "*", " "))
    CrI.meta <- paste0("(", round(exp(EM.meta[, 3]), 2), ",", " ", round(exp(EM.meta[, 7]), 2), ")", ifelse(EM.meta[, 3] > 0 | EM.meta[, 7] < 0, "*", " "))
    CrI.meta.nf <- paste0("(", EM.meta.nf[, 2], ",", " ", EM.meta.nf[, 3], ")", ifelse(EM.meta.nf[, 2] > 1 | EM.meta.nf[, 3] < 1, "*", " "))
    CrI.beta <- paste0("(", round(exp(beta[, 3]), 2), ",", " ", round(exp(beta[, 7]), 2), ")", ifelse(beta[, 3] > 0 | beta[, 7] < 0, "*", " "))

    # Tabulate results on comparisons with the reference (both models)
    EM.both.models <- na.omit(data.frame(drug.names.sorted,
                                         round(exp(EM.full[, 1]), 2),
                                         CrI.full,
                                         round(exp(EM.meta[, 1]), 2),
                                         CrI.meta,
                                         EM.meta.nf[, 1],
                                         CrI.meta.nf))
    colnames(EM.both.models) <- c("Comparison", "Mean NMA", "95% CrI NMA",
                                  paste("Mean NMR", covariate.name[2]),
                                  paste("95% CrI NMR", covariate.name[2]),
                                  paste("Mean NMR", covariate.name[3]),
                                  paste("95% CrI NMR", covariate.name[3]))
    rownames(EM.both.models) <- NULL

    # Results on the regression coefficient
    if (is.element(metareg$covar.assumption, c("exchangeable", "independent")) ) {
      reg.coeff <- na.omit(data.frame(drug.names.sorted, round(exp(beta[, 1]), 2), CrI.beta))
      colnames(reg.coeff) <- c("Comparison", "Mean beta", "95% CrI beta")
    } else {
      reg.coeff <- data.frame(round(exp(beta[1]), 2), CrI.beta)
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
  if(!is.factor(covariate) & !is.element(measure, c("Odds ratio", "Ratio of means"))) {
    prepare.EM.metric <- data.frame(rep(length(drug.names):1, 2),
                                    rep(drug.names.sorted, 2),
                                    round(rbind(EM.full[, c(1, 3, 7)],
                                                EM.meta[, c(1, 3, 7)]), 2),
                                    rep(c("Network meta-analysis", "Network meta-regression"), each = length(drug.names)))
    colnames(prepare.EM.metric) <- c("order", "comparison", "mean", "lower", "upper", "analysis")
    rownames(prepare.EM.metric) <- NULL
  } else if(!is.factor(covariate) & is.element(measure, c("Odds ratio", "Ratio of means"))) {
    prepare.EM.metric <- data.frame(rep(length(drug.names):1, 2),
                                    rep(drug.names.sorted, 2),
                                    round(rbind(exp(EM.full[, c(1, 3, 7)]),
                                                exp(EM.meta[, c(1, 3, 7)])), 2),
                                    rep(c("Network meta-analysis", "Network meta-regression"), each = length(drug.names)))
    colnames(prepare.EM.metric) <- c("order", "comparison", "mean", "lower", "upper", "analysis")
    rownames(prepare.EM.metric) <- NULL
  } else if(is.factor(covariate) & !is.element(measure, c("Odds ratio", "Ratio of means"))) {
    prepare.EM.nominal <- data.frame(rep(length(drug.names):1, 3),
                                     rep(drug.names.sorted, 3),
                                     round(rbind(EM.full[, c(1, 3, 7)],
                                                 EM.meta[, c(1, 3, 7)],
                                                 EM.meta.nf), 2),
                                     rep(c("Network meta-analysis",
                                           paste("Network meta-regression at", covariate.name[2]),
                                           paste("Network meta-regression at", covariate.name[3])),
                                         each = length(drug.names)))
    colnames(prepare.EM.nominal) <- c("order", "comparison", "mean", "lower", "upper", "analysis")
    rownames(prepare.EM.nominal) <- NULL
  } else if(is.factor(covariate) & is.element(measure, c("Odds ratio", "Ratio of means"))) {
    prepare.EM.nominal <- data.frame(rep(length(drug.names):1, 3),
                                     rep(drug.names.sorted, 3),
                                     round(rbind(exp(EM.full[, c(1, 3, 7)]),
                                                 exp(EM.meta[, c(1, 3, 7)]),
                                                 EM.meta.nf), 2),
                                     rep(c("Network meta-analysis",
                                           paste("Network meta-regression at", covariate.name[2]),
                                           paste("Network meta-regression at", covariate.name[3])),
                                         each = length(drug.names)))
    colnames(prepare.EM.nominal) <- c("order", "comparison", "mean", "lower", "upper", "analysis")
    rownames(prepare.EM.nominal) <- NULL
  }



  ## Forest plots of reference-comparisons on effect estimate
  p1 <- if(!is.factor(covariate)) {
    ggplot(data = prepare.EM.metric, aes(x = as.factor(order), y = as.numeric(mean), ymin = as.numeric(lower), ymax = as.numeric(upper), colour = analysis, group = analysis)) +
      geom_linerange(size = 2, position = position_dodge(width = 0.5), width = 0.1) +
      geom_hline(yintercept = ifelse(!is.element(measure, c("Odds ratio", "Ratio of means")), 0, 1), lty = 2, size = 1.3, col = "grey53") +
      geom_point(size = 1.5,  colour = "white", stroke = 0.3, position = position_dodge(width = 0.5)) +
      geom_text(aes(x = as.factor(order), y = round(as.numeric(mean), 2), label = paste0(round(as.numeric(mean), 2), " ", "(", round(as.numeric(lower), 2), ",", " ", round(as.numeric(upper), 2), ")"),
                    hjust = 0, vjust = -0.5), color = "black", size = 4.0, check_overlap = F, parse = F, position = position_dodge(width = 0.5), inherit.aes = T) +
      labs(x = "", y = measure, colour = "Analysis") +
      scale_x_discrete(breaks = as.factor(1:length(drug.names.sorted)), labels = drug.names.sorted[length(drug.names.sorted):1]) +
      scale_color_manual(breaks = c("Network meta-analysis", "Network meta-regression"), values = c("#009E73", "black")) +
      geom_label(aes(x = unique(as.factor(order)[is.na(mean)]), y = ifelse(!is.element(measure, c("Odds ratio", "Ratio of means")), 0, 1), hjust = 0, vjust = 1, label = "Reference intervention"),
                 fill = "beige", colour = "black", fontface = "plain", size = 4) +
      scale_y_continuous(trans = ifelse(!is.element(measure, c("Odds ratio", "Ratio of means")), "identity", "log10")) +
      coord_flip() +
      theme_classic() +
      theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12),
            axis.title.x = element_text(color = "black", face = "bold", size = 12), legend.position = "none",
            legend.text =  element_text(color = "black", size = 12), legend.title =  element_text(color = "black", face = "bold", size = 12))
  } else {
    ggplot(data = prepare.EM.nominal, aes(x = as.factor(order), y = as.numeric(mean), ymin = as.numeric(lower), ymax = as.numeric(upper), colour = analysis, group = analysis)) +
      geom_linerange(size = 2, position = position_dodge(width = 0.5), width = 0.1) +
      geom_hline(yintercept = ifelse(!is.element(measure, c("Odds ratio", "Ratio of means")), 0, 1), lty = 2, size = 1.3, col = "grey53") +
      geom_point(size = 1.5,  colour = "white", stroke = 0.3, position = position_dodge(width = 0.5)) +
      geom_text(aes(x = as.factor(order), y = round(as.numeric(mean), 2), label = paste0(round(as.numeric(mean), 2), " ", "(", round(as.numeric(lower), 2), ",", " ", round(as.numeric(upper), 2), ")"),
                    hjust = 0, vjust = -0.5), color = "black", size = 4.0, check_overlap = F, parse = F, position = position_dodge(width = 0.5), inherit.aes = T) +
      labs(x = "", y = measure, colour = "Analysis") +
      scale_x_discrete(breaks = as.factor(1:length(drug.names.sorted)), labels = drug.names.sorted[length(drug.names.sorted):1]) +
      scale_color_manual(breaks = c("Network meta-analysis", paste("Network meta-regression at", covariate.name[3]), paste("Network meta-regression at", covariate.name[2])),
                         values = c("#009E73", "#D55E00", "#56B4E9")) +
      geom_label(aes(x = unique(as.factor(order)[is.na(mean)]), y = ifelse(!is.element(measure, c("Odds ratio", "Ratio of means")), 0, 1), hjust = 0, vjust = 1, label = "Reference intervention"),
                 fill = "beige", colour = "black", fontface = "plain", size = 4) +
      scale_y_continuous(trans = ifelse(!is.element(measure, c("Odds ratio", "Ratio of means")), "identity", "log10")) +
      coord_flip() +
      theme_classic() +
      theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12),
            axis.title.x = element_text(color = "black", face = "bold", size = 12), legend.position = "none",
            legend.text =  element_text(color = "black", size = 12), legend.title =  element_text(color = "black", face = "bold", size = 12))
  }



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
  forest.plots <- ggarrange(p1, p2, nrow = 1, ncol = 2, labels = c("A)", "B)"), common.legend = F, legend = "bottom")



  ## Prepare data-frame for scatter-plots/ error-bars (comparisons with the reference)
  if (!is.element(measure, c("Odds ratio", "Ratio of means"))) {
    if (is.element(metareg$covar.assumption, c("exchangeable", "independent")) ) {
      EM.meta.mean <- na.omit(rep(EM.meta[, 1], each = length(covariate)) + (rep(beta[, 1], each = length(covariate))*covariate))
      EM.meta.lower <- na.omit(rep(EM.meta[, 3], each = length(covariate)) + (rep(beta[, 3], each = length(covariate))*covariate))
      EM.meta.upper <- na.omit(rep(EM.meta[, 7], each = length(covariate)) + (rep(beta[, 7], each = length(covariate))*covariate))
      slope <- paste0("slope, (95% CrI):", " ", beta[-length(drug.names), 1], ", (", beta[-length(drug.names), 3], ",", beta[-length(drug.names), 7], ")")
    } else {
      EM.meta.mean <- na.omit(rep(EM.meta[, 1], each = length(covariate)) + (rep(beta[1], each = length(covariate))*covariate))
      EM.meta.lower <- na.omit(rep(EM.meta[, 3], each = length(covariate)) + (rep(beta[3], each = length(covariate))*covariate))
      EM.meta.upper <- na.omit(rep(EM.meta[, 7], each = length(covariate)) + (rep(beta[7], each = length(covariate))*covariate))
      slope <- paste0("slope, (95% CrI):", " ", beta[1], ", (", beta[3], ",", beta[7], ")")
    }
    intercept <- na.omit(rep(EM.meta[, 1], each = length(covariate)))
  } else {
    if (is.element(metareg$covar.assumption, c("exchangeable", "independent")) ) {
      EM.meta.mean <- round(exp(na.omit(rep(EM.meta[, 1], each = length(covariate)) + (rep(beta[, 1], each = length(covariate))*covariate))), 2)
      EM.meta.lower <- round(exp(na.omit(rep(EM.meta[, 3], each = length(covariate)) + (rep(beta[, 3], each = length(covariate))*covariate))), 2)
      EM.meta.upper <- round(exp(na.omit(rep(EM.meta[, 7], each = length(covariate)) + (rep(beta[, 7], each = length(covariate))*covariate))), 2)
      slope <- round(exp(rep(beta[-length(drug.names), 1], each = length(covariate))), 2)
    } else {
      EM.meta.mean <- round(exp(na.omit(rep(EM.meta[, 1], each = length(covariate)) + (rep(beta[1], each = length(covariate))*covariate))), 2)
      EM.meta.lower <- round(exp(na.omit(rep(EM.meta[, 3], each = length(covariate)) + (rep(beta[3], each = length(covariate))*covariate))), 2)
      EM.meta.upper <- round(exp(na.omit(rep(EM.meta[, 7], each = length(covariate)) + (rep(beta[7], each = length(covariate))*covariate))), 2)
      slope <- round(exp(rep(beta[1], length(covariate))), 2)
    }
    intercept <- round(exp(na.omit(rep(EM.meta[, 1], each = length(covariate)))), 2)
  }

  prepare <- data.frame(paste(rep(drug.names.sorted[-length(drug.names)], each = length(covariate)), "versus", na.omit(drug.names.sorted[is.na(EM.full)])),
                        rep(covariate, length(drug.names) - 1),
                        EM.meta.mean,
                        EM.meta.lower,
                        EM.meta.upper,
                        slope,
                        intercept)
  colnames(prepare) <- c("comparison", "covariate", "mean", "lower", "upper", "slope", "intercept")


  ## Keep nodes with mean SUCRA > 0.20
  sucra.thres <- rep(sucra.meta[-1, 1][order(sucra.full[-1, 1], decreasing = T)], length(covariate))
  selection <- subset(prepare, sucra.thres > 0.40)



  ## Create the panel of scatterplots on the regression coefficient of the basic parameters
  p3 <- if(!is.factor(covariate) & length(drug.names) <= 16) {
    ggplot(data = prepare, aes(x = covariate, y = mean, ymin = lower, ymax = upper) ) +
            geom_errorbar(colour = "#D55E00", size = 1, position = position_dodge(width = 0.5), width = 0.1) +
            geom_point(size = 1.5, colour = "black") +
            geom_line(aes(y = mean), size = 1, color = "#009E73") +
            geom_hline(yintercept = ifelse(!is.element(measure, c("Odds ratio", "Ratio of means")), EM.full[, 1], exp(EM.full[, 1])), lty = 1, size = 1, col = "black") +
            geom_hline(yintercept = ifelse(!is.element(measure, c("Odds ratio", "Ratio of means")), 0, 1), lty = 2, size = 1, col = "grey") +
            #geom_label(aes(x = min(covariate), y = 0, hjust = 0, vjust = 1, label = paste(intercept, ifelse(slope > 0, "+", "-"), abs(slope), "*covariate")), fill = "beige", colour = "black", fontface = "plain", size = 3.1) +
            facet_wrap(vars(factor(comparison, levels = unique(prepare$comparison))), scales = "free_y") +
            labs(x = covariate.name, y = measure) +
            scale_y_continuous(trans = ifelse(!is.element(measure, c("Odds ratio", "Ratio of means")), "identity", "log10")) +
            theme_classic() +
            theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12), legend.position = "none",
                  axis.title.x = element_text(color = "black", face = "bold", size = 12), strip.text = element_text(size = 11))
  } else if(!is.factor(covariate) & length(drug.names) > 16) {
    ggplot(data = selection, aes(x = covariate, y = round(mean, 2), ymin = lower, ymax = upper) ) +
            geom_errorbar(colour = "#D55E00", size = 1, position = position_dodge(width = 0.5), width = 0.1) +
            geom_point(size = 1.5, colour = "black") +
            geom_line(aes(y = mean), size = 1, color = "#009E73") +
            geom_hline(yintercept = ifelse(!is.element(measure, c("Odds ratio", "Ratio of means")), EM.full[, 1], exp(EM.full[, 1])), lty = 1, size = 1, col = "black") +
            geom_hline(yintercept = ifelse(!is.element(measure, c("Odds ratio", "Ratio of means")), 0, 1), lty = 2, size = 1.5, col = "grey") +
           # geom_label(aes(x = min(covariate), y = 0, hjust = 0, vjust = 1, label = paste(intercept, ifelse(slope > 0, "+", "-"), abs(slope), "*covariate")), fill = "beige", colour = "black", fontface = "plain", size = 3.1) +
            facet_wrap(vars(factor(comparison, levels = unique(prepare$comparison))), scales = "free_y") +
            labs(x = covariate.name, y = measure) +
            scale_y_continuous(trans = ifelse(!is.element(measure, c("Odds ratio", "Ratio of means")), "identity", "log10")) +
            theme_classic() +
            theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12), legend.position = "none",
                  axis.title.x = element_text(color = "black", face = "bold", size = 12), strip.text = element_text(size = 11))
  } else if(is.factor(covariate) & length(drug.names) <= 16) {
    ggplot(data = prepare, aes(x = as.factor(covariate), y = mean, ymin = lower, ymax = upper, group = comparison)) +
            geom_errorbar(size = 2, position = position_dodge(width = 0.1), width = 0.1, colour = "#009E73") +
            geom_point(size = 1.5, colour = "white", stroke = 0.3, position = position_dodge(width = 0.1)) +
            geom_line(size = 1, color = "#D55E00") +
            geom_hline(yintercept = ifelse(!is.element(measure, c("Odds ratio", "Ratio of means")), EM.full[, 1], exp(EM.full[, 1])), lty = 1, size = 1, col = "black") +
            geom_hline(yintercept = ifelse(!is.element(measure, c("Odds ratio", "Ratio of means")), 0, 1), lty = 2, size = 1, col = "grey") +
            geom_text(aes(x = as.factor(covariate), y = round(as.numeric(mean), 2), label = round(as.numeric(mean), 2)),
                      color = "black", hjust = -0.2, vjust = -0.5, size = 4.0, check_overlap = F, parse = F, position = position_dodge(width = 0.8), inherit.aes = T) +
            facet_wrap(vars(factor(comparison, levels = unique(prepare$comparison))), scales = "free_y") +
            labs(x = covariate.name, y = measure) +
            scale_x_discrete(breaks = as.factor(c(0, 1)), labels = covariate) +
            scale_y_continuous(trans = ifelse(!is.element(measure, c("Odds ratio", "Ratio of means")), "identity", "log10")) +
            theme_classic() +
            theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12), legend.position = "none",
                  axis.title.x = element_text(color = "black", face = "bold", size = 12), strip.text = element_text(size = 11))
  } else if(is.factor(covariate) & length(drug.names) > 16) {
    ggplot(data = selection, aes(x = as.factor(covariate), y = mean, ymin = lower, ymax = upper, group = comparison)) +
            geom_errorbar(size = 2, position = position_dodge(width = 0.5), width = 0.1, colour = "#009E73") +
            geom_point(size = 1.5, colour = "white", stroke = 0.3, position = position_dodge(width = 0.5)) +
            geom_line( size = 1, color = "#D55E00") +
            geom_hline(yintercept = ifelse(!is.element(measure, c("Odds ratio", "Ratio of means")), EM.full[, 1], exp(EM.full[, 1])), lty = 1, size = 1, col = "black") +
            geom_hline(yintercept = ifelse(!is.element(measure, c("Odds ratio", "Ratio of means")), 0, 1), lty = 2, size = 1, col = "grey") +
            geom_text(aes(x = as.factor(covariate), y = round(as.numeric(mean), 2), label = round(as.numeric(mean), 2)),
                      color = "black", hjust = -0.2, vjust = -0.5, size = 4.0, check_overlap = F, parse = F, position = position_dodge(width = 0.8), inherit.aes = T) +
            facet_wrap(vars(factor(comparison, levels = unique(prepare$comparison))), scales = "free_y") +
            labs(x = covariate.name, y = measure) +
            scale_x_discrete(breaks = as.factor(c(0, 1)), labels = covariate) +
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
              interval.plots = forest.plots,
              Scatterplots = p3))

}

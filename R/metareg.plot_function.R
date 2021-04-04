#' Plot the results from the meta-regression analysis
#'
#' @export
metareg.plot <- function(full, metareg, covariate, covar.values, drug.names, expon) {


  options(warn = -1)

  ## The results on the following parameters will be used:
  # Posterior results on the SUCRA value under NMA
  sucra.full <- round(full$SUCRA, 2)

  # Posterior results on the SUCRA value under meta-regression
  sucra.meta <- round(metareg$SUCRA, 2)

  # Posterior results on the effect estimates of comparisons with the reference under NMA (sorted by posterior mean NMA-SUCRA in decreasing order)
  EM.full <- rbind(rep(NA, 9), round(full$EM.ref, 2))[order(sucra.full[, 1], decreasing = T), ]

  # Posterior results on the effect estimates of comparisons with the reference under meta-regression (sorted by posterior mean NMA-SUCRA in decreasing order)
  EM.meta <- rbind(rep(NA, 9), round(metareg$EM.re, 2))[order(sucra.full[, 1], decreasing = T), ]

  # Posterior mean of regression coefficients for the basic parameters (sorted by posterior mean NMA-SUCRA in decreasing order)
  beta <- round(metareg$beta, 2)[order(sucra.full[, 1], decreasing = T), ]

  # Posterior results on between-trial standard deviation under NMA
  tau.full <- round(full$tau, 2)

  # Posterior results on between-trial standard deviation under meta-regression
  tau.meta <- round(metareg$tau, 2)

  # Posterior mean of model assessment measures under NMA
  model.assess.NMA <- round(full$model.assessment, 2)

  # Posterior mean of model assessment measures under meta-regression
  model.assess.meta <- round(metareg$model.assessment, 2)

  # Sort the drugs by their SUCRA in decreasing order
  drug.names.sorted <- drug.names[order(sucra.full[, 1], decreasing = T)]

  # Covariate
  covariate <- na.omit(unlist(covariate))


  ## A data-frame with the effect estimates and regression coefficients of reference-comparisons from both analyses (Sort by NMA-SUCRA in decreasing order)
  if (!is.character(covar.values[[1]]) & (missing(expon) || expon == F)) {


    ## Calculate the effect estimate at two other values of the metric covariate: one smaller and one larger than the mean
    EM.meta.s <- round(EM.meta[, c(1, 3, 7)] + beta[, c(1, 3, 7)]*(covar.values[[1]] - mean(covariate)), 2)
    EM.meta.l <- round(EM.meta[, c(1, 3, 7)] + beta[, c(1, 3, 7)]*(covar.values[[2]] - mean(covariate)), 2)


    ## Isolate the 95% CrI and indicate whether the evidence is string (*) or weak (no '*')
    CrI.full <- paste0("(", EM.full[, 3], ",", " ", EM.full[, 7], ")", ifelse(EM.full[, 3] > 0 | EM.full[, 7] < 0, "*", " "))
    CrI.meta <- paste0("(", EM.meta[, 3], ",", " ", EM.meta[, 7], ")", ifelse(EM.meta[, 3] > 0 | EM.meta[, 7] < 0, "*", " "))
    CrI.meta.s <- paste0("(", EM.meta.s[, 2], ",", " ", EM.meta.s[, 3], ")", ifelse(EM.meta.s[, 2] > 0 | EM.meta.s[, 3] < 0, "*", " "))
    CrI.meta.l <- paste0("(", EM.meta.l[, 2], ",", " ", EM.meta.l[, 3], ")", ifelse(EM.meta.l[, 2] > 0 | EM.meta.l[, 3] < 0, "*", " "))
    CrI.beta <- paste0("(", beta[, 3], ",", " ", beta[, 7], ")", ifelse(beta[, 3] > 0 | beta[, 7] < 0, "*", " "))


    EM.both.models <- na.omit(data.frame(drug.names.sorted, EM.full[, 1], CrI.full, EM.meta.s[, 1], CrI.meta.s, EM.meta[, 1], CrI.meta, EM.meta.l[, 1], CrI.meta.l, beta[, 1], CrI.beta))
    colnames(EM.both.models) <- c("Comparison", "Poster. mean NMA", "95% CrI NMA", "Poster. mean NMReg (< mean)", "95% CrI NMReg (< mean)",
                                  "Poster. mean NMReg (at mean)", "95% CrI NMReg (at mean)", "Poster. mean NMReg (> mean)", "95% CrI NMReg (> mean)",
                                  "Poster. mean beta", "95% CrI beta")
    rownames(EM.both.models) <- NULL

  } else if (!is.character(covar.values[[1]]) & expon == T) {

    ## Calculate the effect estimate at two other values of the metric covariate: one smaller and one larger than the mean
    EM.meta.s <- round(exp(EM.meta[, c(1, 3, 7)] + beta[, c(1, 3, 7)]*(covar.values[[1]] - mean(covariate))), 2)
    EM.meta.l <- round(exp(EM.meta[, c(1, 3, 7)] + beta[, c(1, 3, 7)]*(covar.values[[2]] - mean(covariate))), 2)


    ## Isolate the 95% CrI and indicate whether the evidence is string (*) or weak (no '*')
    CrI.full <- paste0("(", round(exp(EM.full[, 3]), 2), ",", " ", round(exp(EM.full[, 7]), 2), ")", ifelse(EM.full[, 3] > 0 | EM.full[, 7] < 0, "*", " "))
    CrI.meta <- paste0("(", round(exp(EM.meta[, 3]), 2), ",", " ", round(exp(EM.meta[, 7]), 2), ")", ifelse(EM.meta[, 3] > 0 | EM.meta[, 7] < 0, "*", " "))
    CrI.meta.s <- paste0("(", EM.meta.s[, 2], ",", " ", EM.meta.s[, 3], ")", ifelse(EM.meta.s[, 2] > 1 | EM.meta.s[, 3] < 1, "*", " "))
    CrI.meta.l <- paste0("(", EM.meta.l[, 2], ",", " ", EM.meta.l[, 3], ")", ifelse(EM.meta.l[, 2] > 1 | EM.meta.l[, 3] < 1, "*", " "))
    CrI.beta <- paste0("(", round(exp(beta[, 3]), 2), ",", " ", round(exp(beta[, 7]), 2), ")", ifelse(beta[, 3] > 0 | beta[, 7] < 0, "*", " "))

    EM.both.models <- na.omit(data.frame(drug.names.sorted, round(exp(EM.full[, 1]), 2), CrI.full, round(exp(EM.meta.s[, 1]), 2), CrI.meta.s, round(exp(EM.meta[, 1]), 2), CrI.meta,
                                        round(exp(EM.meta.l[, 1]), 2), CrI.meta.l, round(exp(beta[, 1]), 2), CrI.beta))
    colnames(EM.both.models) <- c("Comparison", "Poster. mean NMA", "95% CrI NMA", "Poster. mean NMReg (< mean)", "95% CrI NMReg (< mean)",
                                  "Poster. mean NMReg (at mean)", "95% CrI NMReg (at mean)", "Poster. mean NMReg (> mean)", "95% CrI NMReg (> mean)",
                                  "Poster. mean beta", "95% CrI beta")
    rownames(EM.both.models) <- NULL

  } else if (is.character(covar.values[[1]]) & (missing(expon) || expon == F)) {


    ## Calculate the effect estimate at the non-reference level of the binary covariate
    EM.meta.nf <- EM.meta[, c(1, 3, 7)] + beta[, c(1, 3, 7)]*1


    ## Isolate the 95% CrI and indicate whether the evidence is string (*) or weak (no '*')
    CrI.full <- paste0("(", EM.full[, 3], ",", " ", EM.full[, 7], ")", ifelse(EM.full[, 3] > 0 | EM.full[, 7] < 0, "*", " "))
    CrI.meta <- paste0("(", EM.meta[, 3], ",", " ", EM.meta[, 7], ")", ifelse(EM.meta[, 3] > 0 | EM.meta[, 7] < 0, "*", " "))
    CrI.meta.nf <- paste0("(", EM.meta.nf[, 2], ",", " ", EM.meta.nf[, 3], ")", ifelse(EM.meta.nf[, 2] > 0 | EM.meta.nf[, 3] < 0, "*", " "))
    CrI.beta <- paste0("(", beta[, 3], ",", " ", beta[, 7], ")", ifelse(beta[, 3] > 0 | beta[, 7] < 0, "*", " "))

    EM.both.models <- na.omit(data.frame(drug.names.sorted, round(EM.full[, 1], 2), CrI.full, round(EM.meta[, 1], 2), CrI.meta,
                                         round(EM.meta.nf[, 1], 2), CrI.meta.nf, round(beta[, 1], 2), CrI.beta))
    colnames(EM.both.models) <- c("Comparison", "Poster. mean NMA", "95% CrI NMA", "Poster. mean NMReg (ref)", "95% CrI NMReg (ref)",
                                  "Poster. mean NMReg (non-ref)", "95% CrI NMReg (non-ref)", "Poster. mean beta", "95% CrI beta")
    rownames(EM.both.models) <- NULL


  } else if (is.character(covar.values[[1]]) & expon == T) {


    ## Calculate the effect estimate at the non-reference level of the binary covariate
    EM.meta.nf <- round(exp(EM.meta[, c(1, 3, 7)] + beta[, c(1, 3, 7)]*1), 2)


    ## Isolate the 95% CrI and indicate whether the evidence is string (*) or weak (no '*')
    CrI.full <- paste0("(", round(exp(EM.full[, 3]), 2), ",", " ", round(exp(EM.full[, 7]), 2), ")", ifelse(EM.full[, 3] > 0 | EM.full[, 7] < 0, "*", " "))
    CrI.meta <- paste0("(", round(exp(EM.meta[, 3]), 2), ",", " ", round(exp(EM.meta[, 7]), 2), ")", ifelse(EM.meta[, 3] > 0 | EM.meta[, 7] < 0, "*", " "))
    CrI.meta.nf <- paste0("(", EM.meta.nf[, 2], ",", " ", EM.meta.nf[, 3], ")", ifelse(EM.meta.nf[, 2] > 1 | EM.meta.nf[, 3] < 1, "*", " "))
    CrI.beta <- paste0("(", round(exp(beta[, 3]), 2), ",", " ", round(exp(beta[, 7]), 2), ")", ifelse(beta[, 3] > 0 | beta[, 7] < 0, "*", " "))

    EM.both.models <- na.omit(data.frame(drug.names.sorted, round(exp(EM.full[, 1]), 2), CrI.full, round(exp(EM.meta[, 1]), 2), CrI.meta,
                                         EM.meta.nf[, 1], CrI.meta.nf, round(exp(beta[, 1]), 2), CrI.beta))
    colnames(EM.both.models) <- c("Comparison", "Poster. mean NMA", "95% CrI NMA", "Poster. mean NMReg (ref)", "95% CrI NMReg (ref)",
                                  "Poster. mean NMReg (non-ref)", "95% CrI NMReg (non-ref)", "Poster. mean beta", "95% CrI beta")
    rownames(EM.both.models) <- NULL


  }



  ## The 95% CrIs of the between-trial standard deviation under NMA and meta-regression
  CrI.tau.full <- paste0("(", tau.full[, 3], ",", " ", tau.full[, 7], ")")
  CrI.tau.meta <- paste0("(", tau.meta[, 3], ",", " ", tau.meta[, 7], ")")



  ## A data-frame on the model assessment results and between-trial standard deviation under NMA & meta-regression
  table.model.assess <- data.frame(c("Network meta-analysis", "Meta-regression"), rbind(model.assess.NMA[c(1, 3, 2)], model.assess.meta[c(1, 3, 2)]), rbind(cbind(tau.full[, 1], tau.full[, 2], CrI.tau.full), cbind(tau.meta[, 1], tau.meta[, 2], CrI.tau.meta)))
  colnames(table.model.assess) <- c("Analysis", "DIC", "Post. mean Dev.", "pD", "Post. median tau", "Post, SD tau", "95% CrI tau")



  ## Create a data-frame with the SUCRA values
  prepare.sucra <- data.frame(rep(length(drug.names):1, 2), rep(drug.names.sorted, 2), round(rbind(sucra.full[, c(1, 3, 7)], sucra.meta[, c(1, 3, 7)]), 2), rep(c("Unadjusted", "Meta-regression"), each = length(drug.names)))
  colnames(prepare.sucra) <- c("order", "intervention", "mean", "lower", "upper", "analysis")
  rownames(prepare.sucra) <- NULL



  ## Prepare data for ggplot2 (forest-plot)
  if(!is.character(covar.values[[1]]) & (missing(expon) || expon == F)) {

    ## Create a data-frame with effect estimates of basic parameters before and after adjustment
    prepare.EM.metric <- data.frame(rep(length(drug.names):1, 4), rep(drug.names.sorted, 4), round(rbind(EM.full[, c(1, 3, 7)], EM.meta[, c(1, 3, 7)], EM.meta.s, EM.meta.l), 2),
                                    rep(c("Unadjusted", "At the mean value", paste("At", round(covar.values[[1]], 2), "(below the mean value)"), paste("At", round(covar.values[[2]], 2), "(above the mean value)")), each = length(drug.names)))
    colnames(prepare.EM.metric) <- c("order", "comparison", "mean", "lower", "upper", "analysis")
    rownames(prepare.EM.metric) <- NULL


    ## Arguments for ggplot2 that depend on 'expon'
    intercept.value <- 0; trans.value <- "identity"


  } else if(!is.character(covar.values[[1]]) & expon == T) {


    ## Create a data-frame with effect estimates of basic parameters before and after adjustment
    prepare.EM.metric <- data.frame(rep(length(drug.names):1, 4), rep(drug.names.sorted, 4), round(rbind(exp(EM.full[, c(1, 3, 7)]), exp(EM.meta[, c(1, 3, 7)]), EM.meta.s, EM.meta.l), 2),
                                    rep(c("Unadjusted", "At the mean value", paste("At", round(covar.values[[1]], 2), "(below the mean value)"), paste("At", round(covar.values[[2]], 2), "(above the mean value)")), each = length(drug.names)))
    colnames(prepare.EM.metric) <- c("order", "comparison", "mean", "lower", "upper", "analysis")
    rownames(prepare.EM.metric) <- NULL


    ## Arguments for ggplot2 that depend on 'expon'
    intercept.value <- 1; trans.value <- "log10"

  } else if(is.character(covar.values[[1]]) & (missing(expon) || expon == F)) {


    ## Create a data-frame with effect estimates of basic parameters before and after adjustment
    prepare.EM.nominal <- data.frame(rep(length(drug.names):1, 3), rep(drug.names.sorted, 3), round(rbind(EM.full[, c(1, 3, 7)], EM.meta[, c(1, 3, 7)], EM.meta.nf), 2),
                                     rep(c("Unadjusted", paste(covar.values[[1]], "(reference level)"), paste(covar.values[[2]], "(non-reference level)")), each = length(drug.names)))
    colnames(prepare.EM.nominal) <- c("order", "comparison", "mean", "lower", "upper", "analysis")
    rownames(prepare.EM.nominal) <- NULL


    ## Arguments for ggplot2 that depend on 'expon'
    intercept.value <- 0; trans.value <- "identity"


  } else if(is.character(covar.values[[1]]) & expon == T) {


    ## Create a data-frame with effect estimates of basic parameters before and after adjustment
    prepare.EM.nominal <- data.frame(rep(length(drug.names):1, 3), rep(drug.names.sorted, 3), round(rbind(exp(EM.full[, c(1, 3, 7)]), exp(EM.meta[, c(1, 3, 7)]), EM.meta.nf), 2),
                                     rep(c("Unadjusted", paste(covar.values[[1]], "(reference level)"), paste(covar.values[[2]], "(non-reference level)")), each = length(drug.names)))
    colnames(prepare.EM.nominal) <- c("order", "comparison", "mean", "lower", "upper", "analysis")
    rownames(prepare.EM.nominal) <- NULL


    ## Arguments for ggplot2 that depend on 'expon'
    intercept.value <- 1; trans.value <- "log10"

  }



  ## Forest plots of reference-comparisons on effect estimate
  if(!is.character(covar.values[[1]])) {

    p1 <- ggplot(data = prepare.EM.metric, aes(x = as.factor(order), y = as.numeric(mean), ymin = as.numeric(lower), ymax = as.numeric(upper), colour = analysis, group = analysis)) +
            geom_linerange(size = 2, position = position_dodge(width = 0.5)) +
            geom_hline(yintercept = intercept.value, lty = 2, size = 1.3, col = "grey53") +
            geom_point(size = 1.5,  colour = "white", stroke = 0.3, position = position_dodge(width = 0.5)) +
            geom_text(aes(x = as.factor(order), y = round(as.numeric(mean), 2), label = round(as.numeric(mean), 2), hjust = 0, vjust = -0.4),
                     color = "black", size = 4.0, check_overlap = F, parse = F, position = position_dodge(width = 0.5), inherit.aes = T) +
            labs(x = "", y = "Effect estimate", colour = "Analysis") +
            scale_x_discrete(breaks = as.factor(1:length(drug.names.sorted)), labels = drug.names.sorted[length(drug.names.sorted):1]) +
            scale_color_manual(breaks = c("Unadjusted", "At the mean value", paste("At", round(covar.values[[1]], 2), "(below the mean value)"),
                                          paste("At", round(covar.values[[2]], 2), "(above the mean value)")), values = c("#009E73", "#D55E00", "#56B4E9", "#CC79A7")) +
            geom_label(aes(x = unique(as.factor(order)[is.na(mean)]), y = intercept.value, hjust = 0, vjust = 1, label = "Reference intervention"),
                       fill = "beige", colour = "black", fontface = "plain", size = 4) +
            scale_y_continuous(trans = trans.value) +
            coord_flip() +
            theme_classic() +
            theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12),
                  axis.title.x = element_text(color = "black", face = "bold", size = 12), legend.position = "none",
                  legend.text =  element_text(color = "black", size = 12), legend.title =  element_text(color = "black", face = "bold", size = 12))

  } else {

    p1 <- ggplot(data = prepare.EM.nominal, aes(x = as.factor(order), y = as.numeric(mean), ymin = as.numeric(lower), ymax = as.numeric(upper), colour = analysis, group = analysis)) +
            geom_linerange(size = 2, position = position_dodge(width = 0.5)) +
            geom_hline(yintercept = intercept.value, lty = 2, size = 1.3, col = "grey53") +
            geom_point(size = 1.5,  colour = "white", stroke = 0.3, position = position_dodge(width = 0.5)) +
            geom_text(aes(x = as.factor(order), y = round(as.numeric(mean), 2), label = round(as.numeric(mean), 2), hjust = 0, vjust = -0.4),
                      color = "black", size = 4.0, check_overlap = F, parse = F, position = position_dodge(width = 0.5), inherit.aes = T) +
            labs(x = "", y = "Effect estimate", colour = "Analysis") +
            scale_x_discrete(breaks = as.factor(1:length(drug.names.sorted)), labels = drug.names.sorted[length(drug.names.sorted):1]) +
            scale_color_manual(breaks = c("Unadjusted", paste(covar.values[[1]], "(reference level)"), paste(covar.values[[2]], "(non-reference level)")),
                                          values = c("#009E73", "#D55E00", "#56B4E9")) +
            geom_label(aes(x = unique(as.factor(order)[is.na(mean)]), y = intercept.value, hjust = 0, vjust = 1, label = "Reference intervention"),
                       fill = "beige", colour = "black", fontface = "plain", size = 4) +
            scale_y_continuous(trans = trans.value) +
            coord_flip() +
            theme_classic() +
            theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12),
                  axis.title.x = element_text(color = "black", face = "bold", size = 12), legend.position = "none",
                  legend.text =  element_text(color = "black", size = 12), legend.title =  element_text(color = "black", face = "bold", size = 12))

  }



  ## Forest plots of SUCRA under NMA and meta-regression
  p2 <-  ggplot(data = prepare.sucra, aes(x = as.factor(order), y = as.numeric(mean), ymin = as.numeric(lower), ymax = as.numeric(upper), colour = analysis, group = analysis)) +
           geom_linerange(size = 2, position = position_dodge(width = 0.5)) +
           geom_point(size = 1.5,  colour = "white", stroke = 0.3, position = position_dodge(width = 0.5)) +
           geom_text(aes(x = as.factor(order), y = round(as.numeric(mean), 2), label = round(as.numeric(mean)*100, 0), hjust = 0, vjust = -0.4),
                     color = "black", size = 4.0, check_overlap = F, parse = F, position = position_dodge(width = 0.5), inherit.aes = T) +
           labs(x = "", y = "Surface under the cumulative ranking curve value", colour = "Analysis") +
           scale_x_discrete(breaks = as.factor(1:length(drug.names.sorted)), labels = drug.names.sorted[length(drug.names.sorted):1]) +
           scale_color_manual(breaks = c("Unadjusted", "Meta-regression"), values = c("#009E73", "black")) +
           scale_y_continuous(labels = scales::percent) +
           coord_flip(clip = "off") +
           theme_classic() +
           theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12),
                 axis.title.x = element_text(color = "black", face = "bold", size = 12), legend.position = "none",
                 legend.text =  element_text(color = "black", size = 12), legend.title =  element_text(color = "black", face = "bold", size = 12))



  ## Bring together both forest-plots
  forest.plots <- ggarrange(p1, p2, nrow = 1, ncol = 2, labels = c("A)", "B)"), common.legend = F, legend = "bottom")



  ## Prepare data-frame for scatter-plots/ error-bars (ggplot2)
  if (missing(expon) || expon == F) {

    cov.values <- rep(covariate, each = length(drug.names))
    EM.meta.mean <- EM.meta[, 1] + beta[, 1] %o% covariate
    EM.meta.lower <- EM.meta[, 3] + beta[, 3] %o% covariate
    EM.meta.upper <- EM.meta[, 7] + beta[, 7] %o% covariate
    slope <- rep(beta[, 1], length(covariate))
    intercept <- rep(EM.meta[, 1], length(covariate))
    prepare <- na.omit(data.frame(paste(rep(drug.names.sorted, length(covariate)), "versus", na.omit(drug.names.sorted[is.na(EM.full)])), cov.values, melt(EM.meta.mean)[, 3],  melt(EM.meta.lower)[, 3],  melt(EM.meta.upper)[, 3], slope, intercept))
    colnames(prepare) <- c("comparison", "covariate", "mean", "lower", "upper", "slope", "intercept")

  } else {

    cov.values <- rep(covariate, each = length(drug.names))
    EM.meta.mean <- round(exp(EM.meta[, 1] + beta[, 1] %o% covariate), 2)
    EM.meta.lower <- round(exp(EM.meta[, 3] + beta[, 3] %o% covariate), 2)
    EM.meta.upper <- round(exp(EM.meta[, 7] + beta[, 7] %o% covariate), 2)
    slope <- exp(rep(beta[, 1], length(covariate)))
    intercept <- exp(rep(EM.meta[, 1], length(covariate)))
    prepare <- na.omit(data.frame(paste(rep(drug.names.sorted, length(covariate)), "versus", na.omit(drug.names.sorted[is.na(EM.full)])), cov.values, melt(EM.meta.mean)[, 3],  melt(EM.meta.lower)[, 3],  melt(EM.meta.upper)[, 3], slope, intercept))
    colnames(prepare) <- c("comparison", "covariate", "mean", "lower", "upper", "slope", "intercept")

  }


  # Keep nodes with statistically significant inconsistency OR with inconsistent sign in the direct and indirect estimate
  sucra.thres <- rep(sucra.meta[-1, 1][order(sucra.full[-1, 1], decreasing = T)], length(covariate))
  selection <- subset(prepare, sucra.thres > 0.20)



  ## Create the panel of scatterplots on the regression coefficient of the basic parameters
  if(!is.character(covar.values[[1]]) & length(drug.names) <= 16) {

    p3 <- ggplot(data = prepare, aes(x = covariate, y = round(mean, 2), ymin = lower, ymax = upper) ) +
            geom_errorbar(colour = "#D55E00", size = 1, position = position_dodge(width = 0.5), width = 0.1) +
            geom_point(size = 1.5,  colour = "black") +
            geom_line(aes(y = mean),  size = 1, color = "#009E73") +
            geom_hline(yintercept = intercept.value, lty = 2, size = 1, col = "grey") +
            #geom_label(aes(x = min(covariate), y = 0, hjust = 0, vjust = 1, label = paste(intercept, ifelse(slope > 0, "+", "-"), abs(slope), "*covariate")), fill = "beige", colour = "black", fontface = "plain", size = 3.1) +
            facet_wrap(vars(factor(comparison, levels = unique(prepare$comparison))), scales = "free_y") +
            labs(x = covar.values[[3]], y = "") +
            scale_y_continuous(trans = trans.value) +
            theme_classic() +
            theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12), legend.position = "none",
                  axis.title.x = element_text(color = "black", face = "bold", size = 12), strip.text = element_text(size = 11))


  } else if(!is.character(covar.values[[1]]) & length(drug.names) > 16) {


    p3 <- ggplot(data = selection, aes(x = covariate, y = round(mean, 2), ymin = lower, ymax = upper) ) +
            geom_errorbar(colour = "#D55E00", size = 1, position = position_dodge(width = 0.5), width = 0.1) +
            geom_point(size = 1.5,  colour = "black") +
            geom_line(aes(y = mean), size = 1, color = "#009E73") +
            geom_hline(yintercept = intercept.value, lty = 2, size = 1.5, col = "grey") +
           # geom_label(aes(x = min(covariate), y = 0, hjust = 0, vjust = 1, label = paste(intercept, ifelse(slope > 0, "+", "-"), abs(slope), "*covariate")), fill = "beige", colour = "black", fontface = "plain", size = 3.1) +
            facet_wrap(vars(factor(comparison, levels = unique(prepare$comparison))), scales = "free_y") +
            labs(x = covar.values[[3]], y = "") +
            scale_y_continuous(trans = trans.value) +
            theme_classic() +
            theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12), legend.position = "none",
                  axis.title.x = element_text(color = "black", face = "bold", size = 12), strip.text = element_text(size = 11))

  } else if(is.character(covar.values[[1]]) & length(drug.names) <= 16) {

    p3 <- ggplot(data = prepare, aes(x = as.factor(covariate), y = mean, ymin = lower, ymax = upper, group = comparison)) +
            geom_errorbar(size = 2, position = position_dodge(width = 0.1), width = 0.1, colour = "#009E73") +
            geom_point(size = 1.5,  colour = "white", stroke = 0.3, position = position_dodge(width = 0.1)) +
            geom_line(size = 1, color = "#D55E00") +
            geom_hline(yintercept = intercept.value, lty = 2, size = 1, col = "grey") +
            geom_text(aes(x = as.factor(covariate), y = round(as.numeric(mean), 2), label = round(as.numeric(mean), 2)),
                      color = "black", hjust = -0.2, vjust = -0.5, size = 4.0, check_overlap = F, parse = F, position = position_dodge(width = 0.8), inherit.aes = T) +
            facet_wrap(vars(factor(comparison, levels = unique(prepare$comparison))), scales = "free_y") +
            labs(x = covar.values[[3]], y = "") +
            scale_x_discrete(breaks = as.factor(c(0, 1)), labels = covar.values) +
            scale_y_continuous(trans = trans.value) +
            theme_classic() +
            theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12), legend.position = "none",
                  axis.title.x = element_text(color = "black", face = "bold", size = 12), strip.text = element_text(size = 11))



  } else if(is.character(covar.values[[1]]) & length(drug.names) > 16) {

    p3 <- ggplot(data = selection, aes(x = as.factor(covariate), y = mean, ymin = lower, ymax = upper, group = comparison)) +
            geom_errorbar(size = 2, position = position_dodge(width = 0.5), width = 0.1, colour = "#009E73") +
            geom_point(size = 1.5,  colour = "white", stroke = 0.3, position = position_dodge(width = 0.5)) +
            geom_line( size = 1, color = "#D55E00") +
            geom_hline(yintercept = intercept.value, lty = 2, size = 1, col = "grey") +
            geom_text(aes(x = as.factor(covariate), y = round(as.numeric(mean), 2), label = round(as.numeric(mean), 2)),
                      color = "black", hjust = -0.2, vjust = -0.5, size = 4.0, check_overlap = F, parse = F, position = position_dodge(width = 0.8), inherit.aes = T) +
            facet_wrap(vars(factor(comparison, levels = unique(prepare$comparison))), scales = "free_y") +
            labs(x = covar.values[[3]], y = "") +
            scale_x_discrete(breaks = as.factor(c(0, 1)), labels = covar.values) +
            scale_y_continuous(trans = trans.value) +
            theme_classic() +
            theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12), legend.position = "none",
                  axis.title.x = element_text(color = "black", face = "bold", size = 12), strip.text = element_text(size = 11))

  }



  ## Write the table as .xlsx
  write_xlsx(EM.both.models, paste0(getwd(),"Table Meta-regression results.xlsx"))
  write_xlsx(table.model.assess, paste0(getwd(),"Table Meta-regression assessment.xlsx"))



  return(list(EM.both.models = EM.both.models, table.model.assess = table.model.assess, forest.plots = forest.plots, scatterplots = p3))

}

#' A panel of forestplots on the posterior summaries of the observed comparisons in the network under the compared models
#'
#' @param nma A vector with the posterior mean of the treatment effects of the observed comparisons obtained via the \code{UME.plot} function.
#' @param ume A vector with the posterior standard deviation of the treatment effects of the observed comparisons obtained via the \code{UME.plot} function.
#' @param effect.size Text label on the effect measurd considered that appears on the x-axis.
#' @param expon Logical indicating whether the results should appear on the exponential scale. The default which is \code{expon = T}. It is relevant for
#' the odds ratio and ratio of means.
#'
#' @return A panel of forestplots on the posterior summaries of the observed comparisons in the network under the NMA, original and refined UME models.
#' The results refer to the posterior mean and 95% credible interval of the log OR. Grey panels refer to the omitted comparisons.
#' Red and green colours indicate weak and strong evidence, respectively. Namely, the corresponding 95% credible interval includes
#' and excludes the null value, respectively.

forestplot.panel.UME <- function(full, ume, drug.names, effect.size, expon) {



  # Obtain all unique pairwise comparisons using the 'combn' functions
  poss.pair.comp <- data.frame(t(combn(1:length(drug.names), 2))[, 2], t(combn(1:length(drug.names), 2))[, 1])
  colnames(poss.pair.comp) <- c("treat1", "treat2")
  poss.pair.comp$comp <- paste0(poss.pair.comp[, 1], "vs", poss.pair.comp[, 2])
  poss.pair.comp.clean <- poss.pair.comp[is.element(poss.pair.comp$comp, ume$obs.comp), ]


  ## Replace intervention id with their original name
  # For treat1 (non-baseline arm)
  for(i in sort(unique(unlist(poss.pair.comp.clean[, 1])))) {
    poss.pair.comp.clean[poss.pair.comp.clean$treat1 == i, 1] <- drug.names[i]
  }

  # For treat2 (baseline arm)
  for(i in sort(unique(unlist(poss.pair.comp.clean[, 2])))) {
    poss.pair.comp.clean[poss.pair.comp.clean$treat2 == i, 2] <- drug.names[i]
  }


  # Create a vector with the comparisons (non-baseline versus beseline) using the 'paste' function
  comparison <- paste(poss.pair.comp.clean[, 1], "vs", poss.pair.comp.clean[, 2])


  # Indicate frail comparisons
  #frail.comp.ind <- ifelse(is.element(comparison, ume$frail.comp), "yes", "no")
  frail.comp.ind <- ifelse(is.element(poss.pair.comp.clean$comp, ume$frail.comp), "yes", "no")


  # Obtain posterior mean and posterior standard deviation of observed comparisons from NMA and UME model
  EM.full <- full$EM
  ume.mean <- ume$EM[, 1]
  ume.sd <- ume$EM[, 2]


  # Calculate bounds of 95% CrI
  ume.lower <- ume.mean - 1.96*ume.sd
  ume.upper <- ume.mean + 1.96*ume.sd


  ## Keep only the effect estimates according to the 'poss.pair.comp.clean' - Consistency model
  nma.mean <- round(EM.full[is.element(poss.pair.comp$comp, ume$obs.comp), 1], 2)
  nma.lower <- round(EM.full[is.element(poss.pair.comp$comp, ume$obs.comp), 3], 2)
  nma.upper <- round(EM.full[is.element(poss.pair.comp$comp, ume$obs.comp), 7], 2)


  # Indicate statistical significance
  nma.stat.signif <- ifelse(nma.lower > 0 | nma.upper < 0, "strong", "weak")
  ume.stat.signif <- ifelse(ume.lower > 0 | ume.upper < 0, "strong", "weak")

  # Create the dataframe
  data.set <- data.frame(c(nma.mean, ume.mean), c(nma.lower, ume.lower), c(nma.upper, ume.upper), c(nma.stat.signif, ume.stat.signif),
                         rep(comparison, 2), rep(c("NMA model", "UME model"), each = length(ume$obs.comp)), frail.comp.ind)
  colnames(data.set) <- c("mean", "lower", "upper", "stat.sign", "comp", "analysis", "frail")


  ## Obtain forestplot
  if (missing(expon) || expon == F) {
    ggplot(data = data.set, aes(x = as.factor(analysis), y = mean, ymin = lower, ymax = upper, colour = stat.sign)) +
      geom_rect(aes(fill = frail),xmin = -Inf,xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.2) +
      geom_linerange(size = 2, position = position_dodge(width = 0.5)) +
      geom_hline(yintercept = 0, lty = 2, size = 1, col = "black") +
      geom_point(size = 1.5,  colour = "white", stroke = 0.3, position = position_dodge(width = 0.5)) +
      geom_text(aes(x = as.factor(analysis), y = mean, label = mean), color = "black", hjust = 0.35, vjust = -0.25, size = 3.3, check_overlap = F, parse = F,
                position = position_dodge(width = 0.8), inherit.aes = T) +
      facet_wrap(vars(factor(comp, levels = unique(data.set$comp))), scales = "free_x") +
      scale_fill_manual(breaks = c("yes", "no"), values = c("grey53", "white")) +
      scale_color_manual(breaks = c("strong", "weak"), values = c("#009E73", "#D55E00")) +
      labs(x = "", y = effect.size, colour = "Evidence", fill = "") +
      coord_flip() +
      theme_classic() +
      guides(fill = FALSE) +
      theme(axis.text.x = element_text(color = "black", size = 11), axis.text.y = element_text(color = "black", size = 11),
            axis.title.x = element_text(color = "black", face = "bold", size = 11), strip.text = element_text(size = 11),
            legend.position = "bottom", legend.text = element_text(color = "black", size = 11),
            legend.title = element_text(color = "black", face = "bold", size = 11))
  } else {
    ggplot(data = data.set, aes(x = as.factor(analysis), y = round(exp(mean), 2), ymin = round(exp(lower), 2), ymax = round(exp(upper), 2), colour = stat.sign)) +
      geom_rect(aes(fill = frail),xmin = -Inf,xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.2) +
      geom_linerange(size = 2, position = position_dodge(width = 0.5)) +
      geom_hline(yintercept = 1, lty = 2, size = 1, col = "black") +
      geom_point(size = 1.5,  colour = "white", stroke = 0.3, position = position_dodge(width = 0.5)) +
      geom_text(aes(x = as.factor(analysis), y = round(exp(mean), 2), label = round(exp(mean), 2)), color = "black", hjust = 0.35, vjust = -0.25, size = 3.3, check_overlap = F, parse = F,
                position = position_dodge(width = 0.8), inherit.aes = T) +
      facet_wrap(vars(factor(comp, levels = unique(data.set$comp))), scales = "free_x") +
      scale_fill_manual(breaks = c("yes", "no"), values = c("grey53", "white")) +
      scale_color_manual(breaks = c("strong", "weak"), values = c("#009E73", "#D55E00")) +
      scale_y_continuous(trans = "log10") +
      labs(x = "", y = effect.size, colour = "Evidence", fill = "") +
      coord_flip() +
      theme_classic() +
      guides(fill = FALSE) +
      theme(axis.text.x = element_text(color = "black", size = 11), axis.text.y = element_text(color = "black", size = 11),
            axis.title.x = element_text(color = "black", face = "bold", size = 11), strip.text = element_text(size = 11),
            legend.position = "bottom", legend.text = element_text(color = "black", size = 11),
            legend.title = element_text(color = "black", face = "bold", size = 11))
  }


}

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

forestplot.panel.UME <- function(full, ume, drug.names) {


  EM.full <- full$EM
  EM.ume <- ume$EM
  obs.comp <- ume$obs.comp
  frail.comp <- ume$frail.comp
  measure <- effect.measure.name(full$measure)

  ## Possible and observed comparisons
  possible.comp <- possible.observed.frail.comparisons(drug.names, obs.comp)


  # Indicate frail comparisons
  frail.comp.ind <- ifelse(is.element(possible.comp$obs.comp[, 3], frail.comp), "yes", "no")


  # Calculate bounds of 95% CrI
  if (measure != "OR" & measure != "ROM") {
    ume.lower <- EM.ume[, 1] - 1.96*EM.ume[, 2]
    ume.upper <- EM.ume[, 1] + 1.96*EM.ume[, 2]
  } else {
    ume.lower <- exp(EM.ume[, 1] - 1.96*EM.ume[, 2])
    ume.upper <- exp(EM.ume[, 1] + 1.96*EM.ume[, 2])
  }



  ## Keep only the effect estimates according to the 'poss.pair.comp.clean' - Consistency model
  if (measure != "OR" & measure != "ROM") {
    ume.mean <- round(EM.ume[, 1], 2)
    nma.mean <- round(EM.full[is.element(possible.comp$poss.comp[, 4], obs.comp), 1], 2)
    nma.lower <- round(EM.full[is.element(possible.comp$poss.comp[, 4], obs.comp), 3], 2)
    nma.upper <- round(EM.full[is.element(possible.comp$poss.comp[, 4], obs.comp), 7], 2)
  } else {
    ume.mean <- round(exp(EM.ume[, 1]), 2)
    nma.mean <- round(exp(EM.full[is.element(possible.comp$poss.comp[, 4], obs.comp), 1]), 2)
    nma.lower <- round(exp(EM.full[is.element(possible.comp$poss.comp[, 4], obs.comp), 3]), 2)
    nma.upper <- round(exp(EM.full[is.element(possible.comp$poss.comp[, 4], obs.comp), 7]), 2)
  }


  # Indicate statistical significance
  if (measure != "OR" & measure != "ROM") {
    nma.stat.signif <- ifelse(nma.lower > 0 | nma.upper < 0, "strong", "weak")
    ume.stat.signif <- ifelse(ume.lower > 0 | ume.upper < 0, "strong", "weak")
  } else {
    nma.stat.signif <- ifelse(nma.lower > 1 | nma.upper < 1, "strong", "weak")
    ume.stat.signif <- ifelse(ume.lower > 1 | ume.upper < 1, "strong", "weak")
  }


  # Create the dataframe
  data.set <- data.frame(c(nma.mean, ume.mean), c(nma.lower, ume.lower), c(nma.upper, ume.upper), c(nma.stat.signif, ume.stat.signif),
                        rep(possible.comp$obs.comp[, 4], 2), rep(c("NMA model", "UME model"), each = length(obs.comp)), frail.comp.ind)
  colnames(data.set) <- c("mean", "lower", "upper", "stat.sign", "comp", "analysis", "frail")


  ## Obtain forestplot
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
    scale_y_continuous(trans = ifelse(measure != "OR" & measure != "ROM", "identity", "log10")) +
    labs(x = "", y = measure, colour = "Evidence", fill = "") +
    coord_flip() +
    theme_classic() +
    guides(fill = FALSE) +
    theme(axis.text.x = element_text(color = "black", size = 11), axis.text.y = element_text(color = "black", size = 11),
          axis.title.x = element_text(color = "black", face = "bold", size = 11), strip.text = element_text(size = 11),
          legend.position = "bottom", legend.text = element_text(color = "black", size = 11),
          legend.title = element_text(color = "black", face = "bold", size = 11))

}
